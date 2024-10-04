# RUN EPIDEMIC TREE ALGORITHM
# Run using Serial Interval (SI) and Distance Kernel (DK) distributions
# Use convolved distributions for wildlife/ dogs of unknown origin 

# Packages
library(treerabid) 
library(data.table)
library(lubridate)
library(dplyr)
library(lubridate)
library(magrittr)
library(foreach)
library(iterators)
library(doRNG)
library(igraph)
library(glue)

## Distance kernel and serial interval parameters
DK_params <- read.csv("output/DK_params.csv")
SI_params <-  read.csv("output/SI_params.csv")

# Update important SI and DK parameters from Mancy paper (note they haven't
# changed a huge amount on adding data to end of 2022)
params_treerabid$SI_meanlog <- SI_params$SI_ml
params_treerabid$SI_sdlog <- SI_params$SI_sdlog
params_treerabid$DK_shape_weibull <- filter(DK_params,kernel_type=="DK" & dist=="weibull")$par1est
params_treerabid$DK_scale_weibull <- filter(DK_params,kernel_type=="DK" & dist=="weibull")$par2est
params_treerabid$DK2_shape_weibull <- filter(DK_params,kernel_type=="DK2" & dist=="weibull")$par1est
params_treerabid$DK2_scale_weibull <- filter(DK_params,kernel_type=="DK2" & dist=="weibull")$par2est

# clean up (no cases with NA location or time & filter to start/end dates)
case_dt <- readRDS(file = "output/clean_bite_data.rda")

case_dt %<>%
  dplyr::filter(!is.na(Symptoms.started), 
                !is.na(UTM.Easting), 
                !is.na(UTM.Northing), 
                Symptoms.started >= ymd("2002-01-01"),
                Symptoms.started >= "2002-01-01", 
                Symptoms.started <= ymd("2022-12-31")) %>%
  # get uncertainty in days
  mutate(days_uncertain = case_when(Symptoms.started.accuracy == "+/- 14 days" ~ 14L, 
                                    Symptoms.started.accuracy == "+/- 7 days" ~ 7L,
                                    Symptoms.started.accuracy == "+/- 28 days" ~ 28L, 
                                    Symptoms.started.accuracy == "0" ~ 0L, 
                                    TRUE ~ 0L), 
         owned = ifelse(Owner %in% "Known", TRUE, FALSE)) 

# Table 1: agreement to CT data (use known tree for each type of cutoff situation)
ct_data <- data.table(id_case = case_dt$ID,
                      id_biter = case_dt$Biter.ID, 
                      x_coord = case_dt$UTM.Easting,
                      y_coord = case_dt$UTM.Northing,
                      owned = case_dt$owned, 
                      date_symptoms = case_dt$Symptoms.started,
                      days_uncertain = case_dt$days_uncertain)
ct_data <- ct_data[ct_data[, .I[which.max(date_symptoms)],  by = "id_case"]$V1] # pick one case only per id
fwrite(ct_data, "Output/tree_ct_data.csv")

# Use the `best` dists/cutoffs & known source to generate trees + incs ----
pars_selected <-
  tidyr::expand_grid(si_pdist = "lnorm", 
                     dist_pdist = "weibull",
                     convolve = "mixed",
                     prune = TRUE, 
                     cutoff = c(0.95, 0.975, 0.99), 
                     use_known = TRUE, 
                     nsim = 1000)
pars_selected$seed <- 45 * 1:nrow(pars_selected)

comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY = FALSE, fill = TRUE)
}

cl <- parallel::makeCluster(parallel::detectCores() - 1)
doParallel::registerDoParallel(cl)

system.time({
  consensus_known <-
    foreach(i = iter(pars_selected, by = "row"), .combine = comb) %do% {
      
      # ct data with dates w/out uncertainty & only one record per ID
      case_dates <- data.table(id_case = ct_data$id_case, 
                               symptoms_started = ct_data$date_symptoms)
      
      # Distribution functions from treerabid
      if(i$convolve == "mixed") {
        si_pdist <- get(paste0("si_", i$si_pdist, "1"))
        dist_pdist <- get(paste0("dist_", i$dist_pdist, "_mixed"))
      }
      
      if(i$convolve == "convolved") {
        si_pdist <- get(paste0("si_", i$si_pdist, "2"))
        dist_pdist <- get(paste0("dist_", i$dist_pdist, "2"))
      }
      
      if(i$convolve == "baseline") {
        si_pdist <- get(paste0("si_", i$si_pdist, "1"))
        dist_pdist <- get(paste0("dist_", i$dist_pdist, "1"))
      }
      
      ttrees <- 
        boot_trees(id_case = case_dt$ID,
                   id_biter = case_dt$Biter.ID, 
                   x_coord = case_dt$UTM.Easting,
                   y_coord = case_dt$UTM.Northing,
                   owned = case_dt$owned, 
                   date_symptoms = case_dt$Symptoms.started,
                   days_uncertain = case_dt$days_uncertain,
                   use_known_source = i$use_known,
                   prune = i$prune,
                   si_fun = si_pdist,
                   dist_fun = dist_pdist, 
                   params = params_treerabid, 
                   cutoff = i$cutoff,
                   N = i$nsim, 
                   seed = i$seed) 
      
      # Summarize the trees
      # do this outside of function to get min t_diff as well
      links_all <- ttrees[, .(links = .N,
                              t_diff_min_days = min(t_diff),
                              t_diff_median_days = median(t_diff),
                              dist_diff_meters = median(dist_diff)),
                          by = c("id_case", "id_progen")][, prob := links/i$nsim]
      links_consensus <- build_consensus_links(links_all, case_dates)
      tree_ids <- c(mcc = 
                      build_consensus_tree(links_consensus, ttrees, links_all,
                                           type = "mcc", output = "sim"), 
                    majority = 
                      build_consensus_tree(links_consensus, ttrees, links_all,
                                           type = "majority", output = "sim"))
      ttrees$mcc <- ifelse(ttrees$sim %in% tree_ids["mcc"], 1, 0)
      ttrees$majority <- ifelse(ttrees$sim %in% tree_ids["majority"], 1, 0)
      set.seed(i$seed)
      out_trees <- ttrees[sim %in% c(sample((1:i$nsim)[-tree_ids], 100), tree_ids)]
      
      list(links = cbind(links_consensus, i), 
           ttrees_all = data.table(out_trees, cutoff = i$cutoff))
    }
})

parallel::stopCluster(cl)

# Write out files
consensus_links <- consensus_known$links
lapply(split(consensus_links, interaction(consensus_links$cutoff, consensus_links$convolve)), check_loops)
# fwrite(consensus_known$ttrees_all, "output/trees_sampled_best.gz")
# fwrite(consensus_links, "output/consensus_links_best.csv")

# Write out incursions 
incs_best <- consensus_links[is.na(id_progen)]
incursions_95 <- filter(ct_data,id_case %in% filter(incs_best,cutoff==0.95)$id_case)
incursions_975 <- filter(ct_data,id_case %in% filter(incs_best,cutoff==0.975)$id_case)
incursions_99 <- filter(ct_data,id_case %in% filter(incs_best,cutoff==0.99)$id_case)
write.csv(incursions_95,file=paste0("output/serengeti_incursions_treerabid_prune_95.csv"))
write.csv(incursions_975,file=paste0("output/serengeti_incursions_treerabid_prune_975.csv"))
write.csv(incursions_99,file=paste0("output/serengeti_incursions_treerabid_prune_99.csv"))

