rm(list=ls())

set.seed(0)



# Read in data & models
#___________________

library(brms)

model <- readRDS("output/stan_models/incidence_from_vax_model_village_stan.rds")
model_wo_cases <- readRDS("output/stan_models/incidence_from_vax_model_village_stan_wo_cases.rds")

data_vill <- read.csv("Output/incidence_coverage_model_data_village.csv")
data_dist <- read.csv("Output/incidence_coverage_model_data_district.csv")
data_vill$susc_last2monthMean <- 1-data_vill$vax_last2monthMean
vax_vill <- as.matrix(read.csv("Output/vaccinationCoverageByVillageMonth_Jan2002_Dec2022.csv",header = F,row.names = 1)) 
data_vill_case_rate_adjust_dist <- 0.5*min(data_vill$case_rate_last2monthMean_dist[which(data_vill$case_rate_last2monthMean_dist>0)])


load("Output/neighbour_notNeighbour_susceptibilities.Rdata")
load("output/village_borders.Rdata")


source("R/Functions/Lmeans.R")

case_cap <- max(data_vill$incidence)

simulate <- F


if(simulate==T){
  
  n_samples <-200
  
  
  
  # Sample parameters for simulation
  #___________________
  
  samples_pars <- posterior_samples(model, pars = c("Intercept","beta","p","phi"))
  samples <- sample(1:nrow(samples_pars),n_samples)
  samples_pars <- samples_pars[samples,]
  samples_reffs <- posterior_samples(model, pars = c("gamma_t"))
  samples_reffs <- samples_reffs[samples,]
  samples_pars <- samples_pars[,-which(colnames(samples_pars)=="beta[1]")]
  samples_sigma_village <- posterior_samples(model, pars = c("sigma_village"))
  samples_sigma_village <- samples_sigma_village[samples,]
  
  samples_pars_wo_cases <- posterior_samples(model_wo_cases, pars = c("Intercept","beta","p","phi"))
  samples_pars_wo_cases <- samples_pars_wo_cases[samples,]
  samples_pars_wo_cases <- samples_pars_wo_cases[,-which(colnames(samples_pars_wo_cases)=="beta[1]")]
  samples_reffs_wo_cases <- posterior_samples(model_wo_cases, pars = c("gamma_t"))
  samples_reffs_wo_cases <- samples_reffs_wo_cases[samples,]
  samples_sigma_village_wo_cases <- posterior_samples(model_wo_cases, pars = c("sigma_village"))
  samples_sigma_village_wo_cases <- samples_sigma_village_wo_cases[samples,]
  
  rm(model,model_wo_cases)
  

    
  # Simulate cases using case data for only the first two months 
  #___________________
  
  # Matrix for results
  sim_mat <-  sim_mat_cap <-  matrix(NA,nrow=nrow(data_vill),ncol=nrow(samples_pars))
  
  # Design matrices
  X <- as.matrix(cbind(1,data_vill[,c("susc_last2monthMean","log_case_rate_last2monthMean","log_case_rate_neighbours_last2monthMean","log_case_rate_notNeighbours_last2monthMean","log_dog_density","HDR")],
                       "power_mean_neighbours_last2MonthMean"=NA,"power_mean_notNeighbours_last2MonthMean"=NA,"village_re"=NA,log(data_vill$dogs)))
  power_mean_neighbours <- power_mean_notNeighbours <- matrix(NA, nrow=nrow(vax_vill),ncol=ncol(vax_vill))
  for(s in 1:nrow(samples_pars)){
    pars_s <- samples_pars[s,]
    re_s <- samples_reffs[s,]

    for(v in 1:nrow(vax_vill)){
      power_mean_neighbours[v,] <- sapply(1:ncol(vax_vill),function(x) powMean(x=S_n[[v]][,x],p=as.numeric(pars_s["p"]),wts = W_n[[v]][,x]))
      power_mean_notNeighbours[v,] <- sapply(1:ncol(vax_vill),function(x) powMean(x=S_nn[[v]][,x],p=as.numeric(pars_s["p"]),wts = W_nn[[v]][,x]))
    }
    X[,"power_mean_neighbours_last2MonthMean"] <- (c(cbind(NA,power_mean_neighbours[,-ncol(vax_vill)]))+c(cbind(NA,NA,power_mean_neighbours[,-((ncol(vax_vill)-1):ncol(vax_vill))])))/2
    X[,"power_mean_notNeighbours_last2MonthMean"] <- (c(cbind(NA,power_mean_notNeighbours[,-ncol(vax_vill)]))+c(cbind(NA,NA,power_mean_notNeighbours[,-((ncol(vax_vill)-1):ncol(vax_vill))])))/2
    X[,"village_re"] <- rep(as.numeric(re_s),ncol(vax_vill))
    
    for(m in 3:max(data_vill$month)){
      
      # Monthly design matrices
      m_indices <- which(data_vill$month==m)
      last_m_indices <- which(data_vill$month==(m-1))
      two_m_ago_indices <- which(data_vill$month==(m-2))
      X_m <- X_cap_m <- X[m_indices,]
      if(m>4){
        X_m[,"log_case_rate_last2monthMean"] <- log(((sim_mat[last_m_indices,s]/data_vill$dogs[last_m_indices]+sim_mat[two_m_ago_indices,s]/data_vill$dogs[two_m_ago_indices])/2) + data_vill_case_rate_adjust_dist)
        X_m[,"log_case_rate_notNeighbours_last2monthMean"] <- log(((rowSums(not_bordering*matrix(sim_mat[two_m_ago_indices,s],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering)))/rowSums(not_bordering*matrix(data_vill$dogs[two_m_ago_indices],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering))) +
                                                                        rowSums(not_bordering*matrix(sim_mat[last_m_indices,s],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering)))/rowSums(not_bordering*matrix(data_vill$dogs[last_m_indices],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering))))/2) + data_vill_case_rate_adjust_dist)
        X_cap_m[,"log_case_rate_last2monthMean"] <- log(((sim_mat_cap[last_m_indices,s]/data_vill$dogs[last_m_indices]+sim_mat_cap[two_m_ago_indices,s]/data_vill$dogs[two_m_ago_indices])/2) + data_vill_case_rate_adjust_dist)
        X_cap_m[,"log_case_rate_notNeighbours_last2monthMean"] <- log(((rowSums(not_bordering*matrix(sim_mat_cap[two_m_ago_indices,s],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering)))/rowSums(not_bordering*matrix(data_vill$dogs[two_m_ago_indices],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering))) +
                                                                      rowSums(not_bordering*matrix(sim_mat_cap[last_m_indices,s],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering)))/rowSums(not_bordering*matrix(data_vill$dogs[last_m_indices],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering))))/2) + data_vill_case_rate_adjust_dist)
        for(v in 1:nrow(vax_vill)){
          neighbours <- which(bordering[v,])
          case_rates <- cbind(sim_mat[two_m_ago_indices,s],sim_mat[last_m_indices,s])/cbind(data_vill$dogs[two_m_ago_indices],data_vill$dogs[last_m_indices])
          case_rates_w_park_mara <- rbind(case_rates,colMeans(case_rates[neighbours[which(!neighbours%in%(nrow(X_m)+1:2))],,drop=F]),0.01/12) 
          neighbour_case_rates_last2months <-  case_rates_w_park_mara[neighbours,]
          weights_neighbours <-  matrix(rep(borders_prop[v,neighbours],2),ncol=2) #weight by border length 
          X_m[v,"log_case_rate_neighbours_last2monthMean"] <- log(mean(colSums(neighbour_case_rates_last2months*weights_neighbours)) + data_vill_case_rate_adjust_dist)
          
          case_rates_cap <- cbind(sim_mat_cap[two_m_ago_indices,s],sim_mat_cap[last_m_indices,s])/cbind(data_vill$dogs[two_m_ago_indices],data_vill$dogs[last_m_indices])
          case_rates_w_park_mara_cap <- rbind(case_rates_cap,colMeans(case_rates_cap[neighbours[which(!neighbours%in%(nrow(X_m)+1:2))],,drop=F]),0.01/12) 
          neighbour_case_rates_last2months_cap <-  case_rates_w_park_mara_cap[neighbours,]
          X_cap_m[v,"log_case_rate_neighbours_last2monthMean"] <- log(mean(colSums(neighbour_case_rates_last2months_cap*weights_neighbours)) + data_vill_case_rate_adjust_dist)
          
        }
      }
      
      # Simulate cases
      sim_mat[m_indices,s] <- rnbinom(nrow(X_m),mu=exp(colSums(t(X_m)*c(as.numeric(pars_s[c(1:(ncol(pars_s)-2))]),1,1))),size=as.numeric(pars_s["phi"]))
      sim_mat_cap[m_indices,s] <- pmin(rnbinom(nrow(X_m),mu=exp(colSums(t(X_m)*c(as.numeric(pars_s[c(1:(ncol(pars_s)-2))]),1,1))),size=as.numeric(pars_s["phi"])),round(case_cap*data_vill$dogs[m_indices]))
      
    }
  }
  
  sim_mat_dist <- rowsum((sim_mat),data_vill$month)
  sim_median <- apply(sim_mat_dist,1,quantile,0.5,na.rm=T)[-c(1:2)]
  sim_lower <- apply(sim_mat_dist,1,quantile,0.025,na.rm=T)[-c(1:2)]
  sim_upper <- apply(sim_mat_dist,1,quantile,0.975,na.rm=T)[-c(1:2)]
  sim_mat_dist_cap <- rowsum((sim_mat_cap),data_vill$month)
  sim_median_cap <- apply(sim_mat_dist_cap,1,quantile,0.5,na.rm=T)[-c(1:2)]
  sim_lower_cap <- apply(sim_mat_dist_cap,1,quantile,0.025,na.rm=T)[-c(1:2)]
  sim_upper_cap <- apply(sim_mat_dist_cap,1,quantile,0.975,na.rm=T)[-c(1:2)]
  
  save(sim_mat_dist, sim_lower, sim_upper, sim_median,
       sim_mat_dist_cap, sim_lower_cap, sim_upper_cap, sim_median_cap,
       file="Output/simulated_cases_use_initial_cases.Rdata")
  
  
  
  # Simulate cases using vaccination data only
  #___________________
  
  
  # Matrix for results
  sim_mat <-  sim_mat_cap <-  matrix(NA,nrow=nrow(data_vill),ncol=nrow(samples_pars))
  
  # Design matrices
  X <- as.matrix(cbind(1,data_vill[,c("susc_last2monthMean","log_case_rate_last2monthMean","log_case_rate_neighbours_last2monthMean","log_case_rate_notNeighbours_last2monthMean","log_dog_density","HDR")],
                       "power_mean_neighbours_last2MonthMean"=NA,"power_mean_notNeighbours_last2MonthMean"=NA,"village_re"=NA,log(data_vill$dogs)))
  X_wo_cases <- X[,which(!colnames(X)%in%c("log_case_rate_last2monthMean","log_case_rate_neighbours_last2monthMean","log_case_rate_notNeighbours_last2monthMean"))]
  power_mean_neighbours <- power_mean_notNeighbours <- power_mean_neighbours_wo_cases <- power_mean_notNeighbours_wo_cases <- matrix(NA, nrow=nrow(vax_vill),ncol=ncol(vax_vill))
  for(s in 1:nrow(samples_pars)){
    pars_s <- samples_pars[s,]
    pars_wo_cases_s <- samples_pars_wo_cases[s,]
    re_s <- samples_reffs[s,]
    re_wo_cases_s <- samples_reffs_wo_cases[s,]
    
    for(v in 1:nrow(vax_vill)){
      power_mean_neighbours[v,] <- sapply(1:ncol(vax_vill),function(x) powMean(x=S_n[[v]][,x],p=as.numeric(pars_s["p"]),wts = W_n[[v]][,x]))
      power_mean_notNeighbours[v,] <- sapply(1:ncol(vax_vill),function(x) powMean(x=S_nn[[v]][,x],p=as.numeric(pars_s["p"]),wts = W_nn[[v]][,x]))
      power_mean_neighbours_wo_cases[v,] <- sapply(1:ncol(vax_vill),function(x) powMean(x=S_n[[v]][,x],p=as.numeric(pars_wo_cases_s["p"]),wts = W_n[[v]][,x]))
      power_mean_notNeighbours_wo_cases[v,] <- sapply(1:ncol(vax_vill),function(x) powMean(x=S_nn[[v]][,x],p=as.numeric(pars_wo_cases_s["p"]),wts = W_nn[[v]][,x]))
    }
    X[,"power_mean_neighbours_last2MonthMean"] <- (c(cbind(NA,power_mean_neighbours[,-ncol(vax_vill)]))+c(cbind(NA,NA,power_mean_neighbours[,-((ncol(vax_vill)-1):ncol(vax_vill))])))/2
    X[,"power_mean_notNeighbours_last2MonthMean"] <- (c(cbind(NA,power_mean_notNeighbours[,-ncol(vax_vill)]))+c(cbind(NA,NA,power_mean_notNeighbours[,-((ncol(vax_vill)-1):ncol(vax_vill))])))/2
    X[,"village_re"] <- rep(as.numeric(re_s),ncol(vax_vill))
    X_wo_cases[,"power_mean_neighbours_last2MonthMean"] <- (c(cbind(NA,power_mean_neighbours_wo_cases[,-ncol(vax_vill)]))+c(cbind(NA,NA,power_mean_neighbours_wo_cases[,-((ncol(vax_vill)-1):ncol(vax_vill))])))/2
    X_wo_cases[,"power_mean_notNeighbours_last2MonthMean"] <- (c(cbind(NA,power_mean_notNeighbours_wo_cases[,-ncol(vax_vill)]))+c(cbind(NA,NA,power_mean_notNeighbours_wo_cases[,-((ncol(vax_vill)-1):ncol(vax_vill))])))/2
    X_wo_cases[,"village_re"] <- rep(as.numeric(re_wo_cases_s),ncol(vax_vill))
    
    
    for(m in 3:max(data_vill$month)){
      
      # Monthly design matrices
      m_indices <- which(data_vill$month==m)
      last_m_indices <- which(data_vill$month==(m-1))
      two_m_ago_indices <- which(data_vill$month==(m-2))
      X_m <- X_cap_m <- X[m_indices,]
      X_wo_cases_m <- X_wo_cases_cap_m <- X_wo_cases[m_indices,]
      if(m>4){
        X_m[,"log_case_rate_last2monthMean"] <- log(((sim_mat[last_m_indices,s]/data_vill$dogs[last_m_indices]+sim_mat[two_m_ago_indices,s]/data_vill$dogs[two_m_ago_indices])/2) + data_vill_case_rate_adjust_dist)
        X_m[,"log_case_rate_notNeighbours_last2monthMean"] <- log(((rowSums(not_bordering*matrix(sim_mat[two_m_ago_indices,s],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering)))/rowSums(not_bordering*matrix(data_vill$dogs[two_m_ago_indices],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering))) +
                                                                      rowSums(not_bordering*matrix(sim_mat[last_m_indices,s],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering)))/rowSums(not_bordering*matrix(data_vill$dogs[last_m_indices],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering))))/2) + data_vill_case_rate_adjust_dist)
        X_cap_m[,"log_case_rate_last2monthMean"] <- log(((sim_mat_cap[last_m_indices,s]/data_vill$dogs[last_m_indices]+sim_mat_cap[two_m_ago_indices,s]/data_vill$dogs[two_m_ago_indices])/2) + data_vill_case_rate_adjust_dist)
        X_cap_m[,"log_case_rate_notNeighbours_last2monthMean"] <- log(((rowSums(not_bordering*matrix(sim_mat_cap[two_m_ago_indices,s],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering)))/rowSums(not_bordering*matrix(data_vill$dogs[two_m_ago_indices],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering))) +
                                                                      rowSums(not_bordering*matrix(sim_mat_cap[last_m_indices,s],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering)))/rowSums(not_bordering*matrix(data_vill$dogs[last_m_indices],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering))))/2) + data_vill_case_rate_adjust_dist)
        for(v in 1:nrow(vax_vill)){
          neighbours <- which(bordering[v,])
          case_rates <- cbind(sim_mat[two_m_ago_indices,s],sim_mat[last_m_indices,s])/cbind(data_vill$dogs[two_m_ago_indices],data_vill$dogs[last_m_indices])
          case_rates_w_park_mara <- rbind(case_rates,colMeans(case_rates[neighbours[which(!neighbours%in%(nrow(X_m)+1:2))],,drop=F]),0.01/12) 
          neighbour_case_rates_last2months <-  case_rates_w_park_mara[neighbours,]
          weights_neighbours <-  matrix(rep(borders_prop[v,neighbours],2),ncol=2) #weight by border length 
          X_m[v,"log_case_rate_neighbours_last2monthMean"] <- log(mean(colSums(neighbour_case_rates_last2months*weights_neighbours)) + data_vill_case_rate_adjust_dist)

          case_rates_cap <- cbind(sim_mat_cap[two_m_ago_indices,s],sim_mat_cap[last_m_indices,s])/cbind(data_vill$dogs[two_m_ago_indices],data_vill$dogs[last_m_indices])
          case_rates_w_park_mara_cap <- rbind(case_rates_cap,colMeans(case_rates_cap[neighbours[which(!neighbours%in%(nrow(X_m)+1:2))],,drop=F]),0.01/12) 
          neighbour_case_rates_last2months_cap <-  case_rates_w_park_mara_cap[neighbours,]
          X_cap_m[v,"log_case_rate_neighbours_last2monthMean"] <- log(mean(colSums(neighbour_case_rates_last2months_cap*weights_neighbours)) + data_vill_case_rate_adjust_dist)
          
        }
        
        # Simulate cases
        sim_mat[m_indices,s] <- rnbinom(nrow(X_m),mu=exp(colSums(t(X_m)*c(as.numeric(pars_s[c(1:(ncol(pars_s)-2))]),1,1))),size=as.numeric(pars_s["phi"]))
        sim_mat_cap[m_indices,s] <- pmin(rnbinom(nrow(X_cap_m),mu=exp(colSums(t(X_cap_m)*c(as.numeric(pars_s[c(1:(ncol(pars_s)-2))]),1,1))),size=as.numeric(pars_s["phi"])),round(case_cap*data_vill$dogs[m_indices]))
        
      }else{ # use vaccination-only model
        
        # Simulate cases
        sim_mat[m_indices,s] <- rnbinom(nrow(X_wo_cases_m),mu=exp(colSums(t(X_wo_cases_m)*c(as.numeric(pars_wo_cases_s[c(1:(ncol(pars_wo_cases_s)-2))]),1,1))),size=as.numeric(pars_wo_cases_s["phi"]))
        sim_mat_cap[m_indices,s] <- pmin(rnbinom(nrow(X_wo_cases_m),mu=exp(colSums(t(X_wo_cases_m)*c(as.numeric(pars_wo_cases_s[c(1:(ncol(pars_wo_cases_s)-2))]),1,1))),size=as.numeric(pars_wo_cases_s["phi"])),round(case_cap*data_vill$dogs[m_indices]))
        
      }
      
      
    }
  }
  
  sim_mat_dist <- rowsum((sim_mat),data_vill$month)
  sim_median <- apply(sim_mat_dist,1,quantile,0.5,na.rm=T)[-c(1:2)]
  sim_lower <- apply(sim_mat_dist,1,quantile,0.025,na.rm=T)[-c(1:2)]
  sim_upper <- apply(sim_mat_dist,1,quantile,0.975,na.rm=T)[-c(1:2)]
  sim_mat_dist_cap <- rowsum((sim_mat_cap),data_vill$month)
  sim_median_cap <- apply(sim_mat_dist_cap,1,quantile,0.5,na.rm=T)[-c(1:2)]
  sim_lower_cap <- apply(sim_mat_dist_cap,1,quantile,0.025,na.rm=T)[-c(1:2)]
  sim_upper_cap <- apply(sim_mat_dist_cap,1,quantile,0.975,na.rm=T)[-c(1:2)]
  
  save(sim_mat_dist, sim_lower, sim_upper, sim_median, sim_mat_dist_cap, sim_lower_cap, sim_upper_cap, sim_median_cap,
       file="Output/simulated_cases_vaccination_only.Rdata")
  
  
  
  # Simulate cases with unknown random effects (I think this, with the cap on incidence, is the approach that should be used for simulating Mara reigon)
  #___________________
  
  # Matrix for results
  sim_mat <-  sim_mat_cap <-  matrix(NA,nrow=nrow(data_vill),ncol=nrow(samples_pars))
  
  # Design matrices
  X <- as.matrix(cbind(1,data_vill[,c("susc_last2monthMean","log_case_rate_last2monthMean","log_case_rate_neighbours_last2monthMean","log_case_rate_notNeighbours_last2monthMean","log_dog_density","HDR")],
                       "power_mean_neighbours_last2MonthMean"=NA,"power_mean_notNeighbours_last2MonthMean"=NA,"village_re"=NA,log(data_vill$dogs)))
  X_wo_cases <- X[,which(!colnames(X)%in%c("log_case_rate_last2monthMean","log_case_rate_neighbours_last2monthMean","log_case_rate_notNeighbours_last2monthMean"))]
  power_mean_neighbours <- power_mean_notNeighbours <- power_mean_neighbours_wo_cases <- power_mean_notNeighbours_wo_cases <- matrix(NA, nrow=nrow(vax_vill),ncol=ncol(vax_vill))
  for(s in 1:nrow(samples_pars)){
    pars_s <- samples_pars[s,]
    pars_wo_cases_s <- samples_pars_wo_cases[s,]
    sigma_village_s <- samples_sigma_village[s]
    sigma_village_wo_cases_s <- samples_sigma_village_wo_cases[s]
    re_s <- rnorm(nrow(vax_vill),0,sigma_village_s)
    re_wo_cases_s <- rnorm(nrow(vax_vill),0,sigma_village_wo_cases_s)
    
    for(v in 1:nrow(vax_vill)){
      power_mean_neighbours[v,] <- sapply(1:ncol(vax_vill),function(x) powMean(x=S_n[[v]][,x],p=as.numeric(pars_s["p"]),wts = W_n[[v]][,x]))
      power_mean_notNeighbours[v,] <- sapply(1:ncol(vax_vill),function(x) powMean(x=S_nn[[v]][,x],p=as.numeric(pars_s["p"]),wts = W_nn[[v]][,x]))
      power_mean_neighbours_wo_cases[v,] <- sapply(1:ncol(vax_vill),function(x) powMean(x=S_n[[v]][,x],p=as.numeric(pars_wo_cases_s["p"]),wts = W_n[[v]][,x]))
      power_mean_notNeighbours_wo_cases[v,] <- sapply(1:ncol(vax_vill),function(x) powMean(x=S_nn[[v]][,x],p=as.numeric(pars_wo_cases_s["p"]),wts = W_nn[[v]][,x]))
    }
    X[,"power_mean_neighbours_last2MonthMean"] <- (c(cbind(NA,power_mean_neighbours[,-ncol(vax_vill)]))+c(cbind(NA,NA,power_mean_neighbours[,-((ncol(vax_vill)-1):ncol(vax_vill))])))/2
    X[,"power_mean_notNeighbours_last2MonthMean"] <- (c(cbind(NA,power_mean_notNeighbours[,-ncol(vax_vill)]))+c(cbind(NA,NA,power_mean_notNeighbours[,-((ncol(vax_vill)-1):ncol(vax_vill))])))/2
    X[,"village_re"] <- rep(as.numeric(re_s),ncol(vax_vill))
    X_wo_cases[,"power_mean_neighbours_last2MonthMean"] <- (c(cbind(NA,power_mean_neighbours_wo_cases[,-ncol(vax_vill)]))+c(cbind(NA,NA,power_mean_neighbours_wo_cases[,-((ncol(vax_vill)-1):ncol(vax_vill))])))/2
    X_wo_cases[,"power_mean_notNeighbours_last2MonthMean"] <- (c(cbind(NA,power_mean_notNeighbours_wo_cases[,-ncol(vax_vill)]))+c(cbind(NA,NA,power_mean_notNeighbours_wo_cases[,-((ncol(vax_vill)-1):ncol(vax_vill))])))/2
    X_wo_cases[,"village_re"] <- rep(as.numeric(re_wo_cases_s),ncol(vax_vill))
    
    
    for(m in 3:max(data_vill$month)){
      
      # Monthly design matrices
      m_indices <- which(data_vill$month==m)
      last_m_indices <- which(data_vill$month==(m-1))
      two_m_ago_indices <- which(data_vill$month==(m-2))
      X_m <- X_cap_m <- X[m_indices,]
      X_wo_cases_m <- X_wo_cases_cap_m <- X_wo_cases[m_indices,]
      if(m>4){
        X_m[,"log_case_rate_last2monthMean"] <- log(((sim_mat[last_m_indices,s]/data_vill$dogs[last_m_indices]+sim_mat[two_m_ago_indices,s]/data_vill$dogs[two_m_ago_indices])/2) + data_vill_case_rate_adjust_dist)
        X_m[,"log_case_rate_notNeighbours_last2monthMean"] <- log(((rowSums(not_bordering*matrix(sim_mat[two_m_ago_indices,s],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering)))/rowSums(not_bordering*matrix(data_vill$dogs[two_m_ago_indices],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering))) +
                                                                      rowSums(not_bordering*matrix(sim_mat[last_m_indices,s],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering)))/rowSums(not_bordering*matrix(data_vill$dogs[last_m_indices],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering))))/2) + data_vill_case_rate_adjust_dist)
        X_cap_m[,"log_case_rate_last2monthMean"] <- log(((sim_mat_cap[last_m_indices,s]/data_vill$dogs[last_m_indices]+sim_mat_cap[two_m_ago_indices,s]/data_vill$dogs[two_m_ago_indices])/2) + data_vill_case_rate_adjust_dist)
        X_cap_m[,"log_case_rate_notNeighbours_last2monthMean"] <- log(((rowSums(not_bordering*matrix(sim_mat_cap[two_m_ago_indices,s],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering)))/rowSums(not_bordering*matrix(data_vill$dogs[two_m_ago_indices],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering))) +
                                                                          rowSums(not_bordering*matrix(sim_mat_cap[last_m_indices,s],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering)))/rowSums(not_bordering*matrix(data_vill$dogs[last_m_indices],byrow=T,ncol=ncol(not_bordering),nrow=nrow(not_bordering))))/2) + data_vill_case_rate_adjust_dist)
        for(v in 1:nrow(vax_vill)){
          neighbours <- which(bordering[v,])
          case_rates <- cbind(sim_mat[two_m_ago_indices,s],sim_mat[last_m_indices,s])/cbind(data_vill$dogs[two_m_ago_indices],data_vill$dogs[last_m_indices])
          case_rates_w_park_mara <- rbind(case_rates,colMeans(case_rates[neighbours[which(!neighbours%in%(nrow(X_m)+1:2))],,drop=F]),0.01/12) 
          neighbour_case_rates_last2months <-  case_rates_w_park_mara[neighbours,]
          weights_neighbours <-  matrix(rep(borders_prop[v,neighbours],2),ncol=2) #weight by border length 
          X_m[v,"log_case_rate_neighbours_last2monthMean"] <- log(mean(colSums(neighbour_case_rates_last2months*weights_neighbours)) + data_vill_case_rate_adjust_dist)
          
          case_rates_cap <- cbind(sim_mat_cap[two_m_ago_indices,s],sim_mat_cap[last_m_indices,s])/cbind(data_vill$dogs[two_m_ago_indices],data_vill$dogs[last_m_indices])
          case_rates_w_park_mara_cap <- rbind(case_rates_cap,colMeans(case_rates_cap[neighbours[which(!neighbours%in%(nrow(X_m)+1:2))],,drop=F]),0.01/12) 
          neighbour_case_rates_last2months_cap <-  case_rates_w_park_mara_cap[neighbours,]
          X_cap_m[v,"log_case_rate_neighbours_last2monthMean"] <- log(mean(colSums(neighbour_case_rates_last2months_cap*weights_neighbours)) + data_vill_case_rate_adjust_dist)
          
        }
        
        # Simulate cases
        sim_mat[m_indices,s] <- rnbinom(nrow(X_m),mu=exp(colSums(t(X_m)*c(as.numeric(pars_s[c(1:(ncol(pars_s)-2))]),1,1))),size=as.numeric(pars_s["phi"]))
        sim_mat_cap[m_indices,s] <- pmin(rnbinom(nrow(X_cap_m),mu=exp(colSums(t(X_cap_m)*c(as.numeric(pars_s[c(1:(ncol(pars_s)-2))]),1,1))),size=as.numeric(pars_s["phi"])),round(case_cap*data_vill$dogs[m_indices]))
        
      }else{ # use vaccination-only model
        
        # Simulate cases
        sim_mat[m_indices,s] <- rnbinom(nrow(X_wo_cases_m),mu=exp(colSums(t(X_wo_cases_m)*c(as.numeric(pars_wo_cases_s[c(1:(ncol(pars_wo_cases_s)-2))]),1,1))),size=as.numeric(pars_wo_cases_s["phi"]))
        sim_mat_cap[m_indices,s] <- pmin(rnbinom(nrow(X_wo_cases_m),mu=exp(colSums(t(X_wo_cases_m)*c(as.numeric(pars_wo_cases_s[c(1:(ncol(pars_wo_cases_s)-2))]),1,1))),size=as.numeric(pars_wo_cases_s["phi"])),round(case_cap*data_vill$dogs[m_indices]))
        
      }
      
      
    }
  }
  
  sim_mat_dist <- rowsum((sim_mat),data_vill$month)
  sim_median <- apply(sim_mat_dist,1,quantile,0.5,na.rm=T)[-c(1:2)]
  sim_lower <- apply(sim_mat_dist,1,quantile,0.025,na.rm=T)[-c(1:2)]
  sim_upper <- apply(sim_mat_dist,1,quantile,0.975,na.rm=T)[-c(1:2)]
  sim_mat_dist_cap <- rowsum((sim_mat_cap),data_vill$month)
  sim_median_cap <- apply(sim_mat_dist_cap,1,quantile,0.5,na.rm=T)[-c(1:2)]
  sim_lower_cap <- apply(sim_mat_dist_cap,1,quantile,0.025,na.rm=T)[-c(1:2)]
  sim_upper_cap <- apply(sim_mat_dist_cap,1,quantile,0.975,na.rm=T)[-c(1:2)]
  
  save(sim_mat_dist, sim_lower, sim_upper, sim_median, sim_mat_dist_cap, sim_lower_cap, sim_upper_cap, sim_median_cap,
       file="Output/simulated_cases_vaccination_only_draw_reffs.Rdata")
  

  
}




pdf("Figs/simulate_ts_trial.pdf",width=7.5, height=12.5)
cex.axis <- 0.7
cex.lab <- 0.8
cex.pt <- 0.5
cex.main <- 0.8


# Plot fit to data
#___________________


# Using all data to predict a month ahead
#-----------

# Plot predictions and data over time
load("Output/preds_PI_stan_village_model.Rdata")
par(mar=c(2.5,2.5,1,1))
par(fig=c(0,1,0.75,1))
ylim <- c(0,max(preds_upper,data_vill$cases))
plot(NA,ylim=ylim,xlim=c(1,max(data_vill$month)),bty="l",
     ylab="",xlab="",cex.lab=cex.lab,axes=F,main="Simulate 1 month ahead with known case data",cex.main=cex.main)
axis(2,cex.axis=cex.axis,padj=1,at=seq(0,100,50))
axis(1,at=seq(1,length(2002:2022)*12,24),labels=paste(seq(2002,2022,2)),cex.axis=cex.axis,padj=-1.5)
mtext("Dog cases in district",side=2,line=1.5,cex=cex.lab)
mtext("Date",side=1,line=1.5,cex=cex.lab)
box(bty="l")
months_plot <- 3:max(data_vill$month)
polygon(c(months_plot,rev(months_plot)),c(preds_lower,rev(preds_upper)),col=scales::alpha("dodgerblue",0.4),border=NA)
out_PI <- which(data_dist$cases[months_plot]<preds_lower|data_dist$cases[months_plot]>preds_upper)+2
points(data_dist$cases~data_dist$month,col="navy",pch=20,cex=cex.pt)
points(data_dist$cases[out_PI]~data_dist$month[out_PI],col="red",pch=20,cex=cex.pt)
legend(130,115,c("model 95% PI","data within 95% PI","data outside 95% PI"),col=c("skyblue","navy","red"),pch=c(15,20,20),cex=0.75,bty="n",pt.cex = c(1.5,cex.pt,cex.pt))


# Using first two months of data and then simulated data to project
#-----------

# Plot predictions and data over time
load("Output/simulated_cases_use_initial_cases.Rdata")
par(mar=c(2.5,2.5,1,1))
par(fig=c(0,0.5,0.5,0.75),new=T)
ylim <- c(0,max(sim_upper,data_vill$cases))
ylim <- c(0,400)
plot(NA,ylim=ylim,xlim=c(1,max(data_vill$month)),bty="l",cex.main=cex.main,
     ylab="",xlab="",cex.lab=cex.lab,axes=F,main="Simulate full TS after 1st 2 months")
axis(2,cex.axis=cex.axis,padj=1)
axis(1,at=seq(1,length(2002:2022)*12,24),labels=paste(seq(2002,2022,2)),cex.axis=cex.axis,padj=-1.5)
mtext("Dog cases in district",side=2,line=1.5,cex=cex.lab)
mtext("Date",side=1,line=1.5,cex=cex.lab)
box(bty="l")
months_plot <- 3:max(data_vill$month)
polygon(c(months_plot,rev(months_plot)),c(sim_lower,rev(sim_upper)),col=scales::alpha("dodgerblue",0.4),border=NA)
out_PI <- which(data_dist$cases[months_plot]<sim_lower|data_dist$cases[months_plot]>sim_upper)+2
points(data_dist$cases~data_dist$month,col="navy",pch=20,cex=cex.pt)
points(data_dist$cases[out_PI]~data_dist$month[out_PI],col="red",pch=20,cex=cex.pt)
par(fig=c(0.5,1,0,0.57),new=T)

par(fig=c(0.5,1,0.5,0.75),new=T)
ylim <- c(0,max(sim_upper_cap,data_vill$cases))
ylim <- c(0,400)
plot(NA,ylim=ylim,xlim=c(1,max(data_vill$month)),bty="l",cex.main=cex.main,
     ylab="",xlab="",cex.lab=cex.lab,axes=F,main="Simulate full TS after 1st 2 months with cap")
axis(2,cex.axis=cex.axis,padj=1)
axis(1,at=seq(1,length(2002:2022)*12,24),labels=paste(seq(2002,2022,2)),cex.axis=cex.axis,padj=-1.5)
mtext("Dog cases in district",side=2,line=1.5,cex=cex.lab)
mtext("Date",side=1,line=1.5,cex=cex.lab)
box(bty="l")
months_plot <- 3:max(data_vill$month)
polygon(c(months_plot,rev(months_plot)),c(sim_lower_cap,rev(sim_upper_cap)),col=scales::alpha("dodgerblue",0.4),border=NA)
out_PI <- which(data_dist$cases[months_plot]<sim_lower_cap|data_dist$cases[months_plot]>sim_upper_cap)+2
points(data_dist$cases~data_dist$month,col="navy",pch=20,cex=cex.pt)
points(data_dist$cases[out_PI]~data_dist$month[out_PI],col="red",pch=20,cex=cex.pt)


# Using only vaccination data
#-----------

# Plot predictions and data over time
load("Output/simulated_cases_vaccination_only.Rdata")
par(mar=c(2.5,2.5,1,1))
par(fig=c(0,0.5,0.25,0.5),new=T)
ylim <- c(0,max(sim_upper,data_vill$cases))
ylim <- c(0,400)
plot(NA,ylim=ylim,xlim=c(1,max(data_vill$month)),bty="l",cex.main=cex.main,
     ylab="",xlab="",cex.lab=cex.lab,axes=F,main="Simulate full TS")
axis(2,cex.axis=cex.axis,padj=1)
axis(1,at=seq(1,length(2002:2022)*12,24),labels=paste(seq(2002,2022,2)),cex.axis=cex.axis,padj=-1.5)
mtext("Dog cases in district",side=2,line=1.5,cex=cex.lab)
mtext("Date",side=1,line=1.5,cex=cex.lab)
box(bty="l")
months_plot <- 3:max(data_vill$month)
polygon(c(months_plot,rev(months_plot)),c(sim_lower,rev(sim_upper)),col=scales::alpha("dodgerblue",0.4),border=NA)
out_PI <- which(data_dist$cases[months_plot]<sim_lower|data_dist$cases[months_plot]>sim_upper)+2
points(data_dist$cases~data_dist$month,col="navy",pch=20,cex=cex.pt)
points(data_dist$cases[out_PI]~data_dist$month[out_PI],col="red",pch=20,cex=cex.pt)

par(fig=c(0.5,1,0.25,0.5),new=T)
ylim <- c(0,max(sim_upper_cap,data_vill$cases))
ylim <- c(0,400)
plot(NA,ylim=ylim,xlim=c(1,max(data_vill$month)),bty="l",cex.main=cex.main,
     ylab="",xlab="",cex.lab=cex.lab,axes=F,main="Simulate full TS with cap")
axis(2,cex.axis=cex.axis,padj=1)
axis(1,at=seq(1,length(2002:2022)*12,24),labels=paste(seq(2002,2022,2)),cex.axis=cex.axis,padj=-1.5)
mtext("Dog cases in district",side=2,line=1.5,cex=cex.lab)
mtext("Date",side=1,line=1.5,cex=cex.lab)
box(bty="l")
months_plot <- 3:max(data_vill$month)
polygon(c(months_plot,rev(months_plot)),c(sim_lower_cap,rev(sim_upper_cap)),col=scales::alpha("dodgerblue",0.4),border=NA)
out_PI <- which(data_dist$cases[months_plot]<sim_lower_cap|data_dist$cases[months_plot]>sim_upper_cap)+2
points(data_dist$cases~data_dist$month,col="navy",pch=20,cex=cex.pt)
points(data_dist$cases[out_PI]~data_dist$month[out_PI],col="red",pch=20,cex=cex.pt)


# Using only vaccination data and randomly drawn random effect
#-----------

# Plot predictions and data over time
load("Output/simulated_cases_vaccination_only_draw_reffs.Rdata")
par(mar=c(2.5,2.5,1,1))
par(fig=c(0,0.5,0,0.25),new=T)
ylim <- c(0,max(sim_upper,data_vill$cases))
ylim <- c(0,400)
plot(NA,ylim=ylim,xlim=c(1,max(data_vill$month)),bty="l",cex.main=cex.main,
     ylab="",xlab="",cex.lab=cex.lab,axes=F,main="Simulate full TS with unknown REs")
axis(2,cex.axis=cex.axis,padj=1)
axis(1,at=seq(1,length(2002:2022)*12,24),labels=paste(seq(2002,2022,2)),cex.axis=cex.axis,padj=-1.5)
mtext("Dog cases in district",side=2,line=1.5,cex=cex.lab)
mtext("Date",side=1,line=1.5,cex=cex.lab)
box(bty="l")
months_plot <- 3:max(data_vill$month)
polygon(c(months_plot,rev(months_plot)),c(sim_lower,rev(sim_upper)),col=scales::alpha("dodgerblue",0.4),border=NA)
out_PI <- which(data_dist$cases[months_plot]<sim_lower|data_dist$cases[months_plot]>sim_upper)+2
points(data_dist$cases~data_dist$month,col="navy",pch=20,cex=cex.pt)
points(data_dist$cases[out_PI]~data_dist$month[out_PI],col="red",pch=20,cex=cex.pt)

par(fig=c(0.5,1,0,0.25),new=T)
ylim <- c(0,max(sim_upper_cap,data_vill$cases))
ylim <- c(0,400)
plot(NA,ylim=ylim,xlim=c(1,max(data_vill$month)),bty="l",cex.main=cex.main,
     ylab="",xlab="",cex.lab=cex.lab,axes=F,main="Simulate full TS with unknown REs and cap")
axis(2,cex.axis=cex.axis,padj=1)
axis(1,at=seq(1,length(2002:2022)*12,24),labels=paste(seq(2002,2022,2)),cex.axis=cex.axis,padj=-1.5)
mtext("Dog cases in district",side=2,line=1.5,cex=cex.lab)
mtext("Date",side=1,line=1.5,cex=cex.lab)
box(bty="l")
months_plot <- 3:max(data_vill$month)
polygon(c(months_plot,rev(months_plot)),c(sim_lower_cap,rev(sim_upper_cap)),col=scales::alpha("dodgerblue",0.4),border=NA)
out_PI <- which(data_dist$cases[months_plot]<sim_lower_cap|data_dist$cases[months_plot]>sim_upper_cap)+2
points(data_dist$cases~data_dist$month,col="navy",pch=20,cex=cex.pt)
points(data_dist$cases[out_PI]~data_dist$month[out_PI],col="red",pch=20,cex=cex.pt)

dev.off()

