
# Create 1km grid of the study area with cells assigned to district, ward and
# village
source("R/Create_SD_grid.R")

# Estimate human and dog populations monthly by village and district
source("R/SerengetiDogsVillage.R")
source("R/SerengetiDogsGrid.R")

# Produce Fig. 1
source("R/Fig.1.R")

# Look at parameter options for bounding vaccination coverage in (0,1). And
# produce Fig. S11
source("R/Bound_coverage.R")

# Obtain vaccination coverage by village and district from various sources of
# vaccination data
source("R/CreateVaccination.R")

# Process CT data for downstream scripts
source("R/process_CT_data.R")
