### ==================================================================== ###
### Example Code: Using PSO to generate K-stage clinical trial design
### ==================================================================== ###

### Preliminaries
pkg_need <- c("here", "Rcpp", "RcppArmadillo", "globpso")
for (i in 1:length(pkg_need)) {
  if (!(pkg_need[i] %in% rownames(installed.packages()))) {
    install.packages(pkg_need[i])
  }
}

### Import the main local R file "genKStageDesign.R"
setwd(here::here())
source("genKStageDesign.R")

### Set required number of stages
nStage <- 3
### Set requirements of the clinical trial
cliRequirement <- list(
  p0 = 0.2,    # response rate in the null hypothesis
  p1 = 0.4,    # response rate in the alternative hypothesis
  alpha = 0.1, # upper bound of type I error
  beta = 0.1   # upper bound of type II error
)

### -------------------------------------------------------- ###
### Find Optimal Design 
### -------------------------------------------------------- ###
### Run PSO for optimal design
optimRes <- genKStageDesign(designType = "optimal", # set "optimal" for searching for optimal design
                            nStage = nStage, # the required number of stages
                            cliRequirement = cliRequirement, # the clinical trial requirements
                            OmegaL = 30, # the lower bound of the total sample size
                            OmegaU = 70, # the upper bound of the total sample size
                            m1 = 10,     # the minimal sample size at the first stage
                            mk = 1,      # the minimal sample size at stages after the first stage
                            qk = 0,      # the minimum incremental critical value at each stage
                            psoSetting = list(nSwarm = 128, maxIter = 400), # the PSO configurations
                            seed = 1,    # set Random seed for reproducibility
                            verbose = TRUE)

### Show the resulting optimal design
# print(optimRes)

### -------------------------------------------------------- ###
### Find Minimax Design 
### -------------------------------------------------------- ###
### Run PSO for minimax design
minMaxRes <- genKStageDesign(designType = "minimax", # set "optimal" for searching for optimal design
                             nStage = nStage, # the required number of stages
                             cliRequirement = cliRequirement, # the clinical trial requirements
                             OmegaL = 30, # the lower bound of the total sample size
                             OmegaU = 70, # the upper bound of the total sample size
                             m1 = 10,     # the minimal sample size at the first stage
                             mk = 1,      # the minimal sample size at stages after the first stage
                             qk = 0,      # the minimum incremental critical value at each stage
                             psoSetting = list(nSwarm = 128, maxIter = 400), # the PSO configurations
                             seed = 1,    # set Random seed for reproducibility
                             verbose = TRUE)

### Show the resulting minimax design
# print(minMaxRes)
