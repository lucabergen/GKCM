
# Install necessary packages
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

pak::pkg_install(c("batchtools@0.9.18", "here@1.0.2", "comets@0.2-2", 
                   "drf@1.1.0", "momentchi2@0.1.5"))

library(batchtools)
library(here)

# Load registry with computed experiments
reg_dir <- here("registry")
reg <- loadRegistry(reg_dir, writeable = T)
reg$work.dir <- reg_dir

################################################################################

# All experiments are already computed. To run the simulation yourself 
# (caution: this will take quite long!) uncomment and execute the following: 

# Delete existing directory
# unlink(reg_dir, recursive = TRUE)

# Create empty registry
# reg <- makeExperimentRegistry(file.dir = reg_dir, seed = 42)

# No. of iterations
k <- 100

# Sample sizes
n_vec <- c(500, 1000, 1500, 2000)

# Number of cores used
n_cores <- parallel::detectCores() - 1

# Set cluster functions depending on OS
if (.Platform$OS.type == "windows") {
  reg$cluster.functions <- makeClusterFunctionsSocket(n_cores)
} else {
  reg$cluster.functions <- makeClusterFunctionsMulticore(n_cores)
}

# Load data generating processes and tests
source(here("files/problems.R"))
addProblem(name = "null_1", fun = dgp_null1, seed = 100)
addProblem(name = "null_2", fun = dgp_null2, seed = 200)
addProblem(name = "null_3", fun = dgp_null3, seed = 300)
addProblem(name = "null_4", fun = dgp_null4, seed = 400)
addProblem(name = "alt_1",  fun = dgp_alt1,  seed = 500)
addProblem(name = "alt_2",  fun = dgp_alt2,  seed = 600)
addProblem(name = "alt_3",  fun = dgp_alt3,  seed = 700)

source(here("files/algorithms.R"))
addAlgorithm(name = "PCM",      fun = PCM_wrapper)
addAlgorithm(name = "GCM",      fun = GCM_wrapper)
addAlgorithm(name = "wGCM",     fun = wGCM_wrapper)
addAlgorithm(name = "KCIT",     fun = KCIT_wrapper)
addAlgorithm(name = "RCIT",     fun = RCIT_wrapper)
addAlgorithm(name = "RCoT",     fun = RCoT_wrapper)
addAlgorithm(name = "GKCM_RF",  fun = GKCM_RF_wrapper)
addAlgorithm(name = "GKCM_KRR", fun = GKCM_KRR_wrapper)

# Define parameter grids for problems and algorithms
prob_df <- expand.grid(n = n_vec)
prob_dsgn <- list(
  null_1 = prob_df, 
  null_2 = prob_df,
  null_3 = prob_df,
  null_4 = prob_df,
  alt_1  = prob_df,
  alt_2  = prob_df,
  alt_3  = prob_df
)

algo_df <- expand.grid()
algo_dsgn <- list(
  PCM   = algo_df,
  GCM   = algo_df,
  wGCM  = algo_df,
  KCIT  = algo_df,
  RCIT  = algo_df,
  RCoT  = algo_df,
  GKCM_RF = algo_df,
  GKCM_KRR = algo_df
)

# Add experiments (problem x problem parameters x 
# algorithm x hyperparameters x iterations)
addExperiments(prob_dsgn, algo_dsgn, repls = k)

# Get an overview of jobs and test them
summarizeExperiments()
testJob(id = 1)

# Submit
submitJobs(resources = list(name = reg_dir, memory = 3000, walltime = 1800))

# Check results
getStatus()

# In case there were errors: Inspect them
# err_ids <- findErrors()
# getErrorMessages(err_ids, missing.as.error = TRUE)
