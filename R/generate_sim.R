library(tidyverse)
library(scales)
library(preseqR)
library(foreach)
library(doParallel)
source("R/generate_sim_functions.R")
source("R/richness.R")

run_ids <- c("PV38", "PPV294")

pwgs <- tibble()
for(run_id in run_ids)
  pwgs <- pwgs %>%
    bind_rows(combine_pwgs_design_reports(paste0("Data/",run_id,"/")))

sim_runs <- list(
  c('run_id'='PPV294', 'trn'='PPV294-1', 'tst'='PPV294-2'),
  c('run_id'='PPV294', 'trn'='PPV294-2', 'tst'='PPV294-1'),
  c('run_id'='PV38',   'trn'='PV38-1',   'tst'='PV38-2'),
  c('run_id'='PV38',   'trn'='PV38-2',   'tst'='PV38-1')
)
n_sim <- 1000
res <- numeric()
my.cluster <- parallel::makeCluster(parallel::detectCores(), type="FORK")
doParallel::registerDoParallel(cl = my.cluster)
for(sim_run in sim_runs) {
  cat(paste0(date(),"\n"))
  cat(sprintf("Training on %s to predict in %s\n", sim_run['trn'], sim_run['tst']))
  test_data <- pwgs %>% filter(SampleName==sim_run['tst'])
  truth <- downsample_and_analyze_richness(
    target_coverage=median(test_data$TotalSpanningReadCount),
    data=test_data
  )
  train_data <- pwgs %>% filter(SampleName==sim_run['trn'])
  target_coverage <- rep(
    c(
      100,
      250,
      500,
      1000*seq(
        1,
        floor(min(median(train_data$TotalSpanningReadCount),median(test_data$TotalSpanningReadCount))/1000)
      )
    ),
    each=n_sim
  )
  sim <- simulate_results(train_data, target_coverage, truth$table$EditFracEst)
  res <- c(res, list(list(
    'tst' = sim_run['tst'],
    'trn' = sim_run['trn'],
    'truth' = truth,
    'sim' = sim
  )))
}
cat("Saving...\n")
save.image(file="Data/rdata/richness_sim.RData")
