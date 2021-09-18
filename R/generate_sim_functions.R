read_pwgs_csv <- function(csv_files) {
  required_cols <- c(
    "SampleId",
    "SampleName",
    "OnyxDesignDnaId",
    "DesignId",
    "DesignCategory",
    "OnyxEditReadCount",
    "FracOnyxEditReadCount",
    "TotalSpanningReadCount"
  )
  optional_cols <- c(
    "DesignSetId",
    "DesignPlasmidReferenceCoverage"
  )
  csv <- tibble()
  for(csv_file in csv_files) {
    this_csv <- read.csv(csv_file, header=TRUE, as.is=TRUE, na.strings=c("NA", "#NA", "nan"))
    required <- this_csv %>% select(all_of(required_cols))
    optional <- this_csv %>% select(any_of(optional_cols))
    for(opt_col in optional_cols) {
      if(!is.element(opt_col, colnames(optional)))
        optional <- optional %>% mutate(!!opt_col := NA)
    }
    csv <- csv %>% bind_rows(
      bind_cols(required, optional) %>%
        select(all_of(c(required_cols, optional_cols))) %>%
        arrange(SampleId, DesignId)
    )
  }
  csv$DesignSetId[is.na(csv$DesignSetId)] <- csv$OnyxDesignDnaId[is.na(csv$DesignSetId)]
  return(csv)
}

combine_pwgs_design_reports <- function(prefix) {
  pwgs <- tibble()
  for(csv in Sys.glob(sprintf("%s*_PooledWGSDesignReport.csv", prefix))) {
    sample_name <- gsub("_PooledWGSDesignReport.csv", "", gsub(prefix, '', csv))
    pwgs <- pwgs %>%
      bind_rows(
        x <- read_pwgs_csv(csv)
        %>%
          filter(DesignCategory=="edit") %>%
          mutate(
            SampleName = rep(sample_name, length(SampleName)),
            EditRepresentation = OnyxEditReadCount/TotalSpanningReadCount*length(OnyxEditReadCount),
            DesignRepresentation = DesignPlasmidReferenceCoverage/sum(DesignPlasmidReferenceCoverage)*length(DesignPlasmidReferenceCoverage)
          ) %>%
          add_column(RunId=run_id) %>%
          arrange(OnyxDesignDnaId, DesignId)
      )
  }
  pwgs
}

preseq_models <- function(counts) {
  # Transform counts into the format that preseqR functions want
  counts <- counts[counts > 0 & !is.na(counts)]
  if(all(counts < 2)) {
    # preseq methods don't work when there are no duplicate counts
    null_function <- function(x){ rep(NA, length(x)) }
    return(list('ztnb' = null_function, 'ds' = null_function))
  } else {
    count_tab <- table(counts)
    count_tab <- cbind(as.numeric(names(count_tab)), count_tab)
    # the preseqR functions generate some deprecation warnings
    suppressWarnings(preseq_ztnb <- ztnb.rSAC(count_tab))
    suppressWarnings(preseq_ds <- ds.rSAC(count_tab))
    return(list('ztnb' = preseq_ztnb, 'ds' = preseq_ds))
  }
}

frac_edit_estimate <- function(k, n, alpha=0.05) {
  frac_measured <- mean(n>0)
  p_est <- k[n>0]/n[n>0]
  res <- list(
    est = sum(p_est)/frac_measured,
    sd = sqrt(sum(p_est*(1-p_est)/n))/frac_measured
  )
  res$lo <- res$est - qnorm(1-alpha/2)*res$sd
  res$hi <- res$est + qnorm(1-alpha/2)*res$sd
  return(res)
}

analyze_richness <- function(pwgs, edit_frac=NA) {
  plasmid_count <- pwgs %>% pull(DesignPlasmidReferenceCoverage)
  plasmid_depth <- sum(plasmid_count)
  plasmid_repr <- nrow(pwgs)*plasmid_count/plasmid_depth

  total_depth <- pwgs %>% pull(TotalSpanningReadCount)
  depth_outlier <- abs((total_depth-mean(total_depth))/sd(total_depth)) > 4
  if(is.na(edit_frac)) {
    # estimate edit fraction from data
    edit_fraction <- frac_edit_estimate(
      pwgs$OnyxEditReadCount[!depth_outlier],
      pwgs$TotalSpanningReadCount[!depth_outlier],
    )
  } else {
    # known edit fraction has been supplied
    edit_fraction <- list(
      est = edit_frac,
      sd = NA,
      lo = NA,
      hi = NA
    )
  }
  edit_count <- pwgs %>% pull(OnyxEditReadCount)
  edit_depth <- round(edit_fraction$est * pwgs %>% pull(TotalSpanningReadCount))
  edit_count[depth_outlier] <- NA
  edit_depth[depth_outlier] <- NA
  edit_repr <- nrow(pwgs)*edit_count/edit_depth

  # Plots and beta-binomial fits for plasmids and edits
  repr_range <- range(c(plasmid_repr[plasmid_repr > 0], edit_repr[edit_repr > 0]), na.rm=TRUE)
  if(!is.na(plasmid_depth) & plasmid_depth > 0) {
    res_plasmid <- repr_hist_and_qq_plot(plasmid_count, plasmid_depth, repr_range=repr_range, var_type="Design")
  } else {
    res_plasmid <- list(
      frac_obs=NA,
      sampling_depth=NA,
      cv_sample_est=NA,
      cv_bbinom_est=NA,
      cv_bbinom_cv=NA,
      cv_bbinom_lo=NA,
      cv_bbinom_hi=NA,
      repr_range=NA
    )
  }
  res_edit <- repr_hist_and_qq_plot(edit_count, edit_depth, repr_range=repr_range, var_type="Edit")
  res_edit$edit_fraction <- edit_fraction

  # Add preseq fits
  res_edit$preseq <- preseq_models(edit_count)
  res_plasmid$preseq <- preseq_models(plasmid_count)

  return(list(
    'edit' = res_edit,
    'plasmid' = res_plasmid
  ))
}

downsample <- function(counts, downsample_factor) {
  keep_probability <- 1/downsample_factor
  rbinom(length(counts), round(counts), keep_probability)
}

downsample_pwgs <- function(pwgs_in, downsample_factor=NA, target_coverage=NA) {
  if(is.na(downsample_factor)) {
    if(is.na(target_coverage))
      stop("Must specify either downsample_factor or target_coverage")
    actual_coverage = median(pwgs_in$TotalSpanningReadCount)
    if(target_coverage > actual_coverage)
      stop(sprintf("target_coverage is greater than actual_coverage which is %f", actual_coverage))
    downsample_factor = actual_coverage/target_coverage
  } else if(downsample_factor < 1) {
      stop("downsample_factor must be greater than 1")
  }
    
  edit_reads <- downsample(pwgs_in$OnyxEditReadCount, downsample_factor)
  ref_reads <- downsample(pwgs_in$TotalSpanningReadCount - pwgs_in$OnyxEditReadCount, downsample_factor)
  plasmid_reads <- downsample(pwgs_in$DesignPlasmidReferenceCoverage, downsample_factor)
  pwgs_in %>%
    mutate(
      OnyxEditReadCount = edit_reads,
      TotalSpanningReadCount = ref_reads + edit_reads,
      DesignPlasmidReferenceCoverage = plasmid_reads,
      FracOnyxEditReadCount = OnyxEditReadCount / TotalSpanningReadCount
    ) %>%
    mutate(
      EditRepresentation = OnyxEditReadCount/TotalSpanningReadCount*length(OnyxEditReadCount),
      DesignRepresentation = DesignPlasmidReferenceCoverage/sum(DesignPlasmidReferenceCoverage)*length(DesignPlasmidReferenceCoverage)
    )
}

downsample_and_analyze_richness <- function(target_coverage, seed=NA, data, edit_frac=NA, verbose=TRUE) {
  if(!is.na(seed))
    set.seed(seed)
  full_target_coverage <- median(data$TotalSpanningReadCount)
  pwgs_downsampled <- downsample_pwgs(data, target_coverage=target_coverage)
  temp <- analyze_richness(pwgs_downsampled)
  res_preseq_func <- rbind(c(
    'EditPreseqDs'     = temp$edit$preseq$ds,
    'EditPreseqZtnb'   = temp$edit$preseq$ztnb,
    'DesignPreseqDs'   = temp$plasmid$preseq$ds,
    'DesignPreseqZtnb' = temp$plasmid$preseq$ztnb
  ))
  res_table <- rbind(c(
    'EditFracEstimated'  = 1,
    'LibSize'            = nrow(pwgs_downsampled),
    'TargetCoverage'     = target_coverage,
    'FullTargetCoverage' = full_target_coverage,
    'EditFracEst'   = temp$edit$edit_fraction$est,
    'EditFracLo'    = temp$edit$edit_fraction$lo,
    'EditFracHi'    = temp$edit$edit_fraction$hi,
    'EditCvEst'     = temp$edit$cv_bbinom_est,
    'EditCvLo'      = temp$edit$cv_bbinom_lo,
    'EditCvHi'      = temp$edit$cv_bbinom_hi,
    'EditFracObs'   = temp$edit$frac_obs,
    'DesignDepth'   = temp$plasmid$sampling_depth,
    'DesignCvEst'   = temp$plasmid$cv_bbinom_est,
    'DesignCvLo'    = temp$plasmid$cv_bbinom_lo,
    'DesignCvHi'    = temp$plasmid$cv_bbinom_hi,
    'DesignFracObs' = temp$plasmid$frac_obs
  ))
  if(!is.na(edit_frac)) {
    res_preseq_func <- res_preseq_func %>% rbind(res_preseq_func)
    temp <- analyze_richness(pwgs_downsampled, edit_frac=edit_frac)
    res_table <- res_table %>% rbind(c(
      'EditFracEstimated'  = 0,
      'LibSize'            = nrow(pwgs_downsampled),
      'TargetCoverage'     = target_coverage,
      'FullTargetCoverage' = full_target_coverage,
      'EditFracEst'   = temp$edit$edit_fraction$est,
      'EditFracLo'    = temp$edit$edit_fraction$lo,
      'EditFracHi'    = temp$edit$edit_fraction$hi,
      'EditCvEst'     = temp$edit$cv_bbinom_est,
      'EditCvLo'      = temp$edit$cv_bbinom_lo,
      'EditCvHi'      = temp$edit$cv_bbinom_hi,
      'EditFracObs'   = temp$edit$frac_obs,
      'DesignDepth'   = temp$plasmid$sampling_depth,
      'DesignCvEst'   = temp$plasmid$cv_bbinom_est,
      'DesignCvLo'    = temp$plasmid$cv_bbinom_lo,
      'DesignCvHi'    = temp$plasmid$cv_bbinom_hi,
      'DesignFracObs' = temp$plasmid$frac_obs
    ))
  }
  return(list('table'=as_tibble(res_table), 'preseq_func'=as_tibble(res_preseq_func)))
}

simulate_results <- function(data, target_coverage, edit_frac_est) {
  seed <- 1:length(target_coverage)
  sim <- foreach(tc=target_coverage, seed=seed) %dopar%
    downsample_and_analyze_richness(tc, seed=seed, data=data, edit_frac=edit_frac_est)

  cat("Aggregating results into table\n")
  sim_table <- sim %>%
    lapply(function(x){t(x$table)}) %>%
    unlist() %>%
    matrix(ncol=ncol(sim[[1]]$table), byrow=TRUE, dimnames=list(NULL, colnames(sim[[1]]$table))) %>%
    as_tibble() %>%
    mutate(EditFracEstimated=as.logical(EditFracEstimated))

  cat("Aggregating results into preseq functions\n")
  sim_preseq <- sim %>%
    lapply(function(x){t(x$preseq_func)}) %>%
    unlist() %>%
    matrix(ncol=ncol(sim[[1]]$preseq_func), byrow=TRUE, dimnames=list(NULL, colnames(sim[[1]]$preseq_func))) %>%
    as_tibble()

  sim <- bind_cols(sim_table, sim_preseq)
  cat("Returning results\n")
  return(sim)
}
