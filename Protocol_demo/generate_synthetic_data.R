#!/usr/bin/env Rscript
# FULLY SYNTHETIC DATA: no real participants or fitted UKB parameters are used.
# Usage: Rscript generate_synthetic_data.R output_directory [n] [p] [seed]
# Defaults: n=5000, p=20, seed=20260801. Requires base R only.
# This companion uses the same explicit PRNG/draw order as the tested Python
# generator. Numerical results may differ at floating-point rounding precision.
# IMPORTANT: this R companion has not been executed in the file-creation environment.
# See README.md for test status and run_synthetic_example.R for the analysis driver.
# This file only generates/preprocesses input; Cox results are computed by the driver.

make_rng <- function(seed) {
  if (length(seed) != 1L || !is.finite(seed) || seed < 1 || seed >= 2147483647 || seed != floor(seed)) {
    stop("seed must be an integer between 1 and 2147483646")
  }
  state <- as.double(seed)
  uniform <- function(n) {
    ans <- numeric(n)
    for (i in seq_len(n)) {
      state <<- (16807 * state) %% 2147483647
      ans[i] <- state / 2147483647
    }
    ans
  }
  normal <- function(n) {
    k <- ceiling(n / 2)
    u1 <- uniform(k)
    u2 <- uniform(k)
    radius <- sqrt(-2 * log(u1))
    as.vector(rbind(radius * cos(2*pi*u2), radius * sin(2*pi*u2)))[seq_len(n)]
  }
  categorical <- function(n, probabilities) {
    if (abs(sum(probabilities) - 1) > 1e-12) stop("Probabilities must sum to 1")
    cuts <- cumsum(probabilities)
    vapply(uniform(n), function(u) as.integer(sum(u >= cuts)), integer(1))
  }
  list(uniform=uniform, normal=normal, categorical=categorical)
}

preprocess_metabolites <- function(raw) {
  z <- matrix(NA_real_, nrow(raw), ncol(raw), dimnames=dimnames(raw))
  audit <- vector("list", ncol(raw))
  for (j in seq_len(ncol(raw))) {
    x <- raw[,j]
    if (any(is.infinite(x))) stop(colnames(raw)[j], ": infinite value")
    observed <- x[is.finite(x)]
    if (!length(observed)) stop(colnames(raw)[j], ": all values are missing")
    missing <- !is.finite(x)
    med <- median(observed)
    n_nonpositive <- sum(is.finite(x) & x <= 0)
    x[missing] <- med
    positive <- x[x > 0]
    if (!length(positive)) stop(colnames(raw)[j], ": no positive values")
    replacement <- min(positive)/2
    x[x <= 0] <- replacement
    logged <- log(x)
    mn <- mean(logged)
    s <- sd(logged)
    if (!is.finite(s) || s <= 0) stop(colnames(raw)[j], ": invalid log-scale SD")
    z[,j] <- (logged-mn)/s
    audit[[j]] <- data.frame(
      Metabolite_raw=colnames(raw)[j], n=length(x),
      n_missing_raw=sum(missing), missing_fraction_raw=mean(missing),
      n_nonpositive_raw=n_nonpositive, imputation_median=med,
      nonpositive_replacement=replacement, log_mean=mn, log_sample_sd=s,
      n_missing_preprocessed=sum(is.na(z[,j])),
      preprocessed_mean=mean(z[,j]), preprocessed_sample_sd=sd(z[,j]))
  }
  list(values=z, audit=do.call(rbind,audit))
}

generate_synthetic_data <- function(out_dir="synthetic_data", n=5000L,
                                    n_metabolites=20L, seed=20260801L) {
  if (length(n)!=1L || !is.finite(n) || n<1000 || n!=floor(n)) stop("Use integer n >= 1000")
  if (length(n_metabolites)!=1L || !is.finite(n_metabolites) || n_metabolites<10 || n_metabolites!=floor(n_metabolites)) stop("Use integer n_metabolites >= 10")
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  rng <- make_rng(seed)
  clamp <- function(x, lo, hi) pmin(pmax(x,lo),hi)
  round_half_up <- function(x, scale) floor(x*scale+0.5)/scale
  # Arbitrary simulation distributions, NOT empirical UKB estimates.
  age <- round_half_up(clamp(56+8*rng$normal(n),40,75),10)
  sex_code <- as.integer(rng$uniform(n)<0.47)
  sex <- ifelse(sex_code==1L,"Male","Female")
  ethnicity_code <- rng$categorical(n,c(0.70,0.20,0.10))
  ethnicity <- c("Group_A","Group_B","Group_C")[ethnicity_code+1L]
  education_code <- rng$categorical(n,c(0.30,0.40,0.30))
  education <- c("Low","Medium","High")[education_code+1L]
  deprivation <- 2.5*rng$normal(n)
  bmi <- round_half_up(clamp(27+4*rng$normal(n)+0.6*sex_code,18,43),100)
  smoking_code <- rng$categorical(n,c(0.55,0.30,0.15))
  smoking <- c("Never","Former","Current")[smoking_code+1L]
  alcohol_code <- rng$categorical(n,c(0.25,0.35,0.40))
  alcohol <- c("Never","Occasional","Regular")[alcohol_code+1L]
  tv_time <- round_half_up(clamp(2.8+1.2*rng$normal(n),0,8),100)
  computer_time <- round_half_up(clamp(2+1.1*rng$normal(n),0,8),100)
  sleep <- round_half_up(clamp(7+0.8*rng$normal(n),4,10),100)
  diabetes_prob <- 1/(1+exp(-(-2.7+0.05*(bmi-27)+0.025*(age-56))))
  diabetes <- as.integer(rng$uniform(n)<diabetes_prob)
  cvd_prob <- 1/(1+exp(-(-2.5+0.03*(age-56)+0.25*sex_code+0.30*diabetes)))
  cvd <- as.integer(rng$uniform(n)<cvd_prob)
  baseline <- as.Date("2010-01-01")+floor(1096*rng$uniform(n))
  admin_end <- as.Date("2026-01-01")
  dropout <- rng$uniform(n)<0.10
  dropout_days <- floor(365.25*(2+10*rng$uniform(n)))
  observation_end <- rep(admin_end,n)
  observation_end[dropout] <- pmin(admin_end,baseline[dropout]+dropout_days[dropout])
  observation_days <- as.numeric(observation_end-baseline)
  factors <- do.call(cbind,lapply(seq_len(4L),function(j) rng$normal(n)))
  log_risk <- (0.65*factors[,1]-0.55*factors[,2]+0.40*factors[,3]+
               0.06*(age-56)/8+0.08*sex_code+0.05*(bmi-27)/4+0.12*diabetes+0.08*cvd)
  latent_event_years <- -log(rng$uniform(n))/(0.0028*exp(log_risk))
  latent_event_days <- pmax(1,floor(latent_event_years*365.25+0.5))
  status <- as.integer(latent_event_days<=observation_days)
  diagnosis <- rep(as.Date(NA_character_),n)
  diagnosis[status==1L] <- baseline[status==1L]+latent_event_days[status==1L]
  years_to_diagnosis <- ifelse(status==1L,latent_event_days/365.25,NA_real_)
  time_at_risk_years <- ifelse(status==1L,latent_event_days,observation_days)/365.25
  raw <- matrix(NA_real_,n,n_metabolites,dimnames=list(NULL,paste0("Met",seq_len(n_metabolites))))
  dictionary <- vector("list",n_metabolites)
  block_coefficients <- c(0.65,-0.55,0.40,0)
  for (j in seq_len(n_metabolites)) {
    j0 <- j-1L
    fraction <- (j0+0.5)/n_metabolites
    block <- 1L+sum(fraction>=c(0.30,0.60,0.80))
    loading <- 0.72+0.04*(j0%%3)
    noise_sd <- 0.43+0.02*(j0%%4)
    mu <- 1+0.08*(j0%%5)
    nuisance <- (0.06*(age-56)/8+0.08*(bmi-27)/4+0.06*sex_code+
                 0.06*diabetes+0.04*(smoking_code==2L))
    raw[,j] <- exp(mu+loading*factors[,block]+nuisance+noise_sd*rng$normal(n))
    flags <- rng$uniform(n)
    missing_prob <- 0.01+0.005*(j0%%5)
    raw[flags<missing_prob,j] <- NA_real_
    raw[flags>=missing_prob & flags<missing_prob+0.001,j] <- 0
    raw[flags>=missing_prob+0.001 & flags<missing_prob+0.0015,j] <- -0.01
    dictionary[[j]] <- data.frame(
      Metabolite_raw=paste0("Met",j),Metabolite=paste0("Meta_",j),
      Metabolite_name=sprintf("Synthetic metabolite %02d",j),
      simulation_block=paste0("Block_",block),raw_unit="arbitrary simulation units",
      analysis_unit="SD of log-transformed imputed values",
      block_loading=loading,residual_log_sd=noise_sd,
      latent_block_event_log_hazard_coefficient=block_coefficients[block],
      missing_probability=missing_prob)
  }
  pp <- preprocess_metabolites(raw)
  common <- data.frame(participant_id=sprintf("SYN%06d",seq_len(n)),baseline_date=baseline,
    diagnosis_date=diagnosis,followup_end=observation_end,status=status,
    years_to_diagnosis=years_to_diagnosis,followup_years=observation_days/365.25,
    time_at_risk_years=time_at_risk_years,age=age,sex=sex,ethnicity=ethnicity,
    education=education,deprivation=deprivation,bmi=bmi,smoking=smoking,alcohol=alcohol,
    tv_time=tv_time,computer_time=computer_time,sleep=sleep,diabetes=diabetes,cvd=cvd,
    stringsAsFactors=FALSE)
  cohort_raw <- cbind(common,as.data.frame(raw))
  cohort <- cbind(common,as.data.frame(pp$values))
  dictionary <- do.call(rbind,dictionary)
  write.csv(cohort_raw,file.path(out_dir,"cohort_raw.csv"),row.names=FALSE,na="")
  write.csv(cohort,file.path(out_dir,"cohort.csv"),row.names=FALSE,na="")
  write.csv(dictionary,file.path(out_dir,"metabolite_dictionary.csv"),row.names=FALSE,na="")
  write.csv(pp$audit,file.path(out_dir,"metabolite_preprocessing_qc.csv"),row.names=FALSE,na="")
  saveRDS(cohort,file.path(out_dir,"cohort.rds"))
  saveRDS(cohort_raw,file.path(out_dir,"cohort_raw.rds"))
  saveRDS(list(n=n,n_metabolites=n_metabolites,seed=seed,
               data_type="FULLY_SYNTHETIC",baseline_rate=0.0028,
               latent_block_log_hazard=block_coefficients),
          file.path(out_dir,"simulation_parameters.rds"))
  sufficient <- all(vapply(c("Female","Male"),function(s)
    sum(sex==s & status==0L)>=10*sum(sex==s & status==1L),logical(1)))
  checks <- c(unique_participant_ids=!anyDuplicated(cohort$participant_id),
    valid_observation_intervals=all(observation_days>0),
    valid_event_dates=all(diagnosis[status==1L]>baseline[status==1L] &
                         diagnosis[status==1L]<=observation_end[status==1L]),
    enough_controls_per_sex_for_1_to_10=sufficient,
    at_least_50_cases=sum(status)>=50,
    at_least_10_distinct_case_times=length(unique(years_to_diagnosis[status==1L]))>=10,
    finite_processed_values=all(is.finite(pp$values)),
    zero_processed_means=all(abs(colMeans(pp$values))<1e-10),
    unit_processed_sds=all(abs(apply(pp$values,2,sd)-1)<1e-10))
  write.csv(data.frame(check=names(checks),result=ifelse(checks,"PASS","FAIL")),
            file.path(out_dir,"validation_checks_R.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo_R.txt"))
  if (!all(checks)) stop("A simulation check failed; inspect validation_checks_R.csv")
  message("SYNTHETIC ONLY: n=",n,"; metabolites=",n_metabolites,
          "; incident cases=",sum(status),"; noncases=",sum(status==0L))
  invisible(list(raw=cohort_raw,cohort=cohort,dictionary=dictionary,qc=pp$audit))
}

if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  out <- if(length(args)>=1L) args[1] else "synthetic_data"
  n <- if(length(args)>=2L) as.integer(args[2]) else 5000L
  p <- if(length(args)>=3L) as.integer(args[3]) else 20L
  seed <- if(length(args)>=4L) as.integer(args[4]) else 20260801L
  generate_synthetic_data(out,n,p,seed)
}
