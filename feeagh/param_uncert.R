rm(list=ls())
graphics.off()
cat("\14")

# Load libraries
library(gotmtools)
library(LakeEnsemblR)
library(ggplot2)
library(ggpubr)
library(FME)
library(reshape2)


setwd("feeagh") # Change working directory to example folder



##---------------- calibrate and run models with daily forcing ----------------
# Set config file & models for daily forcing data
config_file <- 'LakeEnsemblR_pars_uncert.yaml'
model <-  c("FLake", "GLM", "GOTM", "Simstrat", "MyLake") 

run_ensemble(config_file = config_file, model = model)

ncdf <- "output/ensemble_output_pars.nc"

out <- load_var(ncdf, var = "temp", return = "array")

res <- analyse_ncdf(ncdf, model)

res$strat


model <- c("FLake", "GOTM", "Simstrat", "GLM")
input_yaml(config_file, "output", "file", paste0("output_", paste0(model, collapse = "_")))
cal_fold <- "mcmc_v1"
file_end <- "_MCMC_202102251516.csv"
# load MCMC pars
m <- model[4]
fils <- list.files(cal_fold, pattern = m, full.names = TRUE)
res <- NA
for( i in fils ) {
  cal_res <- read.csv(i)
  cal_res$run <- 1:nrow(cal_res)
  cal_res$chain <- which(fils == i)
  if(!is.data.frame(res)) {
    res <- cal_res
  } else {
    res <- rbind(res, cal_res)
  }
}

# summary(res)
res <- res[res$qual < 5e5, ]
res$chain <- factor(res$chain)
mlt <- reshape2::melt(res, id.vars = c("qual", "run", "chain"))
ggplot(mlt, aes(run, value, color = chain)) +
  # geom_point(aes(run, value)) +
  geom_line() +
  geom_smooth() +
  facet_wrap(~variable, scales = "free_y") +
  scale_y_log10()

ggplot(mlt, aes(value, qual, color = chain)) +
  geom_point() +
  # geom_line(aes(run, value)) +
  geom_smooth() +
  facet_wrap(~variable, scales = "free_x") +
  scale_x_log10()


# Select chain 3
res2 <- res[res$chain == 3, ]
burn_in <- 1000
res2 <- res2[res2$run > burn_in, ]

# Create distributions
wind_mn <- mean(res2$wind_speed)
wind_sd <- sd(res2$wind_speed)

swr_mn <- mean(res2$swr)
swr_sd <- sd(res2$swr)

mpar_mn <- mean(res2[, 3])
mpar_sd <- sd(res2[, 3])

cov_mat <- cov(res2[, 1:3])
mu <- colMeans(res2[, 1:3])

par_run <- MASS::mvrnorm(100, mu, Sigma = cov_mat)
par_run[par_run[, 1] > 1.5, 1] <- 1.5

par_run_list <- setNames(lapply(model, function(m) {
  fils <- list.files(cal_fold, pattern = m, full.names = TRUE)[1]
  res <- NA
  for( i in fils ) {
    cal_res <- read.csv(i)
    cal_res$run <- 1:nrow(cal_res)
    cal_res$chain <- which(fils == i)
    if(!is.data.frame(res)) {
      res <- cal_res
    } else {
      res <- rbind(res, cal_res)
    }
  }

  print(dim(res))
  # summary(res)
  # res <- res[res$qual < 1e5, ]
  
  # Select chain 3
  res2 <- res
  burn_in <- 500
  res2 <- res2[res2$run > burn_in, ]
  
  # Create distributions
  wind_mn <- mean(res2$wind_speed)
  wind_sd <- sd(res2$wind_speed)
  
  swr_mn <- mean(res2$swr)
  swr_sd <- sd(res2$swr)
  
  mpar_mn <- mean(res2[, 3])
  mpar_sd <- sd(res2[, 3])
  
  print(head(res2))
  cov_mat <- cov(res2[, 1:3])
  mu <- colMeans(res2[, 1:3])
  
  par_run <- MASS::mvrnorm(100, mu, Sigma = cov_mat)
  par_run[par_run[, 1] > 1.5, 1] <- 1.5
  return(par_run)
}), model)

# Run ensembles
lapply(model, function(m) {
  out_file <- paste0("output_", paste0(m, collapse = "_"))
  ncdf <- file.path("output", paste0(out_file, ".nc"))
  # unlink(ncdf, recursive = TRUE)
  input_yaml(config_file, "output", "file", out_file)
  
  max_members <- 100
  input_yaml(config_file, "output", "max_members", max_members)
  
  # export configuration
  export_config(config_file, m)
  
  for(j in seq_len((nrow(par_run_list[[m]])))) { #
    input_yaml_multiple(file = config_file, key1 = "scaling_factors", key2 = m,
                        key3 = "wind_speed", value = par_run_list[[m]][j, 1])
    
    input_yaml_multiple(file = config_file, key1 = "scaling_factors", key2 = m,
                        key3 = "swr", value = par_run_list[[m]][j, 2])
    if(m == "FLake") {
      input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "FLake",
                          key3 = "c_relax_C", value = par_run_list[[m]][j, 3])
    }
    if(m == "GOTM") {
      input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "GOTM",
                          key3 = "k_min", value = par_run_list[[m]][j, 3])
    }
    if(m == "Simstrat") {
      input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "Simstrat",
                          key3 = "a_seiche", value = par_run_list[[m]][j, 3])
    }
    
    
    export_config(config_file, m, dirs = F, time = F, location = F, output_settings = F, init_cond = F, meteo = TRUE, extinction = F, 
                  inflow = F, model_parameters = TRUE)
    
    # run ensemble
    run_ensemble(config_file = config_file, parallel = F,
                 model = m, return_list = FALSE,
                 add = ifelse(j == 1, FALSE, TRUE))
  }
})

out_res <- setNames(lapply(model, function(m) {
  out_file <- paste0("output_", paste0(m, collapse = "_"))
  ncdf <- file.path("output", paste0(out_file, ".nc"))
  analyse_ncdf(ncdf, m, dim = "member")
}), model)

strat_res <- setNames(lapply(out_res, function(x) x$strat[, 1:17]), model)

strat_mlt <- reshape2::melt(strat_res, id.vars = names(strat_res$FLake))
strat_mlt <- strat_mlt[strat_mlt$year == 2013, ]
strat_mlt <- strat_mlt[-1, ]
strat_mlt$L1[strat_mlt$model == "obs"] <- "Obs"
obs_idx <- which(strat_mlt$model == "obs")


obs_dat <- data.frame(xintercept = strat_mlt$StratStart[obs_idx[1]], variable = "Obs")
ggplot(strat_mlt[-obs_idx, ]) +
  geom_density(aes(StratStart, color = L1)) +
  geom_vline(data = obs_dat, aes(xintercept = xintercept, linetype = variable))

obs_dat <- data.frame(xintercept = strat_mlt$TotStratDur[obs_idx[1]], variable = "Obs")
ggplot(strat_mlt[-obs_idx, ]) +
  geom_density(aes(TotStratDur, color = L1)) +
  geom_vline(data = obs_dat, aes(xintercept = xintercept, linetype = variable))

obs_dat <- data.frame(xintercept = strat_mlt$StratEnd[obs_idx[1]], variable = "Obs")
ggplot(strat_mlt[-obs_idx, ]) +
  geom_density(aes(StratEnd, fill = L1), alpha = 0.3, color = NA) +
  geom_vline(data = obs_dat, aes(xintercept = xintercept, linetype = variable))


# Calculate Schmidt Stability
bathy <- read.csv("LakeEnsemblR_bathymetry_standard.csv")
colnames(bathy) <- c("depths", "areas")

max_ss <- setNames(lapply(model, function(m) {
  out_file <- paste0("output_", paste0(m, collapse = "_"))
  ncdf <- file.path("output", paste0(out_file, ".nc"))
  out <- load_var(ncdf, var = "temp", return = "list", dim = "member")
  if(m == "FLake") {
    bath2 <- data.frame(depths = c(0, 16.05), areas = c(rep(max(bathy$areas), 2)))
  } else {
    bath2 <- bathy
  }
  
  max.ss <- lapply(out, function(x) {
    max(ts.schmidt.stability(x, bath2, na.rm = TRUE)[, 2], na.rm = TRUE)
  })
  return(unlist(max.ss))
}), model)

fil <- "output/output_GOTM.nc"
obs <- load_var(fil, var = "temp", return = "list", dim = "model")$Obs
idx <- apply(obs, 2, function(x)sum(is.na(x)))
obs <- obs[, idx != 366]
max_obs_ss <- max(ts.schmidt.stability(obs, bathy, na.rm = TRUE)[, 2], na.rm = TRUE)


max_ss_mlt <- reshape2::melt(max_ss)

obs_dat <- data.frame(xintercept = max_obs_ss, variable = "Obs")
ggplot(max_ss_mlt) +
  geom_density(aes(value, fill = L1), alpha = 0.3, color = NA) +
  geom_vline(data = obs_dat, aes(xintercept = xintercept, linetype = variable))


# Calculate Thermocline depth
med_td <- setNames(lapply(model, function(m) {
  out_file <- paste0("output_", paste0(m, collapse = "_"))
  ncdf <- file.path("output", paste0(out_file, ".nc"))
  out <- load_var(ncdf, var = "temp", return = "list", dim = "member")
  med.td <- lapply(out, function(x) {
    median(ts.thermo.depth(x, na.rm = TRUE)[, 2], na.rm = TRUE)
  })
  return(unlist(med.td))
}), model)

fil <- "output/output_GOTM.nc"
obs <- load_var(fil, var = "temp", return = "list", dim = "model")$Obs
idx <- apply(obs, 2, function(x)sum(is.na(x)))
obs <- obs[, idx != 366]
med_obs_td <- median(ts.thermo.depth(obs, na.rm = TRUE)[, 2], na.rm = TRUE)


med_td_mlt <- reshape2::melt(med_td)

obs_dat <- data.frame(xintercept = med_obs_td, variable = "Obs")
ggplot(med_td_mlt) +
  geom_density(aes(value, fill = L1), alpha = 0.3, color = NA) +
  geom_vline(data = obs_dat, aes(xintercept = xintercept, linetype = variable))



