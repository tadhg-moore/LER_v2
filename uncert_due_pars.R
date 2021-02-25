####################################################################################################
#                                                                                                  #
#                   estimate uncretainty from parameter for lake Feeagh                            #
#                                                                                                  #
####################################################################################################
# Johannes Feldbauer                                                                               #
# created: 23.04.2021                                                                              #
# last edited: 23.02.2021                                                                          #
####################################################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#install.packages("remotes")
remotes::install_github("GLEON/rLakeAnalyzer")
remotes::install_github("aemon-j/GLM3r", ref = "v3.1.1")
remotes::install_github("USGS-R/glmtools", ref = "ggplot_overhaul")
remotes::install_github("aemon-j/FLakeR", ref = "inflow")
remotes::install_github("aemon-j/GOTMr")
remotes::install_github("aemon-j/gotmtools")
remotes::install_github("aemon-j/SimstratR")
remotes::install_github("aemon-j/MyLakeR")

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

# Format met file
# met <- read.csv("LakeEnsemblR_meteo_standard.csv")
# met$Rainfall_millimeterPerDay <- met$Precipitation_meterPerSecond * 8.64e7
# met$Snowfall_millimeterPerDay <- met$Snowfall_meterPerDay * 1000
# summary(met)
# rem_cols <- which(colnames(met) %in% c("Precipitation_meterPerSecond", "Snowfall_meterPerDay"))
# met <- met[, -rem_cols]
# write.csv(met, "LakeEnsemblR_meteo_standardv2.csv", row.names = FALSE, quote = FALSE)
# input_yaml(config_file, "meteo", "file", "LakeEnsemblR_meteo_standardv2.csv")
############

export_config(config_file = config_file, model = model, dirs = F, time = F, location = F, output_settings = T,
              meteo = F, init_cond = F, extinction = F, inflow = F, model_parameters = F)

run_ensemble(config_file = config_file, model = model)

ncdf <- "output/ensemble_output_pars.nc"
calc_fit(ncdf, model = model)

# MCMC method
# model <-  c("GOTM") 
model <-  c("GLM", "GOTM", "Simstrat", "FLake") 
cali_ensemble(config_file = config_file, num = 5, cmethod = "MCMC",
              model = model, parallel = F, out_f = "testv1")

library(parallel)
ncores <- detectCores()
clust <- makeCluster(ncores)
clusterExport(clust, varlist = list("config_file"),
              envir = environment())
clusterEvalQ(clust, library(LakeEnsemblR))

resMCMC <- parLapply(clust, model, function(m) {
  cali_ensemble(config_file = config_file, num = 5000, cmethod = "MCMC",
                model = m, out_f = "mcmc_v1")
})
# stop cluster
stopCluster(clust)


# 


# Calibration results -----
res <- read.csv("mcmc_v1/FLake_MCMC_202102250222.csv")
res$run <- 1:nrow(res)
mlt <- reshape2::melt(res, id.vars = c("qual", "run"))
ggplot(mlt) +
  # geom_point(aes(run, value)) +
  geom_line(aes(run, value)) +
  geom_smooth(aes(run, value)) +
  facet_wrap(~variable, scales = "free_y") +
  scale_y_log10()

ggplot(mlt) +
  geom_point(aes(value, qual)) +
  # geom_line(aes(run, value)) +
  geom_smooth(aes(value, qual)) +
  facet_wrap(~variable, scales = "free_x") +
  scale_x_log10()

idx <- which.min(res$qual)
res[idx,]
##############

cali_files <- list.files("lhc_v1/", full.names = T)
lhc_res <- load_LHC_results(config_file, model = model, res_files = cali_files)
str(lhc_res)

glm <- lhc_res$GLM
ggplot(glm) +
  geom_point(aes(rmse, LL)) +
  scale_y_log10()

mlt <- reshape2::melt(glm[, 7:10], id.vars = "LL")
ggplot(mlt, aes(value, LL)) +
  geom_point() +
  facet_wrap(~variable, scales = "free") +
  geom_smooth() +
  scale_y_log10()

# number of parameter sets to run this for
n_par <- 20

# names of the calibration results
cali <- file.path("cali",
                  list.files("cali")[grepl("202009040748", list.files("cali"))])
# read in results
d_cali <- load_LHC_results(config_file, model, cali)

# find optimum parameters
opt_par <- sapply(model, function(m)
  d_cali[[m]][which.min(d_cali[[m]]$rmse), c(2, (ncol(d_cali[[m]])-2):ncol(d_cali[[m]]))] )

# find range of parameters with good fit
par_rang <- setNames(lapply(model, function(m)
  d_cali[[m]][d_cali[[m]]$rmse <= quantile(d_cali[[m]]$rmse, 0.1), 8:10]), model)
par_rang <- setNames(lapply(model, function(m)
  apply(par_rang[[m]],2, range)), model)

# draw from a lhc to generate sets of parameters
par_run <- setNames(lapply(model, function(m)
  Latinhyper(t(par_rang[[m]]), n_par)), model)


# loop to change parameters and run LER
for (i in seq_len(n_par)) {
  
  for (m in model) {
    
    input_yaml_multiple(file = config_file, key1 = "scaling_factors", key2 = m,
                        key3 = "wind_speed", value = par_run[[m]][i, 1])
    
    input_yaml_multiple(file = config_file, key1 = "scaling_factors", key2 = m,
                        key3 = "swr", value = par_run[[m]][i, 2])
  }
  input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "GLM",
                      key3 = "mixing/coef_mix_hyp", value = par_run$GLM[i, 3])
  input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "GOTM",
                      key3 = "k_min", value = par_run$GOTM[i, 3])
  input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "Simstrat",
                      key3 = "a_seiche", value = par_run$Simstrat[i, 3])
  input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "FLake",
                      key3 = "c_relaxe_C", value = par_run$FLake[i, 3])
  input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "MyLake",
                      key3 = "Phys.par/C_shelter", value = par_run$MyLake[i, 3])
  # export configuration
  export_config(config_file, model)
  
  # run ensemble
  run_ensemble(config_file = config_file, parallel = TRUE,
               model = model, return_list = FALSE,
               add = ifelse(i == 1, FALSE, TRUE))
  
  
}

temp_09 <- list()
temp_32 <- list()
for (i in 1:6) {
  temp_tmp <- load_var("output/ensemble_output_pars.nc", "temp", dim = "member", dim_index = i)
  temp_tmp <- melt(temp_tmp, id.vars = c("datetime"))
  temp_tmp_09 <- subset(temp_tmp, temp_tmp$variable == "wtr_0.9")
  temp_tmp_32 <- subset(temp_tmp, temp_tmp$variable == "wtr_32")
  temp_09[[c(model, "Observations")[i]]] <- temp_tmp_09
  temp_32[[c(model, "Observations")[i]]] <- temp_tmp_32
}
rm(temp_tmp, temp_tmp_09, temp_tmp_32)


res_09 <- lapply(model, function(m){
  mean <-  aggregate(list(temp = temp_09[[m]]$value), by = list(datetime = temp_09[[m]]$datetime),
            mean, na.rm = TRUE)
  q5 <-  aggregate(list(temp = temp_09[[m]]$value), by = list(datetime = temp_09[[m]]$datetime),
                     quantile, 0.05, na.rm = TRUE)
  q95 <-  aggregate(list(temp = temp_09[[m]]$value), by = list(datetime = temp_09[[m]]$datetime),
                     quantile, 0.95, na.rm = TRUE)
  mint <-  aggregate(list(temp = temp_09[[m]]$value), by = list(datetime = temp_09[[m]]$datetime),
                   min, na.rm = TRUE)
  maxt <-  aggregate(list(temp = temp_09[[m]]$value), by = list(datetime = temp_09[[m]]$datetime),
                    max, na.rm = TRUE)
  datetime <- unique(temp_09[[m]]$datetime)
  return(data.frame(mean = mean$temp,
                    q5 = q5$temp,
                    q95 = q95$temp,
                    min = mint$temp,
                    max = maxt$temp,
                    datetime = datetime,
                    model = m,
                    depth = 0.9))
})
res_32 <- lapply(model, function(m){
  mean <-  aggregate(list(temp = temp_32[[m]]$value), by = list(datetime = temp_32[[m]]$datetime),
                     mean, na.rm = TRUE)
  q5 <-  aggregate(list(temp = temp_32[[m]]$value), by = list(datetime = temp_32[[m]]$datetime),
                   quantile, 0.05, na.rm = TRUE)
  q95 <-  aggregate(list(temp = temp_32[[m]]$value), by = list(datetime = temp_32[[m]]$datetime),
                    quantile, 0.95, na.rm = TRUE)
  mint <-  aggregate(list(temp = temp_32[[m]]$value), by = list(datetime = temp_32[[m]]$datetime),
                     min, na.rm = TRUE)
  maxt <-  aggregate(list(temp = temp_32[[m]]$value), by = list(datetime = temp_32[[m]]$datetime),
                     max, na.rm = TRUE)
  datetime <- unique(temp_32[[m]]$datetime)
  return(data.frame(mean = mean$temp,
                    q5 = q5$temp,
                    q95 = q95$temp,
                    min = mint$temp,
                    max = maxt$temp,
                    datetime = datetime,
                    model = m,
                    depth = 32))
})
res_all <- rbind(melt(res_09, id.vars = c("mean", "q5", "q95", "min", "max",
                                          "datetime", "model", "depth")),
             melt(res_32, id.vars = c("mean", "q5", "q95", "min", "max",
                                      "datetime", "model", "depth")))
res_all$min[is.infinite(res_all$min)] <- NA
res_all$max[is.infinite(res_all$max)] <- NA

colfunc <- colorRampPalette(RColorBrewer::brewer.pal(length(model), "Set2"))

p_res <- ggplot(res_all) +
  geom_ribbon(aes(x = datetime, ymin = min, ymax = max, fill = model), alpha = 0.45) +
  geom_line(aes(x = datetime, y = mean, col = model), alpha = 0.95, lwd = 1.15) +
  theme_classic(base_size = 21) + xlab("Date") +
  ylab("Water temperature (°C)") + facet_grid(rows = vars(depth)) +
  scale_color_manual("Model", values = colfunc(5)[c(1:3,5,4)]) +
  scale_fill_manual("Model", values = colfunc(5)[c(1:3,5,4)]) +
  scale_x_datetime(date_labels = "%Y-%m-%d")

ggsave("output/plot_paper_pars.png", p_res, device = "png", width = 13, height = 9)
