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
library(readr)
library(lubridate)


setwd("feeagh") # Change working directory to example folder

wtemp <- read_csv("LakeEnsemblR_wtemp_profile_standard.csv")
start <- "2013-06-12 00:00:00"
ncdf <- file.path("output", paste0("output_params", ".nc"))

config_file <- 'LakeEnsemblR_params_uncert.yaml'
model <-  c("FLake", "GLM", "GOTM", "Simstrat", "MyLake") 
n_mem <- 100

export_config(config_file, model)

set.seed(123)
samp_fact <- data.frame(c_relax_C = rnorm(n_mem, 0.0009, 0.01^2),
                        mixing.coef_mix_hyp = rnorm(n_mem, 1.4, 0.5^2),
                        turb_param.k_min = rnorm(n_mem, 9e-5, 0.00001),
                        a_seiche = rnorm(n_mem, 0.001, 0.0001),
                        Phys.par.C_shelter = rnorm(n_mem, 0.25, 0.1))

for(i in seq_len(n_mem)) {
  
  suppressMessages({
    input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "FLake",
                        key3 = "c_relax_C", value = samp_fact$c_relax_C[i])
    input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "GLM",
                        key3 = "mixing/coef_mix_hyp", value = samp_fact$mixing.coef_mix_hyp[i])
    input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "GOTM",
                        key3 = "k_min", value = samp_fact$turb_param.k_min[i])
    input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "Simstrat",
                        key3 = "a_seiche", value = samp_fact$a_seiche[i])
    input_yaml_multiple(file = config_file, key1 = "model_parameters", key2 = "MyLake",
                        key3 = "Phys.par/C_shelter", value = samp_fact$Phys.par.C_shelter[i])
    
    export_config(config_file = config_file, model = model, dirs = F, time = F, location = F, output_settings = F,
                  meteo = F, init_cond = F, extinction = F, inflow = F, model_parameters = T)
  })
  
  run_ensemble(config_file, model, add = (i != 1), parallel = TRUE)
  print(i)
  
}

plot_heatmap(ncdf)

out <- load_var(ncdf, "temp", return = "array")

out_list <- lapply(1:5, function(x) {
  apply(out[, x, , 3], 2, function(y) quantile(y, c(0.025, 0.5, 0.975), na.rm = TRUE))
})
names(out_list) <- dimnames(out)[2][[1]][1:5]
obs <- out[1, 6, , 3]
obs_df <- data.frame(datetime = as.POSIXct(names(obs), tz = "UTC"),
                     temp = obs)

mlt1 <- reshape2::melt(out_list)
wid <- tidyr::pivot_wider(mlt1, values_from = 3, names_from = 1)
colnames(wid) <- c("datetime", "model", "c025", "c50", "c975")
wid$datetime <- as.POSIXct(as.character(wid$datetime))

ggplot(wid) +
  geom_ribbon(aes(datetime, ymin = c025, ymax = c975, fill = model), alpha = 0.3)

# Surface temp
out_var <- lapply(1:5, function(x) {
  data.frame(datetime = as.POSIXct(dimnames(out)[3][[1]], tz = "UTC"),
             vari = apply(out[, x, , 3], 2, function(y) var(y, na.rm = TRUE)))
})
names(out_var) <- dimnames(out)[2][[1]][1:5]
mlt2 <- reshape2::melt(out_var, id.vars = "datetime")
ggplot(mlt2) +
  geom_line(aes(datetime, value, color = L1))

# Bottom temp
out_var <- lapply(1:5, function(x) {
  data.frame(datetime = as.POSIXct(dimnames(out)[3][[1]], tz = "UTC"),
             vari = apply(out[, x, , 86], 2, function(y) var(y, na.rm = TRUE)))
})
names(out_var) <- dimnames(out)[2][[1]][1:5]
mlt2 <- reshape2::melt(out_var, id.vars = "datetime")
ggplot(mlt2) +
  geom_line(aes(datetime, value, color = L1))
