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
ncdf <- file.path("output", paste0("output_bc", ".nc"))

config_file <- 'LakeEnsemblR_bc_uncert.yaml'
model <-  c("FLake", "GLM", "GOTM", "Simstrat", "MyLake") 
n_mem <- 100

met_vari <- 0.44 

export_config(config_file, model)

set.seed(123)
samp_fact <- data.frame(wind = rnorm(n_mem, 1, met_vari^2),
                        swr = rnorm(n_mem, 1, met_vari^2))

for(i in seq_len(n_mem)) {
  
  for (m in model) {
    input_yaml_multiple(file = config_file, key1 = "scaling_factors", key2 = m,
                        key3 = "wind_speed", value = samp_fact$wind[i])
    
    input_yaml_multiple(file = config_file, key1 = "scaling_factors", key2 = m,
                        key3 = "swr", value = samp_fact$swr[i])
  }
  suppressMessages(export_config(config_file = config_file, model = model, dirs = F, time = F, location = F, output_settings = F,
                meteo = T, init_cond = F, extinction = F, inflow = F, model_parameters = F))
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


out_var <- lapply(1:5, function(x) {
  data.frame(datetime = as.POSIXct(dimnames(out)[3][[1]], tz = "UTC"),
             vari = apply(out[, x, , 3], 2, function(y) var(y, na.rm = TRUE)))
})
names(out_var) <- dimnames(out)[2][[1]][1:5]
mlt2 <- reshape2::melt(out_var, id.vars = "datetime")
ggplot(mlt2) +
  geom_line(aes(datetime, value, color = L1))

out_var <- lapply(1:5, function(x) {
  data.frame(datetime = as.POSIXct(dimnames(out)[3][[1]], tz = "UTC"),
             vari = apply(out[, x, , 86], 2, function(y) var(y, na.rm = TRUE)))
})
names(out_var) <- dimnames(out)[2][[1]][1:5]
mlt2 <- reshape2::melt(out_var, id.vars = "datetime")
ggplot(mlt2) +
  geom_line(aes(datetime, value, color = L1))
