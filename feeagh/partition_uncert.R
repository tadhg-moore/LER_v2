# variance plots

ucert <- c("Initial condtions", "Boundary conditions", "Parameters")
ncdf_fils <- c("output_ic", "output_bc", "output_params")

for( i in seq_len(length(ucert))) {
  ncdf <- file.path("output", paste0(ncdf_fils[i], ".nc"))
  out <- load_var(ncdf, "temp", return = "array")
  # Surface temp
  out_var <- lapply(1:5, function(x) {
    data.frame(datetime = as.POSIXct(dimnames(out)[3][[1]], tz = "UTC"),
               vari = apply(out[, x, , 3], 2, function(y) var(y, na.rm = TRUE)))
  })
  names(out_var) <- dimnames(out)[2][[1]][1:5]
  mlt1 <- reshape2::melt(out_var, id.vars = "datetime")
  mlt1$uncert <- ucert[i]
  if( i == 1) {
    df2 <- mlt1
  } else {
    df2 <- rbind(df2, mlt1)
  }
}

idx <- which(df2$datetime < "2013-06-28")
p1 <- ggplot(df2[idx, ]) +
  geom_line(aes(datetime, value, color = L1), size = 1.1) +
  facet_wrap(~uncert, nrow = 1, scales = "free_y") +
  scale_color_brewer(palette = "Set1", drop = FALSE) +
  ylab("Variance") +
  xlab("") +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 1, shape = 0))) +
  labs(color = "Model") +
  theme_classic(base_size = 26)


# Bottom temperature
for( i in seq_len(length(ucert))) {
  ncdf <- file.path("output", paste0(ncdf_fils[i], ".nc"))
  out <- load_var(ncdf, "temp", return = "array")
  # Surface temp
  out_var <- lapply(1:5, function(x) {
    data.frame(datetime = as.POSIXct(dimnames(out)[3][[1]], tz = "UTC"),
               vari = apply(out[, x, , 86], 2, function(y) var(y, na.rm = TRUE)))
  })
  names(out_var) <- dimnames(out)[2][[1]][1:5]
  mlt1 <- reshape2::melt(out_var, id.vars = "datetime")
  mlt1$uncert <- ucert[i]
  if( i == 1) {
    df3 <- mlt1
  } else {
    df3 <- rbind(df3, mlt1)
  }
}


idx <- which(df3$datetime < "2013-06-28")
p2 <- ggplot(df3[idx, ]) +
  geom_line(aes(datetime, value, color = L1), size = 1.1) +
  facet_wrap(~uncert, nrow = 1, scales = "free_y") +
  scale_color_brewer(palette = "Set1", drop = FALSE) +
  ylab("Variance") +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 1, shape = 0))) +
  labs(color = "Model") +
  xlab("Date") +
  theme_classic(base_size = 26)

library(ggpubr)

g1 <- ggarrange(p1, p2, align = "v", ncol = 1, common.legend = TRUE, labels = "AUTO")
ggsave("variance_plots.png", dpi = 120,width = 384,height = 280, units = 'mm')
