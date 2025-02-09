library(plFA)
data("LSAT")
mod <- "eta =~ y1 + y2 + y3 + y4 + y5"
fit <- cfa(mod, LSAT, std.lv = TRUE, estimator = "PML", test = "mean.var.adjusted")
