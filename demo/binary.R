library(plFA)
data("LSAT")
mod <- "eta =~ y1 + y2 + y3 + y4 + y5"
fit <- cfa(mod, LSAT, std.lv = TRUE, estimator = "PML")


lavaan:::lav_lavaan_step15_baseline(lavoptions = fit@Options,
                                       lavsamplestats = fit@SampleStats,
                                       lavdata = fit@Data,
                                       lavcache = fit@Cache,
                                       lavh1 = fit@h1,
                                       lavpartable = fit@ParTable)
