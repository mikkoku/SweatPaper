setwd("c:/stuff/a/julia/SweatPaper/")
library(spatstat)
xy <- read.csv("data/glands.csv")
meta <- read.csv("data/meta.csv")

subject <- 20
xy <- subset(xy, subjectid==subject)
window <- subset(meta, subjectid==subject)
pp <- ppp(xy$x, xy$y, window=owin(c(window$x0, window$x1), c(window$y0, window$y1)))

plot(pp)
plot(pcf(pp))
