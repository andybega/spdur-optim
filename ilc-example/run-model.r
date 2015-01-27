#
#   Example code from ILC 2014 papers
#

setwd("~/Work/spduration-optim/ilc-example")

library(corpcor)
library(lineprof)
library(microbenchmark)

# Load core estimation functions
core_funcs <- list.files("../core-functions", pattern="*.R")
sapply(paste0("../core-functions/", core_funcs), source, .GlobalEnv)

# Load example data
load("irc-data-mod.rda")

# Theme 1 with all data
model1 <- function() {
  spdur(duration ~ 1 + log10(i.matl.conf.DIStGOV.l1+1) + 
          log10(i.matl.coop.GOVtGOV.l1+1),
        atrisk ~ 1 + ldr.irregular + ldr.foreign + log10(mths.in.power+1),
        data=irc.data, silent=TRUE)
}

# Run time; use times>1 to eliminate overhead variance
microbenchmark(model1(), times=1, unit="s")

# Line by line profiling; might give different results each run
m1_l1 <- lineprof({
  spdur(duration ~ 1 + log10(i.matl.conf.DIStGOV.l1+1) + 
          log10(i.matl.coop.GOVtGOV.l1+1),
        atrisk ~ 1 + ldr.irregular + ldr.foreign + log10(mths.in.power+1),
        data=irc.data, silent=TRUE)
}, interval=0.1)

shine(m1_l1)
