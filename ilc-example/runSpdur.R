################################################################
# Optimizing spweibull in rcpp
if(Sys.info()['user']=='janus829' | Sys.info()['user']=='s7m'){
  setwd("~/Research/WardProjects/spduration-optim/ilc-example")}

# Load libraries
loadPkg=function(toLoad){
  for(lib in toLoad){
  if(! lib %in% installed.packages()[,1])
    { install.packages(lib, repos='http://cran.rstudio.com/') }
  suppressMessages( library(lib, character.only=TRUE) ) }
}

packs=c("corpcor", 'microbenchmark', 'Rcpp', 'RcppArmadillo', 'ggplot2')
loadPkg(packs)

# Load core estimation functions
coreFuncs = list.files("../core-functions", pattern="*.R")
for(func in coreFuncs){source(paste0('../core-functions/', func))}

cppFuncs = list.files("../core-functions", pattern="*.cpp")
for(func in cppFuncs){sourceCpp(paste0('../core-functions/', func))}
################################################################


################################################################
# Run glm - spdur - spdurCpp
load("irc-data-mod.rda")

mdl1Spdur = function() { 
  spdur(duration ~ 1 + log10(i.matl.conf.DIStGOV.l1+1) + 
    log10(i.matl.coop.GOVtGOV.l1+1),
    atrisk ~ 1 + ldr.irregular + ldr.foreign + log10(mths.in.power+1),
    data=irc.data, silent=TRUE) 
}

mdl1SpdurCpp = function() { 
  spdurCpp(duration ~ 1 + log10(i.matl.conf.DIStGOV.l1+1) + 
    log10(i.matl.coop.GOVtGOV.l1+1),
    atrisk ~ 1 + ldr.irregular + ldr.foreign + log10(mths.in.power+1),
    data=irc.data, silent=TRUE) 
}

mdl1Glm = function() {
  glm(failure ~ log10(i.matl.conf.DIStGOV.l1+1) + 
    log10(i.matl.coop.GOVtGOV.l1+1) + ldr.irregular + ldr.foreign + 
    log10(mths.in.power+1) + duration + I(duration^2) + I(duration^3),
    data=irc.data, family="binomial")
}

# Check to make sure results between the two are equivalent
spdurRes=mdl1Spdur()
spdurcppRes=mdl1SpdurCpp()
cbind(spdurRes$'coefficients', spdurRes$'se')
cbind(spdurcppRes$'coefficients', spdurcppRes$'se')
################################################################

################################################################
# Time trials
if(!'spdurTimes.rda' %in% list.files()){
  times = microbenchmark(
    mdl1Spdur(), mdl1SpdurCpp(), mdl1Glm(), 
    times=100, unit="s")
  save(times, file='spdurTimes.rda')
}  

if(!'times' %in% ls()){ load('spdurTimes.rda') }
print(times)
################################################################