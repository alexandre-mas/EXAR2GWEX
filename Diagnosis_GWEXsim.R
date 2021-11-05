###############################################################
########     DIAGNOSIS OF GWEX SIMULATED SERIES      ##########
###############################################################

## Alexandre MAS - 10/2021 - Project EXAR2
## Relies (mostly) on computational and ploting tools developed for EXAR 1
## Objectivs : comparing observed and simulated series at gages and basin level with several statistics, on various time scales
## For rainfall : 
## - reproduction of wet/dry days probability and wet/dry spells occurence and length
## - reproduction of 6H, 1D, 3D mean and max precipitation
## - reproduction of quantile for various accumulation time 
## For temperature : 
## 

## Source R scripts containing functions
source("C:/EXAR/SVN/GWEX/computeStat_lib.r")
source("C:/EXAR/SVN/GWEX/plotStat_lib.r")

## 18/10/2021
## Diagnosis of first series generated without accounting for weather type on the gages of group 3 (catchment Ticino)
## The 3 parameters of EGPD at each gage were estimated by Abubakar with the ROI method
## Figures and results saved in "C:/EXAR_2/ETUDES/Diagnostic_GWEX/Firstsim_noWT/"

## Yearly prec (mean sim)
# data
year.sum.obs = colMeans(rollapply(myObs@obs, width=365, by=365, mean, na.rm=T, by.column=T), na.rm = T)*365
year.sum.sim = colMeans(do.call(rbind, lapply(1:10, function(x) colMeans(rollapply(mySim@sim[,,x], width=365, by=365, sum, na.rm=F, by.column=T), na.rm = T))))
# plot
layout(matrix(c(1,2), nrow=1), widths=c(6,1))
plot(year.sum.obs, type='l', xaxt="n", ylab = "Mean yearly prec (mm)", xlab = NA)
axis(1, at=1:45, labels = names(year.sum.obs), las=2)
lines(year.sum.sim, col='red')
legend("topleft", bty="n", col=c("black", "red"), lty=1, legend = c("Obs", "Sim"))
boxplot((year.sum.sim-year.sum.obs)/year.sum.obs*100, ylab = "Relative error (%)")

## Yearly prec (10 sim boxplot)
# data
year.sum.sim.df.v0 = as.data.frame(do.call(rbind, lapply(1:10, function(x) colMeans(rollapply(mySim@sim[,,x], width=365, by=365, sum, na.rm=F, by.column=T), na.rm = T))))
RE.full.v0= apply(year.sum.sim.df.v0,1,function(x) (x - year.sum.obs)/year.sum.obs*100)

# plot
layout(matrix(c(1,2), nrow=1), widths=c(5,1))
par(mar=c(5,4,4,0))
boxplot(year.sum.sim.regioROI.th01lag4.df, xaxt="n", xlab="",ylab = "Mean yearly prec (mm)")
axis(1, at = 1:45, labels = names(year.sum.obs), las=2)
points(year.sum.obs, col='red')
legend("topright", bty="n",legend=c("sim.regio.th01.lag4","obs"), pch=c(12,1), col=c("black",'red'), pt.cex = 1.4)
title(main="Group 3")
boxplot(list(RE.full.th02, RE.full.th01), ylab = "Relative error (%)", names = c("th = 0.2","th = 0.1"))

## Season prec
sai.obs = c("WIN","WIN","SPR","SPR","SPR","SUM","SUM","SUM","AUT","AUT","AUT","WIN")[as.numeric(format(myObs@date, "%m"))]
for (sai in unique(vec.sai)){
season.sum.obs = colMeans(rollapply(myObs@obs[sai.obs==sai,], width=91, by=91, mean, na.rm=T, by.column=T), na.rm = T)*91
season.sum.sim = colMeans(do.call(rbind, lapply(1:10, function(x) colMeans(rollapply(mySim@sim[vec.sai==sai,,x], width=91, by=91, sum, na.rm=F, by.column=T), na.rm = T))))
layout(matrix(c(1,2), nrow=1), widths=c(6,1))
plot(season.sum.obs, type='l', xaxt="n", ylab = "Mean season prec (mm)", xlab = NA, main = sai)
axis(1, at=1:45, labels = names(season.sum.obs), las=2)
lines(season.sum.sim, col='red')
legend("topleft", bty="n", col=c("black", "red"), lty=1, legend = c("Obs", "Sim"))
boxplot((season.sum.sim-season.sum.obs)/season.sum.obs*100, ylab = "Relative error (%)")
delta_seasonprec[[sai]] = season.sum.obs - season.sum.sim
}

## ecart pour les 4 saisons sur 1 graph
layout(matrix(c(1,2), nrow=1), widths=c(6,1))
plot(NA, xlim=c(0,46), ylim=c(-32,145), ylab = "Mean season diff (mm)", xaxt="n", xlab = "")
axis(1, at=1:45, labels = names(season.sum.obs), las=2)
points(delta_seasonprec$WIN)
points(delta_seasonprec$SPR, col='red')
points(delta_seasonprec$SUM, col='blue')
points(delta_seasonprec$AUT, col='green')
boxplot(delta_seasonprec, border= c("black","red","blue","green"))

##
plot.daily.stats.by.month(obs = myObs@obs, sim = mySim@sim, obs.date = myObs@date, sim.date = mySim@date, type = "WDF", obs.th = 0.5, sim.th = 0)


## 03/11 Diagnosis of the 1x10000 years series on group 3 and 4
library(GWEXWT)
gr=3

# load data
load(paste0("C:/EXAR_2/APPLICATION/GWex/Fit.GWEX/Prec_noWT/SmoothedROI_byseasons/Final/fit.gwex.gr",gr,".th01lag2.Rdata"))
load(paste0("C:/EXAR_2/APPLICATION/GWex/Sim.GWEX/SmoothedROI_byseasons/Final/sim.gwex.gr",gr,"_10000years.1sim.RData"))
# transform 1x10000 years to 10x1000 years
df = which(mySim@date %in% as.Date(c('0999-12-31', '1999-12-31', '2999-12-31', '3999-12-31', '4999-12-31', '5999-12-31', '6999-12-31', '7999-12-31', '8999-12-31', '9999-12-31')))
dd = c(1,df[1:9]+1)
df[c(1,3,5,7,9)] = df[c(1,3,5,7,9)]-1
array.sim = array(NA, dim = c(365242,ncol(myObs@obs),10))
for (y in 1:10){
array.sim[,,y] = mySim@sim[dd[y]:df[y],]
}

# plot
par(mfrow=c(2,3))
plot.scatter.fun(obs = myObs@obs, sim = array.sim, obs.date =  myObs@date, sim.date =  mySim@date[1:365242], n.day.cumul = 1, lab = colnames(myObs@obs), n.fig.page = 6, main = "Annual 1day maxima")
plot.scatter.fun(obs = myObs@obs, sim = array.sim, obs.date =  myObs@date, sim.date =  mySim@date[1:365242], n.day.cumul = 3, lab = colnames(myObs@obs), n.fig.page = 6, main = "Annual 1day maxima")
plot.cdf.spell.length(obs = myObs@obs, sim = array.sim, obs.date =  myObs@date, sim.date =  mySim@date[1:365242], type="wet", obs.th = 0.1, sim.th = 0, lab = colnames(myObs@obs))






