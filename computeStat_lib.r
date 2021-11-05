###===============================###===============================###
### Guillaume Evin
### 29/02/2016, Grenoble
###  LTHE
### guillaume_evin@yahoo.fr
###
### The following functions compute standard statistics for precipitation
### time series
###===============================###===============================###
library(zoo)

#==============================================================================
# Dry day frequency DDF
#==============================================================================
dry.day.frequency = function(mat.prec, # matrix of precipitation
                             th # th above which we consider that a day is wet
){
  # dry day frequency (computeStat_lib)
  return(apply(mat.prec,2,function(x,th) mean(x<=th,na.rm=T),th))
}


#==============================================================================
# Wet day frequency WDF
#==============================================================================
wet.day.frequency = function(mat.prec, # matrix of precipitation
                             th # th above which we consider that a day is wet
){
  # wet day frequency (computeStat_lib)
  return(apply(mat.prec,2,function(x,th) mean(x>th,na.rm=T),th))
}


#==============================================================================
# mean wet day intensity MWDF
#==============================================================================
mean.wet.day.intensity = function(mat.prec, # matrix of precipitation
                                  th # th above which we consider that a day is wet
){
  # mean wet day intensity (computeStat_lib)
  return(apply(mat.prec,2,function(x,th) mean(x[x>th],na.rm=T),th))
}


#==============================================================================
# standard deviation of wet day intensity SDWDF
#==============================================================================
sd.wet.day.intensity = function(mat.prec, # matrix of precipitation
                                th # th above which we consider that a day is wet
){
  # standard deviation of wet day intensity (computeStat_lib) 
  return(apply(mat.prec,2,function(x,th) sd(x[x>th],na.rm=T),th))
}


#==============================================================================
# transition probabilities on several consecutive days
#==============================================================================
lag.trans.proba.vector = function(vec.prec, # vector of precipitation
                                  vec.dates, #vector of dates
                                  th, # threshold
                                  nlag, # number of lag days
                                  dayScale # time resolution
){
  # lag.trans.proba.vector (computeStat_lib): number of consecutive days on 
  # which we compute the probabilities
  
  ndays = nlag+1
  
  # filter consecutive dates
  has.all.days = diff(vec.dates,lag=nlag)==(nlag*dayScale)
  ind.has.all.days = which(has.all.days)
  
  # is.wet.lag gives the boolean values indicating if the time step 0, 1, ... ,0+nlag
  # are wet or not
  is.wet.lag = matrix(nrow=length(ind.has.all.days),ncol=ndays)
  for(i.lag in 0:nlag) is.wet.lag[,(i.lag+1)] = vec.prec[ind.has.all.days+i.lag]>th
  
  # no nan
  nz = apply(is.wet.lag,1,function(x) !any(is.na(x)))
  x.nz = is.wet.lag[nz,]
  
  # matrix of joints probabilities
  comb.pr = expand.grid(lapply(numeric(ndays), function(x) c(F,T)))
  names(comb.pr) = paste0('t',(-nlag:0))
  comb.grid = as.matrix(comb.pr)
  comb.pr$P = NA
  
  for(i.comb in 1:nrow(comb.grid)){
    comb = comb.grid[i.comb,]
    comb.pr$P[i.comb] = mean(apply(x.nz,1, function(x) all(x==comb)))
  }
  
  # matrix of conditional probabilities Pr(i0=wet|i-1=.,i-2=.,...)
  nPr = nrow(comb.pr)
  prF = comb.pr$P[1:(nPr/2)]
  prT = comb.pr$P[(nPr/2+1):nPr]
  PrCond = prT/(prF+prT)
  comb.prT = comb.pr[(nPr/2+1):nPr,1:nlag,drop=F]
  comb.prT$P = PrCond
  
  # return  wconditional probabilities
  return(comb.prT)
}

#==============================================================================
lag.trans.proba.proba = function(mat.prec, # matrix of precipitation
                                 vec.dates, # vector of dates
                                 th, # threshold above which we consider that a day is wet
                                 nlag, # number of lag days
                                 dayScale=1 # time resolution
){
  # lag.trans.proba.proba (computeStat_lib)
  return(apply(mat.prec,2,lag.trans.proba.vector,vec.dates,th,nlag,dayScale))
}


#==============================================================================
# wet-wet transition probability WWTP
#==============================================================================
wet.wet.trans.proba.vector = function(vec.prec, # vector of precipitation
                                      vec.dates, # vector of dates
                                      th, # threshold above which we consider that a day is wet
                                      dayScale # time resolution
){
  # wet.wet.trans.proba.vector (computeStat_lib): wet-wet transition probability WWTP
  
  # filter consecutive dates
  has.next.day = diff(vec.dates,lag=1)==dayScale
  
  # time t and lag time t+1
  t1 = which(has.next.day)
  t2 = t1+1
  
  # no nan
  nz = (!is.na(vec.prec[t1])) & (!is.na(vec.prec[t2]))
  t1nz = t1[nz]
  t2nz = t2[nz]
  
  # return  wet.wet transition probability
  return(mean(vec.prec[t1nz]>th&vec.prec[t2nz]>th)/mean(vec.prec[t1nz]>th))
}

#==============================================================================
wet.wet.trans.proba = function(mat.prec, # matrix of precipitation
                               vec.dates, # vector of dates
                               th, # threshold above which we consider that a day is wet
                               dayScale=1 # time resolution
){
  # wet.wet.trans.proba (computeStat_lib): wet.wet transition probability
  return(apply(mat.prec,2,wet.wet.trans.proba.vector,vec.dates,th,dayScale))
}


#==============================================================================
# dry-wet transition probability DWTP
#==============================================================================
dry.wet.trans.proba.vector = function(vec.prec, # vector of precipitation
                                      vec.dates, # vector of dates
                                      th, # threshold above which we consider that a day is wet
                                      dayScale # time resolution
){
  # dry.wet.trans.proba.vector(computeStat_lib): dry-wet transition probability DWTP
  
  # filter consecutive dates
  has.next.day = diff(vec.dates,lag=1)==dayScale
  
  # time t and lag time t+1
  t1 = which(has.next.day)
  t2 = t1+1
  
  # no nan
  nz = (!is.na(vec.prec[t1])) & (!is.na(vec.prec[t2]))
  t1nz = t1[nz]
  t2nz = t2[nz]
  
  # return  wet.wet transition probability
  return(mean(vec.prec[t1nz]<=th&vec.prec[t2nz]>th)/mean(vec.prec[t1nz]<=th))
}

#==============================================================================
dry.wet.trans.proba = function(mat.prec, # matrix of precipitation
                               vec.dates, # vector of dates
                               th, # th above which we consider that a day is wet
                               dayScale=1 # time resolution
){
  # dry.wet.trans.proba (computeStat_lib): wet.wet transition probability
  return(apply(mat.prec,2,dry.wet.trans.proba.vector,vec.dates,th,dayScale))
}



#==============================================================================
# wrapper to get several statistics for daily precipitation
#==============================================================================
get.daily.stat = function(mat, # matrix of obs/sim
                          vec.dates, # vector of dates
                          type, # one of DDF WDF MWDF SDWDF WWTP DWTP
                          th # threshold above which we consider that a day is wet
){
  # get.daily.stat (computeStat_lib):wrapper to get several statistics for daily variables
  
  # if type is a function, we apply the function directly
  if(is.function(type)){
    return(apply(mat,2,type))
  }else{
    # otherwise, we have different options (for precipitation variables)
    switch(type,
           DDF = dry.day.frequency(mat,th),
           WDF = wet.day.frequency(mat,th),
           MWDF = mean.wet.day.intensity(mat,th),
           SDWDF = sd.wet.day.intensity(mat,th),
           WWTP = wet.wet.trans.proba(mat,vec.dates,th),
           DWTP = dry.wet.trans.proba(mat,vec.dates,th))
  }
}

  


#==============================================================================
# aggregation at a monthly scale
#==============================================================================
get.monthly.agg = function(mat, # matrix of obs/sim
                           vec.dates, # vector of dates
                           fun.agg # function to apply
){
  # get.monthly.TS (computeStat_lib): aggregation at a monthly scale
  # number of stations
  
  # number of spatial entities (gauges, basins)
  p = ncol(mat)
  
  # year/month
  short.date = strftime(vec.dates, "%Y/%m")
  n.Dates = length(unique(short.date))
  
  # aggregation at a monthly scale
  mo.agg = matrix(nrow=n.Dates,ncol=p)
  for(i.st in 1:p){
    R = mat[,i.st]
    R.agg = aggregate(R ~ short.date,FUN = fun.agg, na.action=na.pass)
    mo.agg[,i.st] <- R.agg$R
  }
  
  # month and year
  dates.agg = R.agg$short.date
  mo = as.numeric(substr(dates.agg,start=6,stop=7))
  
  # return TS-like list 
  return(list(mat.agg=mo.agg,month=mo))
}


#==============================================================================
# First letter of months
#==============================================================================
get.letter1.month = function(){
  # get.letter1.month (computeStat_lib): First letter of months
  return(c('J','F','M','A','M','J','J','A','S','O','N','D'))
}


#==============================================================================
# Season label as a function of month indices
#==============================================================================
get.season.label = function(vec.ind){
  # get.season.label (computeStat_lib): season label as a function of month indices
  all.mo = get.letter1.month()
  letter.mo = all.mo[unique(vec.ind)]
  return(paste0(letter.mo,collapse=""))
}


#==============================================================================
# counts lengths of sequences of a certain value in a series
#==============================================================================
rle.length <- function(series, val) {
  # rle.length (computeStat_lib): counts lengths of sequences of a certain value in a series
  rle(series)$lengths[rle(series)$val==val] 
}


#==============================================================================
# add na to a time series when the series is interrupted
#==============================================================================
add.na.series.interrupted = function(x, # vector of values
                                     vec.dates # vector of dates
){
  # add.na.series.interrupted (computeStat_lib): add na to a time series when the series is interrupted
  
  # filter consecutive dates
  cut.seq = diff(vec.dates,lag=1)>1
  
  # where to add NAs
  index.cut.seq = which(cut.seq)
  n.cut = length(index.cut.seq)
  
  # add nas
  if(n.cut>0){
    for(i in 1:n.cut) x = append(x,NA,index.cut.seq[i]+i-1)
  }
  
  return(x)
}

#==============================================================================
test.add.na.series.interrupted = function(){
  x = c(0,0,0,0,1,0,1,0,0,0,1,1,0,1,0,1,0,0,1,0,1,1,1,1,0)
  dates <- c("01/01/90", "02/01/90", "03/01/90", "04/01/90", "05/01/90",
             "06/01/90", "10/01/90", "11/01/90", "12/01/90", "13/01/90",
             "14/01/90", "15/01/90", "16/01/90", "17/01/90", "18/01/90",
             "01/02/90", "02/02/90", "03/02/90", "04/02/90", "05/02/90",
             "06/02/90", "07/02/90", "08/02/90", "11/02/90", "12/02/90")
  vec.dates = as.Date(dates, "%d/%m/%y")
  return(add.na.series.interrupted(x,vec.dates))
}


#==============================================================================
# counts lengths of dry and wet periods in a series
#==============================================================================
length.dry.wet.periods = function(x, # vector of values
                                  vec.dates, # vector of dates
                                  th # threshold
){
  # length.dry.wet.periods (computeStat_lib): counts lengths of dry and wet periods in a series
  
  # vector x with NAs where there is an interruption in the series
  x.interr = add.na.series.interrupted(x,vec.dates)
  
  # transform x in boolean indicating dry (=0) and wet states
  x.states = x.interr>th
  
  # lengths of dry and wet periods
  dry.lengths = rle.length(x.states,0)
  wet.lengths = rle.length(x.states,1)
  
  return(list(dry.lengths=dry.lengths,wet.lengths=wet.lengths))
}



#==============================================================================
# length of dry/wet spells in a series
#==============================================================================
get.length.spell = function(x,vec.dates,th,max.spell){
  # get.length.spell (computeStat_lib.r): return length of dry and wet spells
  #
  # INPUTS:
  # - x: vector of precipitation values
  # - vec.dates: vector of related dates
  # - th: threshold
  # - maxspell: maximum durations we are interested in: vector of length 2 (e.g. c(20,15))
  
  library(zoo)
  
  # vector x with NAs where there is an interruption in the series
  x.interr = add.na.series.interrupted(x,vec.dates)
  
  # transform x in boolean indicating dry (=0) and wet states
  x.states = x.interr>th
  
  # run length encoding
  rle.x = rle(x.states)
  
  # spell lengths
  len.dry = 1:max.spell[["dry"]]
  len.wet = 1:max.spell[["wet"]]
  
  # dry spells
  is.dry = rle.x$values==F
  dry.len = rle.x$lengths[is.dry]
  dry.freq = table(dry.len)
  len.obs = as.numeric(names(dry.freq))
  i.match = match(len.dry,len.obs)
  freq.obs = as.numeric(dry.freq)
  freq.obs = freq.obs/sum(freq.obs)
  vec.dry = freq.obs[i.match]
  
  # wet spells
  is.wet = rle.x$values==T
  wet.len = rle.x$lengths[is.wet]
  wet.freq = table(wet.len)
  len.obs = as.numeric(names(wet.freq))
  i.match = match(len.wet,len.obs)
  freq.obs = as.numeric(wet.freq)
  freq.obs = freq.obs/sum(freq.obs)
  vec.wet = freq.obs[i.match]
  
  return(list(vec.dry=vec.dry,vec.wet=vec.wet))
}


#==============================================================================
# number of consecutive wet days in a series
#==============================================================================
rollapply.wet.spell = function(x,vec.dates,nb.days,th){
  # rollapply.wet.spell (computeStat_lib.r): number of consecutive wet days in a series
  #
  # INPUTS:
  # - x: vector of precipitation values
  # - vec.dates: vector of related dates
  # - nb.days: number of consecutive wet days
  # get.length.spell (computeStat_lib.r)
  
  library(zoo)
  
  # vector x with NAs where there is an interruption in the series
  x.interr = add.na.series.interrupted(x,vec.dates)
  
  # transform x in boolean indicating dry (=0) and wet states
  x.states = x.interr>th
  
  # number of consecutive wet days in a series
  return(rollapply(x.states, width = nb.days, FUN = sum, align = "left"))
}


#==============================================================================
# cumulative amounts of consecutive wet days in a series
#==============================================================================
cumsum.wet.spell = function(x,vec.dates,nb.days,th){
  # cumsum.wet.spell (computeStat_lib): cumulative amounts of consecutive wet days in a series
  
  # number of consecutive wet days over a period of "nb.days" days 
  nb.wet.days=rollapply.wet.spell(x,vec.dates,nb.days)
  # starting index of the periods with "nb.days" wet days 
  ind.per.wet = which(nb.wet.days==nb.days)
  # difference between successive indices
  diff.ind.per.wet = diff(ind.per.wet)
  # number of periods
  nb.per = sum(diff.ind.per.wet>1)
  # initialize matrix of periods starts and ends
  index.per = matrix(nrow=nb.per,ncol=2)
  # loop over the periods
  cmt.per=0
  ind.deb.per = ind.per.wet[1]
  for(i in 1:length(diff.ind.per.wet)){
    # if the difference between the indices is greater than 1, that's the end of a period
    if(diff.ind.per.wet[i]>1){
      cmt.per = cmt.per+1 
      index.per[cmt.per,1]=ind.deb.per # begins period
      index.per[cmt.per,2]=ind.per.wet[i]+nb.days-1 # ends period
      ind.deb.per = ind.per.wet[i+1] # start of the next period
    }
  }
  # compute cumulated amounts
  int.cum = vector(length=nb.per)
  for(i in 1:nb.per) int.cum[i] = sum(x[index.per[i,1]:index.per[i,2]])
  
  return(int.cum)
}


#==============================================================================
# proportion of spatial entities (gauges, basins) which are wet
#==============================================================================
get.p1.prec = function(mat.prec,th){
  # get.p1.prec (computeStat_Lib): roportion of spatial entities (gauges, basins) which are wet
  
  # day which are wet
  is.wet = mat.prec>th
  # proportion of the stations which are wet
  p1.prec = apply(is.wet,1,mean)
  
  return(p1.prec)
}


#==============================================================================
# get yearly statistic of a vector of weather data
#==============================================================================
get.annual.yearly.Stat = function(vec.var,
                              vec.dates,
                              vec.y,
                              n.day.cumul,
                              funYear,
                              funCumul,
                              na.action=na.pass
){
  # get.annual.yearly.Stat (plotStat_lib): get yearly statistic of a vector of weather data
  
  if(n.day.cumul>1){
    cum = rollapply(zoo(vec.var,vec.dates),width=n.day.cumul,FUN=funCumul,align="right")
    vec.y.trim = vec.y[(n.day.cumul):length(vec.var)]
    x = aggregate(cum~vec.y.trim, FUN=funYear, na.action=na.action)$cum  
  }else{
    x = aggregate(vec.var~vec.y, FUN=funYear, na.action=na.action)$vec.var  
  }
  
  return(x)
}


#============================================================================================
# precipitation criteria
#============================================================================================
get.prec.Stats <- function(vec.prec, # precipitation series
                           vec.dates,  # dates
                           pdt = 24,
                           pct.na = 0.1 # pourcentage de valeurs manquantes autorise
){
  # get.prec.Stats (computeStat_lib.r): precipitation criteria
  
  # dates
  Y = as.numeric(format(vec.dates, "%Y"))
  vec.year = min(Y):max(Y)
  n.year = length(vec.year)
  
  ###################
  #   traitement    #
  ###################
  # max et total de pluies par an
  vecmax = vecsum = vecna = rep(NA,n.year)
  j=0
  nb.ans = 0
  for(y in vec.year){
    j=j+1
    is.y = Y==y
    # pluie de l'annee
    Py = vec.prec[is.y]
    # pluies non manquantes
    Pyz = Py[!is.na(Py)] 
    if(length(Pyz)>365.25*(24/pdt)*(1-pct.na)){
      # nombre annees pseudo-completes
      nb.ans = nb.ans+1
      vecsum[j]=mean(Pyz)*365.25*(24/pdt)
      vecmax[j]=max(Pyz)
      vecna[j]=mean(is.na(Py))
    }
  }
  
  # calcul Indice
  PX = mean(vecmax,na.rm=T) # moyenne maxima annuels
  PA = mean(vecsum,na.rm=T) # moyenne totaux annuels
  IC = PX/PA
  return(list(IC=IC,TOT=PA,MAX=PX,NBANS=nb.ans,VECSUM=vecsum,VECMAX=vecmax,VECNA=vecna))
}

#========================================================================
# FitEMP: ajuste distribution empirique
#========================================================================
FitEMP = function(P){
  e = .Machine$double.eps
  x = sort(P)
  n = length(x)
  if(n>0){
    Femp.func = ecdf(x)
    Femp = Femp.func(x)*n/(n+1)
  }else{
    Femp = 1
  }
  return(Femp)
}


#========================================================================
# renvoie quantile Gumbel -ln(-ln(x))
#========================================================================
Q_Gumbel = function(CDF,pow=1){
  # calculate quantile on a Gumbel scale (-log(-log(x)))
  e = .Machine$double.xmin
  MIN_U=-log(-log(e)) # minimum quantile allowed
  MAX_U=-log(-log(1-e)) # maximum quantile allowed
  U = array(dim=length(CDF)) # return
  
  # quantile on a Gumbel scale
  toolow = CDF<=0
  toohigh = CDF>=1
  U[toolow] = MIN_U
  U[toohigh] = MAX_U
  U[!(toohigh|toolow)] = -log(-pow*log(CDF[!(toohigh|toolow)]))
  
  return(U)
}



#========================================================================
# CDF Empirique pour 2 stations
#========================================================================
get.emp.cdf.pair = function(P1,P2,th){
  # station 1
  is.wet1 = P1>=th&!is.na(P1)
  
  # station 2
  is.wet2 = P2>=th&!is.na(P2)
  
  # trouve successions pluies
  is.joint = is.wet1&is.wet2
  P1.wet = P1[is.joint]
  P2.wet = P2[is.joint] 
  
  # cdf empiriques
  cdf.P1 = ecdf(P1.wet)(P1.wet)
  cdf.P2 = ecdf(P2.wet)(P2.wet)
  
  return(list(cdf.P1=cdf.P1,cdf.P2=cdf.P2))
}

#========================================================================
# CDF Empirique pour une station aux temps t et t+1
#========================================================================
get.emp.cdf.lag1 = function(P,th){
  # jours pluvieux
  is.wet = P>=th&!is.na(P)
  
  # trouve successions pluies
  n = length(P)
  is.joint = is.wet[1:(n-1)]&is.wet[2:n]
  Pt1 = P[1:(n-1)][is.joint]
  Pt2 = P[2:n][is.joint] 
  
  # cdf empiriques
  cdf.Pt1 = ecdf(Pt1)(Pt1)
  cdf.Pt2 = ecdf(Pt2)(Pt2)
  
  return(list(cdf.Pt1=cdf.Pt1,cdf.Pt2=cdf.Pt2))
}


#========================================================================
# CDF Copule Empirique pour une station aux temps t et t+1
#========================================================================
get.biv.cdf.lag1 = function(P,th){
  library(copula)
  
  # taille vecteur
  n = length(P)
  
  # jours pluvieux
  is.wet = P>=th&!is.na(P)
  is.joint = is.wet[1:(n-1)]&is.wet[2:n]
  Pt1 = P[1:(n-1)][is.joint]
  Pt2 = P[2:n][is.joint] 
  x = cbind(Pt1,Pt2)
  
  # copule empirique
  cdf.lag1 = F.n(x,x)
  
  # results
  cdf.out = rep(NA,(n-1))
  cdf.out[is.joint] = cdf.lag1
  
  return(cdf.out)
}


#========================================================================
# Copule Empirique
#========================================================================
get.emp.copula = function(u,v,max.cop.dens=3){
  library(fCopulae)
  empCop = dempiricalCopula(u, v, N=30)
  d = expand.grid(empCop$x,empCop$y)
  colnames(d) = c("x","y")
  pdf.emp = as.vector(empCop$z)*nrow(d)
  pdf.emp[pdf.emp<0]=0
  pdf.emp[pdf.emp>max.cop.dens]=max.cop.dens
  d$emp = pdf.emp
  
  return(d)
}


#========================================================================
# Copule Gumbel
#========================================================================
get.Gumbel.copula = function(u,v,x,y,max.cop.dens=3){
  # Gumbel copula
  library(copula)
  cop.obj = gumbelCopula()
  fit.cop <- fitCopula(cop.obj, cbind(u, v), method="itau")
  cdf.cop = dCopula(cbind(x,y), slot(fit.cop,"copula"))
  cdf.cop[cdf.cop<0]=0
  cdf.cop[cdf.cop>max.cop.dens]=max.cop.dens
  d$x = x
  d$y = y
  d$cop = cdf.cop
  
  return(d)
}


#========================================================================
# Krige on Switzerland
#========================================================================
krige.Swiss.stat = function(x,y,Stat,is.pos=T){
  source('../SVN/GWEX/maps_lib.r')
  library(gstat)
  proj='CH1903_LV03'
  
  # Grid
  library(sp)
  n.grid = 200
  coord.min = c(min(x),min(y))
  cellsize = round(c(diff(range(x))/n.grid,diff(range(y))/n.grid))
  grd.sp = GridTopology(coord.min,cellsize,c(n.grid,n.grid))
  grd = as.data.frame(coordinates(grd.sp))
  names(grd) = c("x","y")
  
  # variogram
  a <- data.frame(x=x/1000,y=y/1000,z=Stat)
  g <- gstat(id="X", formula = z~1, locations = ~x+y, data = a)
  v.fit <- fit.variogram(variogram(g), vgm(0,"Sph",100,0))
  v.fit[2,3] = v.fit[2,3]*1000
  
  # We want now to use kriging to predict H at each point on the grid
  a <- data.frame(x=x,y=y,z=Stat)
  krig.DF <- krige(id="Stat",z~1, locations=~x+y, model=v.fit, data=a, newdata=grd)
  if(is.pos) krig.DF$Stat.pred[krig.DF$Stat.pred<0]=0
  
  return(krig.DF)
}


#####################################################
# get.cor.matrix: robust estimation of a correlation
# matrix
#####################################################
get.cor.matrix = function(P){
  # get.cor.matrix (GWex_lib): robust estimation of a correlation
  # matrix
  #
  # INPUT: P (n x p) matrix of precipitation
  
  # pairwise estimation of a correlation matrix with the Kendall tau
  Sig.Raw = cor.Pearson.robust(P)
  
  # modify correlation matrix to obtain a positive definite matrix
  Sig = modify.cor.matrix(Sig.Raw)
  
  return(Sig)
}



#####################################################
# Compute pearson coefficient from the Kendall Tau
#####################################################
cor.Pearson.robust = function(MAT){
  # cor.Pearson.robust (GWex_lib.r): Compute pearson coefficient from the kendall tau
  # This relationship is valid for both Gaussian and Student copula dependences
  # which is not the case of the Spearman rho (see McNeil et al., 2005, p.235)
  #
  # INPUT: MAT, matrix of observations
  
  # kendall correlation coefficients
  cor.kendall = cor(MAT, method="kendall", use="pairwise.complete.obs")
  
  # filter correlations computed on a very low number of obs.
  MAT.nz = !is.na(MAT)
  nPair <- t(MAT.nz) %*% MAT.nz
  cor.kendall[nPair<5] = 0
  
  # relation between kendall and pearson correlations
  cor.pearson = sin(pi*cor.kendall/2)
  
  return(cor.pearson)
}


#####################################################
# Modify a non-positive definite correlation matrix
# in order to have a positive definite matrix
#####################################################
modify.cor.matrix = function(cor.matrix){
  # modify.cor.matrix (GWex_lib.r): Modify a non-positive definite correlation matrix
  # in order to have a positive definite matrix
  # Rousseeuw, P. J. and G. Molenberghs. 1993. Transformation of non positive semidefinite
  # correlation matrices. Communications in Statistics: Theory and Methods 22(4):965-984.
  # Rebonato, R., & Jackel, P. (2000). The most general methodology
  # to create a valid correlation matrix for risk management and
  # option pricing purposes. J. Risk, 2(2), 17-26.
  #
  # INPUT: cor.matrix: positive definite correlation matrix
  #
  # OUTPUT: modified correlation matrix
  
  # eigen values and vectors
  eigen.out = eigen(cor.matrix,symmetric=T)
  eig.val.diag = eigen.out$values
  eig.vec = eigen.out$vectors
  
  # is there negative eigen values?
  is.neg.eigen.val = eig.val.diag<0
  
  # if yes, replace this negative eigen values by small positive values
  if(any(is.neg.eigen.val)){
    eig.val.diag[is.neg.eigen.val] = 10^-10
    eig.val = diag(eig.val.diag)
    # recontruct correlation matrix
    cor.matrix.r = eig.vec %*% eig.val %*% solve(eig.vec)
    # return normalized correlation matrix
    cor.matrix.out = cor.matrix.r / (diag(cor.matrix.r) %*% t(diag(cor.matrix.r)))
  }else{
    # if there are no negative eigen values, return the original correlation matrix
    cor.matrix.out = cor.matrix 
  }
  
  # diagonal values to be equal to 1, it appears that there is no tolerance for some functions
  # (generation of Student variates) and these values differ from 1 by an order of 10^-16
  diag(cor.matrix.out)=1
  return(cor.matrix.out)
}


#####################################################
# get.df.Student: cette fonction estime le paramètre
# nu (degrees of freedom) de loi Student multivariée
#####################################################
get.df.Student = function(P,Sig,max.df=20){
  # get.df.Student (GWex_lib.r): estimate nu parameter (degrees of freedom) from a multivariate t distribution
  # when the correlation matrix Sig is given (McNeil et al. (2005) "Quantitative Risk Management")
  
  # transformation PIT for the margins
  U = get.emp.cdf.matrix(P)
  
  # Compute likelihood for every df value, from 1 to 20, for more than 20, we can approximate by a Gaussian dist.
  vec.lk = rep(NA,max.df)
  for(df in 1:max.df){
    t.data = apply(U, 2, qt, df = df)
    lk = dmvt(x=t.data,sigma=Sig,df=df,log=T) - apply(dt(x=t.data,df=df,log=T),1,sum)
    nz = !is.na(lk)
    vec.lk[df] = sum(lk[nz])
  }
  
  # We take the df value which maximizes the likelihood
  df.hat = which.max(vec.lk)
  
  # if the likelihood return finite values
  if(!is.infinite(vec.lk[df.hat])){
    return(df.hat)
  }else{
    return(max.df+10)
  }
}



#==============================================================================
# Fonction get.emp.cdf.matrix: Empirical cdf values for all columns of a matrix
#==============================================================================
get.emp.cdf.matrix = function(X,th=NULL){
  # get.emp.cdf.matrix (processData_lib): Empirical cdf values for all columns of a matrix
  #
  # X: matrix n x p
  # th (optional): threshold
  library(lmomco)
  
  # number of columns
  p = ncol(X)
  
  # prepare output
  Y = matrix(NA,nrow=nrow(X),ncol=p)
  
  # 
  for(i in 1:p){
    # all data
    X.i = X[,i]
    # filter
    nz = !is.na(X.i)
    # cdf values (Gringorten prob (optimized for gumbel dist))
    X.cdf = lmomco::pp(X.i[nz],a=0.44,sort=F)
    # gaussian variates
    Y[which(nz),i] = X.cdf
  }
  
  return(Y)
}


#========================================================================
# renvoie quantile Gumbel -ln(-ln(x))
#========================================================================
Q_Gumbel = function(CDF){
  return(-log(-log(CDF)))
}

#========================================================================
# renvoie quantile Gumbel -ln(-ln(x))
#========================================================================
qq.emp = function(x,p){
  # sort x values
  x.srt = sort(x)
  # associated probabilites with the Gringorten prob (optimized for gumbel dist)
  p.srt = lmomco::pp(x.srt, a=0.44)
  # Gumbel transformation
  u.srt = Q_Gumbel(p.srt)
  # Gumbel transformation of the expected prob
  u.out = Q_Gumbel(p)
  # return interpolation
  x.out = approx(x=u.srt,y=x.srt,xout=u.out,rule=2)$y
  
  return(x.out)
}