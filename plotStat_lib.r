###===============================###===============================###
### Guillaume Evin
### 12/09/2016, Grenoble
###  LTHE
### guillaume_evin@yahoo.fr
###
### The following functions plot standard statistics for Precipitation
### time series
###===============================###===============================###
library(zoo)


#==============================================================================
# add bar of scores to an existing plot
#==============================================================================
add.bar.score = function(vec.scores,mywidth=0.05){
  # plotStat_lib.r
  #
  # vec.score: vector indicating a score = 1,2,3 or 4 (see Benett et al., 2017)
  #
  # compute proportion score
  prop.score = vector(length=3)
  prop.score[1] = mean(vec.scores==1)
  prop.score[2] = mean(vec.scores%in%c(2,3))
  prop.score[3] = mean(vec.scores==4)
  
  # add bar
  par(new=TRUE,usr=c(-0.2,1.05,-0.05,1),mgp=c(3,0,0))
  barplot(matrix(data = prop.score,ncol = 1), horiz=T, col=c('#1b9e77','gold1','#DB1702'),width=mywidth, add=T,  xaxt='n')
  axis(side = 1,at = c(0,0.25,0.5,0.75,1), pos = 0.01,labels = c('0','25','50','75','100'), tcl=-0.1, cex.axis=0.7)
  # add scores as a percentage
  if(prop.score[1]>0.12) text(prop.score[1]/2+0.008,0.035,labels=paste0(round(prop.score[1]*100),'%'))
  if(prop.score[2]>0.12) text(prop.score[1]+prop.score[2]/2+0.008,0.035,labels=paste0(round(prop.score[2]*100),'%'))
  if(prop.score[3]>0.12) text(prop.score[1]+prop.score[2]+prop.score[3]/2+0.008,0.035,col='white',
                             labels=paste0(round(prop.score[3]*100),'%'))
}

#==============================================================================
compute.score = function(stat.sim,stat.obs){
  # plotStat_lib.r
  #
  
  # cases with missing metrics
  if(any(is.na(stat.sim))) return(0)
  if(is.na(stat.obs)) return(0)
  
  # quantiles, mean and sd
  q.05 = quantile(stat.sim,probs=0.05,na.rm=T)
  q.95 = quantile(stat.sim,probs=0.95,na.rm=T)
  sim.mean= mean(stat.sim,na.rm=T)
  sim.sd= sd(stat.sim,na.rm=T)
   
  # score
  if(stat.obs>=q.05&stat.obs<=q.95){
    score = 1
  }else if(stat.obs>=(sim.mean-3*sim.sd)&stat.obs<=(sim.mean+3*sim.sd)){
    score = 2
  }else if(abs(stat.obs-sim.mean)/stat.obs<=0.05){
    score = 3
  }else{
    score = 4
  }
  
  return(score)
}



#==============================================================================
# add confidence intervals with polygons
#==============================================================================
add.CI.polygon = function(
  X, # matrix of data n x p: each line is a repetition (simulation)
  plot.x = NULL, # values for the x axis against which the polygons are added
  CI.vec = c(0.5,0.9,0.99), # levels of the CI
  col.min = 0.6, # level of the darkest gray
  col.max = 0.95 # level of the lightest gray
){
  # add.CI.polygon (plotStat_lib): add confidence intervals (CI) on a plot using polygons
  # with different shades of gray
  
  # number of CI
  CI.n = length(CI.vec)
  
  # length of the simulations
  p = ncol(X)
  
  # compute quantiles
  CI.X.inf = CI.X.sup = matrix(nrow=p,ncol=CI.n)
  for(i.CI in 1:CI.n){
    CI.X.inf[,i.CI] = apply(X,2,quantile,probs=0.5-CI.vec[i.CI]/2)
    CI.X.sup[,i.CI] = apply(X,2,quantile,probs=0.5+CI.vec[i.CI]/2)
  }
  
  # col
  col.vec = seq(from=col.min,to=col.max,length.out = CI.n)
  
  # values for the x axis against which the polygons are added
  if(is.null(plot.x)) plot.x = 1:p
  
  # plot
  for(i.CI in CI.n:1){
    col.CI = col.vec[i.CI]
    polygon(c(plot.x,rev(plot.x)),c(CI.X.inf[,i.CI],rev(CI.X.sup[,i.CI])), 
            border = NA, col=rgb(col.CI,col.CI,col.CI))
  }
}


#==============================================================================
# get label of standard Precipitation statistics
#==============================================================================
get.daily.label.stat = function(type){
  # get.daily.label.stat (plotStat_lib): get label of standard Precipitation statistics
  switch(type,
         DDF = expression(p[0]),
         WDF = expression(p[1]),
         MWDF = 'Mean wdi [mm/d]',
         SDWDF = 'Std Dev. wdi [mm/d]',
         WWTP = expression(p[11]),
         DWTP = expression(p['01']))
}


#==============================================================================
# plot standard daily statistics at a monthly scale
#==============================================================================
plot.daily.stats.by.month = function(obs, # matrix of observations nTime x nStation
                                     sim, # array of simulations nTime x nStation x nScenarios
                                     obs.date,
                                     sim.date,
                                     type,
                                     obs.th=NULL,
                                     sim.th=NULL,
                                     lab=NULL,
                                     n.fig.page=NULL,
                                     ylab = NULL,
                                     xlab = 'Month',
                                     main=NULL,
                                     fillPage=T,
                                     col.box='black',
                                     ...
){
  # plot.daily.stats.by.month (plotStat_lib): plot standard daily statistics at a monthly scale
  
  #### nombre de stations
  p = ncol(obs)
  
  #### ylab
  if(is.null(ylab)) ylab = get.daily.label.stat(type)
  
  #stats observes
  stat.obs = matrix(nrow=p,ncol=12)
  for(m in 1:12){
    is.m = as.numeric(format(obs.date,'%m'))==m
    stat.obs[,m] = get.daily.stat(mat=obs[is.m,,drop=FALSE], vec.dates=obs.date[is.m], type=type, th=obs.th)
  }
  
  #stats simulees
  n.SIM = dim(sim)[3]
  stat.sim = array(dim=c(p,12,n.SIM))
  for(m in 1:12){
    # filter this month
    is.m = as.numeric(format(sim.date,'%m'))==m
    for(i.s in 1:n.SIM){
      stat.sim[,m,i.s] = get.daily.stat(mat=sim[is.m,,i.s,drop=FALSE], vec.dates=sim.date[is.m], type=type, th=sim.th)
    }
  }
  
  # figures
  my.lim = range(c(stat.sim,stat.obs))
  
  
  for(i.st in 1:p){
    # sim
    sim.i = t(stat.sim[i.st,,])
    
    # obs
    obs.i = stat.obs[i.st,]
    
    # plot
    plot(obs.i,ylim=my.lim,xaxt='n',xlab="",ylab="",...)
    box(col=col.box,lwd=1.3)
    add.CI.polygon(sim.i)
    lines(obs.i,col='black',type="b",pch=20,cex=1.5)
    
    # add x-axis
    axis(1,at = 1:12,labels=get.letter1.month())
    
    # station name
    if(!is.null(lab)) legend('topleft',legend = lab[i.st],bty='n',cex = 1.3)
    # ajoute titres generaux
    if(!is.null(n.fig.page)){
      if(i.st %% n.fig.page == 0){
        mtext(text=xlab,side=1,line=1,cex=1.5,outer=TRUE)
        mtext(text=ylab,side=2,line=2,cex=1.5,outer=TRUE)
        if(!is.null(main)){mtext(text=main,side=3,line=1,cex=2.5,outer=TRUE)}
      }
    }
  }
  
  mtext(text=xlab,side=1,line=1,cex=1.5,outer=TRUE)
  mtext(text=ylab,side=2,line=2,cex=1.5,outer=TRUE)
  if(!is.null(main)){mtext(text=main,side=3,line=1,cex=2.5,outer=TRUE)}
  if(fillPage){
    while(!par('page')) plot.new()
  }
}

#==============================================================================
# plot monthly statistics
#==============================================================================
plot.monthly.stats.by.month = function(obs,
                                       obs.date,
                                       sim,
                                       sim.date,
                                       type, # any standard statistic (ex 'mean' or 'sd')
                                       type.agg = sum,
                                       lab=NULL,
                                       ylab = NULL,
                                       xlab = 'Month',
                                       main=NULL,
                                       n.fig.page=NULL,
                                       fillPage=T,
                                       col.box='black',
                                       ...){
  # plot.monthly.stats.by.month (plotStat_lib): plot monthly statistics
  
  #### number of spatial entities
  p = ncol(obs)
  
  #### ylab
  # if(is.null(ylab)) ylab = paste0('Monthly ',type)
  
  #### observed statistics
  stat.obs = matrix(nrow=p,ncol=12)
  TS.mo = get.monthly.agg(obs,obs.date,fun.agg=type.agg)
  for(m in 1:12){
    is.m = TS.mo$month==m
    R.mo = TS.mo$mat.agg[is.m,]
    stat.obs[,m] = apply(matrix(R.mo,nr=sum(is.m),nc=p),2,type)
  }
  
  #### simulated statistics
  n.rep = dim(sim)[3]
  stat.sim = array(dim=c(n.rep,p,ncol=12))
  for(i in 1:n.rep){
    # for this repetition
    TS.mo = get.monthly.agg(matrix(sim[,,i],nr=nrow(sim),nc=p),sim.date,fun.agg=type.agg)
    
    for(m in 1:12){
      is.m = TS.mo$month==m
      R.mo = TS.mo$mat.agg[is.m,]
      stat.sim[i,,m] = apply(matrix(R.mo,nr=sum(is.m),nc=p),2,type)
    }
  }
  
  #### figures
  my.lim = range(c(stat.sim,stat.obs),na.rm=T)
  
  # for each spatial entity
  for(i.st in 1:p){
    # sim
    sim.i = stat.sim[,i.st,]
    
    # obs
    obs.i = stat.obs[i.st,]
    
    # plot
    plot(obs.i,ylim=my.lim,xaxt='n',xlab="",ylab="")
    box(col=col.box,lwd=1.3)
    add.CI.polygon(sim.i)
    lines(obs.i,col='black',type="b",pch=20,cex=1.5)
    
    # add x-axis
    axis(1,at = 1:12,labels=get.letter1.month())
    
    # station name
    if(!is.null(lab)) legend('topleft',legend = lab[i.st],bty='n',cex = 1.3)
    # ajoute titres generaux
    if(!is.null(n.fig.page)){
      if(i.st %% n.fig.page == 0){
        mtext(text=xlab,side=1,line=1,cex=1.5,outer=TRUE)
        mtext(text=ylab,side=2,line=2,cex=1.5,outer=TRUE)
        if(!is.null(main)){mtext(text=main,side=3,line=1,cex=2.5,outer=TRUE)}
      }
    }
  }
  
  mtext(text=xlab,side=1,line=1,cex=1.5,outer=TRUE)
  mtext(text=ylab,side=2,line=2,cex=1.5,outer=TRUE)
  if(!is.null(main)){mtext(text=main,side=3,line=1,cex=2.5,outer=TRUE)}
  if(fillPage){
    while(!par('page')) plot.new()
  }
}


#==============================================================================
# plot CDFs of the observed and simulated dry or wet spell length frequencies
#==============================================================================
plot.cdf.spell.length = function(obs, # matrix of observations nTime x nStation
                                 sim, # array of simulations nTime x nStation x nScenarios
                                 obs.date,
                                 sim.date,
                                 type, # 'dry' or 'wet'
                                 obs.th,
                                 sim.th,
                                 ind.season=NULL,
                                 lab=NULL,
                                 ylab = NULL,
                                 xlab = NULL,
                                 main=NULL,
                                 n.fig.page=NULL,
                                 fillPage=T,
                                 add.legend=T,
                                 col.box='black',
                                 ...){
  #### number of spatial entities (gauges, stations)
  p = ncol(obs)
  
  #### nombre de simulations
  n.rep = dim(sim)[3]
  
  #### max nb of days plotted by default
  max.spell = list(dry=30,wet=20)
  
  ##### xlab
  if(is.null(xlab)){
    xlab = ifelse(type=='dry','Dry spell length [days]','Wet spell length [days]')
  }
  
  #### legend
  col.poly = rgb(0.9,0.9,0.9)
  
  for(i.st in 1:p){
    spell.len.x = 1:max.spell[[type]]
    
    #### spell length obs
    # obs prec for this station and dates
    R = obs[,i.st]
    D = obs.date
    # filter if for a season only
    if(!is.null(ind.season)){
      dates.mo = as.numeric(format(D, "%m"))
      is.season = dates.mo%in%ind.season
      R = R[is.season]
      D = D[is.season]
    }
    # length of dry/wet spells
    DW.len = length.dry.wet.periods(R,D,obs.th)
    if(type=="dry"){
      len=DW.len$dry.lengths
    }else{
      len=DW.len$wet.lengths 
    }
    len.ecdf = ecdf(len)
    spell.len.obs = 1-len.ecdf(spell.len.x)
    spell.len.obs[spell.len.obs==0] = 10^-9
    
    #### spell length sim
    spell.len.sim = matrix(nrow=n.rep,ncol=length(spell.len.x))
    for(i.sim in 1:n.rep){
      # obs prec for this station and dates
      R = sim[,i.st,i.sim]
      D = sim.date
      # filter if for a season only
      if(!is.null(ind.season)){
        dates.mo = as.numeric(format(D, "%m"))
        is.season = dates.mo%in%ind.season
        R = R[is.season]
        D = D[is.season]
      }
      # length of dry/wet spells
      DW.len = length.dry.wet.periods(R,D,sim.th)
      if(type=="dry"){
        len=DW.len$dry.lengths
      }else{
        len=DW.len$wet.lengths
      }
      len.ecdf = ecdf(len)
      spell.len.sim[i.sim,] = len.ecdf(spell.len.x)
    }
    # CI 95%
    qInf = apply(1-spell.len.sim,2,quantile,0.005)
    qInf[qInf==0] = 10^-9
    QSup = apply(1-spell.len.sim,2,quantile,0.995)
    QSup[QSup==0] = 10^-9
    
    #### figure
    plot(spell.len.x,spell.len.obs,type="b",pch=20,cex=1.5,ylim=c(0.001,1),log="xy",xlab='',ylab='',...)
    box(col=col.box,lwd=1.3)
    polygon(c(spell.len.x,rev(spell.len.x)),
            c(qInf,rev(QSup)),col=col.poly)
    lines(spell.len.x,spell.len.obs,type="b",pch=20,cex=1.5)
    
    # legend
    if(add.legend){
      legend('topright',legend=c('Obs.','Sim.'),bty='n',
             col = c("black",col.poly),lty = c(1,0), lwd = c(1.5, 0),
             pch = c(20,22),pt.bg = c(NA,col.poly),pt.cex = 2, cex=1.3) 
    }
    
    # main label
    if(!is.null(lab)) legend('bottomleft',legend = lab[i.st], cex=1.3, bty='n') 
    
    # ajoute titres generaux
    if(!is.null(n.fig.page)){
      if(i.st%%n.fig.page==0){
        mtext(text=xlab,side=1,line=1,cex=1.5,outer=TRUE)
        mtext(text=ylab,side=2,line=2,cex=1.5,outer=TRUE)
        if(!is.null(main)){mtext(text=main,side=3,line=1,cex=2.5,outer=TRUE)}
      }
    }
  }
  
  mtext(text=xlab,side=1,line=1,cex=1.5,outer=TRUE)
  mtext(text=ylab,side=2,line=2,cex=1.5,outer=TRUE)
  if(!is.null(main)){mtext(text=main,side=3,line=1,cex=2.5,outer=TRUE)}
  if(fillPage){
    while(!par('page')) plot.new()
  }
}



#==============================================================================
# plot CDFs of the observed and simulated dry or wet spell length frequencies
#==============================================================================
plot.freq.spell.length = function(OBS,
                                  SIM,
                                  th,
                                  lab=NULL,
                                  xlab1='Dry spell length [days]',
                                  xlab2='Wet spell length [days]',
                                  ylab='Frequency',
                                  main=NULL,
                                  n.fig.page=NULL,
                                  fillPage=T){
  # plot.freq.spell.length (plotStat_lib): plot CDFs of the observed and simulated 
  # dry or wet spell length frequencies
  #### nombre de stations
  
  p = ncol(OBS$obs)
  
  ### nombre de simulations
  n.sim = dim(SIM$sim)[1]
  
  ### durees etudiees
  max.spell = list(dry=30,wet=20)
  len.dry = 1:max.spell[["dry"]]
  len.wet = 1:max.spell[["wet"]]
  
  for(i.st in 1:p){
    #### spell length obs
    # obs prec for this station and dates
    R = OBS$obs[,i.st]
    D = OBS$dates
    
    #### length of dry/wet spells
    DW.len.obs = get.length.spell(R,D,th,max.spell)
    vec.dry.obs = DW.len.obs$vec.dry
    vec.wet.obs = DW.len.obs$vec.wet
    
    #### spell length sim
    mat.dry.sim = matrix(nrow=max.spell[["dry"]],ncol=n.sim)
    mat.wet.sim = matrix(nrow=max.spell[["wet"]],ncol=n.sim)
    for(i.sim in 1:n.sim){
      # obs prec for this station and dates
      R = SIM$sim[i.sim,,i.st]
      # length of dry/wet spells
      DW.len.sim = get.length.spell(R,D,0,max.spell)
      mat.dry.sim[,i.sim] = DW.len.sim$vec.dry
      mat.wet.sim[,i.sim] = DW.len.sim$vec.wet 
    }
    # mean
    mean.dry.sim = apply(mat.dry.sim,1,median)
    mean.wet.sim = apply(mat.wet.sim,1,median)
    
    #### figure dry
    barplot(rbind(vec.dry.obs,mean.dry.sim), beside=T, names.arg=len.dry,
            xlab=xlab1,ylab=ylab,cex.lab=1.3)
    legend('topright',legend=c('Obs.','Sim.'),fill=c('black','gray'), bty='n',cex=1.3)
    # main label
    if(!is.null(lab)) title(lab[i.st],cex=1.3) 
    
    #### figure wet
    barplot(rbind(vec.wet.obs,mean.wet.sim), beside=T, names.arg=len.wet,
            xlab=xlab2,ylab=ylab,cex.lab=1.3)
    legend('topright',legend=c('Obs.','Sim.'),fill=c('black','gray'), bty='n',cex=1.3)
    # main label
    if(!is.null(lab)) title(lab[i.st],cex=1.3) 
  }
}





#==============================================================================
# plot scatterplot of observed vs simulated max, for a fixed period
#==============================================================================
plot.scatter.fun = function(obs,
                            sim,
                            obs.date,
                            sim.date,
                            n.day.cumul=1,
                            lab=NULL,
                            n.fig.page=NULL,
                            funYear=max,
                            funCumul=sum,
                            xlab='Return period (-ln(-ln(F)))',
                            ylab='Maxima [mm]',
                            main=NULL,
                            fillPage=T,
                            col.box='black',
                            na.action=na.action,
                            ...)
{
  library(zoo)
  library(lmomco)
  
  #### nombre de stations
  p = ncol(obs)
  
  #### years observed
  y.obs = strftime(obs.date, "%Y")
  
  #### years simulated
  n.SIM = dim(sim)[3]
  y.sim = strftime(sim.date, "%Y")
  
  for(i.st in 1:p){
    
    #### max observed
    # get max obs
    x.obs = get.annual.yearly.Stat(obs[,i.st],obs.date,y.obs,n.day.cumul,funYear,funCumul,na.action=na.action)
    nz = !is.na(x.obs)
    x.obs.sort = sort(x.obs[nz]) 
    U.obs = -log(-log(lmomco::pp(x.obs.sort, a=0.44)))
    
    #### max simulated
    # get max sim
    mat.x.sim = matrix(nrow=n.SIM,ncol=length(unique(y.sim)))
    for(i.sim in 1:n.SIM){
      mat.x.sim[i.sim,] = get.annual.yearly.Stat(vec.var = sim[,i.st,i.sim],vec.dates = sim.date,
                                                 vec.y = y.sim,n.day.cumul = n.day.cumul,
                                                 funYear,funCumul,na.action=na.action)
    }
    x.sim.sort = t(apply(mat.x.sim,1,sort))
    U.sim = -log(-log(lmomco::pp(x.sim.sort[1,], a=0.44)))
    
    #### plot
    myxlim = range(c(U.obs,U.sim))
    myylim=range(c(x.obs.sort,x.sim.sort))
    plot(-100,-100,xlim=myxlim,ylim=myylim,xlab='',ylab='',...)
    box(col=col.box,lwd=1.3)
    grid()
    add.CI.polygon(X=x.sim.sort,plot.x=U.sim)
    lines(U.obs,x.obs.sort,lwd=2)
    
    # station name
    if(!is.null(lab)) legend('topleft',legend=lab[i.st],bty='n',cex=1.3)
    
    # ajoute titres generaux
    if(!is.null(n.fig.page)){
      if(i.st %% n.fig.page == 0){
        mtext(text=xlab,side=1,line=1,cex=1.5,outer=TRUE)
        mtext(text=ylab,side=2,line=2,cex=1.5,outer=TRUE)
        if(!is.null(main)){mtext(text=main,side=3,line=1,cex=2.5,outer=TRUE)}
      }
    }
  }
  
  mtext(text=xlab,side=1,line=1,cex=1.5,outer=TRUE)
  mtext(text=ylab,side=2,line=2,cex=1.5,outer=TRUE)
  if(!is.null(main)){mtext(text=main,side=3,line=1,cex=2.5,outer=TRUE)}
  if(fillPage){
    while(!par('page')) plot.new()
  }
}



#==============================================================================
# plot.NHRH.fit represents the fitted NHRH distribution on top of the obervations
# on an histogram and a cdf, for all the stations
plot.NHRH.fit = function(P.mat,fit.out,th=0.2){
  # plotStat_lib
  #
  # INPUTS:
  # - P.mat: matrix of observations
  # - fit.out: output from fit.GWex
  # - th: threshold
  
  # nombre de stations
  n.Stat = fit.out$n.Stat
  
  # Dates
  vec.Dates = as.Date(rownames(P.mat))
  
  # noms des stations
  lab = colnames(P.mat)
  
  # For each period
  for(per in fit.out$periods){
    # filter
    is.per = get.vec.dates.filter(vec.Dates,per)
    
    # pour chaque station
    for(i.st in 1:n.Stat){
      # pluies de la station
      P = P.mat[is.per,i.st]
      nz = is.prec(P,th)
      all.P = P[nz]
      
      # parametres de la loi LGPD
      par.fit = fit.out$NHRH.par[[per]][i.st,]
      plot.dist.fit(all.P,par.fit,type="LGPD",th=th,mylab="Precipitation [mm]")
      mtext(text = paste0(per,", ",lab[i.st]), side=3, line=-2, outer=T, cex = 1.5)
    }
  }
}


#==============================================================================
# plot.dist.fit represents the fitted distribution on top of the obervations
# on an histogram and a cdf
plot.dist.fit = function(obs,par.fit,type,th,mylab="x",mycex=1.3){
  # plotStat_lib
  #
  # INPUTS:
  # - obs: vector of observations
  # - par.fit: vector of parameters
  # - type: name of the distribution
  # - th: threshold
  
  # prepare figure
  par(mfrow=c(1,2))
  
  #____________ ajustement avec la densite
  # prepare figure
  max.obs = quantile(obs,probs = 0.95)
  x.obs = seq(from=th,to=max.obs,by=0.1)
  # density of the fitted distribution
  if(type=="LGPD"){
    dfit = dNHRH.GI(x.obs,par.fit[1],par.fit[2],par.fit[3])
  }else if(type=="gamma"){
    dfit = dgamma(x.obs,par.fit[1],par.fit[2])
  }else if(type=="kde"){
    dfit = dkde(x.obs,par.fit)
  }
  # histogramme + fit
  hist(obs[obs<max.obs],breaks=seq(from=th,to=max.obs+2,length.out=30),cex.lab=mycex,
       cex.axis=mycex,freq=F,xlim=c(th,max.obs),main="",xlab=mylab,ylab="Density f(x)")
  lines(x.obs,dfit,lty=1,lwd=2,col="red")
  
  #____________ ajustement avec le cdf
  # prepare figure
  max.obs = max(obs,na.rm=T)
  plot(-10,-10,xlim=c(-2,10),ylim=c(0,max.obs),cex.lab=mycex,cex.axis=mycex,
       xlab="U=-ln{-ln[F(x)]}",ylab=mylab)
  grid()
  # distribution empirique
  Femp = FitEMP(obs)
  U_emp = Q_Gumbel(Femp)
  points(U_emp,sort(obs),col='black',pch=20,cex=1.5)
  # cdf of the fitted distribution
  x.obs = seq(from=th,to=max.obs,by=1)
  if(type=="LGPD"){
    pfit = pNHRH.GI(x.obs,par.fit[1],par.fit[2],par.fit[3])
  }else if(type=="gamma"){
    pfit = pgamma(x.obs,par.fit[1],par.fit[2])
  }else if(type=="kde"){
    pfit = pkde(x.obs,par.fit)
  }
  U_fit = Q_Gumbel(pfit)
  lines(U_fit,x.obs,col='red',lwd=2)
}


#==============================================================================
plot.copula.amount.PairOfSites = function(OBS,
                                          SIM,
                                          th,
                                          lab){
  #### nombre de stations
  n.Stat = ncol(OBS$obs)
  list.comb = combn(n.Stat,2)
  
  ####  parametres figure
  ncut = 50
  max.cop.dens = 3
  my.at = c(seq(from=0,to=max.cop.dens,length.out=ncut),max.cop.dens+1)
  col.l <- colorRampPalette(c('blue', 'cyan', 'green', 'yellow', 'orange', 
                              'red','purple'))
  
  for(i.comb in 1:ncol(list.comb)){
    i1 = list.comb[1,i.comb]
    i2 = list.comb[2,i.comb]
    
    # label
    xlab = paste0('Precipitation Day at station ',lab[i1],' [mm]') 
    ylab = paste0('Precipitation Day at station ',lab[i2],' [mm]')
    
    # observed
    cdf.obs = get.emp.cdf.pair(OBS$obs[,i1],OBS$obs[,i2])
    cop.obs = get.emp.copula(cdf.obs$cdf.P1, cdf.obs$cdf.P2,max.cop.dens)
    
    # simulated
    cdf.sim = get.emp.cdf.pair(SIM$sim[1,,i1],SIM$sim[1,,i2],th)
    cop.sim = get.emp.copula(cdf.sim$cdf.P1, cdf.sim$cdf.P2,max.cop.dens)
    
    # figure
    library(latticeExtra)
    library(gridExtra)
    
    plot1 <- levelplot(emp ~ x * y, cop.obs, cuts=ncut, col.regions=col.l,
                       xlab=xlab, ylab=ylab, at=my.at, main="Observed",
                       panel=panel.2dsmoother, args=list(span=0.3))
    
    plot2 <- levelplot(emp ~ x * y, cop.sim, cuts=ncut, col.regions=col.l,
                       xlab=xlab, ylab=ylab, at=my.at, main="Simulated",
                       panel=panel.2dsmoother, args=list(span=0.3))
    
    grid.arrange(plot1,plot2, ncol=2)
  }
}


#==============================================================================
plot.copula.amount.lag1 = function(OBS,
                                   SIM,
                                   th,
                                   lab){
  #### nombre de stations
  n.Stat = ncol(OBS$obs)
  list.comb = combn(n.Stat,2)
  
  ####  parametres figure
  ncut = 50
  max.cop.dens = 2
  my.at = c(seq(from=0,to=max.cop.dens,length.out=ncut),max.cop.dens+0.1)
  col.l <- colorRampPalette(c('blue', 'cyan', 'green', 'yellow', 'orange', 
                              'red','purple'))
  
  for(i.st in 1:n.Stat){
    # label
    xlab = 'Precipitation Day t [mm]'
    ylab = 'Precipitation Day t+1 [mm]'
    
    # observed
    cdf.obs = get.emp.cdf.lag1(OBS$obs[,i.st],th)
    cop.obs = get.emp.copula(cdf.obs$cdf.Pt1, cdf.obs$cdf.Pt2,max.cop.dens)
    
    # simulated
    cdf.sim = get.emp.cdf.lag1(SIM$sim[1,,i.st],th)
    cop.sim = get.emp.copula(cdf.sim$cdf.Pt1, cdf.sim$cdf.Pt2,max.cop.dens)
    
    # figure
    library(latticeExtra)
    library(gridExtra)
    
    plot1 <- levelplot(emp ~ x * y, cop.obs, cuts=ncut, col.regions=col.l,
                       xlab=xlab, ylab=ylab, at=my.at, main="Observed",
                       panel=panel.2dsmoother, args=list(span=0.3))
    
    plot2 <- levelplot(emp ~ x * y, cop.sim, cuts=ncut, col.regions=col.l,
                       xlab=xlab, ylab=ylab, at=my.at, main="Simulated",
                       panel=panel.2dsmoother, args=list(span=0.3))
    
    grid.arrange(plot1,plot2, ncol=2)
  }
}


#==============================================================================
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
pairs.panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs.add.abline <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1,...)
  grid()
}

my.pairs.df = function(df){
  pairs(df, lower.panel = pairs.add.abline, upper.panel = pairs.panel.cor)
}


#====================================
plot.density.cop =function(u, v, n=20, plot.points=T, zlim=NULL){
  # bivariate density estimation
  library(MASS)
  bivn.kde <- kde2d(u, v, n = n)
  
  # plot density and contour lines
  par(mar=c(3,3,4,1))
  if(is.null(zlim)){
    image(bivn.kde,col=rev(heat.colors(101)))  
  }else{
    image(bivn.kde,col=rev(heat.colors(101)),zlim=zlim)  
  }
  contour(bivn.kde,col="black",add=TRUE)
  if(plot.points) points(cbind(u,v),cex=.2)
}



#==============================================================================
pdf.plot.fun = function(filename,plot.fun,args){
  ### pdf.plot.scatter.fun: observed versus simulated statistics
  
  # number of rows
  ncol.fig = 4
  nrow.fig = ncol.fig*1.25
  n.fig.page = nrow.fig*ncol.fig
  args[['n.fig.page']] = n.fig.page
  
  # start pdf
  pdf(file=filename, width = 8.27, height = 11.69)
  layout(matrix(c(1:n.fig.page),nrow.fig,ncol.fig))
  par(oma=c(2.5,2,4,0),mar=c(2.5,1.5,0.1,0.1),mgp=c(2,0.5,0),mfrow=c(nrow.fig,ncol.fig))
  do.call(plot.fun,args)
  dev.off()
}

#==============================================================================
pdf.systematic.diagnostic = function(f.main, # file that will contain the daignostic (pdf file)
                                     varType, # type of variable: 'Prec' or 'Temp'
                                     obs.mat, # matrix of observations
                                     obs.date, # vector of date for the observations
                                     obs.th, # threshold of the dry/wet days for obs
                                     sim.array, # array of simulated scenarios
                                     sim.date, # vector of date for the scenarios
                                     sim.th, # threshold of the dry/wet days for sim
                                     MetaData, # Metadata for the stations
                                     modelType='GWEX', # Type of model: 'GWEX' or 'SCAMP'
                                     my.cex.axis=0.9, # expansion for axis ticks
                                     eval.6H = F, # evaluation at a 6-hour scale
                                     f.arealEstimates = NULL, # RData file containing the areal estimates
                                     f.6h.Estimates = NULL, # RData file containing the areal estimates at a 6h-scale
                                     StationPrec_SEL = c('VAR','COY','MUR','LTB','ANT','GLA'), # 6 representative precipitation stations
                                     Basin15_SEL = c('AabLen','XNEXSE','LinWee','XTHXPP','YEMYYY','XVWXNE'), # 6 representative basins
                                     StationTemp_SEL = c('NEU','BAS','SMA','JUN','ANT','GLA')  # 6 representative temperature stations
                                     
){
  #### SNV folder
  dir.SVN = Sys.getenv("SVN_EXAR")
  
  #### list of spatial aggregation scale
  list.sal = c('BASIN1','BASIN5','BASIN15')
  
  ##================= AREAL ESTIMATES ==================###
  if(!is.null(f.arealEstimates)){
    load(f.arealEstimates)
  }else{
    f.arealEstimates = paste0(f.main,'arealEstimates.RData')
    ####  mean elevation of the catchments (basin.mean.elevation)
    load(paste0(dir.SVN,'DATA/GEODATA/CATCHMENTS/Catchments_EXAR_mean_elevation.RData'))
    
    ## observed areal estimates
    list.arealObs = list()
    for(sal in list.sal){
      # number of basins this dissection
      load(paste0(dir.SVN,'DATA/OBS/DAILY/',sal,'/',varType,'/ArealObs_',varType,'.RData'))
      # selection period
      vec.Dates.arealobs = as.Date(rownames(arealObs))
      # areal estimates
      list.arealObs[[sal]] = arealObs[vec.Dates.arealobs%in%obs.date,]
    }
    
    if (modelType == 'GWEX') {
      ### simulated areal estimates
      list.arealSim  = list()
      for(sal in list.sal){
        # Thiessen weights
        TW = get.TW.stations(sal,MetaData$xCoord,MetaData$yCoord)
        # areal estimates
        list.arealSim[[sal]] = get.areal.sim(var.Gauge=sim.array,thiessen.weights=as.matrix(TW),varType=varType,
                                             station.Elevation=MetaData$Elevation,basin.Elevation=basin.mean.elevation[[sal]])
      }
    } else if (modelType == 'SCAMP') {
      load(paste0(f.main,'SCAMP_',varType,'.RData'))
    }
    
    ###  save areal estimates
    save(list.arealObs,list.arealSim,file=f.arealEstimates)
  }
  
  ##===== OBSERVED PSEUDO-OBS AND AREAL ESTIMATES AT 6-HOURS ====###
  if(eval.6H){
    if(!is.null(f.6h.Estimates)){
      load(f.6h.Estimates)
    }else{
      f.6h.Estimates = paste0(f.main,'6hEstimates.RData')
      ####  mean elevation of the catchments (basin.mean.elevation)
      load(paste0(dir.SVN,'DATA/GEODATA/CATCHMENTS/Catchments_EXAR_mean_elevation.RData'))
      
      ### Observed pseudo-observations for the recent periods (no meteorological analog): mat.DISAG
      load(paste0(dir.SVN,'DATA/OBS/PSEUDO_HOURLY/GAUGES/',varType,'/DisagObsGauge_Recent.RData'))
      if(any(colnames(mat.DISAG)!=colnames(obs.mat))) stop('mismatch station names in 6h disag')
      # aggregate
      obs.mat.6h = agg.matrix(mat.DISAG,6,average=(varType=='Temp'))
      obs.date.6h = as.Date(rownames(obs.mat.6h))
      
      ### observed areal estimates
      list.arealObs.6h = list()
      for(sal in list.sal){
        # areal obs. for this dissection: arealObs
        load(paste0(dir.SVN,'DATA/OBS/PSEUDO_HOURLY/',sal,'/',varType,'/ArealObs_',varType,'.RData'))
        # aggregate
        list.arealObs.6h[[sal]] = agg.matrix(arealObs,6,average=(varType=='Temp'))
      }
      
      ### disaggregate simulations
      # matrices for the disaggregation: list.matDisag
      load(paste0(dir.SVN,'DATA/OBS/PSEUDO_HOURLY/GAUGES/',varType,'/Mat4Disag.RData'))
      # disaggregate each scenario and aggregate to 6h
      nSim6h = dim(sim.array)[1]*4
      nStat = dim(sim.array)[2]
      nScen = dim(sim.array)[3]
      sim.array.6h = array(dim=c(nSim6h,nStat,nScen))
      # for the progress bar
      pb <- txtProgressBar()
      for(iScen in 1:nScen){
        sim.array.1h = getDisagObs1H(varType, list.matDisag$obs,  sim.array[,,iScen], list.matDisag$pi)$disagMAT
        sim.array.6h[,,iScen] = agg.matrix(sim.array.1h,6,average=(varType=='Temp'))
        # progress bar
        setTxtProgressBar(pb, iScen/nScen)
      }
      close(pb)
      
      ### simulated dates
      sim.date.6h = rep(sim.date, each = 4)
      
      ### simulated areal estimates
      list.arealSim.6h = list()
      for(sal in list.sal){
        # Thiessen weights
        TW = get.TW.stations(sal,MetaData$xCoord,MetaData$yCoord)
        # areal estimates
        list.arealSim.6h[[sal]] = get.areal.sim(var.Gauge=sim.array.6h,thiessen.weights=as.matrix(TW),varType=varType,
                                             station.Elevation=MetaData$Elevation,basin.Elevation=basin.mean.elevation[[sal]])
      }
      
      ###  save areal estimates
      save(obs.mat.6h,obs.date.6h,sim.array.6h,sim.date.6h,list.arealObs.6h,list.arealSim.6h,
           file=f.6h.Estimates)
    }
  }
 
  ###======================= START PDF ==========================###
  
  #### INITIALIZE ####
  # number of rows
  ncol.fig = 3
  nrow.fig = 6
  n.fig.page = nrow.fig*ncol.fig
  
  # start pdf figure
  pdf(file=paste0(f.main,'diagnostic.pdf'), width = 8.27, height = 11.69)
  layout(matrix(c(1:n.fig.page),nrow.fig,ncol.fig))
  par(oma=c(2.5,4,4,0),mar=c(2.5,1.5,0.1,0.1),mgp=c(2,0.5,0),mfrow=c(nrow.fig,ncol.fig))
  
  ###===================== PRECIPITATION ========================###
  if(varType=='Prec'){
    
    if(eval.6H){
      ############ 6-Hour MAXIMA ############
      # at the stations
      sel = match(StationPrec_SEL,colnames(obs.mat.6h))
      plot.scatter.fun(obs=obs.mat.6h[,sel], sim=sim.array.6h[,sel,], obs.date=obs.date.6h,sim.date=sim.date.6h, 
                       funYear=max, ylab='Maxima [mm]', lab=StationPrec_SEL, cex.axis=my.cex.axis, xaxt='n',
                       main='6-Hour maxima',fillPage = F,col.box='red')
      mtext('At the stations (6/105)',side=2,at=0.85,outer=T,col='red')
      # at the 15 basins
      sel = match(Basin15_SEL,colnames(list.arealSim.6h[['BASIN15']]))
      plot.scatter.fun(obs=list.arealObs.6h[['BASIN15']][,sel], sim=list.arealSim.6h[['BASIN15']][,sel,], obs.date=obs.date.6h, 
                       sim.date=sim.date.6h, funYear=max, ylab='', lab=Basin15_SEL, cex.axis=my.cex.axis, xaxt='n', 
                       main='',fillPage = F,col.box='blue')
      mtext('At the basins (6/15)',side=2,at=0.5,outer=T,col='blue')
      # at the main basins
      plot.scatter.fun(obs=list.arealObs.6h[['BASIN5']], sim=list.arealSim.6h[['BASIN5']], obs.date=obs.date.6h, sim.date=sim.date.6h, 
                       funYear=max, ylab='', lab=colnames(list.arealSim[['BASIN5']]), cex.axis=my.cex.axis, xaxt='n', 
                       main='',fillPage = F, col.box='green')
      plot.scatter.fun(obs=matrix(list.arealObs.6h[['BASIN1']]), sim=list.arealSim.6h[['BASIN1']], obs.date=obs.date.6h, sim.date=sim.date.6h, 
                       funYear=max, ylab='', lab='Whole catchment', cex.axis=my.cex.axis, xaxt='n', main='',fillPage = F, col.box='orange')
      mtext('At the main basins',side=2,at=0.2,outer=T,col='green')
    }
    
    ############ 1-Day MAXIMA ############
    # at the stations
    sel = match(StationPrec_SEL,colnames(obs.mat))
    plot.scatter.fun(obs=obs.mat[,sel], sim=sim.array[,sel,], obs.date=obs.date, sim.date=sim.date, 
                     funYear=max, ylab='Maxima [mm]', lab=StationPrec_SEL, cex.axis=my.cex.axis, xaxt='n',
                     main='1-Day maxima',fillPage = F,col.box='red')
    mtext('At the stations (6/105)',side=2,at=0.85,outer=T,col='red')
    # at the 15 basins
    sel = match(Basin15_SEL,colnames(list.arealSim[['BASIN15']]))
    plot.scatter.fun(obs=list.arealObs[['BASIN15']][,sel], sim=list.arealSim[['BASIN15']][,sel,], obs.date=obs.date, 
                     sim.date=sim.date, funYear=max, ylab='', lab=Basin15_SEL, cex.axis=my.cex.axis, xaxt='n', 
                     main='',fillPage = F,col.box='blue')
    mtext('At the basins (6/15)',side=2,at=0.5,outer=T,col='blue')
    # at the main basins
    plot.scatter.fun(obs=list.arealObs[['BASIN5']], sim=list.arealSim[['BASIN5']], obs.date=obs.date, sim.date=sim.date, 
                     funYear=max, ylab='', lab=colnames(list.arealSim[['BASIN5']]), cex.axis=my.cex.axis, xaxt='n', 
                     main='',fillPage = F, col.box='green')
    plot.scatter.fun(obs=matrix(list.arealObs[['BASIN1']]), sim=list.arealSim[['BASIN1']], obs.date=obs.date, sim.date=sim.date, 
                     funYear=max, ylab='', lab='Whole catchment', cex.axis=my.cex.axis, xaxt='n', main='',fillPage = F, col.box='orange')
    mtext('At the main basins',side=2,at=0.2,outer=T,col='green')
    
    ############ 3-Day MAXIMA ############
    # at the stations
    sel = match(StationPrec_SEL,colnames(obs.mat))
    plot.scatter.fun(obs=obs.mat[,sel], sim=sim.array[,sel,], obs.date=obs.date, sim.date=sim.date, 
                     funYear=max, n.day.cumul=3, ylab='Maxima [mm]', lab=StationPrec_SEL, cex.axis=my.cex.axis, 
                     xaxt='n', main='3-Day maxima',fillPage = F,col.box='red')
    mtext('At the stations (6/105)',side=2,at=0.85,outer=T,col='red')
    # at the 15 basins
    sel = match(Basin15_SEL,colnames(list.arealSim[['BASIN15']]))
    plot.scatter.fun(obs=list.arealObs[['BASIN15']][,sel], sim=list.arealSim[['BASIN15']][,sel,], obs.date=obs.date, sim.date=sim.date,
                     funYear=max, n.day.cumul=3, ylab='', lab=Basin15_SEL, cex.axis=my.cex.axis, 
                     xaxt='n', main='',fillPage = F,col.box='blue')
    mtext('At the basins (6/15)',side=2,at=0.5,outer=T,col='blue')
    # at the main basins
    plot.scatter.fun(obs=list.arealObs[['BASIN5']], sim=list.arealSim[['BASIN5']], obs.date=obs.date, sim.date=sim.date, 
                     funYear=max, n.day.cumul=3, ylab='', lab=colnames(list.arealSim[['BASIN5']]), cex.axis=my.cex.axis, 
                     xaxt='n', main='',fillPage = F, col.box='green')
    plot.scatter.fun(obs=matrix(list.arealObs[['BASIN1']]), sim=list.arealSim[['BASIN1']], obs.date=obs.date, sim.date=sim.date, 
                     funYear=max, n.day.cumul=3, ylab='', lab='Whole catchment', cex.axis=my.cex.axis, 
                     xaxt='n', main='',fillPage = F, col.box='orange')
    mtext('At the main basins',side=2,at=0.2,outer=T,col='green')
    
    ############ STATS AT THE STATIONS ############
    # dry spell length
    sel = match(StationPrec_SEL,colnames(obs.mat))
    plot.cdf.spell.length(obs=obs.mat[,sel], sim=sim.array[,sel,], obs.date=obs.date, sim.date=sim.date, type='dry', 
                          obs.th=obs.th, sim.th=sim.th, xlab='',ylab='', lab=StationPrec_SEL, cex.axis=my.cex.axis, 
                          main='Statistics at the stations',fillPage = F,add.legend=F,col.box='red')
    mtext('Dry spell length [days]',side=1,line=-56,outer=T,col='red')
    mtext('Exceedance probability',side=2,at=0.85,outer=T,col='red')
    # wet spell length
    sel = match(StationPrec_SEL,colnames(obs.mat))
    plot.cdf.spell.length(obs=obs.mat[,sel], sim=sim.array[,sel,], obs.date=obs.date, sim.date=sim.date, type='wet', 
                          obs.th=obs.th, sim.th=sim.th, xlab='',ylab='', lab=StationPrec_SEL, cex.axis=my.cex.axis, 
                          main='',fillPage = F,add.legend=F,col.box='blue')
    mtext('Wet spell length [days]',side=1,line=-28.5,outer=T,col='blue')
    mtext('Exceedance probability',side=2,at=0.53,outer=T,col='blue')
    # Pr. of staying in a wet state
    sel = match(StationPrec_SEL,colnames(obs.mat))
    plot.daily.stats.by.month(obs=obs.mat[,sel], sim=sim.array[,sel,], obs.date=obs.date, sim.date=sim.date, 
                              obs.th=obs.th, sim.th=sim.th, type='WWTP', xlab='', ylab='', lab=StationPrec_SEL, cex.axis=my.cex.axis, 
                              xaxt='n', main='',fillPage = F,col.box='green')
    mtext('Month',side=1,outer=T,line=-1,col='green')
    mtext('Pr. of staying in a wet state',side=2,at=0.18,outer=T,col='green')
    
    ############ Monthly Mean ############
    myfun = function(x) mean(x,na.rm=T)
    # at the stations
    sel = match(StationPrec_SEL,colnames(obs.mat))
    plot.monthly.stats.by.month(obs=obs.mat[,sel], sim=sim.array[,sel,],type=myfun,
                                obs.date=obs.date, sim.date=sim.date,
                                ylab='Monthly mean [mm]', lab=StationPrec_SEL,
                                main='Monthly mean',
                                fillPage = F,col.box='red',
                                cex.axis=my.cex.axis, 
                                xaxt='n')
    mtext('At the stations (6/105)',side=2,at=0.85,outer=T,col='red')
    # at the 15 basins
    sel = match(Basin15_SEL,colnames(list.arealSim[['BASIN15']]))
    plot.monthly.stats.by.month(obs=list.arealObs[['BASIN15']][,sel], sim=list.arealSim[['BASIN15']][,sel,],type=myfun,
                                obs.date=obs.date, sim.date=sim.date,
                                ylab='', lab=Basin15_SEL,
                                main='',
                                fillPage = F,col.box='blue',
                                cex.axis=my.cex.axis, 
                                xaxt='n')
    mtext('At the basins (6/15)',side=2,at=0.5,outer=T,col='blue')
    # at the main basins
    plot.monthly.stats.by.month(obs=list.arealObs[['BASIN5']], sim=list.arealSim[['BASIN5']],type=myfun,
                                obs.date=obs.date, sim.date=sim.date,
                                ylab='', lab=colnames(list.arealSim[['BASIN5']]),
                                main='',
                                fillPage = F,col.box='green',
                                cex.axis=my.cex.axis, 
                                xaxt='n')
    plot.monthly.stats.by.month(obs=matrix(list.arealObs[['BASIN1']]), sim=list.arealSim[['BASIN1']],type=myfun,
                                obs.date=obs.date, sim.date=sim.date,
                                ylab='', lab='Whole catchment',
                                main='',
                                fillPage = F,col.box='orange',
                                cex.axis=my.cex.axis, 
                                xaxt='n')
    mtext('At the main basins',side=2,at=0.2,outer=T,col='green')
    
    ############ Monthly SD ############
    myfun = function(x) sd(x,na.rm=T)
    # at the stations
    sel = match(StationPrec_SEL,colnames(obs.mat))
    plot.monthly.stats.by.month(obs=obs.mat[,sel], sim=sim.array[,sel,],type=myfun,
                                obs.date=obs.date, sim.date=sim.date,
                                ylab='Monthly Standard deviation [mm]', lab=StationPrec_SEL,
                                main='Monthly Standard deviation',
                                fillPage = F,col.box='red',
                                cex.axis=my.cex.axis, 
                                xaxt='n')
    mtext('At the stations (6/26)',side=2,at=0.85,outer=T,col='red')
    # at the 15 basins
    sel = match(Basin15_SEL,colnames(list.arealSim[['BASIN15']]))
    plot.monthly.stats.by.month(obs=list.arealObs[['BASIN15']][,sel], sim=list.arealSim[['BASIN15']][,sel,],type=myfun,
                                obs.date=obs.date, sim.date=sim.date,
                                ylab='', lab=Basin15_SEL,
                                main='',
                                fillPage = F,col.box='blue',
                                cex.axis=my.cex.axis, 
                                xaxt='n')
    mtext('At the basins (6/15)',side=2,at=0.5,outer=T,col='blue')
    # at the main basins
    plot.monthly.stats.by.month(obs=list.arealObs[['BASIN5']], sim=list.arealSim[['BASIN5']],type=myfun,
                                obs.date=obs.date, sim.date=sim.date,
                                ylab='', lab=colnames(list.arealSim[['BASIN5']]),
                                main='',
                                fillPage = F,col.box='green',
                                cex.axis=my.cex.axis, 
                                xaxt='n')
    plot.monthly.stats.by.month(obs=matrix(list.arealObs[['BASIN1']]), sim=list.arealSim[['BASIN1']],type=myfun,
                                obs.date=obs.date, sim.date=sim.date,
                                ylab='', lab='Whole catchment',
                                main='',
                                fillPage = F,col.box='orange',
                                cex.axis=my.cex.axis, 
                                xaxt='n')
    mtext('At the main basins',side=2,at=0.2,outer=T,col='green')
    
    ###===================== TEMPERATURE ========================###
  }else if(varType=='Temp'){
    if(eval.6H){
      ############ 6-Hour MAXIMA ############
      # at the stations
      sel = match(StationTemp_SEL,colnames(obs.mat.6h))
      plot.scatter.fun(obs=obs.mat.6h[,sel], sim=sim.array.6h[,sel,], obs.date=obs.date.6h,sim.date=sim.date.6h, 
                       funYear=max, ylab='Maxima [mm]', lab=StationTemp_SEL, cex.axis=my.cex.axis, xaxt='n',
                       main='6-Hour maxima',fillPage = F,col.box='red',na.action = na.omit)
      mtext('At the stations (6/26)',side=2,at=0.85,outer=T,col='red')
      # at the 15 basins
      sel = match(Basin15_SEL,colnames(list.arealSim.6h[['BASIN15']]))
      plot.scatter.fun(obs=list.arealObs.6h[['BASIN15']][,sel], sim=list.arealSim.6h[['BASIN15']][,sel,], obs.date=obs.date.6h, 
                       sim.date=sim.date.6h, funYear=max, ylab='', lab=Basin15_SEL, cex.axis=my.cex.axis, xaxt='n', 
                       main='',fillPage = F,col.box='blue',na.action = na.omit)
      mtext('At the basins (6/15)',side=2,at=0.5,outer=T,col='blue')
      # at the main basins
      plot.scatter.fun(obs=list.arealObs.6h[['BASIN5']], sim=list.arealSim.6h[['BASIN5']], obs.date=obs.date.6h, sim.date=sim.date.6h, 
                       funYear=max, ylab='', lab=colnames(list.arealSim[['BASIN5']]), cex.axis=my.cex.axis, xaxt='n', 
                       main='',fillPage = F, col.box='green',na.action = na.omit)
      plot.scatter.fun(obs=matrix(list.arealObs.6h[['BASIN1']]), sim=list.arealSim.6h[['BASIN1']], obs.date=obs.date.6h, sim.date=sim.date.6h, 
                       funYear=max, ylab='', lab='Whole catchment', cex.axis=my.cex.axis, xaxt='n', main='',fillPage = F, col.box='orange',
                       na.action = na.omit)
      mtext('At the main basins',side=2,at=0.2,outer=T,col='green')
    }
    ############ 1-Day MAXIMA ############
    # at the stations
    sel = match(StationTemp_SEL,colnames(obs.mat))
    plot.scatter.fun(obs=obs.mat[,sel], sim=sim.array[,sel,], obs.date=obs.date, sim.date=sim.date, 
                     funYear=max, ylab='Maxima [Â°C]', lab=StationTemp_SEL, cex.axis=my.cex.axis, xaxt='n',
                     main='1-Day maxima',fillPage = F,col.box='red',na.action = na.omit)
    mtext('At the stations (6/26)',side=2,at=0.85,outer=T,col='red')
    # at the 15 basins
    sel = match(Basin15_SEL,colnames(list.arealSim[['BASIN15']]))
    plot.scatter.fun(obs=list.arealObs[['BASIN15']][,sel], sim=list.arealSim[['BASIN15']][,sel,], obs.date=obs.date, 
                     sim.date=sim.date, funYear=max, ylab='', lab=Basin15_SEL, cex.axis=my.cex.axis, xaxt='n', 
                     main='',fillPage = F,col.box='blue',na.action = na.omit)
    mtext('At the basins (6/15)',side=2,at=0.5,outer=T,col='blue')
    # at the main basins
    plot.scatter.fun(obs=list.arealObs[['BASIN5']], sim=list.arealSim[['BASIN5']], obs.date=obs.date, sim.date=sim.date, 
                     funYear=max, ylab='', lab=colnames(list.arealSim[['BASIN5']]), cex.axis=my.cex.axis, xaxt='n', 
                     main='',fillPage = F, col.box='green',na.action = na.omit)
    plot.scatter.fun(obs=matrix(list.arealObs[['BASIN1']]), sim=list.arealSim[['BASIN1']], obs.date=obs.date, sim.date=sim.date, 
                     funYear=max, ylab='', lab='Whole catchment', cex.axis=my.cex.axis, xaxt='n', 
                     main='',fillPage = F, col.box='orange',na.action = na.omit)
    mtext('At the main basins',side=2,at=0.2,outer=T,col='green')
    
    ############ 1-Day MINIMA ############
    # at the stations
    sel = match(StationTemp_SEL,colnames(obs.mat))
    plot.scatter.fun(obs=obs.mat[,sel], sim=sim.array[,sel,], obs.date=obs.date, sim.date=sim.date, 
                     funYear=min, ylab='Minima [Â°C]', lab=StationTemp_SEL, cex.axis=my.cex.axis, xaxt='n',
                     main='1-Day minima',fillPage = F,col.box='red',na.action = na.omit)
    mtext('At the stations (6/26)',side=2,at=0.85,outer=T,col='red')
    # at the 15 basins
    sel = match(Basin15_SEL,colnames(list.arealSim[['BASIN15']]))
    plot.scatter.fun(obs=list.arealObs[['BASIN15']][,sel], sim=list.arealSim[['BASIN15']][,sel,], obs.date=obs.date, 
                     sim.date=sim.date, funYear=min, ylab='', lab=Basin15_SEL, cex.axis=my.cex.axis, xaxt='n', 
                     main='',fillPage = F,col.box='blue',na.action = na.omit)
    mtext('At the basins (6/15)',side=2,at=0.5,outer=T,col='blue')
    # at the main basins
    plot.scatter.fun(obs=list.arealObs[['BASIN5']], sim=list.arealSim[['BASIN5']], obs.date=obs.date, sim.date=sim.date, 
                     funYear=min, ylab='', lab=colnames(list.arealSim[['BASIN5']]), cex.axis=0.9, xaxt='n', 
                     main='',fillPage = F, col.box='green',na.action = na.omit)
    plot.scatter.fun(obs=matrix(list.arealObs[['BASIN1']]), sim=list.arealSim[['BASIN1']], obs.date=obs.date, sim.date=sim.date, 
                     funYear=min, ylab='', lab='Whole catchment', cex.axis=my.cex.axis, xaxt='n', 
                     main='',fillPage = F, col.box='orange',na.action = na.omit)
    mtext('At the main basins',side=2,at=0.2,outer=T,col='green')
    
    ############ MONTHLY MEAN ############
    myfun = function(x) mean(x,na.rm=T)
    # at the stations
    sel = match(StationTemp_SEL,colnames(obs.mat))
    plot.daily.stats.by.month(obs=obs.mat[,sel], sim=sim.array[,sel,], obs.date=obs.date, sim.date=sim.date, 
                              type=myfun, ylab='Mean [Â°C]', lab=StationTemp_SEL, cex.axis=my.cex.axis, xaxt='n',
                              main='Monthly mean',fillPage = F,col.box='red')
    mtext('At the stations (6/26)',side=2,at=0.85,outer=T,col='red')
    # at the 15 basins
    sel = match(Basin15_SEL,colnames(list.arealSim[['BASIN15']]))
    plot.daily.stats.by.month(obs=list.arealObs[['BASIN15']][,sel], sim=list.arealSim[['BASIN15']][,sel,], obs.date=obs.date, 
                              sim.date=sim.date, type=myfun, ylab='', lab=Basin15_SEL, cex.axis=my.cex.axis, xaxt='n', 
                              main='',fillPage = F,col.box='blue')
    mtext('At the basins (6/15)',side=2,at=0.5,outer=T,col='blue')
    # at the main basins
    plot.daily.stats.by.month(obs=list.arealObs[['BASIN5']], sim=list.arealSim[['BASIN5']], obs.date=obs.date, sim.date=sim.date, 
                              type=myfun, ylab='', lab=colnames(list.arealSim[['BASIN5']]), cex.axis=my.cex.axis, xaxt='n', 
                              main='',fillPage = F, col.box='green')
    plot.daily.stats.by.month(obs=matrix(list.arealObs[['BASIN1']]), sim=list.arealSim[['BASIN1']], obs.date=obs.date, sim.date=sim.date, 
                              type=myfun, ylab='', lab='Whole catchment', cex.axis=my.cex.axis, xaxt='n', 
                              main='',fillPage = F, col.box='orange')
    mtext('At the main basins',side=2,at=0.2,outer=T,col='green')
  }
  
  ### CLOSE ###
  dev.off()
}



#==============================================================================
pdf.systematic.diagnostic.stations = function(f.main, # folder that will contain the results
                                     varType, # type of variable: 'Prec' or 'Temp'
                                     obs.mat, # matrix of observations
                                     obs.date, # vector of date for the observations
                                     obs.th, # threshold of the dry/wet days for obs
                                     sim.array, # array of simulated scenarios
                                     sim.date, # vector of date for the scenarios
                                     sim.th, # threshold of the dry/wet days for sim
                                     MetaData, # Metadata for the stations
                                     modelType='GWEX', # Type of model: 'GWEX' or 'SCAMP'
                                     my.cex.axis=0.9 # expansion for axis ticks
){
  #### SNV folder
  dir.SVN = Sys.getenv("SVN_EXAR")
  
  ###======================= START PDF ==========================###
  
  #### INITIALIZE ####
  # number of rows
  ncol.fig = 3
  nrow.fig = 6
  n.fig.page = nrow.fig*ncol.fig
  
  # start pdf figure
  pdf(file=paste0(f.main,'diagnostic.pdf'), width = 8.27, height = 11.69)
  layout(matrix(c(1:n.fig.page),nrow.fig,ncol.fig))
  par(oma=c(2.5,4,4,0),mar=c(2.5,1.5,0.1,0.1),mgp=c(2,0.5,0),mfrow=c(nrow.fig,ncol.fig))
  
  ###===================== PRECIPITATION ========================###
  if(varType=='Prec'){
    
    
    ############ 1-Day MAXIMA ############
    # at the stations
    plot.scatter.fun(obs=obs.mat, sim=sim.array, obs.date=obs.date, sim.date=sim.date, 
                     funYear=max, ylab='Maxima [mm]', lab=colnames(obs.mat), cex.axis=my.cex.axis, xaxt='n',
                     main='1-Day maxima',fillPage = T,col.box='red')
     
    ############ 3-Day MAXIMA ############
    # at the stations
    plot.scatter.fun(obs=obs.mat, sim=sim.array, obs.date=obs.date, sim.date=sim.date, 
                     funYear=max, n.day.cumul=3, ylab='Maxima [mm]', lab=colnames(obs.mat), cex.axis=my.cex.axis, 
                     xaxt='n', main='3-Day maxima',fillPage = T,col.box='red')
    
    
    ############ STATS AT THE STATIONS ############
    # dry spell length
    plot.cdf.spell.length(obs=obs.mat, sim=sim.array, obs.date=obs.date, sim.date=sim.date, type='dry', 
                          obs.th=obs.th, sim.th=sim.th, xlab='',ylab='', lab=colnames(obs.mat), cex.axis=my.cex.axis, 
                          main='Exceedance probability',fillPage = T,add.legend=F,col.box='red')
    mtext('Dry spell length [days]',side=1,line=-1,outer=T)
    # wet spell length
    plot.cdf.spell.length(obs=obs.mat, sim=sim.array, obs.date=obs.date, sim.date=sim.date, type='wet', 
                          obs.th=obs.th, sim.th=sim.th, xlab='',ylab='', lab=colnames(obs.mat), cex.axis=my.cex.axis, 
                          main='Exceedance probability',fillPage = T,add.legend=F,col.box='blue')
    mtext('Wet spell length [days]',side=1,line=-1,outer=T)
    # Pr. of staying in a wet state
    plot.daily.stats.by.month(obs=obs.mat, sim=sim.array, obs.date=obs.date, sim.date=sim.date, 
                              obs.th=obs.th, sim.th=sim.th, type='WWTP', xlab='', ylab='', lab=colnames(obs.mat), cex.axis=my.cex.axis, 
                              xaxt='n', main='Pr. of staying in a wet state',fillPage = T,col.box='green')
    mtext('Month',side=1,outer=T,line=-1)
    
    ###===================== TEMPERATURE ========================###
  }else if(varType=='Temp'){
    ############ 1-Day MAXIMA ############
    # at the stations
    plot.scatter.fun(obs=obs.mat, sim=sim.array, obs.date=obs.date, sim.date=sim.date, 
                     funYear=max, ylab='Maxima [°C]', lab=colnames(obs.mat), cex.axis=my.cex.axis, xaxt='n',
                     main='1-Day maxima',fillPage = T,col.box='red',na.action = na.omit)
   
    
    ############ 1-Day MINIMA ############
    # at the stations
    plot.scatter.fun(obs=obs.mat, sim=sim.array, obs.date=obs.date, sim.date=sim.date, 
                     funYear=min, ylab='Minima [°C]', lab=colnames(obs.mat), cex.axis=my.cex.axis, xaxt='n',
                     main='1-Day minima',fillPage = T,col.box='red',na.action = na.omit)
    
    
    ############ MONTHLY MEAN ############
    myfun = function(x) mean(x,na.rm=T)
    # at the stations
    plot.daily.stats.by.month(obs=obs.mat, sim=sim.array, obs.date=obs.date, sim.date=sim.date, 
                              type=myfun, ylab='Mean [°C]', lab=colnames(obs.mat), cex.axis=my.cex.axis, xaxt='n',
                              main='Monthly mean',fillPage = T,col.box='red')
  }
  
  ### CLOSE ###
  dev.off()
}



