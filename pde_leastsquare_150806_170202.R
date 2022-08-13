# viability and nash data from death curves8 150806 and subpopulation growth from 170202
library(optimx)
library(ggplot2)
library(deSolve)
library(reshape2)
setwd("C:/Siavash/Codes/pde/pde_plots_data")

###############################################################################
# read in the spreadsheet of sample metadata
flaskIDs<-read.csv('formaldehyde_switching_flaskIDs.csv', header=T) 
flaskIDs$treatment<-factor(flaskIDs$treatment)
flaskIDs$transfer<-factor(flaskIDs$transfer)

#######################################################################
# CFU/mL

CFU<-read.csv('formaldehyde_switching_colonycounts.csv', header=T)
CFU$cellcounts<-as.numeric(as.character(CFU$cellcounts))
CFU$HCHO_in_plate<-factor(CFU$HCHO_in_plate)
CFU_mean<-dcast(CFU, time+flask+HCHO_in_plate~., value.var='cellcounts', fun.aggregate=mean)
colnames(CFU_mean)[4]<-'CFU_permL_mean'
CFU_mean$CFU_permL_sd<-dcast(CFU, time+flask+HCHO_in_plate~., value.var='cellcounts', fun.aggregate=sd)[,4]
CFU<-merge(CFU_mean, flaskIDs, by='flask')
colnames(CFU)[3]<-'Population'
#levels(CFU$Population) <-c('all cells','tolerant to 2mM','tolerant to 4mM')
levels(CFU$transfer)<-c('original population','removed from formaldehyde')
CFU<-subset(CFU, !(flask=='G' & time==28.5)) # removing outlier. Something clearly must have gone wrong with dilutions here.

#######################################################################
# Formaldehyde
Nash<-read.csv('HCHO_switching_Nashassay.csv', header=T)
wellIDs<-read.csv('HCHO_switching_Nashassay_wellIDs.csv', header=T)
Nash<-melt(Nash, id.vars='time', variable.name='well', value.name='A412')
Nash<-merge(Nash, wellIDs, by='well')

# deal with standard curves
Nash_std<-subset(Nash, row %in% c('E','F','G') & time %in% c(1.5, 20.5, 48.5, 69.5))
#plot_stcurve<-ggplot() + geom_point(data=Nash_std, aes(x=A412, y=mM.formaldehyde, color=factor(time)))
#plot_stcurve
Nash_stcurve<-lm(Nash_std$mM.formaldehyde~Nash_std$A412)
Nash_std$predict<-predict(Nash_stcurve)
#plot_stcurve + geom_line(data=Nash_std, aes(x=A412, y=predict))
Nash_slope<-Nash_stcurve$coefficients[[2]]
Nash_intercept<-Nash_stcurve$coefficients[[1]]

Nash_sampledata<-subset(Nash, row %in% c('A','B','C') & is.na(A412)==F)
Nash_sampledata$mM.formaldehyde<-Nash_sampledata$A412*Nash_slope+Nash_intercept

# something funny happened at time-0. I think it needs to be blanked separately
Nash_std0<-subset(Nash, row %in% c('E','F','G') & time==0)
Nash_stcurve0<-lm(Nash_std0$mM.formaldehyde~Nash_std0$A412)
Nash_std0$predict<-predict(Nash_stcurve0)
#plot(Nash_std0$mM.formaldehyde~Nash_std0$A412)
#lines(Nash_std0$predict~Nash_std0$A412, col='red')
Nash_slope0<-Nash_stcurve0$coefficients[[2]]
Nash_intercept0<-Nash_stcurve0$coefficients[[1]]
Nash_sampledata$mM.formaldehyde[Nash_sampledata$time==0]<-Nash_sampledata$A412[Nash_sampledata$time==0]*Nash_slope0+Nash_intercept0
Nash_sampledata<-merge(Nash_sampledata, flaskIDs, by='flask')

# average the Nash technical replicates
Nash_sampledata<-dcast(Nash_sampledata, flask+column+treatment+rep+transfer+time~row,value.var='mM.formaldehyde')
Nash_sampledata$HCHO_mean<-rowMeans(Nash_sampledata[,7:9])
Nash_sampledata$HCHO_sd<-apply(Nash_sampledata[,7:9],1,sd)

# merge datasets
Nash_mean<-merge(subset(Nash_sampledata, select=c('flask','time','HCHO_mean','HCHO_sd')), CFU)

## plot, finally
plot_HCHO_CFU<-ggplot() + geom_line(data=Nash_mean, aes(x=time, y=HCHO_mean, group=flask))
plot_HCHO_CFU<-plot_HCHO_CFU + geom_point(data=Nash_mean, size=1.5, aes(x=time, y=HCHO_mean))
plot_HCHO_CFU<-plot_HCHO_CFU + geom_errorbar(data=Nash_mean, size=1, width=1, aes(x=time, ymin=HCHO_mean-HCHO_sd, ymax=HCHO_mean+HCHO_sd))
plot_HCHO_CFU<-plot_HCHO_CFU + geom_line(data=Nash_mean, aes(x=time, y=CFU_permL_mean, group=interaction(flask, Population), color=Population))
plot_HCHO_CFU<-plot_HCHO_CFU + geom_point(data=Nash_mean, size=1.5, aes(x=time, y=CFU_permL_mean, color=Population))
plot_HCHO_CFU<-plot_HCHO_CFU + geom_errorbar(data=Nash_mean, size=1, width=1, 
                                             aes(x=time, ymin=CFU_permL_mean-CFU_permL_sd, ymax=CFU_permL_mean+CFU_permL_sd, color=Population))
plot_HCHO_CFU<-plot_HCHO_CFU + theme_bw() + scale_y_log10(breaks=c(1e-01, 5e-01, 1, 5, 10, 1e+02, 1e+04, 1e+06, 1e+08)) + facet_grid(.~treatment, scale='free_x')
plot_HCHO_CFU<-plot_HCHO_CFU + xlab('Time (hours)') + ylab('Formaldehyde (mM)(black)\nviable cells (CFU/mL/)(color)')
plot_HCHO_CFU<-plot_HCHO_CFU + ggtitle('CM2730 growth on methanol', subtitle='panels = formaldehyde in medium (mM)')
plot_HCHO_CFU
#jpeg('subpopulations_with_Nash_170213.jpg', width=10, height=5, res=200, units='in')
#plot_HCHO_CFU
#dev.off()

setwd("C:/Siavash/Codes/phenmod/plots_data")

death_viability <- read.csv("killingcurve_8615_results_recalc_160923_CFUs_final_recalc160923.csv")
death_nash <- read.csv("RECOUNT_killingcurve_8615_alldata_final_nash_final_data2.csv")

tmp2 <- merge(death_viability,death_nash, all.x = TRUE, all.y = TRUE)
tmp2 <-tmp2[tmp2$strain !="CM3745",]
death_total <- subset(tmp2, select = -c(strain,NashA412,standard_curve))
mean_death <- mean(death_total[death_total$time==0,]$cells)
#death_total$cells[death_total$time==0] <- mean_death
#colnames(growth_total)[5] <- "formaldehyde_nash"
#combined <- rbind(growth_total,death_total)
#death_total <- combined
death_total <- death_total[order(death_total$time),]

names(death_total) <- c("time","treatment", "cells","rep", "nash")
death_total$nash <- as.numeric(as.character(death_total$nash))
#death_total$cells[death_total$cells == 0.0] <- 1
#death_total$cells <- death_total$cells + 1
death_total$time <- death_total$time/60 
# ploting the data 

Treatment <- as.factor(death_total$treatment)
plot_death_n <-ggplot() + geom_line(data=death_total, size=1.2,aes(x=time, y=cells, col=Treatment, group=treatment)) + scale_y_log10() + scale_x_continuous(limits = c(0,3))  
plot_death_n 

Treatment <- as.factor(death_total[!is.na(death_total$nash),]$treatment)
plot_death_fex<-ggplot() + geom_line(data=death_total[!is.na(death_total$nash),], size=1.2,aes(x=time, y=nash, col=Treatment, group=treatment)) + scale_x_continuous(limits = c(0,3))  
plot_death_fex

nash_mean2 <- Nash_mean[,c(2,3,5,6,8,9)]
colnames(nash_mean2)[2:4] <- c("nash","pop","cells") 
death_total$rep <- "A"
death_total$pop <- 0
#nash_mean2 <- subset(nash_mean2,treatment==4)
all_data <- rbind(death_total,nash_mean2)
all_data$treatment <- as.numeric(all_data$treatment)
NClasses <- 100

#maxTol = highest tolerance category
maxTol = 10;

#h = step size (function of maxTol and n)
h <- maxTol/NClasses;

x <- seq(0,10,length.out=NClasses+1)
#DE parameters

f0 <- 4 #formaldehyde at time 0
m0 <- 15 #methanol at time 0
s0 <- 0


pde <- function (time, state, parms) {
  with (as.list(parms), {
    m <- state[1] #methanol concentration
    s <- state[2] #succinate concentration
    f <- state[3] #formaldehyde concentration
    n <- state[4:(N+4)] #N phenotypic states
    
    dm <- -vmax*m/(m+km)*sum(n) 
    ds <- -vmaxs*s/(s+ks)*sum(n) 
    df <-  vmax*(m/(m+km)-gamma*f/(f+kf))*sum(n)
    dn <- n*(f*rf/(f+kf)+s*rs/(s+ks)) - n*alpha*f*(x<f) 
    
    # Defining the Flux P1: Chemotaxis toward the F concentration in space 
    nn <- c(n,0)
    nn[1] <- 0
    FluxP1 <- rep(0,length(dn))
    if (f<0.03){FluxP1 <- V*diff(nn)}
    #FluxP1 <- -V*diff(nn)
    
    FluxP2 <- D * diff(c(n[1], n, n[N]), diff=2) ###note that there should technically be division by dx^2
    ###however, to recover d^2 n/dx^2, you have to then multiply by dx^2, making it pointless to do the division
    
    ## Rate of change = Flux gradient + Biology
    dn <- dn  + FluxP2 + FluxP1
    
    return (list(c(dm,ds, df, dn)))
  })
}

model_multi_levels<-function(model, time_range, baseline_state, parameters, f_init_range){
  out<-data.frame() # create an empty dataframe for putting all of the outputs into
  for (i in f_init_range) { # for-loop to run through all initial external formaldehyde levels
    state_i<-baseline_state # start off with your baseline state variables, of which all except f will remain the same
    state_i[3]<-i # replace the f from the baseline set of state variables with the f value you want to model
    out_i <- as.data.frame(ode(y=state_i, times=time_range, func=model, parms=parameters, method="lsoda")) # run the model
    f_init<-rep(i, length(time_range)) # create a column called "f_init" to append to the data tables, so you know which initial f value yielded these outputs
    out<-rbind(out, cbind(out_i, f_init)) # append the results of this model run, along with a column indicating the initial formaldehyde concentration, to the master output dataframe
  }
  colnames(out)[1:4] <- c("time","methanol","succinate","formaldehyde")
  colnames(out)[5:(NClasses+5)] <- paste("n",seq(0,100),sep="") # Non-cumulative populations 
  colnames(out)[(NClasses+6)] <- "f_init"
  out$"N0" <- rowSums(out[5:(NClasses+5)])
  out$"N2" <- rowSums(out[25:(NClasses+5)])
  out$"N4" <- rowSums(out[45:(NClasses+5)])
  return(out)
}


fixedparms = list(
  alpha = 0.2, #alpha
  km = 0.02, #K[m]
  ks = 0.02,
  kf = 0.2, #K[f]1e-3
  vmax = 2.5e-8, #V[max]1e-9
  vmaxs = 1e-8, 
  #gamma = 2, #gamma
  N = NClasses, #number of states in x
  rf = 0.18, #growth rate
  rs = 0.26,
  V = 0,
  D = 0,
  x=x)

IHS <- function(x, theta){  # IHS transformation
  (1/theta)*asinh(theta * x)
}

all.resids <- NULL
resid.fun.pde <- function(gparms){
  #### Need to solve the DEs for each strain and concentration
  for(conc in unique(all_data$treatment)){
    rep_range <- "A"
    if(conc<5) rep_range <- c("A","B","C")
      for (replica in rep_range){
    dat <- subset(all_data, treatment==conc & rep==replica)
    times <- sort(unique(dat$time))
    N0 <- dat$cells[dat$time==0 & dat$pop==0]*1.11
    n0 <- N0*h*dgamma((0:NClasses)*h,shape = 0.7236445, rate = 2.9015106)
    n0[1] <- n0[2]
    #n0 <- N0*h*dgamma((0:NClasses)*h,shape = 13.03363, rate = 10.02192)
    fitstate <- c(15,0,conc,n0)
    working.parms <- c(fixedparms, 
                       gamma = gparms[1]          
    )
    out <- model_multi_levels(model=pde, time_range=times, baseline_state = fitstate, parameters=working.parms, f_init_range=conc)
    resid_n0 <- sum((IHS(out$N0, 1) - IHS(dat$cells[dat$pop==0], 1))^2,na.rm = T)
    resid_n2 <- 0 
    resid_n4 <- 0
    if(conc==4) {
    resid_n2 <- sum((IHS(out$N2, 1) - IHS(dat$cells[dat$pop==2], 1))^2,na.rm = T)
    resid_n4 <- sum((IHS(out$N4, 1) - IHS(dat$cells[dat$pop==4], 1))^2,na.rm = T)
    }
    #resid_nash <- sum(((out$formaldehyde - dat$nash)/dat$nash)^2, na.rm=T)
    resid_nash <- sum((out$formaldehyde - dat$nash)^2, na.rm=T)
    resids <- sum(resid_n0,resid_n2,resid_n4,resid_nash)
    #resids <- sum(resid_n0,resid_nash)
    
    all.resids <- rbind(all.resids, resids) 
    }
  }
  sum(all.resids) + 10^5*sum(gparms[gparms<0])^2  
}

#### Need to find reasonable starting parameters..
guess.parms.pde <- c(
  2
)
resid.fun.pde(guess.parms.pde)
optimize(resid.fun.pde,interval=c(1,20))
bestparms.pde <- optimx(guess.parms.pde, resid.fun.pde, method="Nelder-Mead", control=list(trace=2))

