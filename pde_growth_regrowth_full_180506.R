# A simulation for both growth, regrowth 
setwd("C:/Siavash/Codes/pde/pde_plots_data")
library(ggplot2)
library(MASS)
library(reshape2)
library(RColorBrewer)
library(deSolve)
library(matrixStats)
library(fitdistrplus)

####################### First growth data
CFUcounts <- read.csv("death_at_4mM_colonycounts_170814.csv")
DIL_K <- 1#90.90909090909 #this is the dilution constant that scanner (?) seems to have multiplied by
CFUcounts$Colonies <- CFUcounts$CFU/DIL_K #convert back to colony counts
#drop formaldehyde levels beyond observed growth
maxF <- max(CFUcounts$mM_HCHO[CFUcounts$Colonies > 0])
CFUcounts <- CFUcounts[CFUcounts$mM_HCHO <= maxF,]
#create xSet data frame for ensuring that each
#flask, [HCHO], and timepoint are present in 
#our data
tps <- unique(CFUcounts$Timepoint)
xSet <- data.frame(Flask = rep(LETTERS[1:3],each = (maxF+1)*length(tps)),mM_HCHO = rep(0:maxF, each=length(tps)), Timepoint = tps )
CFUcounts <- aggregate(CFUcounts$Colonies, by = CFUcounts[,c("Flask","Timepoint","mM_HCHO")], sum)
CFUcounts <- merge(xSet,CFUcounts,all.x=T)
CFUcounts$x[is.na(CFUcounts$x)]<- 0 #force NA to 0; I think this is safe???

CFUcounts <- CFUcounts[do.call(order,CFUcounts[,c("Flask","Timepoint","mM_HCHO")]),] #just make sure things are in order
#impose the assumption on the data and calculate the cells of particular tolerance level
for(f in LETTERS[1:3]){
  for(tt in tps) CFUcounts[CFUcounts$Flask == f & CFUcounts$Timepoint == tt,"x"] <- round(-1 * diff(c(rev(cummax(rev(CFUcounts[CFUcounts$Flask == f & CFUcounts$Timepoint == tt,"x"]))),0)))
}

#division by 9 due to 9 samples per sum (3 flasks, 3 reps)
#rounding to make the table look nicer
summaryTab <- round(xtabs(x ~ mM_HCHO + Timepoint, CFUcounts)/9)
summaryTab

#plot(round(xtabs(x ~ Timepoint, CFUcounts)/9), type = "l", ylab = "Average Population Size")

dfSummary <- as.data.frame(summaryTab)
dfSummary$mM_HCHO<- as.numeric(as.character(dfSummary$mM_HCHO))
dfSummary$Timepoint <- as.numeric(as.character(dfSummary$Timepoint))
dfSummary <- subset(dfSummary,Timepoint>0)
t0.dist <- dfSummary[dfSummary$Timepoint == 2 & dfSummary$Freq > 0,]
counts <- as.numeric(unlist(apply(t0.dist,1, function(x) rep(x[1],x[3]))))

####################### Second regrowth data 
CFUcounts_regrowth <- read.csv('formaldehyde_switching3_tolerance_summary.csv', header=T)
# taking the mean and sd of freq from three replicates 
CFUcounts_regrowth <- CFUcounts_regrowth[,c(2:4,8:9)]

DIL_K <- 1#90.90909090909 #this is the dilution constant that scanner (?) seems to have multiplied by
CFUcounts_regrowth$Colonies <- CFUcounts_regrowth$CFU_mean/DIL_K #convert back to colony counts
#drop formaldehyde levels beyond observed growth
maxF_regrowth <- max(CFUcounts_regrowth$mM_HCHO[CFUcounts_regrowth$Colonies > 0])
CFUcounts_regrowth <- CFUcounts_regrowth[CFUcounts_regrowth$mM_HCHO <= maxF_regrowth,]
#create xSet data frame for ensuring that each
#flask, [HCHO], and timepoint are present in 
#our data
tps_regrowth <- unique(CFUcounts_regrowth$Timepoint)
conc_regrowth <- seq(0,10,2)
xSet_regrowth <- data.frame(Flask_Rep = rep(LETTERS[1:3],each = length(conc_regrowth)*length(tps_regrowth),2),mM_HCHO = rep(conc_regrowth, each=length(tps_regrowth),2), Timepoint = tps_regrowth, Substrate=c(rep("methanol",length(tps_regrowth)*length(conc_regrowth)),rep("succinate",length(tps_regrowth)*length(conc_regrowth)) ))
CFUcounts_regrowth <- aggregate(CFUcounts_regrowth$Colonies, by = CFUcounts_regrowth[,c("Flask_Rep","Timepoint","mM_HCHO","Substrate")], sum)
CFUcounts2_regrowth <- merge(xSet_regrowth,CFUcounts_regrowth,all.x=T)
CFUcounts2_regrowth$x[is.na(CFUcounts2_regrowth$x)]<- 0 #force NA to 0; I think this is safe???

CFUcounts2_regrowth <- CFUcounts2_regrowth[do.call(order,CFUcounts2_regrowth[,c("Flask_Rep","Timepoint","mM_HCHO","Substrate")]),] #just make sure things are in order
#impose the assumption on the data and calculate the cells of particular tolerance level
for (S in c("methanol","succinate")){
  for(f in LETTERS[1:3]){
    for(tt in tps_regrowth) CFUcounts2_regrowth[CFUcounts2_regrowth$Flask_Rep == f & CFUcounts2_regrowth$Timepoint == tt & CFUcounts2_regrowth$Substrate==S,"x"] <- round(-1 * diff(c(rev(cummax(rev(CFUcounts2_regrowth[CFUcounts2_regrowth$Flask_Rep == f & CFUcounts2_regrowth$Timepoint == tt  & CFUcounts2_regrowth$Substrate==S,"x"]))),0)))
  }
}
#division by 9 due to 9 samples per sum (3 flasks, 3 reps)
#rounding to make the table look nicer
summaryTab_regrowth <- round(xtabs(x ~ mM_HCHO + Timepoint + Substrate, CFUcounts2_regrowth)/9)
summaryTab_regrowth

plot_cfu_data_regrowth <- ggplot() + geom_line(data=CFUcounts2_regrowth, size=1.2, 
                                               aes(x=mM_HCHO, y=x, col=factor(Timepoint), group=interaction(Timepoint,Flask_Rep))) 
plot_cfu_data_regrowth <- plot_cfu_data_regrowth + facet_grid(Substrate~Flask_Rep) + scale_y_log10(breaks=c(1e+1,1e+2,1e+3,1e+4,1e+5)) + xlab("Tolerance (mM)") + ylab("CFU") #+ scale_x_discrete(breaks=seq(0,10,2))
plot_cfu_data_regrowth <- plot_cfu_data_regrowth + scale_color_manual(values=rev(brewer.pal(length(unique(CFUcounts2$Timepoint))+2, 'YlOrRd')), name='time (hours)')
#plot_cfu_data
#plot(round(xtabs(x ~ Timepoint, CFUcounts)/9), type = "l", ylab = "Average Population Size")

dfSummary_regrowth <- as.data.frame(summaryTab_regrowth)
dfSummary_regrowth$mM_HCHO <- as.numeric(as.character(dfSummary_regrowth$mM_HCHO))
dfSummary_regrowth$Timepoint <- as.numeric(as.character(dfSummary_regrowth$Timepoint))

colnames(dfSummary_regrowth)[4] <- "cfu"
Time <- as.factor(dfSummary_regrowth$Timepoint)
plot_data_mean_regrowth <- ggplot() + geom_line(data=dfSummary_regrowth, size=1.2, 
                                                aes(x=mM_HCHO, y=cfu, col=Time, group=Timepoint)) 
plot_data_mean_regrowth <- plot_data_mean_regrowth + facet_grid(Substrate~.) + scale_y_log10(breaks=c(1e+1,1e+2,1e+3)) + xlab("Tolerance (mM)") + ylab("CFU")
plot_data_mean_regrowth <- plot_data_mean_regrowth + scale_color_manual(values=rev(brewer.pal(length(unique(dfSummary_regrowth$Timepoint))+2, 'YlOrRd')), name='time (hours)')
#plot_data_mean_regrowth

t0.dist_regrowth <- dfSummary_regrowth[dfSummary_regrowth$Timepoint == 0 & dfSummary_regrowth$cfu > 0,]
t0.dist.succinate <- subset(t0.dist_regrowth,Substrate=="succinate")
t0.dist.succinate <- t0.dist.succinate[,-3]
t0.dist.methanol <- subset(t0.dist_regrowth,Substrate=="methanol")
t0.dist.methanol <- t0.dist.methanol[,-3]
counts.succinate <- as.numeric(unlist(apply(t0.dist.succinate,1, function(x) rep(x[1],x[3]))))
counts.p.succinate <- counts.succinate + 1e-10 #add a small constant to put our data on the support for the gamma distribution
counts.methanol <- as.numeric(unlist(apply(t0.dist.methanol,1, function(x) rep(x[1],x[3]))))
counts.p.methanol <- counts.methanol + 1e-10 #add a small constant to put our data on the support for the gamma distribution
#the fit.gamma function calculates the sums of squares to compare
#the proposed gamma density to observed density
fit.gamma <- function(shape, rate, data){
  qFun <- stepfun(unique(data), c(0,cumsum(tabulate(match(data, unique(data))))/length(data)), right = T)
  qData <- qFun(data) #get quantiles for the data
  qG <- pgamma(data,shape,rate)
  sum((qData - qG)^2)
}
#minimize the sums of squares
#best.gam.succinate <- optim(c(shape = 1, rate = 1), function(x) fit.gamma(x["shape"], x["rate"], counts.p.succinate), lower = 0, method = "L-BFGS-B")
#best.gam.succinate$par
#best.gam.methanol <- optim(c(shape = 1, rate = 1), function(x) fit.gamma(x["shape"], x["rate"], counts.p.methanol), lower = 0, method = "L-BFGS-B")
#best.gam.methanol$par

#maxTol = highest tolerance category
maxTol = 10;
#n = number of classes to model 
NClasses <- 100;
#h = step size (function of maxTol and n)
h <- maxTol/NClasses;


N0.succinate <- sum(dfSummary_regrowth$cfu[dfSummary_regrowth$Timepoint==0 & dfSummary_regrowth$Substrate=="succinate"])#mean(CFUcounts$x[CFUcounts$mM_HCHO==0 & CFUcounts$Timepoint==0 & CFUcounts$Substrate=="succinate"])
N0.methanol <-  sum(dfSummary_regrowth$cfu[dfSummary_regrowth$Timepoint==0 & dfSummary_regrowth$Substrate=="methanol"])#mean(CFUcounts$x[CFUcounts$mM_HCHO==0 & CFUcounts$Timepoint==0 & CFUcounts$Substrate=="methanol"])
#use gamma distributed initial conditions on tolerance
n0.succinate <- N0.succinate*h*dgamma((0:NClasses)*h,shape=0.9301252, rate = 0.4645393)
n0.succinate[1] <- n0.succinate[2] # The first one is inf
n0.methanol <- N0.methanol*h*dgamma((0:NClasses)*h,shape = 0.9960762, rate = 0.4613910)
n0.methanol[1] <- n0.methanol[2] # The first one is inf

x <- seq(0,10,length.out=NClasses+1)

f0_regrowth <- 0 #formaldehyde at time 0
m0.methanol <- 15 #methanol at time 0
m0.succinate <- 0
s0.methanol <- 0
s0.succinate <- 3.5

state.succinate <- c(m0.succinate ,s0.succinate, f0_regrowth, n0.succinate)
state.methanol <- c(m0.methanol, s0.methanol, f0_regrowth, n0.methanol)
#initial population size

IC <- read.csv("CJM_IC.csv",header = T)
n0 <- IC$Cells
#n0[61:101] <- n0[60]
# vector of phenotypic states

f0 <- 4 #formaldehyde at time 0
m0 <- 15 #methanol at time 0
s0 <- 0

state_growth <- c(m0,s0,f0,n0)

parms = list(
  alpha = 0.18, #alpha0.145 is good if we have cond for death: x<=f
  ks = 0.02, #K[s]
  km = 0.02, #K[m]
  kf = 0.2, #K[f]
  vmax = 2.5e-8, #V[max]4e-8 
  vmaxs = 1e-8,  #V[maxs]
  gamma = 1.23, #gamma
  N = NClasses, #number of states in x
  rf = 0.18, #growth rate in formaldehyde
  rs = 0.26, # growth rate in succinate
  V1 = 2.6881, #advection constant for dn/dx 
  D = 3.591034, #diffusion constant for d^2 n/dx^2 
  V2 = 0,#1.46, # Advection to the left no matter what
  V3 = 0,#1.41, # Advection to the right if methanol is around 
  x=x) 

# model consists of two IDE and one PDE with NClasses discrete states (subpopulations)
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
    if (f<0.03){FluxP1 <- V1*diff(nn)}
    #FluxP1 <- -V*diff(nn)
    
    #FluxP3 <- V2*diff(nn)
    FluxP4 <- rep(0,length(dn))
    if(f>0.1) FluxP4 <- -V3*diff(nn)
    FluxP2 <- D * diff(c(n[1], n, n[N]), diff=2) ###note that there should technically be division by dx^2
    ###however, to recover d^2 n/dx^2, you have to then multiply by dx^2, making it pointless to do the division
    
    ## Rate of change = Flux gradient + Biology
    dn <- dn  + FluxP2 + FluxP1  + FluxP4
    
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
  out$"total_pop" <- rowSums(out[5:(NClasses+5)])
  
  return(out)
}
########################### This part matches the data to the model 
dfSummary2 <- dfSummary
times_growth=unique(dfSummary$Timepoint)
first_growth<-model_multi_levels(model=pde, time_range=times_growth, baseline_state = state_growth, parameters=parms, f_init_range=f0)

#plot(total_pop~time,first_growth,log='y')
#plot(formaldehyde~time,first_growth,log='y',ylim=c(1,4))

long_first <- melt(first_growth, id.vars = c("time","methanol","succinate","formaldehyde","total_pop","f_init"),variable.name = "phenotype",value.name ="cfu")
levels(long_first$phenotype) <- x

long_first$phenotype <- as.numeric(as.character(long_first$phenotype))
plot_first_growth <- ggplot() + geom_point(data=long_first, size=1.2, 
                                           aes(x=phenotype, y=cfu, col=factor(time), group=time))
plot_first_growth <- plot_first_growth + facet_grid(.~time) + scale_color_manual(values=rev(brewer.pal(length(unique(long_first$time))+2, 'YlOrRd')), name='time (hours)')
plot_first_growth <- plot_first_growth + scale_y_log10() + scale_x_continuous(breaks = seq(0,10))
#plot_first_growth

first_growth_binned <- first_growth[,5:(NClasses+5)] # if I want binning freq, this should be <- first_growth_freq
first_growth_binned <- as.matrix(first_growth_binned)
first_aggregate <- data.frame()
for (i in 1:length(times_growth)){
  tmp <- aggregate(first_growth_binned[i,grepl("n",colnames(first_growth_binned))], by=list("phenotype" = floor(seq(0,maxTol, len = NClasses+1))), sum)
  first_aggregate <- rbind(first_aggregate,tmp) 
}
first_aggregate$"time" <- rep(times_growth,each=11)
first_aggregate$type <- "model"

colnames(first_aggregate)[2] <- c("cfu")
colnames(dfSummary2) <- c("phenotype","time","cfu")
dfSummary2$type <- "data"
long_hist <- long_first[,c(1,7,8)]
long_hist$type <- "cont. model"
data_model <- rbind(first_aggregate,dfSummary2)

# plot_first_binned <-ggplot() + geom_point(data=data_model, size=2.2, 
#                                           aes(x=phenotype, y=cfu, col=factor(time), group=interaction(time,type))) + geom_hline(yintercept=1)
# plot_first_binned <- plot_first_binned + theme_bw(base_size = 17) + xlab('tolerance (mM)') + ylab('cfu') + scale_y_log10(limits=c(1e+1,1e+8))
# plot_first_binned <- plot_first_binned + facet_grid(type~time) + scale_color_manual(values=rev(brewer.pal(length(unique(data_model$time))+2, 'YlOrRd')), name='time (hours)')
# plot_first_binned

data_model_hist <- subset(data_model,cfu>1)
hist_first_binned <-ggplot(data=data_model_hist, size=2.2, 
                           aes(x=phenotype, y=cfu, fill=factor(time), group=interaction(time,type))) + geom_bar(stat = "identity") + scale_y_log10()  + theme(text = element_text(size=20))
hist_first_binned <- hist_first_binned + facet_grid(type~time) + ylab("viable cells (CFU/mL)") + xlab("tolerance (mM)") + guides(fill=guide_legend("time (hours)"))
hist_first_binned

long_hist <- subset(long_hist,cfu>=1)
hist_long <-ggplot(data=long_hist, size=2.2, 
                           aes(x=phenotype, y=cfu, fill=factor(time), group=interaction(time,type))) + geom_bar(stat = "identity") + scale_y_log10()  + theme(text = element_text(size=20))
hist_long <- hist_long + facet_grid(type~time) + ylab("viable cells (CFU/mL)") + xlab("tolerance (mM)") + guides(fill=guide_legend("time (hours)"))
hist_long

first_growth_freq <- first_growth[,5:(NClasses+5)]/first_growth$total_pop
first_growth_freq <- as.matrix(first_growth_freq)
first_aggregate_freq <- data.frame()
for (i in 1:length(times_growth)){
  tmp <- aggregate(first_growth_freq[i,grepl("n",colnames(first_growth_freq))], by=list("phenotype" = floor(seq(0,maxTol, len = NClasses+1))), sum)
  first_aggregate_freq <- rbind(first_aggregate_freq,tmp) 
}
first_aggregate_freq$"time" <- rep(times_growth,each=11)
first_aggregate_freq$type <- "model"
colnames(first_aggregate_freq)[2] <- "pdf"

dfSummaryFreq <- data.frame()
for (td in unique(dfSummary2$time)){
  tmp <- subset(dfSummary2,time==td)
  total_cells <- sum(tmp$cfu)
  tmp$pdf <- tmp$cfu/total_cells
  dfSummaryFreq <- rbind(dfSummaryFreq,tmp)
}
dfSummaryFreq <- dfSummaryFreq[,-3]

total_freq <- rbind(first_aggregate_freq,dfSummaryFreq)

plot_first_freq <-ggplot() + geom_line(data=total_freq, size=1.2, 
                                       aes(x=phenotype, y=pdf, col=factor(time), group=interaction(time,type))) 
plot_first_freq <- plot_first_freq + theme_bw(base_size = 17) + xlab('tolerance (mM)') + ylab('frequency') 
plot_first_freq <- plot_first_freq + facet_grid(type~time) + scale_color_manual(values=rev(brewer.pal(length(unique(total_freq$time))+2, 'YlOrRd')), name='time (hours)')
plot_first_freq
########################## This part is to check basic growth/dynamic
gg_color_hue <- function(n) { # this is a cute function I found online, for making evenly-spaced color codes
  hues = seq(15, 375, length=n+1) # n = the number of factors you're plotting
  hcl(h=hues, l=65, c=100)[1:n] 
} 
times_growth <- seq(0,100)
growth_f_range <- c(0,2,3,4)
second_growth<- model_multi_levels(model=pde, time_range=times_growth, baseline_state = c(15,0,0,388500*h*dgamma((0:NClasses)*h,13.03363, rate = 10.02192)), parameters=parms, f_init_range=growth_f_range)
treatment = as.factor(second_growth$f_init)
plot_second_growth <- ggplot() + geom_line(data=second_growth,size=1.2,aes(x=time,y=total_pop,col=treatment))  + theme(text = element_text(size=20))
plot_second_growth <- plot_second_growth + scale_y_log10() + scale_color_manual(values=gg_color_hue(length(growth_f_range)), name='formaldehyde \ntreatment (mM)') + ylab("viable cells (CFU/mL)") + xlab("time (hrs)")
plot_second_growth

plot_second_met <- ggplot(data=second_growth,aes(x=time,y=methanol,col=treatment)) + geom_line(size=1.2) 
plot_second_met <- plot_second_met #+ scale_y_log10() 
#plot_second_met

plot_second_form <- ggplot() + geom_line(data=second_growth,size=1.2,aes(x=time,y=formaldehyde,col=treatment)) 
#plot_second_form

#grid.arrange(plot_second_growth,plot_second_met,plot_second_form)

second_growth_binned <- subset(second_growth,f_init==4)
second_growth_binned <- second_growth_binned[nrow(second_growth_binned),5:(NClasses+5)] # if I want binning freq, this should be <- first_growth_freq
second_growth_binned <- as.matrix(second_growth_binned)
second_aggregate <- aggregate(second_growth_binned[1,grepl("n",colnames(second_growth_binned))], by=list("phenotype" = floor(seq(0,maxTol, len = NClasses+1))), sum)
plot(x~phenotype,second_aggregate,log='y',ylab='viable cells (CFU/mL)',xlab='tolerance (mM)',type='h',ylim=c(1,1e+8))
hist_second_binned <-ggplot(data=second_aggregate, size=2.2, aes(x=phenotype, y=x)) + geom_bar(stat = "identity") + scale_y_log10()  + theme(text = element_text(size=20))
hist_second_binned <- hist_second_binned  + ylab("viable cells (CFU/mL)") + xlab("tolerance (mM)") + guides(fill=guide_legend("time (hours)")) 
hist_second_binned

third_growth<- model_multi_levels(model=pde, time_range=times_growth, baseline_state = c(0,3.5,0,388500*h*dgamma((0:NClasses)*h,13.03363, rate = 10.02192)), parameters=parms, f_init_range=c(3,4,5))
treatment = as.factor(third_growth$f_init)
plot_third_growth <- ggplot() + geom_line(data=third_growth,size=1.2,aes(x=time,y=total_pop,col=treatment)) 
plot_third_growth <- plot_third_growth + scale_y_log10()
#plot_third_growth

death_f_range <- c(5, 7.5, 10, 12.5, 15, 20)
first_death <- model_multi_levels(model=pde, time_range=seq(0,3), baseline_state = c(15,0,0,77166668*h*dgamma((0:NClasses)*h,13.03363, rate = 10.02192)), parameters=parms, f_init_range=death_f_range)
treatment = as.factor(first_death$f_init)
plot_first_death <- ggplot() + geom_line(data=first_death,size=1.2,aes(x=time,y=total_pop,col=treatment)) + scale_y_log10(limits=c(1e+2,2e+8))
plot_first_death <- plot_first_death + theme(text = element_text(size=20)) + scale_color_manual(values=gg_color_hue(length(death_f_range)), name='formaldehyde \ntreatment (mM)') + ylab("viable cells (CFU/mL)") + xlab("time (hrs)")
plot_first_death

plot_first_death_form <- ggplot() + geom_line(data=first_death,size=1.2,aes(x=time,y=formaldehyde,col=treatment)) 
#plot_first_death_form 

#grid.arrange(plot_first_death,plot_first_death_form)

# fourth growth, stationary distribution  
n0 <- rep(0,101)
n0[80] <- 1000

f0 <- 0 #formaldehyde at time 0
m0 <- 0 #methanol at time 0
s0 <- 3.5

state_fourth_growth <- c(m0,s0,f0,n0)
times_growth <- seq(0,300)
growth_f_range <- 0
fourth_growth<- model_multi_levels(model=pde, time_range=times_growth, baseline_state = state_fourth_growth, parameters=parms, f_init_range=growth_f_range)
treatment = as.factor(fourth_growth$f_init)
plot_fourth_growth <- ggplot() + geom_line(data=fourth_growth,size=1.2,aes(x=time,y=total_pop,col=treatment))  + theme(text = element_text(size=20))
plot_fourth_growth <- plot_fourth_growth + scale_y_log10() + scale_color_manual(values=gg_color_hue(length(growth_f_range)), name='formaldehyde \ntreatment (mM)') + ylab("viable cells (CFU/mL)") + xlab("time (hrs)")
plot_fourth_growth

plot_fourth_met <- ggplot(data=fourth_growth,aes(x=time,y=methanol,col=treatment)) + geom_line(size=1.2) 
plot_fourth_met <- plot_fourth_met #+ scale_y_log10() 
#plot_fourth_met

fourth_growth_binned <- fourth_growth
fourth_growth_binned <- fourth_growth_binned[nrow(fourth_growth_binned),5:(NClasses+5)] # if I want binning freq, this should be <- first_growth_freq
fourth_growth_binned <- as.matrix(fourth_growth_binned)
fourth_aggregate <- aggregate(fourth_growth_binned[1,grepl("n",colnames(fourth_growth_binned))], by=list("phenotype" = floor(seq(0,maxTol, len = NClasses+1))), sum)
plot(x~phenotype,fourth_aggregate,ylab='viable cells (CFU/mL)',xlab='tolerance (mM)',type='h',log='y',ylim=c(1,1e+8))
hist4 <- subset(fourth_aggregate,x>1)
hist_first_binned <-ggplot(data=hist4, size=2.2, aes(x=phenotype, y=x)) + geom_bar(stat = "identity") + scale_y_log10()  + theme(text = element_text(size=20))
hist_first_binned <- hist_first_binned  + ylab("viable cells (CFU/mL)") + xlab("tolerance (mM)") + guides(fill=guide_legend("time (hours)"))
hist_first_binned

########################### Regrowth 
dfSummary_regrowth2 <- dfSummary_regrowth
times_regrowth=unique(dfSummary_regrowth2$Timepoint)

regrowth_succinate <- as.data.frame(ode(y=state.succinate, times=times_regrowth, func=pde, parms=parms, method="lsoda"))
colnames(regrowth_succinate) <- c("time","methanol","succinate","formaldehyde",paste("n",seq(0,100),sep = ""))
regrowth_succinate$"media" <- "succinate"
regrowth_methanol <- as.data.frame(ode(y=state.methanol, times=times_regrowth, func=pde, parms=parms, method="lsoda"))
colnames(regrowth_methanol) <- c("time","methanol","succinate","formaldehyde",paste("n",seq(0,100),sep = ""))
regrowth_methanol$"media" <- "methanol"
total_regrowth <- rbind(regrowth_succinate,regrowth_methanol)
long_regrowth <- melt(total_regrowth, id.vars = c("time","methanol","succinate","formaldehyde","media"),variable.name = "phenotype",value.name ="cfu")

levels(long_regrowth$phenotype) <- x

long_regrowth$phenotype <- as.numeric(as.character(long_regrowth$phenotype))
plot_first_growth <- ggplot() + geom_point(data=long_regrowth, size=1.2, 
                                           aes(x=phenotype, y=cfu, col=factor(time), group=time))
plot_first_growth <- plot_first_growth + facet_grid(media~.) + scale_color_manual(values=rev(brewer.pal(length(unique(long_regrowth$time))+2, 'YlOrRd')), name='time (hours)')
plot_first_growth <- plot_first_growth + scale_y_log10(limits=c(1,1e+8)) + scale_x_continuous(breaks = seq(0,10))
#plot_first_growth


# Aggregating succinate media 
regrowth_binned_succinate <- regrowth_succinate[,5:105] # if I want binning freq, this should be <- first_growth_freq
regrowth_binned_succinate <- as.matrix(regrowth_binned_succinate)
succinate_aggregate <- data.frame()
for (i in 1:length(times_regrowth)){
  tmp <- aggregate(regrowth_binned_succinate[i,grepl("n",colnames(regrowth_binned_succinate))], by=list("phenotype"=seq(0, maxTol, len = NClasses+1) %/% 2*2), sum)
  succinate_aggregate <- rbind(succinate_aggregate,tmp) 
}
succinate_aggregate$"time" <- rep(times_regrowth,each=6)
succinate_aggregate$media <- "succinate"

# aggregating methanol media 
regrowth_binned_methanol <- regrowth_methanol[,5:105] # if I want binning freq, this should be <- first_growth_freq
regrowth_binned_methanol <- as.matrix(regrowth_binned_methanol)
methanol_aggregate <- data.frame()
for (i in 1:length(times_regrowth)){
  tmp <- aggregate(regrowth_binned_methanol[i,grepl("n",colnames(regrowth_binned_methanol))], by=list("phenotype" = seq(0, maxTol, len = NClasses+1) %/% 2*2), sum)
  methanol_aggregate <- rbind(methanol_aggregate,tmp) 
}
methanol_aggregate$"time" <- rep(times_regrowth,each=6)
methanol_aggregate$media <- "methanol"

total_aggregate <- rbind(succinate_aggregate,methanol_aggregate)
colnames(total_aggregate)[2] <- c("cfu")
total_aggregate$type <- "model"

colnames(dfSummary_regrowth2)[1:3] <- c("phenotype","time","media")
dfSummary_regrowth2$type <- "data"
data_model_regrowth <- rbind(total_aggregate,dfSummary_regrowth2)

plot_regrowth_binned <-ggplot() + geom_line(data=data_model_regrowth, size=1.2, 
                                            aes(x=phenotype, y=cfu, col=factor(time), group=interaction(time,type))) 
plot_regrowth_binned <- plot_regrowth_binned + theme_bw(base_size = 17) + xlab('tolerance (mM)') + ylab('viable cells (CFU/mL)') + scale_y_log10() + coord_cartesian(ylim =c(1e-1,1e+9))
plot_regrowth_binned <- plot_regrowth_binned + facet_grid(type~media) + scale_color_manual(values=rev(brewer.pal(length(unique(data_model_regrowth$time))+2, 'YlOrRd')), name='time (hours)')
plot_regrowth_binned

