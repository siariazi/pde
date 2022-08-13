library(ReacTran)
library(reshape2)
library(animation)
library(ggplot2)
setwd("C:/Siavash/Codes/pde/pde_plots_data")
#############
# number of classes 
NClasses <- 100
N <- NClasses+1
# setting up the grid
GR <- setup.grid.1D(x.up =-0.05, x.down = 10.05 , N=N) #NB: to get a mid point sequence that matches seq(0,10,by=0.1), you code it like this
x <- round(GR$x.mid,1)

# parameters in the model 
fixedparms = list(
  alpha = 0.189, #alpha0.145 is good if we have cond for death: x<=f
  f= 0,
  m = 0,
  s = 0,
  rm = 0.195, #growth rate in formaldehyde
  rs = 0.267, # growth rate in succinate
  x = x,
  D = 0,#0.02,#abs(0.2),
  V = -0.1
) 

# pde function for diffusion and advection 
pdefn <- function(time, state, parms){
  with (as.list(parms), {
    tran <- tran.1D(C=state, flux.up = 0, flux.down = 0, D = D, v = V, dx = GR)
    Death <- state*alpha*f*(x<f)
    Growth <- state*(rm*(m != 0) + rs*(s != 0))
    return(list(dn = tran$dC + Growth - Death))
  })
}

# pde function for advection only 
pdefn <- function(time, state, parms){
  with (as.list(parms), {
    tran <- advection.1D(C=state, flux.up = 0, flux.down = 0, v = V, dx = GR)
    Death <- state*alpha*f*(x<f)
    Growth <- state*(rm*(m != 0) + rs*(s != 0))
    return(list(dn = tran$dC + Growth - Death))
  })
}

# set a time interval for simulation 
ts <- seq(0,20.5,0.2)

# making a data frame for initial condition of the model 
IC <- data.frame(phenotype=seq(0,10,0.1),cells=rep(0,NClasses+1))

# assigning a gaussian distn. 
IC$cells <- 1e+4*dnorm(seq(0,10,0.1), mean = 5, sd = 0.1)
#IC$cells[41:81] <- exp(-2.505*IC$phenotype+13.285)[1:41] # form pde_reactran_181028.R

# plot the intial codition 
plot(cells~phenotype,IC,log='y')

# simulate the results 
out <- as.data.frame(ode.1D(y=IC$cells, times=ts, func=pdefn, parms=fixedparms, dimens = length(fixedparms$x)))
colnames(out) <- c("time",paste("n",seq(0,NClasses),sep="")) # Non-cumulative populations 
long_out <- melt(out, id.vars = "time",variable.name = "phenotype",value.name ="cfu")
levels(long_out$phenotype) <- x
long_out$phenotype <- as.numeric(as.character(long_out$phenotype))

# this function makes color scale for the code below 
gg_color_hue <- function(n) { 
  hues = seq(15, 375, length=n+1) # n = the number of factors you're plotting
  hcl(h=hues, l=65, c=100)[1:n] 
} 
# plotting smooth version, don't run this when there are too many time points 
# dim: 1000-351
hist_long <- ggplot(data=long_out,aes(x=phenotype,y=cfu,group=time,color=as.factor(time))) + geom_line(size=1.8) 
hist_long <- hist_long + facet_grid(.~time) + scale_y_log10(limits=c(1e-1,1e+3)) + theme_bw(base_size = 20)
hist_long <- hist_long + scale_x_continuous(breaks=seq(0,10,2)) + ylab("viability (CFU/mL)") + xlab('tolerance (mM)')
hist_long <- hist_long + scale_color_manual(values=gg_color_hue(length(ts)), name='time (hrs)')
#hist_long

# animation
ani.options(interval=0.07,ani.width=600, ani.height=600)
saveVideo(expr = {for (i in unique(long_out$time)) {
  out_t <- subset(long_out,time==i)
  # ploting regrowth
  par(pch=16,cex=1.5,cex.axis=1.3,cex.lab=1.3,cex.main=1.3,lwd=2)
  # ploting cfu~phenotype
  #plot(cfu~phenotype,out_t,log='y',type='l',ylim=c(1e-6,1e+8),ylab = "viability (CFU/mL)",xlab="tolerance (mM)")
  plot(cfu~phenotype,out_t,log='y',type='l',ylim=c(1e-1,1e+5),ylab = "viability (CFU/mL)",xlab="tolerance (mM)")
  #title(main = paste("Time: ",i))
  ani.pause()}}, video.name = "pde_reactran_adv.mp4", img.name = "Rplot",
  ffmpeg = ani.options("ffmpeg"))

