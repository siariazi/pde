library(deSolve)
library(rstan)
library(ggplot2)
library(reshape2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Combination of v = 0.1 and d = 0 makes negative cells. 
# In this code I cannot alternate parameters cuz they are hard coded 
#### Genrating a data 
setwd("C:/Siavash/EfgA/subpop/sp_plot_data")
res.freq <- read.csv("formaldehyde_resistance_frequency.csv", header=T)
res.freq$frequency.of.resistant.cells[2] <- 1
res.freq$alt <- 1 - res.freq$frequency.of.resistant.cells
freqVec <- diff(c(0,res.freq$alt,1))
names(freqVec) <- 0:6
plot(0:6, freqVec, xlab = "Phenotype", ylab="Frequency", type="b")
library(fitdistrplus)
sampleVals <- sample(0:6,10000,replace = T,prob = freqVec)
fg<-fitdistr(sampleVals,"gamma")
rZ <- rgamma(shape = fg$estimate["shape"], rate = fg$estimate["rate"],1000)
#hist(rZ)
bins <- seq(1,5,length.out = 21)
binnedZ <- cut(rZ,bins,labels = 1.1 + (0:19)*0.2)
freqVals <- table(binnedZ)/length(rZ)
#plot(seq(0,6,length=100),1- pgamma(shape = fg$estimate["shape"], rate = fg$estimate["rate"],seq(0,6,length=100)),type="l", xlab = "Phenotype (Z)", ylab = "P(z > Z)")
#points(res.freq[,1:2],col="red")

NClasses <- 20
f0 <- 4 #formaldehyde at time 0
m0 <- 15 #methanol at time 0
n0 <- rep(0,NClasses) #starting population sizes for phenotypic classes

x <- as.numeric(names(freqVals))

for (i in seq(1,NClasses)){
  n0[i] <- freqVals[i]*1e6
}
state <- c(f0,m0,n0)

tMax <- 100 #number of time steps for the ode solution

parms = list(
  alpha = 0.28, #alpha
  km = 0.02, #K[m]
  kf = 1e-3, #K[f]1e-5
  vmax = 4e-8, #V[max]1e-9
  gamma = 1.1, #gamma
  N = NClasses, #number of states in x
  r = 0.23, #growth rate
  V= 0.1, #advection constant for dn/dx 0.001
  D= 0, #diffusion constant for d^2 n/dx^2 0.01, If I change D to 0.1 I don't get negative cells
  x=x )

#alpha needs to be a vector of length N
pde <- function (time, state, parms) {
  with (as.list(parms), {
    f <- state[1] #methanol concentration
    m <- state[2] #formaldehyde concentration
    n <- state[3:(N+2)] #N phenotypic states
    
    df <- vmax*(m/(m+km)-gamma*f/(f+kf))*sum(n) 
    dm <- -vmax*m/(m+km)*sum(n) 
    dn <- n*r*f/(f+kf) - n*alpha*(f-x)*(x < f) 
    
    # Defining the Flux P1: Chemotaxis toward the F concentration in space 
    FluxP1 = c(rep(0,NClasses))
    nn <- c(n,0)
    nn[1] <- 0
    P1 <- diff(nn)
    for (i in 1:NClasses){
      if (x[i]<f) {FluxP1[i]=-V*P1[i]}
      else if (x[i]>f) {FluxP1[i]=V*P1[i]}
      else {FluxP1[i]=0}
    }
    # if (x<f){FluxP1 <- -V*diff(nn)}
    # if (x>f){FluxP1 <- V*diff(nn)}
    # if (x==f){FluxP1 <- 0}
    
    FluxP2 <- D * diff(c(n[1], n, n[N]), diff=2) ###note that there should technically be division by dx^2
    ###however, to recover d^2 n/dx^2, you have to then multiply by dx^2, making it pointless to do the division
    
    ## Rate of change = Flux gradient + Biology
    dn <- dn  + FluxP2 + FluxP1
    
    return (list(c(df, dm, dn)))
  })
}
Z.max <- max(x)
Z.min <- min(x)

model_multi_levels<-function(model, time_range, baseline_state, parameters, f_init_range){
  out<-data.frame() # create an empty dataframe for putting all of the outputs into
  for (i in f_init_range) { # for-loop to run through all initial external formaldehyde levels
    state_i<-baseline_state # start off with your baseline state variables, of which all except f will remain the same
    state_i[1]<-i # replace the f from the baseline set of state variables with the f value you want to model
    out_i <- as.data.frame(ode(y=state_i, times=time_range, func=model, parms=parameters, method="lsoda")) # run the model
    f_init<-rep(i, length(time_range)) # create a column called "f_init" to append to the data tables, so you know which initial f value yielded these outputs
    out<-rbind(out, cbind(out_i, f_init)) # append the results of this model run, along with a column indicating the initial formaldehyde concentration, to the master output dataframe
  }
  colnames(out)[1:3] <- c("time","formaldehyde","methanol")
  colnames(out)[4:(NClasses+3)] <- paste(rep('n',NClasses), seq(1,NClasses,1), sep='') 
  colnames(out)[(NClasses+4)] <- "f_init"
  
  out$"ntot" <- rowSums(out[,4:(NClasses+3)])
  gg_names <- paste(rep('ntot',Z.max),c(seq(1,Z.max)),sep='')
  out$ntot1 <- rowSums(out[,(NClasses/Z.max+3):(NClasses+3)])
  out$ntot2 <- rowSums(out[,(2*NClasses/Z.max+3):(NClasses+3)])
  out$ntot3 <- rowSums(out[,(3*NClasses/Z.max+3):(NClasses+3)])
  out$ntot4 <- rowSums(out[,(4*NClasses/Z.max+3):(NClasses+3)])
  out$ntot5 <- as.vector(out[,NClasses+3])
  return(out)
}

times_growth=seq(0,tMax,by=1)
growth_f_range<-seq(0,4) # use one range of low fex levels for growth
growth_out<-model_multi_levels(model=pde, time_range=times_growth, baseline_state=state, parameters=parms, f_init_range=growth_f_range)


for (i in c(growth_f_range)){
  out <- subset(growth_out,f_init==i)
  out[,4:(NClasses+3)] <- abs(out[,4:(NClasses+3)])
  out <- data.matrix(out)
  BlueRed <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                                "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),bias=1)
  
  filled.contour(log10(out[,4:(NClasses+3)]), xlab = "time", ylab = "phenotype", 
                 color.palette = BlueRed, nlevels = NClasses, key.title = title(main="Popn"), 
                 plot.axes = {axis(1,at = c(0,0.2,0.4,0.6,0.8,1), labels = c(0,0.2,0.4,0.6,0.8,1)*tMax); 
                   axis(2,at = c(0,0.2,0.4,0.6,0.8,1), labels = c(0,0.2,0.4,0.6,0.8,1)*(Z.max-Z.min)+Z.min) })
}


plot_growth_model<-ggplot() + geom_line(data=growth_out, size=1.2, 
                                        aes(x=time, y=ntot, col=as.factor(f_init), group=f_init)) # I'm multiplying n by 10^5 so that the y axis looks right
plot_growth_model<-plot_growth_model + theme_bw() + xlab('time (hrs)') + ylab('viable cells (CFU/mL)')#plot_growth_model <- plot_growth_model + facet_wrap(~f_init)
plot_growth_model<-plot_growth_model + scale_y_log10(breaks=c(1e+02,1e+04,1e+06,1e+8)) 
plot_growth_model

# for now I'm not incorporating different cumulative subtypes 
pde_stan <- '
functions {

real[] pde(real t, real[] y, real[] eta,real[] x_r,int[] x_i) {
    //intialize parameters
    int NClasses; // number of parameters 
    real zmin; // minimum phenotype 
    real zmax; // maximum phenotype 
    real x[NClasses]; // vector keeping all the phenotypic states 
    real n[NClasses]; // vector keeping all the subtypes for differencing 
    real nn[NClasses+1]; // vector of n with nn[1] = 0 and nn[Nclasses+1]=0
    real P1[NClasses];   // advection flux wihout sign 
    real FluxP1[NClasses]; // Advection with sign 
    real tmp[NClasses+2];  // a tmp vector for calculating diffusion c(n[1], n, n[NClasses]) in R
    real tmp2[NClasses+1];  // a vector to keep diff of tmp (first derivatives)
    real P2[NClasses];     // vector P2 contians the second derivatives, (diff of tmp2)
    real FluxP2[NClasses]; // diffusion flux with its constant
    real gamma;
    real alpha;
    //vector[NClasses+2] dydt;
    real dydt[NClasses+2];
    real totpop;

    NClasses = 20;
    zmin = 1.1; 
    zmax = 4.9;  

    x[1] = zmin;
    for (i in 2:NClasses) x[i] = x[i-1]+(zmax-zmin)/(NClasses-1);
    gamma = eta[1];
    alpha = eta[2];

    totpop = sum(y[3:(NClasses+2)]);
    //calculate the ode values, the equations are same as model in the r code 
    dydt[1] = (4e-8)*(y[2]/(y[2]+0.02)-gamma*y[1]/(y[1]+(1e-3)))*totpop; // df
    dydt[2] = -(4e-8)*y[2]/(y[2]+0.02)*totpop; // dm 
    for (i in 3:(NClasses+2)){                 // dn  
        for (j in 1:NClasses){dydt[i]=y[i]*0.23*y[1]/(y[1]+(1e-3)) - 
              y[i]*alpha*(y[1]-x[j])*(x[j] < y[1]);}} 

    
    for (i in 3:(NClasses+2)) n[i]=y[i];
    for (i in 1:NClasses){nn[i] = n[i];}
    nn[NClasses+1] = 0;
    nn[1] = 0;
    for (i in 1:(NClasses))P1[i] = nn[i+1]-nn[i];
    for (i in 1:NClasses){
    if (x[i]<y[1]){FluxP1[i] = -0.1*P1[i];}
    else if (x[i]>y[1]){FluxP1[i] = 0.1*P1[i];}
    else {FluxP1[i] = 0;}
    }

    tmp[1] = n[1];
    for (i in 2:(NClasses+1)){tmp[i]=n[i-1];}
    tmp[NClasses+2] = n[NClasses];    
    for (i in 1:(NClasses+1)){tmp2[i] = tmp[i+1]-tmp[i];}
    for (i in 1:NClasses){P2[i] = tmp2[i+1]-tmp2[i];}
    for (k in 1:NClasses) FluxP2[k] = 0.1*P2[k]; 

    // Rate of change = Flux gradient + Biology
    for (i in 3:(NClasses+2)){dydt[i] = dydt[i] + FluxP1[i] + FluxP2[i];}
    return dydt;
  } //end pde
} //end functions

data {
  int<lower=1> T; // time 
  int<lower=1> C; //number of f initial conditions
  real f0[C]; // pass a vector of length C initial fex values
  real y0[(1+20)];     // initial values of state variables; f_0 in f0 now
  real t0;        // initial time
  real ts[T];     // real array that has n time steps
  int y[T,C];       // getting the data, here I pass cell count 
  real rel_tol;
  real abs_tol;
  int max_steps;
}

transformed data {
  real x_r[0];      
  int x_i[0];
}

parameters { 
  real<lower=0,upper=10> gamma; 
  real<lower=0,upper=10> alpha;
} 

model {
  real y_hat[T,20+2];
  real eta[2];
  real s0[2+20]; //this will now have the initial conditions

  // priors, I want to fit gamma and alpha. so they have a high sd, the rest supposed to be constant 
  //gamma ~ lognormal(10,5); 
  //alpha ~ lognormal(10,5); 

  eta[1] = gamma;
  eta[2] = alpha;

  //initial conditions for m, fin, and n are not changing
  //so we define them outside of the loop
  s0[2] = y0[1];
  for (i in 3:(20+2)) s0[i] = y0[i-1];

  for (c in 1:C){
  //at the beginning of each c we create the correct fex initial condition
  s0[1] = f0[c];
  y_hat = integrate_ode_bdf(pde, s0, t0, ts, eta, x_r, x_i, rel_tol, abs_tol, max_steps); //note we pass s0 now, I can enter the values of rel_tol... directly here 
  for (i in 3:(20+2)){y[,c] ~ poisson(y_hat[,i]);}
   //note that I vectorized sampling and made sure to get column c of data
  }
}
'
# compiling the model 
pde_stan.comp<-stan_model(model_code=pde_stan)

############ now Stan
new_data <- subset(growth_out,select = c(time,f_init,ntot))
new_data2 <- dcast(new_data,time~f_init,value.var = "ntot")
y_data <- subset(new_data2,select = -c(time))
y_data <- as.matrix(y_data)
y_data <- round(y_data)
y_data <- y_data[2:101,]

# data: y says a matrix of 2:71 not counting 0, the last column which is N data 
pde_data <- list(T = tMax, C=length(growth_f_range), f0 = growth_f_range, y0 = state[2:22], t0 = 0, ts = 1:tMax, y = y_data, rel_tol = 1e-6, abs_tol = 1e-8,max_steps = 1e6)

# sampling 
fit <- sampling(pde_stan.comp, pde_data, chains=1, iter=100, verbose = T, control = list(adapt_delta = 0.99))
pairs(fit)

