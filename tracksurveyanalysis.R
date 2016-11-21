#Track survey N-mixture model
#Author: Jessie Golding
#11/21/2016

### R code for to simulate data for MDAM development ###

#Create function to simulate carnivore track count survey data for multiple species 
#with both track and genetic id (not every track has a genetic id but all genetic IDs associated with tracks)

sim.fun <- function(n.sites){
  
  ##  Setup the logistics of sampling
  
  #  Number of sites
  n.sites <<- n.sites
  #  Number of visits to each site
  n.reps <- 3
  #  Number of sp
  n.sp <- 4
  #  Number of observations
  n.obs <- n.sites * n.reps *n.sp
  
  
  #  Indices for long format

  
  #Generate site info
  site <- rep(1:n.sites, each = n.reps*n.sp)
  
  #Generate survey replicate information
  reps <- rep(rep(1:n.reps, n.sites), n.sp)
  
  #Generate species information
  sp <-rep(1:n.sp, each = n.sites*n.reps)
  
  
  #  Detection probability of primary observer 
  P <- vector("numeric")

  P[1] <- 0.5
  #  Detection probability of secondary observer
  P[2] <- 0.3
  
  #  Sum of p's should be less than 1, where the remainder represents the
  #  proportion of the sampled population not observed
  stopifnot((P[1] + P[2]) < 1)
  cat("\nProbability of not observing carnivores", 1 - (P[1] + (P[2] * (1 - P[1]))),
      "\n\n")
  
  ##  Biological Parameters
  
  #  Mean abundance across sites, one for each species 
  #  The numbers are meant to be very different so we can see how the model handles them
  lambda <- c(20, 50, 10, 35)
  
  #  Proportion of the population captured at each session
  #p.cap <- (P[1]*(1-P[2])) + (P[2] * P[1])
 
  #  Proportion of population not captured at each session
  #p.nocap <- 1 - p.cap

	# probability you observe something on a track survey
	p.cap <- 1 - ((1-P[1])*(1-P[2]))
  
  ##  Simulation proper
  
  #  Initialize matrices to hold values of abundance corrected for availability
  N <- array(NA, dim = c(n.sites, n.reps, n.sp))
  
  # Initialize matrices to hold values of observations and probability of detection
  #Columns are outcomes of the multinomial
  y <- cp <- matrix(NA, nrow = n.obs, ncol = 2)

  
  #  Initialize matrix to hold values of true abundance
  M <- matrix(NA, n.sites, n.sp)
  for(i in 1:n.sites){
    M[i,] <- rpois(n.sp, lambda)
  }
  
  #  Abundance corrected for availability during each survey rep
  for(i in 1:n.sites){
    for(j in 1:n.reps){
      for(k in 1:n.sp){
        N[i,j,k] <- rbinom(1, M[i,k], p.cap)
      }
    }
  }
  
  # Number observed
  for(i in 1:n.obs){
    cp[i,] <- c((P[1]*(1-P[2])), (P[1]*P[2]))
    y[i,1] <- rbinom(1, N[site[i], reps[i], sp[i]], cp[i,1])
    y[i,2] <- rbinom(1, y[i,1], cp[i,2])	
  }
  
  #  Put the data together in long format
  input <- data.frame(cbind(y[,1:2], site, reps, sp))
  colnames(input)[1:2] <- c("y1", "y2")
  
  ################################################################################
  #  JAGS model to estimate parameters
  sink("model_multinomial_multisp_sim_track.txt")
  cat("
      model{
      #  Priors
      #  Linear predictor on abundance, setup for species variation only,
      #  abundance assumed the same at every site
      for(i in 1:n.sp){
      log.n[i] ~ dnorm(0, 0.001)
      mu.lambda[i] <- exp(log.n[i])
      }
      #  Population size of each species at each site
      for(i in 1:n.sites){
      for(k in 1:n.sp){
      N[i,k] ~ dpois(mu.lambda[k])
      }
      }
      #  Individual observer detection probability, no variation
      for(i in 1:2){
      p[i] ~ dbeta(1, 1)
      }
      
      #  Likelihood
      for(i in 1:n.obs){
      #  Indices always follow site, reps, species order
      #  Capture probabilities
      #  probability of detecting a track and no genetics
      cp[i,1] <- p[1]*(1-p[2])
      #  probability of detecting a track and genetics
      cp[i,2] <- p[1] * p[2]
      #  Seen by some method
      pcap[i] <- 1-((1-p[1])*(1-p[2]))
  
      #  Adust the prob of capture to the prop available
      #  2 is for number of outcomes (probablities for track and genetics)
      #  Realizations
      #  Number captured (ncap) and population size (N)
      ncap[i] ~ dbin(pcap[i], round(N[site[i],sp[i]]))
      y[i,1] ~ dbin(cp[i,1], ncap[i])
      y[i,2] ~ dbin(cp[i,1], y[i,1])
      
      }
      }
      ", fill = T)
  sink()
  ################################################################################
  #Format JAGS data

  data <- list("y" = input[,1:2],
               "n.obs" = nrow(input),
               "n.sites" = length(unique(input$site)),
               "site" = input$site,
               "n.sp"=length(unique(input$sp)),
               "sp"=input$sp)
  
  #R2jags requires the data is in the global environment. Because this is in a 
  #function need to write it to the global environment each time.
  list2env(data, envir=globalenv())
  
  require(R2jags)
  inits <- function(){list(
    log.n = log(lambda),
    p = c(0.3, 0.5),
    N = M*2 )}
  parms <- c("p", "N", "mu.lambda")
  out <- jags.parallel(data=names(data), inits, parms, "model_multinomial_multisp_sim_track.txt", 3, 30000, 1000, 1)
  
  summ <- list("P" = round(cbind(P, out$BUGS$mean$p, 100 * abs(P - out$BUGS$mean$p)/P), 2),
               "N" = round(cbind(M, out$BUGS$mean$N, 100 * abs(M - out$BUGS$mean$N)/M), 2))
  
  coverage<-list("Pcov" = ifelse(P>(quantile(out$BUGS$sims.list$p,.025)) & P<(quantile(out$BUGS$sims.list$p,.975)), 1, 0),
                 "Ncov" = ifelse(M>(quantile(out$BUGS$sims.list$N,.025)) & M<(quantile(out$BUGS$sims.list$N,.975)), 1, 0))
  
  data<-list(summ,coverage)
  
}


#Create function to specify the number of times and for how many sites the sim.fun 
#should run 

sim.fun.rep<-function(n.times, n.sites){
  replicate(n.times, sim.fun(n.sites), simplify = F)
}

sim.out<-sim.fun.rep(1,20)