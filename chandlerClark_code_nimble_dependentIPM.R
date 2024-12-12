library(nimble)

# Define the NIMBLE model
code <- nimbleCode({
  # Priors
  psi ~ dbeta(1, 1)           # inclusion
  phi ~ dbeta(1, 1)           # survival
#  gamma ~ dunif(0, 3)         # per-capita recruitment
  p_0h ~ dbeta(1, 1)            # baseline capture probability of hair traps
  p_0c ~ dbeta(1, 1)            # baseline capture probability of camera traps
  sigmaK ~ dunif(0, 50)       # scale parameter of detection function
  delta ~ dunif(0, 1)         # recruitment prob.
#  sigma <- sigmaK * 1000      # now in meters
  
  # Loop over time steps (for population and recruitment dynamics)
  for(t in 1:T) {
    N[t] <- sum(z[1:M, t])       # Population size at time t
 #   EB[t] <- N[t] * gamma     # Expected births (actually, recruits)
 #  V[t] <- max(M - sum(a[, t]), 0.01)  # Bears available to be recruited
 #   b[t] <- min(EB[t] / V[t], 0.999)    # Probability of being recruited
  }
  
  # Loop over individuals
  for(i in 1:M) {
    z[i, 1] ~ dbern(psi)            # Initial state of each individual
    A[i, 1] <- z[i, 1]              # Recruited yet?
    s[i, 1] ~ dunif(0, x.max)  # Initial x-location
    s[i, 2] ~ dunif(0, y.max)  # Initial y-location
    
    for(j in 1:J) {
      d2[i, j] <- (s[i, 1] - detector.xy[j, 1])^2 + (s[i, 2] - detector.xy[j, 2])^2
      p_h[i, j] <- p_0h * exp(-d2[i, j] / (2 * sigmaK^2))
      p_c[i, j] <- p_0c * exp(-d2[i, j] / (2 * sigmaK^2))
    }
    
    # For subsequent time steps
    for(t in 2:T) {
    #  mu[i, t - 1] <- z[i, t - 1] * phi + (1 - a[i, t - 1]) * b[t - 1] 
      mu[i, t - 1] <- z[i, t-1]*phi+ (1- A[i, t-1])*delta # probability of
      z[i, t] ~ dbern(mu[i, (t - 1)])  # State of individual at time t
      A[i, t] <- max(z[i, 1:(t-1)])
   #   a[i, t] <- max(z[i, 1:t])      # Recruited yet?
    }
    
    # Calculate whether this bear was ever alive
    zi[i] <- (sum(z[i,1: T]) > 0)
 
  
  # Loop over years with SCR data (J x K x 2 detection years)
 
  for(j in 1:J) {
    for(k in 1:K) {
      for(t in 1:2) {
        y_hair[i, j, k, t] ~ dbern(p_h[i, j] * z[i, t])  # SCR data year 1
    #    y_hair[i, j, k, t] ~ dbern(p_h[i, j] * z[i, 2])  # SCR data year 2
 
        y2[i, j, k, t] ~ dbern(p_c[i, j] * z[i, t])  # Detection data year 1
    #    y2[i, j, k, t] ~ dbern(p_c[i, j] * z[i, 2])  # Detection data year 2
 
      }
    }
  }
  }
  
  # Set up trap count and binary outcome for detection data
  cuts[1] <- 0.9
  cuts[2] <- 999
  for(j in 1:J) {
    for(k in 1:K) {
      for(t in 1:2) {
        C[j, k, t] <- sum(y2[1:M, j, k, t])    # Trap count
        y_cam[j, k, t] ~ dinterval(C[j, k, t], cuts[1:2])  # Trap binary data
      }
    }
  }
  
  # Summing over the "zi" to count how many bears were ever alive
  Never <- sum(zi[1:M])
})

constants<- list(M= 20, J= 20, K= 8, T= 2, x.max= 1000, y.max= 1000)
# Define model with NIMBLE
model <- nimbleModel(code, 
                     data = list(y_hair= y_hair, y_cam= y_cam, detector.xy= detector.xy),
                     constants= constants,
                     inits = list(psi = 0.8, 
                                  phi = 0.5,
                              #    gamma = 0.01,
                                  p_0h = 0.6,
                                  p_0c = 0.8,
                                  delta = 0.01,
                                  sigmaK = 0.01))

sample<-nimbleMCMC(model = model,
                   monitors = c("N", "p_0c", "p_0h", "sigmaK", "phi", "psi" ),
                   niter=1000,
                   nburnin=50, 
                   thin= 2, 
                   nchains=2, 
                   progressBar = TRUE ,
                   samplesAsCodaMCMC = T)

summary(sample)

