rm(list = ls())

library(abind)

dir<- "C:\\Users\\mehnazjahid\\Desktop\\UVic\\CRT\\dependency\\R"
#setwd(dir)

# parameters
M <- 150

# how to read an array
ar1<- c(1,2,3,4,5)
ar2<- c("a","b","c")
array <- array(c(ar1, ar2), dim = c(3, 4, 2,2)) # default: bycolumn=T. One vector of objects (say, only ar1) will do.
array
array[1,4,1,1]

#making a fake C-R data
capt<- rep(c(0,0,0,0,0,1,1,0,0,0,0,1,0,0), 100)
y1a<- array(capt, dim=c(72, 50, 8, 2))
class(y1a)
y1a[1,1,1,1]

y1b <- array(0, dim = (c(M -72, 50, 8, 2)))

y1 <- abind(y1a, y1b, along = 1)


# Fake CT detection data
det<- rep(c(1,1,1,1,0,1,1), 80)
O<- array(det, dim=c(50, 8, 2))

#y1<- read.csv("binary_bearGenetic.csv", header=T)

constants<- list(M= M, K= dim(y1)[3], J= dim(y1)[2], T=2)
library(nimble)

code <- nimbleCode({
  psi ~ dbeta(1,1) # M*psi = E[N(1)]
  phi ~ dbeta(1,1) # survival ### MJ- do we need this? We only have one year of data.
  gamma ~ dunif(0, 3) # per-capita recruitment
  p0 ~ dbeta(1,1) # baseline capture probability
  p_cam ~ dbeta(1,1)
  delta ~ dbeta(1,1)
  sigmaK ~ dunif(0, 50) ### MJ- should it be only sigma?
  for(t in 1:T) { #(T-1)) { ### MJ- We have two years of data 
    N[t] <- sum(z[1:M,t]) # Population size at time t
    EB[t] <- N[t]*gamma # Expected births (actually, recruits)
    V[t] <- max(M - sum(a[1:M,t]), 0.01) # Bears available to be recruited
  #  b[t] <- min(EB[t] / V[t], 0.999) # Probability of being recruited #MJ- delta in paper.
  }
  for(i in 1:M) {
    z[i,1] ~ dbern(psi)
    a[i,1] <- z[i,1] # recruited yet?
    s[i,1] ~ dunif(xlims[1], xlims[2])
    s[i,2] ~ dunif(ylims[1], ylims[2])
    for(j in 1:J) {
      d2[i,j] <- (s[i,1]-X[j,1])^2 + (s[i,2]-X[j,2])^2
      p[i,j] <- p0*exp(-d2[i,j] / (2*sigma^2)) # MJ- where are we getting the sigma from?
    }
    for(t in 2:T) {
     # mu[i,t-1] <- z[i,t-1]*phi + (1 - a[i,t-1])*b[t-1]
      mu[i,t-1] <- z[i,t-1]*phi + (1 - a[i,t-1])*delta
      z[i,t] ~ dbern(mu[i,t-1])
      a[i,t] <- max(z[i,1:t]) # recruited yet? Once z(i,t)=1, then a(i,t:T)=1
    }
    for(j in 1:J) {
      for(k in 1:K) {
        # The years with SCR data
        y1[i,j,k,1] ~ dbern(p[i,j]*z[i,1])
        y1[i,j,k,2] ~ dbern(p[i,j]*z[i,2])
     #   y1[i,j,k,3] ~ dbern(p[i,j]*z[i,5])
        
        # The years with detection data
        y2[i,j,k,1] ~ dbern(p_cam*z[i,1])
        y2[i,j,k,2] ~ dbern(p_cam*z[i,2])
    #    y2[i,j,k,3] ~ dbern(p_cam*z[i,6])
      }
    }
    zi[i] <- (sum(z[i,1:T]) > 0) # Was this bear ever alive?
  }
  # JAGS trick to enforce Eq.6
  #cuts[1] <- 0.9 #1 # If C<0.9, then O=0
  #cuts[2] <- 999 #2 # If C>0.9 & C<999,, then O=1
  for(j in 1:J) {
    for(k in 1:K) {
      for(t in 1:2) {
        C[j,k,t] <- sum(y2[1:M,j,k,t]) # Trap count
        #O[j,k,t] ~ dinterval(C[j,k,t], cuts) # Trap binary data
        S[j,k,t] <- nimStep(C[j,k,t] - .5)
        O[j,k,t] ~ dbern(S[j,k,t])
      }
    }
  }
  Never <- sum(zi[1:M]) # Bears ever alive

})

model <- nimbleModel(code, 
                     data = list(y1 = y1, O = O),
                     constants= constants,
                     inits = list(psi = 0.01, 
                                  phi = 0.01,
                                  gamma = 0.01,
                                  p0 = 0.01,
                                  p_cam = 0.01,
                                  delta = 0.01,
                                  sigmaK = 0.01))

niter <- 100
nburnin <- 50
thin <- 1

sample<-nimbleMCMC(model = model,
                   monitors = c("N"),
                   niter=niter,
                   nburnin=nburnin, 
                   thin= thin, 
                   nchains=1, 
                   progressBar = TRUE ,
                   samplesAsCodaMCMC = T)                                  
summary(sample)
