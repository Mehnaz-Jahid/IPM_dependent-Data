
dir <- "C:/Mehnaz/Uvic pc/UVic/CRT/dependency/R"
setwd(dir)
set.seed(100)
# -----------### SECR simulation ###---------------

### State model ###

M<- 1000
J<- 50
K<- 8

x.max<- 1000
y.max<- 1000

sxy<- matrix(, M, 2)

# Model_secr <- nimbleCode({
  ##-----------------------------------------------------------------
 
    ## Activity centre locations
    for (i in 1:M) {
      sxy[i, 1] = runif(1, 0, x.max) 
      sxy[i, 2] = runif(1, 0, y.max)
    }
  ## INDIVIDUAL INCLUSION
  psi = 0.5
  z<- c(NA, M)
  for (i in 1:M) {
    z[i] = rbinom(1, n= 1, prob= psi) 
  }
  N <- sum(z[1:M]) 
  
  ### Observation model ###
  
  sigma = 50
  p0_1 = 0.9
  
  detector.xy<- matrix(, J, 2)
  # Trap locations
  for (i in 1:J) {
    detector.xy[i, 1] = runif(1, 150, 850) 
    detector.xy[i, 2] = runif(1, 150, 850)
  }
  
  d_squared_1 = matrix (, M, J)
#  p_1 = matrix (, M, J)
   p_1<- array(, dim= c(M, J, K))
   
#  y_1 = matrix (, M, J)
  y_1<- array(, dim= c(M, J, K))
  
  for (i in 1:M) {
    for(j in 1:J){
      for(k in 1:K){
    d_squared_1[i, j] <- (sxy[i, 1] - detector.xy[j,1])^2 +
      (sxy[i, 2] - detector.xy[j,2])^2
    p_1[i, j, k] <- p0_1 * exp((-d_squared_1[i,j])/(2*sigma^2))

y_1[i, j, k] = rbinom(1, 1, p_1[i, j, k]* z[i]) 
      }
    }
  }

y_1[1,30,]


# -----------### Open population SECR simulation ###---------------

M<- 1000 # no. of augmented individuals
J<- 50   # no. of traps 
K<- 8    # no. of occasions (secondary sampling occasions)
T<- 2    # no. of years (primary sampling occasions)

x.max<- 1000
y.max<- 1000

sxy<- matrix(, M, 2)

# Model_secr <- nimbleCode({
##-----------------------------------------------------------------
#------------------------------
  ## Activity centre locations
  for (i in 1:M) {
    sxy[i, 1] = runif(1, 0, x.max) 
    sxy[i, 2] = runif(1, 0, y.max)
  }
## INDIVIDUAL INCLUSION
psi = 0.5 
phi = 0.8 # annual survival probability
# gamma = 0.1 # per-capita recruitment rate
delta = 0.2
z<- matrix( , M, T)
N<- c(NA, T)

for (i in 1:M) {
  z[i,1] = rbinom(1, n= 1, prob= psi) 
for (t in 2:T) {

  A = max(z[i, 1:t-1])
  p = z[i, t-1]*phi+ A*delta
  z[i, t] = rbinom(1, n= 1, p)
 }

}

N = apply(z, 2, sum)
### Observation model ###

sigma = 50
p0_1 = 0.9

detector.xy<- matrix(, J, 2)
# Trap locations
for (i in 1:J) {
  detector.xy[i, 1] = runif(1, 150, 850) 
  detector.xy[i, 2] = runif(1, 150, 850)
}

d_squared_1 = matrix (, M, J)
#  p_1 = matrix (, M, J)
p_1<- array(, dim= c(M, J, K))

#  y_1 = matrix (, M, J)
y_1<- array(, dim= c(M, J, K, T))
y_12<- array(, dim= c(M, J, K, T))

for (i in 1:M) {
  for(j in 1:J){
    for(k in 1:K){
      for (t in 1:T){
      d_squared_1[i, j] <- (sxy[i, 1] - detector.xy[j,1])^2 +
        (sxy[i, 2] - detector.xy[j,2])^2
      p_1[i, j, k] <- p0_1 * exp((-d_squared_1[i,j])/(2*sigma^2))
      
      y_1[i, j, k, t] = rbinom(1, 1, p_1[i, j, k]* z[i, t]) 
      y_12[i, j, k, t] = rbinom(1, 1, p_1[i, j, k]* z[i, t]) # latent capture data
      }
    }
  }
}
 y_1[31:50, 3, , 1]
 View(y_1[, 21, , 1])
 View(y_12[, 3, , 1])
z[10,] 

##### Survey data ####

o<- array(, dim= c(J, K, T))
y_2<- array(, dim= c(J, K, T))


for (j in 1:J){
  for (k in 1:K){
    for (t in 1:T){
      o[j, k, t] = sum(y_12[, j, k, t])
      y_2[j, k, t] = ifelse(o[j, k, t]>0, 1, 0)
    }
  }
}

View(y_2[,,1])
View(y_1[, 1, , 1])
View (z[,])


z= as.data.frame(z)
z[z$V1 == 0 & z$V2 == 0,]

library(dplyr)
dplyr::filter(z, V1 %in% 1, V2 %in% 0)

# Case 1: Not captured in hair snares, but detected in camera traps.

# Case 2: Captured in hair snares, not detected in camera traps.

# Case 3: Captured in both hair snares and camera traps. This can be attained by 
# using y_1 to generate the survey data.





