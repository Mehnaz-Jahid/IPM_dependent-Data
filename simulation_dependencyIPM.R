
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

### The survey data simulated was not dependent. 

### Attempt to simulate a dependent data-

p_c<- 0.8 # detection probability of camera traps # detection function from SECR as p0=p_c
p_h<- 0.6 # capture probability of the hair traps # detection function from SECR as p0= p_h
p_ch<- 0.5 # probability that the individual captured in both camera and hair traps
p_capt<- p_c+p_h-p_ch # detection prob. from SECR with p0=p_capt

y_hair<- array(, dim= c(M, J, K, T))
y_cam_indi<- array(, dim= c(M, J, K, T)) # latent camera trap captures
y_cam<- array(, dim= c(J, K, T)) # occupancy camera trap capt. history

for (i in 1:M){
  for (j in 1:J){
    for (k in 1:K){
      for (t in 1:T) {
        
        u_capt<- runif(1, 0, 1)
        u<- runif(1,0,1)
        
        if (u_capt>p_capt || u<p_h) {
          y_hair[i, j, k, t]= 0
          y_cam[j, k, t] = 1 
        } 
        
        if (u_capt>p_capt || u<p_c) {
          y_hair[i, j, k, t]= 1
          y_cam[j, k, t] = 0 
        } 
        
        if (u_capt>p_capt || u>p_h) {
          y_hair[i, j, k, t]= 1
          y_cam[j, k, t] = 1 
        } 
      }
    }
  }
}
  

# u_capt[i]<- runif(1, 0, 1)
# u<- runif(1,0,1)
# ifelse (u_capt[i]>p_capt & u<p_h, y_hair[i, j, k, t]= 0 &  )
# ifelse (u_capt[i]>p_capt & p_h<u<p_c, y_hair[i, j, k, t]= 1 & y_cam[J, K, T]= 0)
# ifelse (u_capt[i]>p_capt & u>p_c, y_hair[i, j, k, t]= 1 & y_cam[J, K, T]= 1)



# Open population spatially explicit and dependent capture-recapture and presence-absence data #

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
# p0_1 = 0.9

p_0c<- 0.8 # detection probability of camera traps at AC# detection function from SECR as p0=p_c
p_0h<- 0.6 # capture probability of the hair traps at AC# detection function from SECR as p0= p_h
p_0ch<- 0.5 # probability that the individual captured in both camera and hair traps at AC
p_0capt<- p_0c+p_0h-p_0ch # Probability that the individual captured by either hair or camera traps at AC. detection prob. from SECR with p0=p_capt

detector.xy<- matrix(, J, 2)
# Trap locations
for (i in 1:J) {
  detector.xy[i, 1] = runif(1, 150, 850) 
  detector.xy[i, 2] = runif(1, 150, 850)
}

d_squared_1 = matrix (, M, J)
#  p_1 = matrix (, M, J)
# p_1<- array(, dim= c(M, J, K))
p_c<- array(, dim= c(M, J, K))
p_h<- array(, dim= c(M, J, K))
p_ch<- array(, dim= c(M, J, K))
p_capt<- array(, dim= c(M, J, K))


#  y_1 = matrix (, M, J)
# y_1<- array(, dim= c(M, J, K, T))

y_hair<- array(, dim= c(M, J, K, T))
y_cam_indi<- array(, dim= c(M, J, K, T)) # latent camera trap captures
y_cam_sum<- array(, dim= c(J, K, T)) 
y_cam<- array(, dim= c(J, K, T)) # camera trap capt. history
# y_12<- array(, dim= c(M, J, K, T))

for (i in 1:M) {
  for(j in 1:J) {
    for(k in 1:K) {
      for (t in 1:T) {
        d_squared_1[i, j] <- (sxy[i, 1] - detector.xy[j,1])^2 +
          (sxy[i, 2] - detector.xy[j,2])^2
       # p_1[i, j, k] <- p0_1 * exp((-d_squared_1[i,j])/(2*sigma^2))
        p_c[i, j, k]<- p_0c * exp((-d_squared_1[i,j])/(2*sigma^2))
        p_h[i, j, k]<- p_0h * exp((-d_squared_1[i,j])/(2*sigma^2))
        p_ch[i, j, k]<- p_0ch * exp((-d_squared_1[i,j])/(2*sigma^2))
        p_capt[i, j, k]<- p_0capt * exp((-d_squared_1[i,j])/(2*sigma^2))
        
        u_capt<- runif(1, 0, 1)
        u<- runif(1,0,1)
        
        if (u_capt>p_capt[i, j, k]*z[i, t] || u<p_h[i, j, k]*z[i, t]) {
          y_hair[i, j, k, t]= 0
          y_cam_indi[i, j, k, t] = 1 
          y_cam_sum[j, k, t] <- sum(y_cam_indi[i, j, k, t])
          y_cam[j, k, t] = ifelse(y_cam_sum[j, k, t]>0, 1, 0)
        } 
        
        if (u_capt>p_capt[i, j, k]*z[i, t] || u<p_c[i, j, k]*z[i, t]) {
          y_hair[i, j, k, t]= 1
          y_cam[j, k, t] = 0 
          y_cam_sum[j, k, t] <- sum(y_cam_indi[i, j, k, t])
          y_cam[j, k, t] = ifelse(y_cam_sum[j, k, t]>0, 1, 0)
        } 
        
        if (u_capt>p_capt[i, j, k]*z[i, t] || u>p_h[i, j, k]*z[i, t]) {
          y_hair[i, j, k, t]= 1
          y_cam_indi[i, j, k, t] = 1 
          y_cam_sum[j, k, t] <- sum(y_cam_indi[i, j, k, t])
          y_cam[j, k, t] = ifelse(y_cam_sum[j, k, t]>0, 1, 0)
        }
  #      y_1[i, j, k, t] = rbinom(1, 1, p_1[i, j, k]* z[i, t]) 
  #      y_12[i, j, k, t] = rbinom(1, 1, p_1[i, j, k]* z[i, t]) # latent capture data
      }
    }
  }
}
   
y_hair[345, , , 1]     
y_cam[, , 2] 

