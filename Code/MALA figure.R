#This script is from (Michelot, 2019). The matérn covariates have been 
#changed to Perlin noise covariates, and lapply is used instead of mclapply,
#the plot has also been changed

library(Rhabit)

# Load packages
library(raster)
library(ggplot2)
library(viridis)
library(parallel)
library(ambient)
set.seed(1)

#######################
## Define covariates ##
#######################
# Generate two random covariates
lim <- c(-1, 1, -1, 1)*100
resol <- 1
ncov <- 2
covlist <- list()
xgrid <- seq(lim[1], lim[2], by = resol)
ygrid <- seq(lim[3], lim[4], by = resol)
coords <- as.matrix(expand.grid(xgrid, ygrid))
for(i in 1:ncov) {
  vals = 3*noise_perlin(c(length(xgrid), length(ygrid)), frequency = 0.05)
  covlist[[i]] = list(x = xgrid, y = ygrid, z = matrix(vals, nrow = length(xgrid)))
}

# Include squared distance to centre of map as covariate
xgrid <- seq(lim[1], lim[2], by=resol)
ygrid <- seq(lim[3], lim[4], by=resol)
xygrid <- expand.grid(xgrid,ygrid)
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/100
covlist[[3]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))

# Compute utilisation distribution
beta <- c(4,2,-0.1)
UD <- getUD(covariates=covlist, beta=beta)
UDrast <- rasterFromXYZ(rasterToGGplot(UD))



###################
## Simulate data ##
###################
nobs <- 1e3
alldt <- c(5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01)
ntrack <- 10
speed <- 5

rates <- rep(NA, length(alldt))
t0 <- Sys.time()
for(iter in 1:length(alldt)) {
    cat("Iteration",iter,"-- dt =",alldt[iter],"\n")
    
    dt <- alldt[iter]
    time <- (1:nobs)*dt
    alltimes <- list()
    for(i in 1:ntrack)
        alltimes[[i]] <- time
    
    # Generate tracks
    alldat <- lapply(alltimes, function(time) {
        Rhabit:::simMALA(beta = beta, gamma2 = speed, times = time, 
                loc0 = runif(2, -50, 50), cov_list = covlist)
    })
    
    rates[iter] <- mean(sapply(alldat, function(dat) dat$acc/(dat$acc+dat$rej)))
    cat("Acceptance rate:",rates[iter],"\n")
    print(Sys.time() - t0)
}

###########################
## Plot acceptance rates ##
###########################
ggplot() +
  geom_path(aes(x=alldt, y=rates)) +
  geom_point(aes(x=alldt, y=rates)) +
  scale_x_continuous(trans = "log", breaks=alldt) +
  labs(x = "Time interval", y = "Acceptance rate")
  

