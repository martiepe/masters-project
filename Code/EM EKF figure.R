library(Rhabit)
library(raster)
library(ggplot2)
library(viridis)
library(reshape2)
library(gridExtra)
library(mvtnorm)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(ambient)
library(mvnfast)
library(optimParallel)
setSessionTimeLimit(123)


#perlin covariates
lim <- c(-1, 1, -1, 1)*150
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
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/(100)
covlist[[3]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))





#compute hessian using bicubic interpolation
hessian <- function(z, covlist, par){
  n = floor(z[1])
  m = floor(z[2])
  x = z[1]
  y = z[2]
  
  Delta = 1
  #lower left corner of square being interpolated upon
  
  f = covlist[[1]]$z[(n+100):(n+103), (m+100):(m+103)]*par[1] + 
    covlist[[2]]$z[(n+100):(n+103), (m+100):(m+103)]*par[2] + 
    covlist[[3]]$z[(n+100):(n+103), (m+100):(m+103)]*par[3]
  
  
  F11 = matrix(c(f[2,2], f[2,3], f[3,2], f[3,3]), nrow = 2, byrow = T)
  
  
  F21 = matrix(c((f[3, 2] - f[1, 2])/(2*Delta),
                 (f[3, 3] - f[1, 3])/(2*Delta),
                 (f[4, 2] - f[2, 2])/(2*Delta),
                 (f[4, 3] - f[2, 3])/(2*Delta)), nrow = 2, byrow = T)
  
  
  F12 = matrix(c((f[2, 3] - f[2, 1])/(2*Delta),
                 (f[2, 4] - f[2, 2])/(2*Delta),
                 (f[3, 3] - f[3, 1])/(2*Delta),
                 (f[3, 4] - f[3, 2])/(2*Delta)), nrow = 2, byrow = T)
  
  
  F22 = matrix(c((f[3,3] - f[3, 1] - f[1, 3] + f[1, 1])/(4*Delta^2), 
                 (f[3,4] - f[3, 2] - f[1, 4] + f[1, 2])/(4*Delta^2),
                 (f[4,3] - f[4, 1] - f[2, 3] + f[2, 1])/(4*Delta^2),
                 (f[4,4] - f[4, 2] - f[2, 4] + f[2, 2])/(4*Delta^2)), nrow = 2, byrow = T)
  
  
  F = cbind(rbind(F11, F21), rbind(F12, F22))
  
  K = matrix(c(1, 0, 0, 0, 0, 0, 1, 0, -3, 3, -2, -1, 2, -2, 1, 1), nrow = 4)
  
  A = t(K) %*% F %*% K
  
  
  #the zero point of the polynomial is in our case the index of the bottom left vertex
  Hxx = c(0, 0, 2, 6*(x-n)) %*% A %*% c(1, (y-m) , (y-m)^2, (y-m)^3)
  
  Hyy = c(1, (x-n), (x-n)^2, (x-n)^3) %*% A %*% c(0, 0, 2, 6*(y-m))
  
  Hxy = c(0, 1, 2*(x-n), 3*(x-n)^2) %*% A %*% c(0, 1, 2*(y-m), 3*(y-m)^2)
  
  H = cbind(rbind(Hxx, Hxy), rbind(Hxy, Hyy))
  
  return(H)
}



#find gradient of covariates at matrix of locations using bilinear interpolation
bilinearGradVec <- function(loc_mat, cov_list) {
  x_grid <- cov_list[[1]]$x
  y_grid <- cov_list[[1]]$y
  n_cov <- length(cov_list)
  n_obs <- nrow(loc_mat)
  
  ix <- findInterval(loc_mat[,1], x_grid)
  iy <- findInterval(loc_mat[,2], y_grid)
  
  valid <- ix > 0 & ix < length(x_grid) & iy > 0 & iy < length(y_grid)
  
  x1 <- x_grid[ix]
  x2 <- x_grid[ix + 1]
  y1 <- y_grid[iy]
  y2 <- y_grid[iy + 1]
  
  dx <- x2 - x1
  dy <- y2 - y1
  lx <- loc_mat[,1]
  ly <- loc_mat[,2]
  
  grad_array <- array(NA_real_, dim = c(n_cov, n_obs, 2))
  
  for (j in seq_len(n_cov)) {
    f11 <- cov_list[[j]]$z[cbind(ix,     iy)]
    f21 <- cov_list[[j]]$z[cbind(ix + 1, iy)]
    f12 <- cov_list[[j]]$z[cbind(ix,     iy + 1)]
    f22 <- cov_list[[j]]$z[cbind(ix + 1, iy + 1)]
    
    dfdx <- ((y2 - ly) * (f21 - f11) + (ly - y1) * (f22 - f12)) / (dy * dx)
    dfdy <- ((x2 - lx) * (f12 - f11) + (lx - x1) * (f22 - f21)) / (dy * dx)
    
    grad_array[j, valid, 1] <- dfdx[valid]
    grad_array[j, valid, 2] <- dfdy[valid]
  }
  
  grad_array
}


#simulate Langevin process
sim_langevin <- function(beta, gammasq, dt, Tmax, loc_0, covlist){
  sigma = sqrt(gammasq*dt)
  n_obs = Tmax/dt
  locs = rbind(matrix(loc_0, nrow=1), mvnfast::rmvn(n_obs-1, c(0,0), sigma = diag(sigma, 2,2), isChol = TRUE)) 
  for (i in 1:(n_obs-1)) {
    locs[i+1, ] = locs[i+1, ] + locs[i, ] + (gammasq*dt/2)*bilinearGrad(locs[i, ], covlist)%*%beta 
  }
  
  return(locs)
}

#likelihood using extended kalman filter
lik3 <- function(par){
  #log-likelihood
  l = 0
  
  for (j in 1:1) {
    #defining transition covariance matrix
    Q = diag(delta*par[4],2,2)
    
    #control matrix
    B = diag(delta*par[4]/2,2,2)
    
    #initial covariance guess
    P = 10*Q
    
    #initial state
    z = X[1, ]
    
    for (i in 2:nrow(X)) {
      #control vector
      u = bilinearGrad(z, covlist) %*% par[1:3]
      
      
      F_k = diag(1,2,2) + (delta*par[4]/2)*hessian(z, covlist, par)
      
      #predicted state estimate
      z_p = z + B %*% u 
      
      #predicted estimate covariance 
      P = F_k %*% P %*% t(F_k) + Q
      
      #innovation covariance
      S = P 
      
      #updated state estimate
      z =  X[i, ]
      
      #updated estimate covariance
      P = diag(0,2,2)
      
      #adding likelihood contribution of i-th state
      l = l - dmvnorm(c(X[i, ] - z_p), mean = c(0,0), sigma = S, log = T)
      
      
      for (k in 1:m) {
        #control vector
        u = bilinearGrad(z, covlist) %*% par[1:3]
        
        
        F_k = diag(1,2,2) + (delta*par[4]/2)*hessian(z, covlist, par) 
        
        #predicted state estimate
        z_p = z + B %*% u
        
        #predicted estimate covariance 
        P = F_k %*% P %*% t(F_k) + Q
        
        #innovation covariance
        S = P 
        
        #updated state estimate
        z = z_p 
        
        #updated estimate covariance
        P = P
        
      }
      
    }
  }
  if(is.infinite(l)){
    return(1e10)
  }
  return(l)
  
  
  
}




#beta1
results = matrix(data = NA, nrow = 200, ncol = 5)


for (ik in 38:100) {
  beta <- c(4,2,-0.1)
  dt = 0.01
  n_obs = 2500
  Tmax = n_obs*100*dt
  
  time <- seq(0, Tmax, by = dt)
  
  # Generate tracks 
  track = sim_langevin(beta, speed, dt, Tmax, c(0,0), covlist)
  
  thin = 10
  delta = dt*thin
    
    
  #thinning tracks
  X = track
  n = nrow(X)
  X = X[(0:(n%/%thin -1))*thin +1, ]
  times = time[(0:(n%/%thin -1))*thin +1]
    
    
  #shotening tracks to 5000 observations
  X = X[1:5000, ]
  times = times[1:5000]
  ID = rep(1, 5000)
  
  #computing gradient of track observations
  gradArray = bilinearGradVec(X, covlist)
  gradArray <- aperm(gradArray, perm = c(2, 3, 1))
  
  
  fit = langevinUD(locs=X, times=times, ID=ID, grad_array=gradArray)
  
  print(c(1, fit$betaHat, fit$gamma2Hat))
  
  results[2*ik-1, 1:3] = fit$betaHat
  results[2*ik-1, 4] = fit$gamma2Hat
  results[2*ik-1, 5] = 1
  
  
  
  m = 10
  delta = dt*thin/(m+1)
  cl <- makeCluster(12)     # set the number of processor cores
  clusterExport(cl, varlist = c("delta", "X", "bilinearGrad", "hessian", "covlist", "dmvnorm", "m", "lik3"))
  setDefaultCluster(cl=cl) # set 'cl' as default cluster
  o = optimParallel(par=c(0,0,0,1), fn=lik3, lower=c(-Inf, -Inf, -Inf, 0))
  setDefaultCluster(NULL)
  stopCluster(cl)
  
  par = o$par
  results[2*ik, ] = c(par, 2)
  print(c(2, par))
  print(ik)
}



df = data.frame(beta1 = results[,1], beta2 = results[,2], beta3 = results[,3], gammasq = results[,4], method = factor(results[,5], levels = c(1, 2),labels = c("EM", "EKF")))
                                                                                                                     
save(df, file="EM_EKF_plot.Rda")                                                                                                               
                

p1 = ggplot(df, aes(x = method, y = beta1)) +
  geom_boxplot() +
  geom_hline(yintercept = 4, linetype = "dashed", color = "red") +
  labs(title = expression(beta[1]), y = expression(beta[1])) +
  theme_bw() 

p2 = ggplot(df, aes(x = method, y = beta2)) +
  geom_boxplot() +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  labs(title = expression(beta[2]), y = expression(beta[2])) +
  theme_bw() 

p3 = ggplot(df, aes(x = method, y = beta3)) +
  geom_boxplot() +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "red") +
  labs(title = expression(beta[3]), y = expression(beta[3])) +
  theme_bw() 

p4 = ggplot(df, aes(x = method, y = gammasq)) +
  geom_boxplot() +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
  labs(title = expression(gamma^2), y = expression(gamma^2)) +
  theme_bw() 

grid.arrange(p1,p2,p3,p4)
























