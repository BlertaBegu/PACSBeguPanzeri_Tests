################################################################################
###################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################### 
################ EXAMPLE 1 - MIXTURE OF 4 GAUSSIAN DISTRIBUTIONS ###############
################################################################################

## Introduction ----------------------------------------------------------------

library(fdaPDE)
library(mvtnorm) # Library to generate data
library(mixtools) # Library to draw ellipses
rm(list=ls())
graphics.off()

source("helper_functions_1.R")

## Creation of a spatial 2D mesh over a squared domain -------------------------
Xbound <- seq(-6, 6, length.out=11)
Ybound <- seq(-6, 6, length.out=11)
grid_XY <- expand.grid(Xbound, Ybound)
Bounds <- grid_XY[(grid_XY$Var1 %in% c(-6, 6)) | (grid_XY$Var2 %in% c(-6, 6)), ]
mesh <- create.mesh.2D(nodes = Bounds, order = 1)
mesh <- refine.mesh.2D(mesh, maximum_area = 0.3, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

x11()
plot(mesh)

## Creation of a temporal 1D mesh over a non-negative interval -----------------
mesh_time <- seq(0,1, length.out=5)

## Spatio-temporal data from a mixture of Gaussian distributions ---------------
set.seed(10)

# Number of distributions
D <- 4 

# Number of observations
N_tot <- 2500

# Data generation
times <- c()
locations <- c()
distribution <- c()
for (n in 1:N_tot) {
  t <- runif(1,0,1)
  p <- get_priors(t)
  d <- sample(c(1, 2, 3, 4), 1, replace=T, prob=p)
  
  distribution <- c(distribution, d)
  parameters <- get_parameters(d, t)
  l <- mvtnorm::rmvnorm(n=1, mean=parameters[[1]], sigma=parameters[[2]], method="svd")
  
  locations <- rbind(locations, l)
  times <- c(times, t)
}
data <- cbind(locations, times)

## Spatio-Temporal Density Estimation ------------------------------------------
lambda <- c(0.1, 0.01)
lambda_time <- c(0.01, 0.001)
solution <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                        mesh_time = mesh_time, lambda = lambda,
                        lambda_time = lambda_time, fvec=NULL, heatStep=0.1,
                        heatIter=10, print=TRUE, nfolds=2, nsimulations=1000,
                        step_method="Fixed_Step", direction_method="Gradient",
                        preprocess_method="RightCV", flagMass=0, flagLumped=0)

FEMfunction = FEM.time(solution$g, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)

## Visualization ---------------------------------------------------------------

# Fine grid
n = 100
X <- seq(-6, 6, length.out = n)
Y <- seq(-6, 6, length.out = n)
grid <- expand.grid(X, Y)

# Time instants during which the solution is evaluated
t <- c(0.2, 0.5, 0.8)

# Plots
for (time_index in 1:length(t)) {
  x11(width=11.5,height=4)
  par(mfrow=c(1,3))
  
  # First plot: Sample
  plot(0,0, xlim=c(-6,6), ylim=c(-6,6), xlab='', ylab='', col='white',
       main='Sample', asp=1, xaxt='n', yaxt='n')
  for (d in 1:D) {
    plot_samples(as.data.frame(data[distribution==d & abs(data[,3]-t[time_index])<(mesh_time[2]-mesh_time[1])/2,1:2]))
  }
  
  # Second plot: True density at t[time_index]
  data_grid <- get_true_data(grid, t[time_index])
  
  image2D(x=X, y=Y, z=matrix(as.matrix(data_grid), n, n), col=heat.colors(100),
          xlab="", ylab="", contour=list(drawlabels = FALSE),
          main=paste("True density at t = ", t[time_index]), zlim=c(0,0.08),
          asp=1, xaxt="n", yaxt="n")

  # Third plot: Estimated density at t[time_index]
  evaluation <- eval.FEM.time(FEM.time = FEMfunction, locations = grid,
                              time.instants = t[time_index])
  image2D(x=X, y=Y, z=matrix(exp(evaluation), n, n), col=heat.colors(100),
         xlab="", ylab="", contour=list(drawlabels = FALSE),
         main=paste("Estimated density at t = ", t[time_index]), zlim=c(0,0.08),
         asp=1, xaxt="n", yaxt="n") 
}
