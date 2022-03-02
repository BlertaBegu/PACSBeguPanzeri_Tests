################################################################################
###################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################### 
################### EXAMPLE 2 - KENT DENSITY OVER 2.5D DOMAIN ##################
################################################################################

## Introduction ----------------------------------------------------------------

library(fdaPDE)
library(Directional) # Library for generating Kent distributions
library(rgl) # Library for titles in 3D plots
rm(list=ls())
graphics.off()

source("helper_functions_2.R")

## Creation of a spatial 2.5D mesh with 606 nodes ------------------------------
vertices <- read.table("sphere.vertices.txt", quote="\"", comment.char="")
triangles <- read.table("sphere.triangles.txt", quote="\"", comment.char="")
mesh <- create.mesh.2.5D(nodes = vertices[,1:3], triangles = triangles[,1:3])
FEMbasis <- create.FEM.basis(mesh)

## Creation of a temporal 1D mesh over a non-negative interval -----------------
mesh_time <- seq(0, 1, length.out=5)

## Spatio-temporal data from Kent density --------------------------------------

# Data generation (5000 observations)
load("data")

## Spatio-temporal Density Estimation ------------------------------------------
lambda = 0.1
lambda_time = 0.01

solution <- DE.FEM.time(data = data[,1:3], data_time = data[,4],
                        FEMbasis = FEMbasis, mesh_time = mesh_time,
                        lambda = lambda, lambda_time = lambda_time,
                        fvec=NULL, heatIter=10, print=TRUE,
                        direction_method = 'Gradient',
                        preprocess_method='NoCrossValidation')

FEMfunction = FEM.time(solution$g, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)

## Visualization ---------------------------------------------------------------

# Fine grid
mesh.eval <- refine.by.splitting.mesh.2.5D(refine.by.splitting.mesh.2.5D(mesh))
FEMbasis.eval <- create.FEM.basis(mesh.eval)

# Time instants during which the solution is evaluated
t <- c(0.2,0.5,0.8)

# Plots
for (time_index in 1:length(t)) {

  # First plot: Sample
  plot(mesh)
  pch3d(data[abs(data[,4]-t[time_index])<(mesh_time[2]-mesh_time[1])/2,1:3], pch=19, cex=0.5, col="red")
  
  # Second plot: true density at t[time_index]
  data_sphere <- dens.func(mesh.eval$nodes, t[time_index])
  plot.FEM(FEM(data_sphere, FEMbasis.eval))
  
  # Third plot: estimated density at t[time_index]
  evaluation <- eval.FEM.time(FEM.time = FEMfunction, locations = mesh.eval$nodes, time.instants = t[time_index])
  plot.FEM(FEM(exp(evaluation), FEMbasis.eval))
  
}  
