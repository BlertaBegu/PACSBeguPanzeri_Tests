# Define the true density function
mixture_true <- function(x, mu1, k1, beta1, gamma11, gamma12, 
                         mu2, k2, beta2, gamma21, gamma22, 
                         mu3, k3, beta3, gamma31, gamma32, 
                         mu4, k4, beta4, gamma41, gamma42, 
                         mu5, k5, beta5, gamma51, gamma52)
{  
  G1 <- cbind(mu1, gamma11, gamma12)
  G2 <- cbind(mu2, gamma21, gamma22)
  G3 <- cbind(mu3, gamma31, gamma32)
  G4 <- cbind(mu4, gamma41, gamma42)
  G5 <- cbind(mu5, gamma51, gamma52)
  
  
  return (dkent(x, G1, param=c(k1, beta1)) / 5 +
            dkent(x, G2, param=c(k2, beta2)) / 5 + 
            dkent(x, G3, param=c(k3, beta3)) / 5 + 
            dkent(x, G4, param=c(k4, beta4)) / 5 +
            dkent(x, G5, param=c(k5, beta5)) / 5 )
}

# Definition of the parameters of the mixture of 5 Kent Distrib
get_mu <- function(time){
  mu1 <- c(-0.5+time, -0.5+time, 0.8+time) 
  mu1 <- mu1 / sqrt( sum(mu1^2) )
  mu2 <- c(-0.3, -0.3, 0.2)
  mu2 <- mu2 / sqrt( sum(mu2^2) )
  mu3 <- c(0.5, -0.5, 0.8)
  mu3 <- mu3 / sqrt( sum(mu3^2) )
  mu4 <- c(0.2, -1, 0)
  mu4 <- mu4 / sqrt( sum(mu4^2) )
  mu5 <- c(0.6+time, -0.5+time, 0.3+time)
  mu5 <- mu5 / sqrt( sum(mu5^2) )
  
  return (list(mu1,mu2,mu3,mu4,mu5))
}

{
  
  gamma11 <- c(-0.7789378, 0.6157424, 0.1188163)
  gamma12 <- c(-0.5695773, -0.6154000, -0.5448528)
  k1=18
  beta1=0
  
  gamma21 <- c(-0.8651146, 0.3803316, -0.3269933)
  gamma22 <- c(0.1482597, -0.4288975, -0.8911038)
  k2=15
  beta2=7
  
  gamma31 <- c(-0.66647307, -0.74323532,-0.05843723)
  gamma32 <- c(0.5753645, -0.4629244, -0.6742824)
  k3=20
  beta3=10
  
  gamma41 <- c( 0.6364658, -0.0303920, -0.7707059)
  gamma42 <-  c(-0.7545879, -0.2314437, -0.6140285)
  k4=20
  beta4=7
  
  gamma51 <- c( 0.6364658, -0.0303920, -0.7707059)
  gamma52 <- c( 0.7545879, -0.2314437, -0.6140285)
  k5=20
  beta5=4
  
}


dens.func <- function(x, time)
{
  param <- get_mu(time)
  mu1 <- param[[1]]
  mu2 <- param[[2]]
  mu3 <- param[[3]]
  mu4 <- param[[4]]
  mu5 <- param[[5]]
  mixture_true(x, mu1, k1, beta1, gamma11, gamma12, 
                  mu2, k2, beta2, gamma21, gamma22, 
                  mu3, k3, beta3, gamma31, gamma32, 
                  mu4, k4, beta4, gamma41, gamma42, 
                  mu5, k5, beta5, gamma51, gamma52)
}

# Plot a FEM object with jet colormap with range [m,M] 
plot.FEM = function(FEM, M=NULL, m=NULL, ...){
  
  if (is.null(m)) { m = min(FEM$coeff)}
  if (is.null(M)) { M = max(FEM$coeff)}
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order=FEM$FEMbasis$mesh$order
  nodes=FEM$FEMbasis$mesh$nodes
  edges=matrix(rep(0,6*ntriangles),ncol=2)
  for(i in 0:(ntriangles-1)){
    edges[3*i+1,]=c(triangles[3*order*i+1],triangles[3*order*i+2])
    edges[3*i+2,]=c(triangles[3*order*i+1],triangles[3*order*i+3])
    edges[3*i+3,]=c(triangles[3*order*i+2],triangles[3*order*i+3])
  }
  edges=edges[!duplicated(edges),]
  edges<-as.vector(t(edges))
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  p=jet.col(n=128,alpha=0.8)
  # p <- colorRampPalette(c("#0E1E44","#3E6DD8","#68D061","#ECAF53", "#EB5F5F","#E11F1C"))(128)
  palette(p)
  
  ncolor=length(p)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
    rgl.pop("lights") 
    light3d(specular="black") 
    
    diffrange = M - m
    
    col = coeff[triangles,isurf]
    col= (col - min(coeff[,isurf]))/diffrange*(ncolor-1)+1
    
    rgl.triangles(x = nodes[triangles ,1], y = nodes[triangles ,2],
                  z=nodes[triangles,3],
                  color = col,...)
    # rgl.lines(x = nodes[edges ,1], y = nodes[edges ,2],
    #           z=nodes[edges,3],
    #           color = "black",...)
    aspect3d("iso")
    
    if (nsurf > 1 && isurf<nsurf)
    {readline("Press a button for the next plot...")}
  }
}

zoom = 1
userMatrix = rbind(c( 0.96563137,  0.1774523, -0.1899119,    0),
                   c( -0.03294301,  0.8083354,  0.5877997,    0),
                   c( 0.25781897, -0.5613416,  0.7863999,    0),
                   c(0.00000000,  0.0000000,  0.0000000 ,   1))
windowRect = c(150,  150, 420,  420)