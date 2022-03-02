## Mixture weights (non-normalized prior probabilities) ------------------------
get_priors <- function(time) {
  var1 <- mean(eigen(get_parameters(1,time)[[2]])$values)
  var2 <- mean(eigen(get_parameters(2,time)[[2]])$values)
  var3 <- mean(eigen(get_parameters(3,time)[[2]])$values)
  var4 <- mean(eigen(get_parameters(4,time)[[2]])$values)
  priors = c(var1, var2, var3, var4)
  return (priors)
}


## Distribution parameters -----------------------------------------------------
get_parameters <- function(d, time) {
  if (d == 1) {
    # First Gaussian distribution
    mu <- c(-2, -1.5-0.5*time) # Mean
    sigma <- matrix(c(0.8, -0.2-0.4*time, -0.2-0.4*time, 0.8), 2, 2) # Covariance matrix
  }
  else if (d == 2) {
    # Second Gaussian distribution
    mu <- c(2+time, -2-time) # Mean
    sigma <- matrix(c(1.5-0.5*time, 0, 0, 1.5-0.5*time), 2, 2) # Covariance matrix
  }
  else if (d == 3) {
    # Third Gaussian distribution
    mu <- c(-2, 1.5+1.5*time) # Mean
    sigma <- matrix(c(0.8+time, 0, 0, 0.8), 2, 2) # Covariance matrix
  }
  else if (d == 4) {
    ## Fourth Gaussian distribution
    mu <- c(2, 2-time) # Mean
    sigma <- matrix(c(1, 0.9-0.3*time, 0.9-0.3*time, 1), 2, 2) # Covariance matrix
  }
  return (list(mu, sigma))
}


## Function to plot data -------------------------------------------------------
plot_samples <- function(samples){
  points(samples, col='black', pch=20, xlim=c(-6,6), ylim=c(-6,6))
  mixtools::ellipse(sapply(samples,mean), cov(samples), alpha=0.05, lty=2, lwd=1, col='black')
  mixtools::ellipse(sapply(samples,mean), cov(samples), alpha=0.5, lty=2, lwd=1, col='black')
  points(mean(samples[,1]), mean(samples[,2]), pch=20, cex=1.5, lwd=2, col='black', xlim=c(-6,6), ylim=c(-6,6))
  grid()
}


## Function to generate data from true density for visualization ---------------
get_true_data <- function(grid, time){
  p <- get_priors(time)
  p <- p/sum(p)
  data_grid <- p[1] * mvtnorm::dmvnorm(grid, mean=get_parameters(1,time)[[1]], sigma=get_parameters(1,time)[[2]]) +
               p[2] * mvtnorm::dmvnorm(grid, mean=get_parameters(2,time)[[1]], sigma=get_parameters(2,time)[[2]]) +
               p[3] * mvtnorm::dmvnorm(grid, mean=get_parameters(3,time)[[1]], sigma=get_parameters(3,time)[[2]]) +
               p[4] * mvtnorm::dmvnorm(grid, mean=get_parameters(4,time)[[1]], sigma=get_parameters(4,time)[[2]])
  
  return (data_grid)
}