library(geigen)
library(bayesreg)
library(mgcv)
library(psre) # TP basis
library(fda) # Fourier basis
library(orthopolynom) # Legendre Polynomials basis
library(rstan)
library(rstanarm) # stan_gamm4
library(bamlss)
library(rjags)
library(stringi) # for stri_paste
library(stringr) # for str_replace_all
library(plotrix) # for plotCI
library(caret) # for confusionMatrix
library(reshape2) # for melt
library(ggplot2)
library(xtable)
library(mvtnorm) # multivariate normal distribution
library(MixAll) # DebTrivedi data set for Poisson regression
library(acid) #confband.kneib
################################################################################
# kappa of horseshoe prior that is U-shaped beta(0.5, 0.5) distributed
x = seq(0, 1, length = 100)
plot(x, dbeta(x, 0.5, 0.5), type='l', xlab = "kappa", 
     ylab = "Density up to proportionality")

# Test functions
h1 <- function(x) {2*x - 1}
h2 <- function(x) {x^2 - 25/12}
h3 <- function(x) {x^3 + 1/4}
h4 <- function(x) {sin(5*x)}
h5 <- function(x) {7*exp(-exp(5*x))}
h6 <- function(x){exp(-x) - 0.4*sinh(2.5)}
h7 <- function(x, a = 2.5) {x = a*x
1.5*(sin(1.25*pi*x + 0.5)/exp(x))}
h8 <- function(x, a = 0.5, b = 1) {x = a*(x + b)
0.25*(20*dnorm(20*(x - 0.15)) + 5*dnorm(5*(x - 0.6)))}
h9 <- function(x, a = 0.5) {x = a*x 
0.1*sin(2*pi*x) + 0.2*cos(2*pi*x) + 0.3*(sin(2*pi*x))^2 + 
  0.4*(cos(2*pi*x))^3 + 0.5*(sin(2*pi*x))^3}

# Plotting functions h1 - h9 ###################################################
par(mfrow = c(3, 3))
for (i in 1:9){
  h <- H[[i]]
  plot(h, ylab = bquote('h'[.(i)](x)) , xlim = c(-1, 1)) # main = bquote('h'[.(i)](x))
}
par(mfrow = c(1, 1))

# Creating functions for error terms, plot design matrix and rmse calculation
norm_case <- function(h, mean = 0, s = 1){
  # adding a normally distributed iid errors to standardized values
  sc <- scale(h) #  standardization
  e <- rnorm(length(h), mean, s)
  e <- e*attr(sc, 'scaled:scale')
  # scale e based on standard deviation of h
  y <- h + e
  return(y)}
create.additive <- function(F, s = 1, logistic = FALSE){ 
  # additive model with normal distributed error or additive logistic model
  
  n <- dim(F)[1]
  f <- rowSums(F)
  
  sc <- scale(f) #  standardization
  
  e <- rnorm(n, sd = s)
  e <- e*attr(sc, 'scaled:scale')
  # scale e based on standard deviation of h
  
  if (!logistic) {
    y <- f + e}
  
  if (logistic) {
    prob <- exp(f)/(exp(f) + 1) # logistic additive model
    y <- rbinom(n, size = 1, prob = prob)}
  
  return(y)
}

plot.design.matrix <- function(x,X){ 
# function to plot design matrix of basis functions
  p <- dim(X)[2]
  ord <- order(x)
  plot(x[ord], X[ord, 1], type = "l", xlab = "x", ylab = "Basis", 
       ylim = c(min(X), max(X)), col = 1)
  for(j in 2:p){
    lines(x[ord], X[ord, j], col = j)
  }
}

RMSE <- function(y_hat, y){ # function to calculate RMSE
  as.numeric(sqrt(mean((y_hat - y)^2)))
}
# Functions to create basis evaluation #########################################
# Truncated Polynomials basis
create.TP.basis <- function(x, degree = 3, nknots = 7, constant = FALSE){
  
  if (is.matrix(x) == FALSE){
    x <- matrix(x, ncol = 1)}
  n = dim(x)[1]
  p = dim(x)[2]
  if (constant) {TP <- cbind(rep(1, n), tpb(x[, 1], degree, nknots))}
  if(!constant) {TP <- tpb(x[, 1], degree, nknots)}
  if (p > 1){
    for (i in (2:p)){
      TP <- cbind(TP, tpb(x[, 1], degree, nknots))}}
  
  return(TP)
}

################################################################################
# B-spline basis

create.B.basis <- function(x, nbasis = 10, constant = FALSE, const_value = 1){
  if (is.matrix(x) == FALSE){
    x <- matrix(x, ncol = 1)}
  n = dim(x)[1]
  p = dim(x)[2]
  fit <- smoothCon(s(x[, 1], bs = "ps", k = nbasis), data = as.data.frame(x[, 1]))
  B <- fit[[1]]$X
  if (constant) {B <- (cbind(rep(const_value, n), B))}
  if(!constant){B <- B}
  if (p > 1){
    for (i in (2:p)){
      fit <- smoothCon(s(x[, i], bs = "ps", k = nbasis), data = as.data.frame(x[, i]))
      B <- cbind(B, fit[[1]]$X)}}
  return(B)}

################################################################################
# Mixed Model basis based on Eigenvalue Problem

create.MM.basis <- function(x, nbasis = 10, rescale = FALSE, constant = FALSE){
  k = nbasis + 1 # constant basis function will be removed later on
  
  MMB <- list()
  if (is.matrix(x) == FALSE){
    x <- matrix(x, ncol = 1)}
  p = dim(x)[2]
  for (i in 1:p){ 
    fit <- smoothCon(s(x[, i], bs = "ps", k = k), data = as.data.frame(x[, i]))
    B <- fit[[1]]$X
    omega <- matrix(unlist(fit[[1]]$S), ncol = k)
    n = length(x)
    
    Pei <- eigen(omega)
    Q <- Pei$vectors
    Pei$values[c(k - 1, k)] <- 1 # replace the last two eigenvalues by 1
    sqrtlambda = diag(sqrt(Pei$values))
    
    if (!rescale) {MM <- B %*% Q}
    if (rescale) {MM <- B %*% Q %*% solve(sqrtlambda)}
    
    MM <- MM[, dim(MM)[2]:1] 
    # reverse the column order so that the linear functions are the first two
    
    MM <- MM[,-c(1,2)] # remove the two linear basis functions
    MMB[[i]] <- cbind(scale(x[, i]), MM)} # attach standardized linear function
  
  
  if (constant) {
    MMfull <- (cbind(rep(1, dim(MMB[[1]])[1]), MMB[[1]]))}
  if(!constant) {
    MMfull <- MMB[[1]]}
  
  if (p > 1){
    for(j in 2:p){
      MMfull <- cbind(MMfull, MMB[[j]])}}
  return(MMfull)}
################################################################################
# Demmler-Reinsch basis based on Generalized Eigenvalue Problem

create.DR.basis <- function(x, nbasis = 10, rescale = FALSE, constant = FALSE) {
  k = nbasis + 1 # constant basis function will be removed later on
  
  DRB <- list()
  if (is.matrix(x) == FALSE){
    x <- matrix(x, ncol = 1)}
  p = dim(x)[2]
  for (i in 1:p){
    fit <- smoothCon(s(x[, i], bs = "ps", k = k), data = as.data.frame(x[, i]))
    B <- fit[[1]]$X
    omega <- matrix(unlist(fit[[1]]$S), ncol = k)
    
    n = length(x[, i])
    G <- t(B) %*% B/n # Gramian Matrix
    Gei <- geigen::geigen(omega, G)
    names(Gei)
    A = Gei$vectors
    Gei$values[1:2] <- 1            
    sqrtgamma = diag(sqrt(Gei$values))
    
    if (!rescale) {DR <- B %*% A}
    if (rescale) {DR <- B %*% A %*% solve(sqrtgamma)}
    
    
    DR <- DR[, -c(1,2)] # remove the two linear basis functions
    DRB[[i]] <- cbind(scale(x[, i]), DR)} # attach standardized linear function
  
  if (constant) {
    DRfull <- (cbind(rep(1, dim(DRB[[1]])[1]), DRB[[1]]))}
  if(!constant) {
    DRfull <- DRB[[1]]}
  if (p > 1){
    for(j in 2:p){
      DRfull <- cbind(DRfull, DRB[[j]])}}
  return(DRfull)}

################################################################################
# Fourier Basis

create.F.basis <- function(x, domain = c(-1, 1), nbasis = 10, constant = FALSE, period = 2){
  
  k = nbasis + 1
  if (is.matrix(x) == FALSE){
    x <- matrix(x, ncol = 1)}
  p <- dim(x)[2]
  fb <- create.fourier.basis(rangeval = domain, nbasis = k, period = period)
  FB <- eval.basis(x[, 1], basisobj = fb)
  if (nbasis %% 2 == 1){FB <- FB[, -(nbasis + 2)]}
  if (constant) {FB[, 1] <- 1}
  if(!constant) {FB <- FB[, -1]}
  if (p > 1){
    for (i in 2:p){
      if (nbasis %% 2 == 1){
        FB <- cbind(FB, eval.basis(x[, i], basisobj = fb)[, c(-1, -(nbasis + 2))])}
      else {FB <- cbind(FB, eval.basis(x[, i], basisobj = fb)[, -1])}}}
  return(FB)
}

################################################################################
# Orthogonal Polynomials basis
create.OP.basis <- function(x, nbasis = 10, constant = FALSE, const_value = 1){
  
  
  if (is.matrix(x) == FALSE){
    x <- matrix(x, ncol = 1)}
  
  p <- dim(x)[2]
  
  OP <- poly(x[, 1], degree = nbasis, raw = FALSE)
  
  if (constant) {OP <- (cbind(rep(const_value, dim(OP)[1]), OP))}
  if(!constant) {OP <- OP}
  
  if (p > 1){
    for (i in 2:p){
      OP <- cbind(OP, poly(x[, i], degree = nbasis, raw = FALSE))}}
  
  return(OP)
}

get.OP.evaluation <- function(x, nbasis = 3, x_new){
  
  if (is.matrix(x) == FALSE){
    x <- matrix(x, ncol = 1)}
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  GP = cbind(rep(1, n), poly(x[, 1], degree = nbasis, raw = TRUE)) 
  # Monomial basis or General Polynomials
  
  qr = qr(GP) ## conduct QR decomposition of moniomial design matrix
  R = qr.R(qr) ## this is the transition matrix
  signs <- sign(diag(R))
  
  GP_new <- cbind(rep(1, dim(x_new)[1]), poly(x_new[, 1], degree = nbasis, raw = TRUE)) 
  OP <- (GP_new %*% solve(R))
  for (i in 1:length(signs)){
    if (signs[i] == -1){
      OP[, i] <- -1*OP[, i]
    }
  }
  
  OP <- OP[, -1]
  
  if (p > 1){
    for (i in 2:p){
      GP = cbind(rep(1, n), poly(x[, i], degree = nbasis, raw = TRUE)) 
      # Monomial basis or General Polynomials
      
      qr = qr(GP) ## conduct QR decomposition of moniomial design matrix
      R = qr.R(qr) ## this is the transition matrix
      signs <- sign(diag(R))
      
      GP_new <- cbind(rep(1, dim(x_new)[1]), poly(x_new[, i], 
                                                  degree = nbasis, raw = TRUE)) 
      OP_new <- (GP_new %*% solve(R))
      for (i in 1:length(signs)){
        if (signs[i] == -1){
          OP_new[, i] <- -1*OP_new[, i]
        }
      }
      OP <- cbind(OP, OP_new[, -1])}}
  
  return(OP)
}

################################################################################
# Legendre Polynomial basis

create.LP.basis <- function(x, nbasis = 10, constant = FALSE, const_value = 1){
  if (is.matrix(x) == FALSE){
    x <- matrix(x, ncol = 1)}
  
  p <- dim(x)[2]
  
  l <- legendre.polynomials(n = nbasis)
  lpl <- polynomial.values(l, x) # list of LP evaluated basis values
  
  LPB <- lpl[[2]]
  
  for (j in 3:(nbasis+1)){
    LPB <- cbind(LPB, lpl[[j]])
  }
  
  
  LP <- LPB[, seq(from = 1, to = p*nbasis, by = p)]
  if (p > 1) {
    for (j in 2:p){ # rearrange the columns so that it has the same structure as 
      # other basis (column 1 to nbasis are for first covariate, 
      # {nbasis + 1} to  2*basis for second covariate and so on
      LP <- cbind(LP, LPB[, seq(from = j, to = p*nbasis, by = nbasis)])
    }}
  
  if (constant){LP <- cbind(rep(const_value, dim(LP)[1]), LP)}
  
  return(LP)
}


basis_functions <- list(create.TP.basis, create.B.basis, create.MM.basis, create.DR.basis, 
                        create.F.basis, create.OP.basis, create.LP.basis)
names(basis_functions) <- c("TP", "B-spline", "MM", "DR", "Fourier", "OP", "LP")
# Plot 9 basis functions of each basis #########################################
# Basis functions bases
domain = c(-1, 1)
n = 1000
set.seed(934)

x = sort(runif(n, min = domain[1], max = domain[2]))
plot.design.matrix(x, create.B.basis(x = x, nbasis = 9))

for (j in 1:7){
  x = sort(runif(n, min = domain[1], max = domain[2]))
  
  if (j == 1){
    D <- data.frame(basis_functions[[j]](x, degree = 3, nknots = 6))
  } else {D <- data.frame(basis_functions[[j]](x, nbasis = 9))}
  colnames(D) <- 1:9
  x <- do.call("rbind", replicate(9, data.frame(x), simplify = FALSE))
  
  D <- cbind(melt(D, measure.vars = colnames(D)), x)
  g <- ggplot(D, aes(x = x, y = value)) + 
    geom_line() +
    facet_wrap( ~ variable, scales = 'fixed', nrow = 3, ncol = 3) +
    theme(plot.title = element_text(hjust = 0.45)) +
    xlab("x") + ylab("")
  print(g)
}

# An example of how B-spline works #############################################
x = sort(runif(1000, min = -1, max = 1))
B <- create.B.basis(x, nbasis = 10)
gamma <- c(1, 2, 1, 3, 6, 1, 3, 2, 4, 1) # additive component of scaled basis functions
f <- gamma %*% t(B)
plot(x, f, type = "l",ylim = c(0, 5), ylab = "B(x) and f(x)")
for (i in 1:10){
  #lines(x, B[, i], col = i)
  lines(x, gamma[i]*B[, i], col = i)
}

# Estimation configuration for Nonparametric Model
domain = c(-1, 1)
n = 100
set.seed(934)
k = 10
uni = FALSE 
# TRUE if the covariates are uniformly distributed, 
#FALSE indicates the use of standard normal distributed covariates
s = 2.5

H <- c(h1, h2, h3, h4, h5, h6, h7, h8, h9) # The test functions

# Create a list to store the results
mu.coefs <- matrix(nrow = 100, ncol = k + 1)
mus <- replicate(n = 7, mu.coefs, simplify = FALSE)
RMSEs <- matrix(nrow = 100, ncol = 7)
l <- list(mus, RMSEs)
names(l) <- c("Mean of coeeficients", "RMSEs")
N <- replicate(n = 9, l, simplify = FALSE)

names(N) <- paste0("h", (1:9))

# Begin for the loop of each configuration for fixed set of n and s
for (g in 1:9){ # Loop over all 9 functions
  h <- H[[g]]
  
  set.seed(934)
  for (i in 1:100) { # 100 replicates simulation of bayesreg
    if (uni) {x <- sort(runif(n, min = domain[1], max = domain[2]))}
    if (!uni) {x = rnorm(n)
    x = 2*(x - min(x))/(max(x) - min(x)) - 1}## transform to [-1,1]
    
    y <- norm_case(h(x), s = s)
    
    for (j in 1:7){ # For each of the 7 basis expansion
      base <- basis_functions[[j]](x, constant = TRUE)
      name <- names(basis_functions)[j]
      
      df <- data.frame(cbind(y, base))
      names(df) <- c("y", paste("V", 1:dim(base)[2], sep = ""))
      
      spline <- bayesreg(y ~. -1, data = df, prior = 'hs')
      betas <- c(spline$mu.beta0, spline$mu.beta)
      N[[g]][[1]][[j]][i, ] <- betas
      
      f <- t(betas) %*% t(base)
      N[[g]][[2]][i, j] <- RMSE(f, h(x))}}}

saveRDS(N, "NNonparametric100n0.4s.RData") # N stands for non-uniform covariates
# Nonparametric RMSEs Boxplot ##################################################
# list from each setting of Nonparametric Model with uniform covariates
N <- list(N1, N2, N3, N4, N5, N6, N7, N8, N9) 
names(N) <- c("1: n = 100, s = 0.4", "2: n = 100, s = 1", "3: n = 100, s = 2.5",
              "4: n = 500, s = 0.4", "5: n = 500, s = 1", "6: n = 500, s = 2.5",
              "7: n = 1000, s = 0.4", "8: n = 1000, s = 1", "9: n = 1000, s = 2.5")
RMSE_N <- list()

mu_RMSE_N <- data.frame(matrix(nrow = 9, ncol = 7))
rownames(mu_RMSE_N) <- c("1: n = 100, s = 0.4", "2: n = 100, s = 1", "3: n = 100, s = 2.5",
                         "4: n = 500, s = 0.4", "5: n = 500, s = 1", "6: n = 500, s = 2.5",
                         "7: n = 1000, s = 0.4", "8: n = 1000, s = 1", "9: n = 1000, s = 2.5")
colnames(mu_RMSE_N) <- names(basis_functions)
mu_RMSE_N <- replicate(n = 9, mu_RMSE_N, simplify = FALSE) 
# replicate one RMSE table for each function

for (j in 1:9){ # each function
  for (i in 1:9){ # each setting
    mu_RMSE_N[[j]][i, ] <- colMeans(N[[i]][[j]][[2]])}} # calculate the mean RMSE

# latex table for mean RMSE of different functions
for (i in 1:9){
  print(paste0("h", i))
  print(xtable(mu_RMSE_N[[i]], type = 'latex', digits = c(0, 4, 4, 4, 4, 4, 4, 4)))}

for (i in 1:9){ # for each function
  RMSEs <- data.frame(N[[1]][[i]][[2]]) # RMSEs of the setting
  colnames(RMSEs) <- names(basis_functions)
  RMSEs <- melt(RMSEs, measure.vars = colnames(RMSEs)) # reshape RMSE for ggplot
  RMSEs <- cbind(config = names(N)[1], RMSEs)
  RMSE_N[[i]] <- RMSEs
}

# For loop to bind the rmse from other settings together
for (i in 1:9){ # function number
  RMSEs <- RMSE_N[[i]]
  for (j in 2:9){ # setting number
    tRMSEs <- data.frame(N[[j]][[i]][[2]])
    colnames(tRMSEs) <- names(basis_functions)
    tRMSEs <- melt(tRMSEs, measure.vars = colnames(tRMSEs))
    tRMSEs <- cbind(config = names(N)[j], tRMSEs)
    RMSEs <- rbind(RMSEs, tRMSEs)
  }
  RMSE_N[[i]] <- RMSEs
} 

# Boxplot for each setting
for (i in 1:9){
  RMSEs <- RMSE_N[[i]]
  
  p <- ggplot(RMSEs, aes(x = variable, y = log(value^2), fill = variable)) + 
    geom_boxplot() + # coord_flip() + outlier.colour="red", outlier.shape = 4, outlier.size = 4) +
    facet_wrap( ~ config, scales = "fixed", nrow = 3, ncol = 3) + 
    ggtitle(paste0("log MSE for estimation of h", i, " with uniform covariates")) +
    theme(plot.title = element_text(hjust = 0.45), axis.text.x = element_blank()) +
    xlab("Basis") + ylab("log MSE") + scale_fill_discrete(name = "Basis")
  print(p)
}
################################################################################
# list from each setting of Nonparametric Model with non-Uniform covariates
NN <- list(NN1, NN2, NN3, NN4, NN5, NN6)
names(NN) <- c("1: n = 500, s = 0.4", "2: n = 500, s = 1", "3: n = 500, s = 2.5",
               "4: n = 1000, s = 0.4", "5: n = 1000, s = 1", "6: n = 1000, s = 2.5")
RMSE_NN <- list()

mu_RMSE_NN <- data.frame(matrix(nrow = 6, ncol = 7))
rownames(mu_RMSE_NN) <- c("1: n = 500, s = 0.4", "2: n = 500, s = 1", "3: n = 500, s = 2.5",
                          "4: n = 1000, s = 0.4", "5: n = 1000, s = 1", "6: n = 1000, s = 2.5")
colnames(mu_RMSE_NN) <- names(basis_functions)
mu_RMSE_NN <- replicate(n = 9, mu_RMSE_NN, simplify = FALSE) # replicate one RMSE table for each function

for (j in 1:9){ # each function
  for (i in 1:6){ # each setting
    mu_RMSE_NN[[j]][i, ] <- colMeans(NN[[i]][[j]][[2]])}}

# latex table for mean RMSE of different functions
for (i in 1:9){
  print(paste0("h", i))
  print(xtable(mu_RMSE_NN[[i]], type = 'latex', digits = c(0, 4, 4, 4, 4, 4, 4, 4)))}

for (i in 1:9){ # for each function
  RMSEs <- data.frame(NN[[1]][[i]][[2]])
  colnames(RMSEs) <- names(basis_functions)
  RMSEs <- melt(RMSEs, measure.vars = colnames(RMSEs))
  RMSEs <- cbind(config = names(NN)[1], RMSEs)
  RMSE_NN[[i]] <- RMSEs
}

for (i in 1:9){ # function number
  RMSEs <- RMSE_N[[i]]
  for (j in 2:6){ # setting number
    tRMSEs <- data.frame(NN[[j]][[i]][[2]])
    colnames(tRMSEs) <- names(basis_functions)
    tRMSEs <- melt(tRMSEs, measure.vars = colnames(tRMSEs))
    tRMSEs <- cbind(config = names(NN)[j], tRMSEs)
    RMSEs <- rbind(RMSEs, tRMSEs)
  }
  RMSE_N[[i]] <- RMSEs
} 

for (i in 1:9){
  RMSEs <- RMSE_N[[i]]
  
  p <- ggplot(RMSEs, aes(x = variable, y = log(value^2), fill = variable)) + 
    geom_boxplot() +
    facet_wrap( ~ config, scales = "fixed", nrow = 2, ncol = 3) + 
    theme(plot.title = element_text(hjust = 0.45), axis.text.x = element_blank(), 
          legend.position="bottom") +
    xlab("Basis") + ylab("log MSE") + scale_fill_discrete(name = "Basis")
  print(p)
}

# Small Simulation Study to compare MM, DR and OP in (G)AM #####################
n = 500
p = 10 # number of covariates
rho = 0.5 # correlation between covariate 1 and 2
k = 10 ## number of basis functions for each additive component
nr = 6 # number of relevant covariates
uni = FALSE # uniform or non-uniform covariates
s = 2.5

RMSEs <- matrix(nrow = 100, ncol = 8) # for all 4 models

H <- c(h7, h1, h2, h8, h3, h4, h9, h5, h6)

# which basis functions and if it is a logistic model
settings <- list(c(3, FALSE), c(4, FALSE), c(6, FALSE), c(7, FALSE),
                 c(3, TRUE), c(4, TRUE), c(6, TRUE), c(7, TRUE))
names(settings) <- c("AM MM", "AM DR", "AM OP", "AM LP", "GAM MM", "GAM DR", "GAM OP", "GAM LP")

# Begin for the loop of each configuration for fixed set of sample size, 
# scale of Var(f) as standard deviation and number of relevant covariates
set.seed(934)
for (i in 1:100){ # 100 replicates of simulation for each model
  if (uni) {X <- matrix(runif(n*p, min = domain[1], max = domain[2]),n,p)}
  
  Sigma = matrix(0, p, p)
  for(g in 1:p){
    for(q in 1:p){
      Sigma[g,q] = rho^(abs(g-q))
    }
  }
  
  X = rmvnorm(n = n, sigma = Sigma)
  for (j in 1:p){
    X[, j] <- 2*(X[, j] - min(X[, j]))/(max(X[, j]) - min(X[, j])) - 1
  }
  
  F <- matrix(0, n, p)
  
  for (j in 1:nr){ # loop to add relevant covariates into the model
    h <- H[[j]]
    F[, j] <- h(X[, j])
  }
  # transform relevant covariates with function h
  
  eta_true <- rowSums(F)
  
  for (si in 1:8){
    bi <- settings[[si]][1]
    bin <- settings[[si]][2]
    
    for (bin in c(TRUE, FALSE)) {
      y <- create.additive(F, s = s, logistic = bin)
      
      if (bi != 1){B <- basis_functions[[bi]](X, nbasis = k, constant = FALSE)}
      if (bi == 1){B <- basis_functions[[bi]](X, degree = 3, 
                                              nknots = k - 3, constant = FALSE)}
      
      df <- as.data.frame(cbind(y, B))
      names(df) <- c("y", paste("V", 1:(p*k),sep = ""))
      B <- df[, -1]
      if (bin == FALSE) {fit <- bayesreg(y ~. , data = df, model = "normal", 
                                         prior = "hs", n.samples = 4e3)
      } else if (bin == TRUE) {
        fit <- bayesreg(as.factor(y) ~., data = df, model = "logistic", 
                        prior = "hs", n.samples = 4e3)}
      
      eta_hat <- predict(fit, B, type = 'linpred')
      
      RMSEs[i, si] <- RMSE(eta_hat, eta_true)
    }
  }
  
  
  if (i %% 5 == 0) {print(paste("Iteration ", i, " complete"))}}
# RMSEs boxplot to compare MM, DR, OP and LP in AM and GAM #####################
RMSEs <- melt(RMSEs, measure.vars = colnames(RMSEs))
RMSEs$variable <- str_replace_all(RMSEs$variable, c("X1" = "MM", "X2" = "DR", 
                                                    "X3" = "OP", "X4" = "LP", 
                                                    "X5" = "MM", "X6" = "DR", 
                                                    "X7" = "OP", "X8" = "LP"))
n <- dim(RMSEs)[1]
model_type <- data.frame(matrix(nrow = n, ncol = 1))
model_type[1:(n/2), ] <- "Additive Model"
model_type[(n/2+1):n, ]<- "Generalized Additive Model"
RMSEs <- cbind(model_type, RMSEs)
colnames(RMSEs) <- c("model", "variable", "value")

ggplot(RMSEs, aes(x = variable, y = log(value^2), fill = variable)) + 
  geom_boxplot() + 
  facet_wrap( ~ model, scales = "fixed", nrow = 1, ncol = 2) + 
  theme(plot.title = element_text(hjust = 0.45)) + scale_fill_discrete(name = "Basis") +
  xlab("Basis") + ylab("RMSE")
# Additive Model and Generalized Additive Model Simulation #####################
# create jag file for jagam
jag.file <- paste(getwd(),"/test.jags",sep = "")

# Configuration for (Generalized) Additive Model
domain = c(-1, 1)
p = 10 # number of covariates
rho = 0.5 # correlation between covariate 1 and 2
k = 10 ## number of basis functions for each additive component
bi = 6 # which basis function

# Setttings of different AM and GAM models
AM_settings <- list(c(TRUE, FALSE, 500, 1, 3), c(TRUE, FALSE, 500, 2.5, 3), c(TRUE, FALSE, 1000, 1, 3), 
                    c(TRUE, FALSE, 1000, 2.5, 3), c(TRUE, FALSE, 500, 1, 6), c(TRUE, FALSE, 500, 2.5, 6), 
                    c(TRUE, FALSE, 1000, 1, 6), c(TRUE, FALSE, 1000, 2.5, 6), c(TRUE, FALSE, 500, 1, 9), 
                    c(TRUE, FALSE, 500, 2.5, 9), c(TRUE, FALSE, 1000, 1, 9) , c(TRUE, FALSE, 1000, 2.5, 9))

NAM_settings <- list(c(FALSE, FALSE, 500, 1, 3), c(FALSE, FALSE, 500, 2.5, 3), c(FALSE, FALSE, 1000, 1, 3), 
                    c(FALSE, FALSE, 1000, 2.5, 3), c(FALSE, FALSE, 500, 1, 6), c(FALSE, FALSE, 500, 2.5, 6), 
                    c(FALSE, FALSE, 1000, 1, 6), c(FALSE, FALSE, 1000, 2.5, 6), c(FALSE, FALSE, 500, 1, 9), 
                    c(FALSE, FALSE, 500, 2.5, 9), c(FALSE, FALSE, 1000, 1, 9) , c(FALSE, FALSE, 1000, 2.5, 9))
# N for non-uniform covariates

GAM_settings <- list(c(TRUE, TRUE, 500, NA, 3), c(TRUE, TRUE, 500, NA, 6), c(TRUE, TRUE, 500, NA, 9), 
                     c(TRUE, TRUE, 1000, NA, 3), c(TRUE, TRUE, 1000, NA, 6), c(TRUE, TRUE, 1000, NA, 9))
# note that s does not affect GAM simulation

NGAM_settings <- list(c(FALSE, TRUE, 500, NA, 3), c(FALSE, TRUE, 500, NA, 6), c(FALSE, TRUE, 500, NA, 9), 
                  c(FALSE, TRUE, 1000, NA, 3), c(FALSE, TRUE, 1000, NA, 6), c(FALSE, TRUE, 1000, NA, 9))
# N for non-uniform covariates

mt <- NGAM_settings # choose a model type

# Rearrange the functions to be the same as described in simulation design
H <- c(h7, h1, h2, h8, h3, h4, h9, h5, h6)

mlist <- list()

for (set in 1:length(mt)){ # For each setting
  setting <- mt[[set]]
  uni = setting[1] # is the covariate uniformly distributed between [-1, 1]
  bin = setting[2] # binary logistic regression or simple additive model
  n = setting[3] # sample size
  s = setting[4] # how many times should the error scales on Var(f)
  nr = setting[5] # number of relevant covariates

  mu.coefs <- matrix(nrow = 50, ncol = p*k + 1) # for bayesreg
  eta_hat <- matrix(nrow = 50, ncol = n)
  eta_hat <- replicate(n = 4, eta_hat, simplify = FALSE) # for each models
  RMSEs <- matrix(nrow = 50, ncol = 4) # for all 4 models
  time <- matrix(nrow = 50, ncol = 4) # for all 4 models
  
  l <- list(mu.coefs, eta_hat, RMSEs, time)
  names(l) <- c("Mean of coeeficients", "Fitted values", "RMSEs", "Run time")
  
  set.seed(934)
  sl <- Sys.time()
  for (i in 1:50){ # 50 replicates of simulation for each model
    
    if (uni) {X <- matrix(runif(n*p, min = domain[1], max = domain[2]),n,p)}
    
    if (!uni) {
      Sigma = matrix(0, p, p)
      for(g in 1:p){
        for(q in 1:p){
          Sigma[g,q] = rho^(abs(g-q))
        }
      }
      
      X = rmvnorm(n = n, sigma = Sigma)}
    for (j in 1:p){
      X[, j] <- 2*(X[, j] - min(X[, j]))/(max(X[, j]) - min(X[, j])) - 1
    }
    F <- matrix(0, n, p)
    
    for (j in 1:nr){ # loop to add relevant covariates into the model
      h <- H[[j]]
      F[, j] <- h(X[, j])
    }
    # transform relevant covariates with function h
    
    eta_true <- rowSums(F)
    
    y <- create.additive(F, s = s, logistic = bin)
    
    if (bi != 1){B <- basis_functions[[bi]](X, nbasis = k, constant = FALSE)}
    if (bi == 1){B <- basis_functions[[bi]](X, degree = 3, nknots = k - 3, constant = FALSE)}
    
    df <- as.data.frame(cbind(y, B))
    names(df) <- c("y", paste("V", 1:(p*k),sep = ""))
    B <- df[, -1]
    
    data <- as.data.frame(cbind(y,X)) # for stan_gamm4, bamlss, jagam
    X <- data[, -1]
    V <- paste0("s(V", 2:(p+1), ")")
    V <- stri_paste(V, collapse = "+")
    formula <- as.formula(paste("y", V, sep = "~"))
    
    for (m in 1:4) { # perform the simulation for bayesreg, stan_gamm4, bamlss and jagam
      start <- Sys.time()
      if (m == 1){
        if (bin == FALSE) {fit <- bayesreg(y ~. , data = df, model = "normal", 
                                           prior = "hs", n.samples = 1e3, burnin = 200, thin = 1)
        } else if (bin == TRUE) {
          fit <- bayesreg(as.factor(y) ~., data = df, model = "logistic", 
                          prior = "hs", n.samples = 1e3, burnin = 200, thin = 1)}
      } else if (m == 2){ # stan_gamm4
        if (bin == FALSE) {fit <- stan_gamm4(formula, data = data, chain = 1, iter = 1200, warmup = 200, thin = 1)
        } else if (bin == TRUE) {
          fit <- stan_gamm4(formula, data = data, family = "binomial", chain = 1, iter = 1200, warmup = 200, thin = 1)}
      } else if (m == 3){ # bamlss
        if (bin == FALSE) {fit <- bamlss(formula, data = data, n.samples = 1e3, burnin = 200, thin = 1)
        } else if (bin == TRUE) {
          fit <- bamlss(formula, data = data, family = "binomial", n.samples = 1e3, burnin = 200, thin = 1)}
      } else if (m == 4){ # jagam
        if (bin == FALSE){
          jag <- jagam(formula, file = jag.file, data = data)
        } else if (bin == TRUE){
          jag <- jagam(formula, family = "binomial", file = jag.file, data = data)
        }
        jm <- jags.model("test.jags", data = jag$jags.data,
                         inits = jag$jags.ini)
        sam <- jags.samples(jm, c("b", "rho"), n.iter = 1e3, burnin = 200, thin = 1)
        jam <- sim2jam(sam, jag$pregam)
      }
      
      if (m == 1){
        l[[1]][i, ] <- c(fit$mu.beta0, fit$mu.beta)
        eta_hat <- predict(fit, B, type = 'linpred')
      }else if (m == 2) {eta_hat <- predict(fit, X, type = 'link')
      } else if (m == 3){if (bin == FALSE) {eta_hat <- predict(fit, data, type = 'link')$mu
      } else if (bin == TRUE) {eta_hat <- predict(fit, X, type = 'link')}
      } else if (m == 4) {eta_hat <- predict(jam, newdata = X)
      }
      
      l[[2]][[m]][i, ] <- eta_hat
      l[[3]][i, m] <- RMSE(eta_hat, eta_true)
      
      end <- Sys.time()
      runtime <- as.numeric(end - start, unit = "secs")
      l[[4]][i, m] <- runtime}
    
    if (i %% 5 == 0) {print(paste("Iteration ", i, " complete"))}
  }
  el <- Sys.time()
  print(el - sl)
  
  mlist[[set]] <- l
  
  print(paste("setting", set, "done"))
}

NGAM <- mlist
# AM and NAM RMSEs Time Boxplots ###############################################
names(AM) <- c("01: n = 500, s = 1, nr = 3", "02: n = 500, s = 2.5, nr = 3",
               "03: n = 1000, s = 1, nr = 3", "04: n = 1000, s = 2.5, nr = 3",
               "05: n = 500, s = 1, nr = 6", "06: n = 500, s = 2.5, nr = 6",
               "07: n = 1000, s = 1, nr = 6", "08: n = 1000, s = 2.5, nr = 6",
               "09: n = 500, s = 1, nr = 9", "10: n = 500, s = 2.5, nr = 9",
               "11: n = 1000, s = 1, nr = 9", "12: n = 1000, s = 2.5, nr = 9")

# create a data frame for mean rmse of each model and settings
mu_RMSE <- data.frame(matrix(nrow = 12, ncol = 4)) 
colnames(mu_RMSE) <- c("bayesreg", "stan_gamm4", "bamlss", "jagam")
rownames(mu_RMSE) <- c("01: n = 500, s = 1, nr = 3", "02: n = 500, s = 2.5, nr = 3",
                       "03: n = 1000, s = 1, nr = 3", "04: n = 1000, s = 2.5, nr = 3",
                       "05: n = 500, s = 1, nr = 6", "06: n = 500, s = 2.5, nr = 6",
                       "07: n = 1000, s = 1, nr = 6", "08: n = 1000, s = 2.5, nr = 6",
                       "09: n = 500, s = 1, nr = 9", "10: n = 500, s = 2.5, nr = 9",
                       "11: n = 1000, s = 1, nr = 9", "12: n = 1000, s = 2.5, nr = 9")
mu_time <- mu_RMSE # replicate the same data frame for fitting time


for (j in 1:12){ # setting number
  tRMSEs <- data.frame(AM[[j]][[3]]) # all rmse for the settings
  mu_RMSE[j, ] <- colMeans(tRMSEs) # calculate its mean
  tTime <- data.frame(AM[[j]][[4]]) # all fitting time for the settings
  mu_time[j, ] <- colMeans(tTime) # calculate its mean
  colnames(tRMSEs) <- c("bayesreg", "stan_gamm4", "bamlss", "jagam")
  colnames(tTime) <- c("bayesreg", "stan_gamm4", "bamlss", "jagam")
  tRMSEs <- melt(tRMSEs, measure.vars = colnames(tRMSEs)) # reshape for ggplot
  tTime <- melt(tTime, measure.vars = colnames(tTime))
  tRMSEs <- cbind(config = names(AM)[j], tRMSEs)
  tTime <- cbind(config = names(AM)[j], tTime)
  if (j == 1) {RMSE_AM <- tRMSEs
  Time_AM <- tTime} 
  else {RMSE_AM <- rbind(RMSE_AM, tRMSEs)
  Time_AM <- rbind(Time_AM, tTime)}
}

mu_RMSE_AM <- mu_RMSE
mu_time_AM <- mu_time

# print latex table of mean rmse and mean fitting time
xtable(mu_time_AM, type = latex, digits = c(0, 2, 2, 2, 2))
xtable(mu_RMSE_AM, type = latex, digits = c(0, 2, 2, 2, 2))

# Boxplots for RMSE and fitting time using ggplot
ggplot(RMSE_AM, aes(x = variable, y = log(value^2), fill = variable)) + 
  geom_boxplot() + 
  facet_wrap( ~ config, scales = "fixed", nrow = 3, ncol = 4) + 
  theme(plot.title = element_text(hjust = 0.45), axis.text.x = element_blank(), legend.position="bottom") +
  xlab("Model") + ylab("log MSE") + scale_fill_discrete(name = "Model")

ggplot(Time_AM, aes(x = variable, y = log(value), fill = variable)) + 
  geom_boxplot() + 
  facet_wrap( ~ config, scales = "fixed", nrow = 3, ncol = 4) + 
  theme(plot.title = element_text(hjust = 0.45), axis.text.x = element_blank(), legend.position="bottom") +
  xlab("Model") + ylab("log Time") + scale_fill_discrete(name = "Model")

# GAM and NGAM RMSEs Time Boxplots #############################################
names(GAM) <- c("01: n = 500, nr = 3", "02: n = 500, nr = 6",
                "03: n = 500, nr = 9",  "04: n = 1000 nr = 3",
                "05: n = 1000, nr = 6", "06: n = 1000, nr = 9")

mu_RMSE_GAM <- data.frame(matrix(nrow = 6, ncol = 4))
colnames(mu_RMSE_GAM) <- c("bayesreg", "stan_gamm4", "bamlss", "jagam")
rownames(mu_RMSE_GAM) <- c("01: n = 500, nr = 3", "02: n = 500, nr = 6", 
                           "03: n = 500, nr = 9",  "04: n = 1000 nr = 3",
                           "05: n = 1000, nr = 6", "06: n = 1000, nr = 9")
mu_time_GAM <- mu_RMSE_GAM

# Same process as AM and NAM
for (j in 1:6){ # setting number
  tRMSEs <- data.frame(GAM[[j]][[3]])
  mu_RMSE_GAM[j, ] <- colMeans(tRMSEs)
  tTime <- data.frame(GAM[[j]][[4]])
  mu_time_GAM[j, ] <- colMeans(tTime)
  colnames(tRMSEs) <- c("bayesreg", "stan_gamm4", "bamlss", "jagam")
  colnames(tTime) <- c("bayesreg", "stan_gamm4", "bamlss", "jagam")
  tRMSEs <- melt(tRMSEs, measure.vars = colnames(tRMSEs))
  tTime <- melt(tTime, measure.vars = colnames(tTime))
  tRMSEs <- cbind(config = names(GAM)[j], tRMSEs)
  tTime <- cbind(config = names(GAM)[j], tTime)
  if (j == 1) {RMSE_GAM <- tRMSEs
  Time_GAM <- tTime} 
  else {RMSE_GAM <- rbind(RMSE_GAM, tRMSEs)
  Time_GAM <- rbind(Time_GAM, tTime)}
}

xtable(mu_time_GAM, type = latex, digits = c(0, 2, 2, 2, 2))
xtable(mu_RMSE_GAM, type = latex, digits = c(0, 2, 2, 2, 2))

# Boxplots for RMSE and fitting time using ggplot
ggplot(RMSE_GAM, aes(x = variable, y = log(value^2), fill = variable)) + 
  geom_boxplot() +
  facet_wrap( ~ config, scales = "fixed", nrow = 2, ncol = 3) + 
  theme(plot.title = element_text(hjust = 0.45), axis.text.x = element_blank(), 
        legend.position="bottom") + xlab("Model") + ylab("log MSE") + 
  scale_fill_discrete(name = "Model")

ggplot(Time_GAM, aes(x = variable, y = log(value), fill = variable)) + 
  geom_boxplot() +
  facet_wrap( ~ config, scales = "fixed", nrow = 2, ncol = 3) + 
  theme(plot.title = element_text(hjust = 0.45), axis.text.x = element_blank(), 
        legend.position="bottom") + xlab("Model") + ylab("log Time") + 
  scale_fill_discrete(name = "Model")
################################################################################
# Plot of estimated effect against true functions

H <- c(h7, h1, h2, h8, h3, h4, h9, h5, h6)
names(H) <- c(7, 1, 2, 8, 3, 4, 9, 5, 6)
l <- GAM[[6]]     # Choose a setting
# l <- NGAM[[6]]  # Uncomment for non-uniform covariate case
n = 1000 # sample size
p = 10 # number of covariate
k = 10 # number of basis function
nr = 9 # number of relevant covariate
bi = 6 # Choose which basis function to use
uni = FALSE # uniform or non-uniform covariate
rho = 0.5 #  correlation with neighboring covariates

coefs <- l[[1]]
mu.coefs <- colMeans(coefs) # Calculate the mean estimates

# Generate covariates
if (uni) {X <- matrix(runif(n*p, min = domain[1], max = domain[2]),n,p)
} else {
  Sigma = matrix(0, p, p)
  for(g in 1:p){
    for(q in 1:p){
      Sigma[g,q] = rho^(abs(g-q))
    }
    
  }
  
  X = rmvnorm(n = n, sigma = Sigma)
  for (j in 1:p){
    X[, j] <- 2*(X[, j] - min(X[, j]))/(max(X[, j]) - min(X[, j])) - 1
  }}

B <- basis_functions[[bi]](X, k, constant = FALSE) # Create basis functions
F_hat <- matrix(0, n, p) # Estimated function
for(j in 1:k){
  F_hat[, j] <- B[, ((j-1)*k + 1):((j-1)*k + k)] %*% mu.coefs[((j-1)*k+2):(j*k+1)]
}


F <- matrix(0, n, p) # Matrix for true function values

for (j in 1:nr){
  h <- H[[j]]
  F[, j] <- h(X[, j])
}

par(mfrow = c(3, 3))
# plot first 9 functions of the covariate excluding covariate 10 
# which is always irrelevant in this simulation
for(j in (1:9)){
  q <- names(H)[j]
  ord <- order(X[, j])
  oX <- X[ord, j]
  plot(oX, F[ord, j] - mean(F[, j]), type = "l", col = "red", # True functions
       xlab = bquote(x[.(j)]), ylab = bquote('h'[.(q)](x[.(j)])))
  lines(oX, F_hat[ord, j] - mean(F_hat[, j]), col = "blue") # estimated functions
}
par(mfrow = c(1, 1))

# Application on real-world data ###############################################

# list of matrices to record the 100 replicates of Munich Rent and Medical Care models
mu.coefs <- matrix(nrow = 100, ncol =  21) # for bayesreg Munich Rent
coefs <- matrix(nrow = 100*5000, ncol = 21) # basis coefficients
perf <- matrix(nrow = 100, ncol = 5)
# to record the rmse and fitting time of each model and each replicate run
estimates <- list(mu.coefs, coefs, perf, perf)
names(estimates) <- c("Mean of coefficients", "Coefficients", "RMSEs", "Run Time") 
case <- replicate(n = 2, estimates, simplify = FALSE)
names(case) <- c("Munich Rent", "Medical Care") 
case[[2]][[1]] <- matrix(nrow = 100, ncol =  26) # for bayesreg Medical Care
case[[2]][[2]] <- matrix(nrow = 100*5000, ncol = 26)

# Munich Rent Index ############################################################
rent <- read.table("http://www.bamlss.org/misc/rent99.raw", header = TRUE)
# import data

set.seed(934)
for (q in 1:100){ # 100 replicates of fitting
  sample <- sample(c(TRUE, FALSE), nrow(rent), replace = TRUE, prob = c(0.8,0.2))
  train_rent  <- rent[sample, c("rent", "area", "yearc")] # train data
  test_rent   <- rent[!sample, c("rent", "area", "yearc")] # test data
  
  y_train_rent <- train_rent[, 1]
  X_train_rent <- train_rent[, -1]
  sX_rent <- scale(X_train_rent) # standardisation
  strain_rent <- data.frame(cbind(y_train_rent, sX_rent))
  colnames(strain_rent) <- c("rent", "area", "yearc")
  # create OP basis
  B_train_rent <- create.OP.basis(as.matrix(sX_rent), nbasis = 10, constant = FALSE)
  
  df_rent <- as.data.frame(cbind(y_train_rent, B_train_rent))
  names(df_rent) <- c("y", paste("V", 1:dim(B_train_rent)[2],sep = ""))
  B_train_rent <- df_rent[, -1]
  
  for (m in 1:5){ # 5 models
    st <- Sys.time()
    if (m == 1){ # Fit regression model using horseshoe prior and OP basis for 5,000 samples
      bayes_rent_ghs <- bayesreg(y ~., data = df_rent, model = "gaussian", 
                                 prior = "horseshoe", n.samples = 5e3, burnin = 1e3, thin = 1)}
    if (m == 2){ # lasso instead of horseshoe prior
      bayes_rent_lasso <- bayesreg(y ~., data = df_rent, model = "gaussian", 
                                   prior = "lasso", n.samples = 5e3, burnin = 1e3, thin = 1)}
    if (m == 3){ # Comparing to bayesreg without basis expansion
      bayes_rent_hslm <- bayesreg(rent ~ ., strain_rent, model = "gaussian", 
                                  prior = "horseshoe", n.samples = 5e3, burnin = 1e3, thin = 1)}
    if (m ==  4){ # BAMLSS
      bamlss_rent <- bamlss(rent ~ s(area, k = 10) + s(yearc, k = 10), data = strain_rent, 
                            family = "gaussian", n.samples = 5e3, burnin = 1e3, thin = 1)}
    if (m == 5){ # t distribution with bayesreg and horseshoe prior
      bayes_rent_ths <- bayesreg(y ~., data = df_rent, model = "t", 
                                 prior = "horseshoe", n.samples = 5e3, burnin = 1e3, thin = 1)}
    et <- Sys.time()
    case[[1]][[4]][q, m] <- as.numeric(et - st, unit = "mins")
  }
  
  y_test_rent <- test_rent[, 1]
  stest_rent <- test_rent[, -1] 
  # standardise test data based on mean and standard deviation from train data
  for (i in 1:(dim(stest_rent)[2])){
    mu <- attr(sX_rent, 'scaled:center')[i]
    s <- attr(sX_rent, 'scaled:scale')[i]
    stest_rent[, i] <- (stest_rent[, i] - mu)/s
  }
  
  # get OP basis evaluation of test data
  B_test_rent <- as.data.frame(get.OP.evaluation(as.matrix(sX_rent), nbasis = 10, 
                                                 as.matrix(stest_rent)))
  colnames(B_test_rent) <- c(paste("V", 1:dim(B_test_rent)[2],sep = ""))
  
  for (m in 1:5){ # prediction for test data
    if (m == 1){y_hat <- predict(bayes_rent_ghs, B_test_rent)}
    if (m == 2){y_hat <- predict(bayes_rent_lasso, B_test_rent)}
    if (m == 3){y_hat <- predict(bayes_rent_hslm, stest_rent)}
    if (m == 4){y_hat <- predict(bamlss_rent, stest_rent)$mu}
    if (m == 5){y_hat <- predict(bayes_rent_ths, B_test_rent)}
    
    case[[1]][[3]][q, m] <- RMSE(y_hat, y_test_rent) # save RMSE
    case[[1]][[1]][q, ] <- c(bayes_rent_ghs$mu.beta0, bayes_rent_ghs$mu.beta)
    # save the estimated coefficients
    case[[1]][[2]][(1+(q-1)*5000):(q*5000), ] <- cbind(t(bayes_rent_ghs$beta0), t(bayes_rent_ghs$beta))
  }
  print(paste0("Iteration ", q, " done."))}

# Plot estimated effect of 1st replicate as example
set.seed(934)
sample <- sample(c(TRUE, FALSE), nrow(rent), replace = TRUE, prob = c(0.8,0.2))
train_rent  <- rent[sample, c("rent", "area", "yearc")] # train data
test_rent   <- rent[!sample, c("rent", "area", "yearc")] # test data

y_train_rent <- train_rent[, 1]
X_train_rent <- train_rent[, -1]
sX_rent <- scale(X_train_rent)
strain_rent <- data.frame(cbind(y_train_rent, sX_rent))
colnames(strain_rent) <- c("rent", "area", "yearc")
B_train_rent <- create.OP.basis(as.matrix(sX_rent), nbasis = 10, constant = FALSE)

gamma_rent_ghs <- case[[1]][[1]][1, -1]
gamma_rent_ghs_all <- case[[1]][[2]][1:5000, -1]
F_hat_rent <- matrix(0, dim(X_train_rent)[1], dim(X_train_rent)[2])

k = 10
for(j in 1:2){
  F_hat_rent[, j] <- as.matrix(B_train_rent[, ((j-1)*k + 1):((j-1)*k + k)]) %*% 
    matrix(gamma_rent_ghs[((j-1)*k+1):(j*k)])
}
F_hat_rent_all <- list(replicate(n = 2, matrix(0, dim(X_train_rent)[1], dim(gamma_rent_ghs_all)[2]), simplify = FALSE))
for(j in 1:2){
  F_hat_rent_all[[j]] <- as.matrix(B_train_rent[, ((j-1)*k + 1):((j-1)*k + k)]) %*% 
    as.matrix(t(gamma_rent_ghs_all)[((j-1)*k+1):(j*k), ])
}

par(mfrow = c(2, 2))
# plot estimated effects for Bayesian regression using horseshoe prior and OP basis
for(j in (1:2)){
  # adding confidence intervals
  conf <- confband.kneib(t(F_hat_rent_all[[j]])) # simultaneous confidence band
  ord <- order(X_train_rent[, j])
  oX <- X_train_rent[ord, j]
  plot(oX, F_hat_rent[ord, j] - mean(F_hat_rent[, j]), type = "l", ylim = c(-250, 650), 
       xlab = colnames(train_rent)[1+j], ylab = "estimated effect for bayesreg")
  lines(oX, conf$upper[ord], lty = 2)
  lines(oX, conf$lower[ord], lty = 2)
}

bamlss_rent <- bamlss(rent ~ s(area, k = 10) + s(yearc, k = 10), data = strain_rent, 
                      family = "gaussian", n.samples = 5e3, burnin = 1e3, thin = 1)
# plot bamlss estimated effect
plot(bamlss_rent, ask = FALSE, ylim = c(-250, 650), ylab = "estimated effect for bamlss")
par(mfrow = c(1, 1))

# Spambase #####################################################################
data(spambase)
sp <- data.frame(matrix(nrow = 5, ncol = 6)) # matrix for spam performance
# time, accuracy, true negative, false negative, false positive, true positive
sp <- replicate(n = 5, sp, simplify = FALSE) # replicate for each model

set.seed(934)
for (q in 1:5){
  sample <- sample(c(TRUE, FALSE), nrow(spambase), replace = TRUE, prob = c(0.8,0.2))
  train_spam  <- spambase[sample, ]
  test_spam   <- spambase[!sample, ]
  
  y_train_spam <- as.factor(train_spam[, 1])
  X_train_spam <- train_spam[, -1]
  sX_spam <- scale(X_train_spam) # standardisation
  strain_spam <- cbind(y_train_spam, data.frame(sX_spam))
  # create OP basis
  B_train_spam <- create.OP.basis(as.matrix(sX_spam), nbasis = 8, constant = FALSE)
  
  df_spam <- cbind(y_train_spam, as.data.frame(B_train_spam))
  names(df_spam) <- c("y", paste("V", 1:dim(B_train_spam)[2],sep = ""))
  B_train_spam <- df_spam[, -1]
  
  # formula for bamlss
  names(strain_spam) <- c("is.spam", paste("V", 1:dim(train_spam[, -1])[2] , sep = ""))
  
  V <- paste0("s(V", 1:dim(train_spam[, -1])[2], ", k = 8)")
  V <- stri_paste(V, collapse = "+")
  formula <- as.formula(paste("is.spam", V, sep = "~"))
  
  for (m in 1:5){
    st <- Sys.time()
    if (m == 1){ # Fit a model using logistic horseshoe for 4,000 samples
    bayes_spam_hs <- bayesreg(as.factor(y) ~., data = df_spam, model = "logistic", 
                              prior = "horseshoe", n.samples = 5e3, burnin = 1e3, thin = 1)}
    if (m == 2){ # lasso instead of horseshoe prior
    bayes_spam_lasso <- bayesreg(as.factor(y) ~., data = df_spam, model = "logistic",
                               prior = "lasso", n.samples = 5e3, burnin = 1e3, thin = 1)}
    if (m == 3){ # bayes glm
    bayes_spam_hsglm <- bayesreg(as.factor(is.spam) ~ ., strain_spam, model = "logistic", 
                                 prior = "horseshoe", n.samples = 5e3, burnin = 1e3, thin = 1)}
    if(m == 4){ ## BAMLSS
      bamlss_spam <- bamlss(formula, data = strain_spam, family = "binomial", 
                          n.samples = 5e3, burnin = 1e3, thin = 1)}
    if (m == 5) { # logistic regression
      glm_spam <- glm(as.factor(is.spam) ~ ., data = strain_spam, family = "binomial")}
    et <- Sys.time()
    rt <- as.numeric(et - st, unit = "mins")
    sp[[m]][q, 1] <- rt
  }
  
  y_test_spam <- as.factor(test_spam[, 1])
  stest_spam <- test_spam[, -1]
  for (i in 1:(dim(stest_spam)[2])){ 
    # standardisation using mean and standard deviation from training data
    mu <- attr(sX_spam, 'scaled:center')[i]
    s <- attr(sX_spam, 'scaled:scale')[i]
    stest_spam[, i] <- (stest_spam[, i] - mu)/s
  }
  names(stest_spam) <- c(paste("V", 1:dim(train_spam[, -1])[2], sep = ""))
  
  # get evaluation of OP basis
  B_test_spam <- as.data.frame(get.OP.evaluation(as.matrix(sX_spam), nbasis = 8, 
                                                 as.matrix(stest_spam)))
  colnames(B_test_spam) <- c(paste("V", 1:dim(B_test_spam)[2],sep = ""))

  for (m in 1:5){ # Spambase prediction
    if (m == 1){prob_test <- predict(bayes_spam_hs, B_test_spam, type = "prob")}
    if (m == 2){prob_test <- predict(bayes_spam_lasso, B_test_spam, type = "prob")}
    if (m == 3){prob_test <- predict(bayes_spam_hsglm, stest_spam, type = "prob")}
    if (m == 4){prob_test <- predict(bamlss_spam, stest_spam, type = "parameter")}
    if (m == 5){prob_test <- predict(glm_spam, stest_spam, type = "response")}
    y_hat <- prob_test
    y_hat[prob_test < 0.5] <- 0 # cut-off value of 0.5
    y_hat[prob_test >= 0.5] <- 1
    y_hat <- as.factor(y_hat)
    
    # Confusion Matrix
    cma <- caret::confusionMatrix(y_hat, y_test_spam, positive = NULL, 
                                  dnn = c("Prediction", "Observation"))
    cm <- cma[2]$table
    sp[[m]][q, 3] <- cm[1, 1] # true negative
    sp[[m]][q, 4] <- cm[1, 2] # false negative
    sp[[m]][q, 5] <- cm[2, 1] # false positive
    sp[[m]][q, 6] <- cm[2, 2] # true positive
    acc <- cma[3]$overall[1] # accuracy
    sp[[m]][q, 2] <- cm[1, 1] <- acc}
  print(paste0("Iteration ", q, " done."))}

# Poisson distributed Medical Care Demand ######################################
data(DebTrivedi)
medcare <- DebTrivedi[, c(1, 6:8, 13, 15, 18)] # subset from Zeileis et.al (2008)
set.seed(934)
for (q in 1:100){
  sample <- sample(c(TRUE, FALSE), nrow(medcare), replace = TRUE, prob = c(0.8,0.2))
  train_med  <- medcare[sample, ]
  test_med   <- medcare[!sample, ]
  
  cont <- c(2, 4, 6) # numeric covariates
  y_train_med <- train_med[, 1]
  sX_med <- scale(train_med[, cont]) # standardisation
  strain_med <- cbind(train_med[, -cont], data.frame(sX_med))
  # create OP basis functions
  B_train_med <- create.OP.basis(as.matrix(sX_med), nbasis = 7, constant = FALSE)
  # All regressors (categorical and basis functions)
  df_med <- as.data.frame(cbind(train_med[, -cont], B_train_med))
  colnames(df_med)[5:dim(df_med)[2]] <- paste("V", 1:dim(B_train_med)[2],sep = "")
  
  for (m in 1:5){
    st <- Sys.time()
    if (m == 1){ # Fit a model using logistic horseshoe for 4,000 samples
      bayes_med_hs <- bayesreg(ofp ~., data = df_med, model = "poisson", 
                               prior = "hs", n.samples = 5e3, burnin = 9e3, thin = 1)}
    if (m == 2){ # lasso instead of horseshoe prior
      bayes_med_lasso <- bayesreg(ofp ~., data = df_med, model = "poisson", 
                                  prior = "lasso", n.samples = 5e3, burnin = 9e3, thin = 1)}
    # Comparing to bayesreg without basis expansion
    if (m == 3){ # bayes glm
      bayes_med_hsglm <- bayesreg(ofp ~., data = strain_med, model = "poisson", 
                                  prior = "horseshoe", n.samples = 5e3, burnin = 9e3, thin = 1)}
    if(m == 4){ # BAMLSS
      bamlss_med <- bamlss(ofp ~ s(hosp, k = 7) + s(numchron, k = 7) + 
                             s(school, k = 7) + health + gender +
                             privins, data = strain_med, family = "poisson", 
                           n.samples = 5e3, burnin = 9e3, thin = 1)}
    if (m == 5) { # Frequentist Poisson Regrerssion
      glm_med <- glm(ofp ~ ., data = strain_med, family = "poisson")}
    et <- Sys.time()
    case[[2]][[4]][q, m] <- as.numeric(et - st, unit = "mins")}
  
  y_test_med <- test_med[, 1]
  stest_med_cont <- test_med[, cont]
  for (i in 1:(dim(stest_med_cont)[2])){
    # Standardisation using mean and standard deviatoon from train data set
    mu <- attr(sX_med, 'scaled:center')[i]
    s <- attr(sX_med, 'scaled:scale')[i]
    stest_med_cont[, i] <- (stest_med_cont[, i] - mu)/s
  }
  stest_med <- cbind(test_med[, -cont], data.frame(stest_med_cont))
  
  # get OP basis evaluation
  B_test_med <- as.data.frame(get.OP.evaluation(as.matrix(sX_med), nbasis = 7, as.matrix(stest_med_cont)))
  # All test regressors (categorical and basis functions)
  XB_test_med <- as.data.frame(cbind(test_med[, -cont], B_test_med))
  colnames(XB_test_med)[5:dim(df_med)[2]] <- paste("V", 1:dim(B_test_med)[2],sep = "")
  
  for (m in 1:5){ # ofp prediction
    if (m == 1){y_hat <- predict(bayes_med_hs, XB_test_med, type = "response")}
    if (m == 2){y_hat <- predict(bayes_med_lasso, XB_test_med, type = "response")}
    if (m == 3){y_hat <- predict(bayes_med_hsglm, stest_med, type = "response")}
    if (m == 4){y_hat <- predict(bamlss_med, stest_med, type = "parameter")}
    if (m == 5){y_hat <- predict(glm_med, stest_med, type = "response")}
    
    RMSE_med <- RMSE(round(y_hat), y_test_med)
    case[[2]][[3]][q, m] <- RMSE_med
    case[[2]][[1]][q, ] <- c(bayes_med_hs$mu.beta0, bayes_med_hs$mu.beta)}
  case[[2]][[2]][(1+(q-1)*5000):(q*5000), ] <- cbind(t(bayes_med_hs$beta0), t(bayes_med_hs$beta))
  print(paste0("Iteration ", q, " done."))
}

# Plot for estimated effects for medical care demand data set
set.seed(934)
sample <- sample(c(TRUE, FALSE), nrow(medcare), replace = TRUE, prob = c(0.8,0.2))
train_med  <- medcare[sample, ]
test_med   <- medcare[!sample, ]

cont <- c(2, 4, 6)
sX_med <- scale(train_med[, cont])
strain_med <- cbind(train_med[, -cont], data.frame(sX_med))
B_train_med <- create.OP.basis(as.matrix(sX_med), nbasis = 7, constant = FALSE)

gamma_med <- case[[2]][[1]][1, -c(1:5)]
gamma_med_all <- case[[2]][[2]][1:5000, -c(1:5)]
F_hat_med <- matrix(0, dim(train_med)[1], length(gamma_med))
F_hat_med_all <- list(replicate(n = 3, matrix(0, dim(train_med)[1], 
                      dim(gamma_med_all)[2]), simplify = FALSE))
k = 7
for(j in 1:3){
  F_hat_med[, j] <- as.matrix(B_train_med[, ((j-1)*k + 1):((j-1)*k + k)]) %*% 
    matrix(gamma_med[((j-1)*k+1):(j*k)])
}

for(j in 1:3){
  F_hat_med_all[[j]] <- as.matrix(B_train_med[, ((j-1)*k + 1):((j-1)*k + k)]) %*% 
    as.matrix(t(gamma_med_all)[((j-1)*k+1):(j*k), ])
}

par(mfrow = c(2, 2))
for(j in (1:3)){
  #lower <- apply(F_hat_med_all[[j]], MARGIN = 1, FUN = quantile, probs = 0.025)
  #upper <- apply(F_hat_med_all[[j]], MARGIN = 1, FUN = quantile, probs = 0.975)
  if (j == 1) {conf <- confband.kneib(t(F_hat_med_all[[j]]))} 
  # simultaneous confidence band did not work for numchron and school
  else {conf <- confband.pw(t(F_hat_med_all[[j]]))} # pointwise method is used
  ord <- order(train_med[, cont][, j])
  oX <- train_med[, cont][ord, j]
  plot(oX, F_hat_med[ord, j] - mean(F_hat_med[, j]), type = "l", ylim = c(-4, 1.5),
       xlab = colnames(train_med[, cont])[j], ylab = "estimated effect on ofp")
  lines(oX, conf$lower[ord], lty = 2)
  lines(oX, conf$upper[ord], lty = 2)
}

# bamlss estimated effect plot
bamlss_med <- bamlss(ofp ~ s(hosp, k = 7) + s(numchron, k = 7) + 
                       s(school, k = 7) + health + gender +
                       privins, data = strain_med, family = "poisson", 
                     n.samples = 5e3, burnin = 9e3, thin = 1)
plot(bamlss_med, ask = FALSE, ylim = c(-3, 1))
par(mfrow = c(1, 1))

# Aggregate result #############################################################
app_performance <- data.frame(matrix(nrow = 5, ncol = 3)) 
# data frame for rmse and accuracy
colnames(app_performance) <- c("Munich Rent", "Spam", "Medical Care")
rownames(app_performance) <- c("bayesreg with horseshoe and OP basis", 
                               "bayesreg with lasso and OP basis", 
                               "bayesreg with  horseshoe and linear predictor", 
                               "bamlss", "alternative")
app_runtime <- app_performance # copy for the same structure and names


mr <- case[[1]] # Munich Rent
mc <- case[[2]] # Medical Care
app_performance[, 1]<- colMeans(mr[[3]]) # Mean RMSE of Munich Rent Index case
app_runtime[, 1] <- colMeans(mr[[4]]) # Mean fitting time of Munich Rent Index case

app_performance[, 3]<- colMeans(mc[[3]]) # Mean RMSE of Medical Care Demand  case
app_runtime[, 3] <- colMeans(mc[[4]]) # MEan fitting time of Medical Care Demand case

cm <- data.frame(matrix(nrow = 2, ncol = 2)) # Average Confusion Matrix
rownames(cm) <- c("y_hat = 0", "y_hat = 1")
colnames(cm) <- c("y = 0", "y = 1")
cm_all <- replicate(n = 5, cm, simplify = FALSE)
ntest <- rowSums(sp[[1]][, 3:6]) # test sample size for each replicates

for (m in 1:5){
  p <- sp[[m]] # Performance of model m
  app_runtime[m, 2] <- mean(p[, 1])
  app_performance[m, 2] <- mean(p[, 2])
  cm <- sp[[m]][, 3:6]/ntest # percentage of true negative, false negative, 
  # false positive, true positive in prediction of test data for confusion matrix
  cm_all[[m]][1, 1] <- mean(cm[, 1])*100
  cm_all[[m]][1, 2] <- mean(cm[, 2])*100
  cm_all[[m]][2, 1] <- mean(cm[, 3])*100
  cm_all[[m]][2, 2] <- mean(cm[, 4])*100
}

models <- c("bayesreg GAM with horseshoe", "bayesreg GAM with lasso", 
            "bayesreg GLM with horseshoe", "bamlss", "GLM")
# latex table for mean performance and fitting time
xtable(app_performance, type = latex, digits = c(0, 4, 4, 4))
xtable(app_runtime, type = latex, digits = c(0, 4, 4, 4))
for (m in 1:5){ # add average confusion matrix entried for each model
  print(paste0("Confusion Matrix in % for ", models[m], ":"))
 print(xtable(cm_all[[m]], type = latex, digits = c(0, 4, 4)))
}


outlier <- which(mc[[3]][, 1]>100) # drop the case where bayesreg with OP has 
# very large prediction to see how it performs without the case
mc_rmse99 <- mc[[3]][-outlier, ]
mc_time99 <- mc[[4]][-outliers, ]
app_performance[, 3]<- colMeans(mc_rmse99) # the average without outlier case
app_runtime[, 3] <- colMeans(mc_time99) # # the average without outlier case

# updated latex table with average from 99 replicates
xtable(app_performance, type = latex, digits = c(0, 4, 4, 4))
xtable(app_time, type = latex, digits = c(0, 4, 4, 4))

# For Boxplot of RMSE and accuracy of application data #########################
mr_rmse <- data.frame(mr[[3]])
colnames(mr_rmse) <- c("bayesreg HOP", "bayesreg LOP", "bayesreg HL", "bamlss", "alternative")
app_perf <- cbind("1. Munich Rent", melt(mr_rmse, measure.vars = colnames(mr_rmse)))
colnames(app_perf) = c("Data", "Model", "Performance")

sp_acc <- data.frame(sp[[1]][, 2])
for (m in 2:5){
  sp_acc <- cbind(sp_acc, sp[[m]][, 2])
}
colnames(sp_acc) <- c("bayesreg HOP", "bayesreg LOP", "bayesreg HL", "bamlss", "alternative")
sp_box <- cbind("2. Spambase", melt(sp_acc, measure.vars = colnames(sp_acc)))
colnames(sp_box) <- c("Data", "Model", "Performance")
app_perf <- rbind(app_perf, sp_box)

mc_rmse99 <- data.frame(mc_rmse99)
colnames(mc_rmse99) <- c("bayesreg HOP", "bayesreg LOP", "bayesreg HL", "bamlss", "alternative")
mc_box <- cbind("3. Medical Care", melt(mc_rmse99, measure.vars = colnames(mc_rmse99)))
colnames(mc_box) <- c("Data", "Model", "Performance")
app_perf <- rbind(app_perf, mc_box)

ggplot(app_perf, aes(x = Model, y = Performance, fill = Model)) + 
  geom_boxplot() + 
  facet_wrap( ~ Data, scales = "free", nrow = 1, ncol = 3) + 
  theme(plot.title = element_text(hjust = 0.45), axis.text.x = element_blank(), 
        legend.position="bottom") + scale_fill_discrete(name = "Model") + 
  xlab("") + ylab("Performance")

# Boxplot for fitting time of application data ##################################
mr_time <- data.frame(mr[[4]])
colnames(mr_time) <- c("bayesreg HOP", "bayesreg LOP", "bayesreg HL", "bamlss", "alternative")
app_time <- cbind("1. Munich Rent", melt(mr_time, measure.vars = colnames(mr_time)))
colnames(app_time) = c("Data", "Model", "Time")

sp_time <- data.frame(sp[[1]][, 1])
for (m in 2:5){
  sp_time <- cbind(sp_time, sp[[m]][, 1])
}
colnames(sp_time) <- c("bayesreg HOP", "bayesreg LOP", "bayesreg HL", "bamlss", "alternative")
sp_box <- cbind("2. Spambase", melt(sp_time, measure.vars = colnames(sp_time)))
colnames(sp_box) <- c("Data", "Model", "Time")
app_time <- rbind(app_time, sp_box)

mc_time99 <- data.frame(mc_time99)
colnames(mc_time99) <- c("bayesreg HOP", "bayesreg LOP", "bayesreg HL", "bamlss", "alternative")
mc_box <- cbind("3. Medical Care", melt(mc_time99, measure.vars = colnames(mc_time99)))
colnames(mc_box) <- c("Data", "Model", "Time")
app_time <- rbind(app_time, mc_box)

ggplot(app_time, aes(x = Model, y = Time, fill = Model)) + 
  geom_boxplot() + 
  facet_wrap( ~ Data, scales = "free", nrow = 1, ncol = 3) + 
  theme(plot.title = element_text(hjust = 0.45), axis.text.x = element_blank(), 
        legend.position="bottom") + scale_fill_discrete(name = "Model") + 
  xlab("") + ylab("Fitting Time")