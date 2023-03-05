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
################################################################################
# Scenarios 
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

norm_case <- function(h, mean = 0, s = 1){
  # adding a normally distributed iid errors to standardized values
  
  sc <- scale(h) #  standardization
  
  e <- rnorm(length(h), mean, s)
  e <- e*attr(sc, 'scaled:scale')
  # scale e based on standard deviation of h
  y <- h + e
  
  return(y)}

plot.design.matrix <- function(x,X){
  p <- dim(X)[2]
  ord <- order(x)
  plot(x[ord], X[ord, 1], type = "l", xlab = "x", ylab = "Basis", 
       ylim = c(min(X), max(X)), col = 1)
  for(j in 2:p){
    lines(x[ord], X[ord, j], col = j)
  }
}

RMSE <- function(y_hat, y){
  as.numeric(sqrt(mean((y_hat - y)^2)))
}
## set up design ###############################################################
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
# B-spline

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
    MMB[[i]] <- cbind(scale(x[, i]), MM)} # attach standardized linear function # 
  
  
  if (constant) {
    MMfull <- (cbind(rep(1, dim(MMB[[1]])[1]), MMB[[1]]))}
  if(!constant) {
    MMfull <- MMB[[1]]}
  
  if (p > 1){
    for(j in 2:p){
      MMfull <- cbind(MMfull, MMB[[j]])}} ## this is the overall design matrix
  return(MMfull)}
################################################################################
# Demmler-Reinsch basis based on Generalized Eigenvalue Problem

create.DR.basis <- function(x, nbasis = 10, rescale = FALSE, constant = FALSE) 
{
  
  k = nbasis + 1 ### => we will remove the constant basis function later on
  
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
    DRB[[i]] <- cbind(scale(x[, i]), DR)} # attach standardized linear function # 
  
  if (constant) {
    DRfull <- (cbind(rep(1, dim(DRB[[1]])[1]), DRB[[1]]))}
  if(!constant) {
    DRfull <- DRB[[1]]}
  if (p > 1){
    for(j in 2:p){
      DRfull <- cbind(DRfull, DRB[[j]])}} ## this is the overall design matrix
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
      if (nbasis %% 2 == 1){FB <- cbind(FB, eval.basis(x[, i], basisobj = fb)[, c(-1, -(nbasis + 2))])}
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
      
      GP_new <- cbind(rep(1, dim(x_new)[1]), poly(x_new[, i], degree = nbasis, raw = TRUE)) 
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
################################################################################
# Estimation configuration for Nonparametric Model
domain = c(-1, 1)
n = 100
set.seed(934)
k = 10
uni = FALSE # if the covariates are uniformly distributed
s = 2.5

# H <- c(h1, h2, h3, h4, h5, h6, h7, h8, h9)
# 
# mu.coefs <- matrix(nrow = 100, ncol = k + 1)
# mus <- replicate(n = 7, mu.coefs, simplify = FALSE)
# RMSEs <- matrix(nrow = 100, ncol = 7)
# 
# l <- list(mus, RMSEs)
# names(l) <- c("Mean of coeeficients", "RMSEs")
# N <- replicate(n = 9, l, simplify = FALSE)
# 
# names(N) <- paste0("h", (1:9))

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

saveRDS(N, "NNonparametric100n0.4s.RData")
################################################################################

create.additive <- function(F, s = 1, logistic = FALSE){ # (logistic) additive model with normal distributed error
  
  n <- dim(F)[1]
  f <- rowSums(F)
  
  sc <- scale(f) #  standardization
  
  e <- rnorm(n, sd = s)
  e <- e*attr(sc, 'scaled:scale')
  # scale e based on standard deviation of h
  
  if (!logistic) {
    y <- f + e}
  
  if (logistic) {
    prob <- exp(f)/(exp(f) + 1) # logistic additive model, do we add error terms here
    y <- rbinom(n, size = 1, prob = prob)}
  
  return(y)
}

# Small Simulation Study to compare MM, DR and OP in (G)AM #####################
n = 500
p = 10 # number of covariates
rho = 0.5 # correlation between covariate 1 and 2
k = 10 ## number of basis functions for each additive component
nr = 6 # number of relevant covariates
uni = FALSE
s = 2.5

RMSEs <- matrix(nrow = 100, ncol = 8) # for all 4 models

H <- c(h7, h1, h2, h8, h3, h4, h9, h5, h6)

settings <- list(c(3, FALSE), c(4, FALSE), c(6, FALSE), c(7, FALSE),
                 c(3, TRUE), c(4, TRUE), c(6, TRUE), c(7, TRUE))
names(settings) <- c("AM MM", "AM DR", "AM OP", "AM LP", "GAM MM", "GAM DR", "GAM OP", "GAM LP")

# Begin for the loop of each configuration for fixed set of sample size, 
# scale of Var(f) as standard deviation and number of relevant covariates
set.seed(934)
sl <- Sys.time()
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
el <- Sys.time()
print(el - sl)

saveRDS(RMSEs, "NRMSEs500n2.5s6nr.RData")
################################################################################
#(Generalized Additive Model)
jag.file <- paste(getwd(),"/test.jags",sep = "")

# Configuration for (Generalized) Additive Model
domain = c(-1, 1)
n = 500
p = 10 # number of covariates
rho = 0.5 # correlation between covariate 1 and 2
k = 10 ## number of basis functions for each additive component
nr = 9 # number of relevant covariates
uni = FALSE # is the covariate uniformly distributed between [-1, 1]
bin = FALSE
s = 1 # note that s does not affect GAM simulation
bi = 6 # which basis function

mu.coefs <- matrix(nrow = 100, ncol = p*k + 1) # for bayesreg
eta_hat <- matrix(nrow = 100, ncol = n)
eta_hat <- replicate(n = 4, eta_hat, simplify = FALSE) # for each models
RMSEs <- matrix(nrow = 100, ncol = 4) # for all 4 models
time <- matrix(nrow = 100, ncol = 4) # for all 4 models

l <- list(mu.coefs, eta_hat, RMSEs, time)
names(l) <- c("Mean of coeeficients", "Fitted values", "RMSEs", "Run time")

H <- c(h7, h1, h2, h8, h3, h4, h9, h5, h6)

# Begin for the loop of each configuration for fixed set of sample size, 
# scale of Var(f) as standard deviation and number of relevant covariates
set.seed(934)
sl <- Sys.time()
for (i in 1:100){ # 100 replicates of simulation for each model
  
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
                                         prior = "hs", n.samples = 4e3)
      } else if (bin == TRUE) {
        fit <- bayesreg(as.factor(y) ~., data = df, model = "logistic", 
                        prior = "hs", n.samples = 4e3)}
    } else if (m == 2){ # stan_gamm4
      if (bin == FALSE) {fit <- stan_gamm4(formula, data = data)
      } else if (bin == TRUE) {
        fit <- stan_gamm4(formula, data = data, family = "binomial")}
    } else if (m == 3){ # bamlss
      if (bin == FALSE) {fit <- bamlss(formula, data = data)
      } else if (bin == TRUE) {
        fit <- bamlss(formula, data = data, family = "binomial", n.samples = 4e3)}
    } else if (m == 4){ # jagam
      if (bin == FALSE){
        jag <- jagam(formula, file = jag.file, data = data)
      } else if (bin == TRUE){
        jag <- jagam(formula, family = "binomial", file = jag.file, data = data)
      }
      jm <- jags.model("test.jags", data = jag$jags.data,
                       inits = jag$jags.ini)
      sam <- jags.samples(jm, c("b", "rho"), n.iter = 4e3)
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
    runtime <- end - start
    l[[4]][i, m] <- runtime}
  
  if (i %% 5 == 0) {print(paste("Iteration ", i, " complete"))}
}
el <- Sys.time()
print(el - sl)

saveRDS(l, "NAMOP5001s9nr.RData")

# Munich Rent ##################################################################
app_performance <- data.frame(matrix(nrow = 5, ncol = 3))
colnames(app_performance) <- c("Munich Rent", "Spam", "Medical Care")
rownames(app_performance) <- c("bayesreg with horseshoe and OP basis", 
                               "bayesreg with lasso and OP basis", 
                               "bayesreg with  horseshoe and linear predictor", 
                               "bamlss", "alternative")
runtime <- data.frame(matrix(nrow = 5, ncol = 3))

rent <- read.table("http://www.bamlss.org/misc/rent99.raw", header = TRUE)
set.seed(934)
sample <- sample(c(TRUE, FALSE), nrow(rent), replace = TRUE, prob = c(0.8,0.2))
train_rent  <- rent[sample, c("rent", "area", "yearc")]
test_rent   <- rent[!sample, c("rent", "area", "yearc")]

y_train_rent <- train_rent[, 1]
X_train_rent <- train_rent[, -1]
sX_rent <- scale(X_train_rent)
strain_rent <- data.frame(cbind(y_train_rent, sX_rent))
colnames(strain_rent) <- c("rent", "area", "yearc")
B_train_rent <- create.OP.basis(as.matrix(sX_rent), nbasis = 10, constant = FALSE)

df_rent <- as.data.frame(cbind(y_train_rent, B_train_rent))
names(df_rent) <- c("y", paste("V", 1:dim(B_train_rent)[2],sep = ""))
B_train_rent <- df_rent[, -1]

for (m in 1:5){
  set.seed(934)
  if (m == 1){ # Fit regression model using horseshoe prior and OP basis for 4,000 samples
    bayes_rent_ghs <- bayesreg(y ~., data = df_rent, model = "gaussian", 
                                 prior = "horseshoe", n.samples = 12e3)}
  if (m == 2){ # lasso instead of horseshoe prior
    bayes_rent_lasso <- bayesreg(y ~., data = df_rent, model = "gaussian", 
                                 prior = "lasso", n.samples = 12e3)}
  # Comparing to bayesreg without basis expansion
  if (m == 3){
    bayes_rent_hslm <- bayesreg(rent ~ ., strain_rent, model = "gaussian", 
                                 prior = "horseshoe", n.samples = 12e3)}
  if (m ==  4){ # Comparing to bayesreg without basis expansion
    bamlss_rent <- bamlss(rent ~ s(area, k = 10) + s(yearc, k = 10), 
                          data = strain_rent, family = "gaussian", n.samples = 12e3)}
  if (m == 5){
    bayes_rent_ths <- bayesreg(y ~., data = df_rent, model = "t", 
                                 prior = "horseshoe", n.samples = 12e3)}
  }

y_test_rent <- test_rent[, 1]
stest_rent <- test_rent[, -1]
for (i in 1:(dim(stest_rent)[2])){
  mu <- attr(sX_rent, 'scaled:center')[i]
  s <- attr(sX_rent, 'scaled:scale')[i]
  stest_rent[, i] <- (stest_rent[, i] - mu)/s
}

B_test_rent <- as.data.frame(get.OP.evaluation(as.matrix(sX_rent), nbasis = 10, 
                                               as.matrix(stest_rent)))
colnames(B_test_rent) <- c(paste("V", 1:dim(B_test_rent)[2],sep = ""))

models <- c("bayesreg AM with horseshoe", "bayesreg AM with lasso", "bayesreg LM", 
            "bamlss", "bayesreg t-distributed AM")
for (m in 1:5){
  if (m == 1){y_hat <- predict(bayes_rent_ghs, B_test_rent)}
  if (m == 2){y_hat <- predict(bayes_rent_lasso, B_test_rent)}
  if (m == 3){y_hat <- predict(bayes_rent_hslm, stest_rent)}
  if (m == 4){y_hat <- predict(bamlss_rent, stest_rent)$mu}
  if (m == 5){y_hat <- predict(bayes_rent_ths, B_test_rent)}
  
  RMSE_rent <- RMSE(y_hat, y_test_rent)
  app_performance[m, 1] <- RMSE_rent
  print(paste0("The RMSE for ", models[m], " is: ", RMSE_rent))
}

gamma_rent_ghs <- bayes_rent_ghs$mu.beta
F_hat_rent <- matrix(0, dim(X_train_rent)[1], dim(X_train_rent)[2])
for(j in 1:2){
  F_hat_rent[, j] <- as.matrix(B_train_rent[, ((j-1)*k + 1):((j-1)*k + k)]) %*% 
    matrix(gamma_rent_ghs[((j-1)*k+1):(j*k)])
  
}


par(mfrow = c(1, 2))
for(j in (1:2)){
  ord <- order(X_train_rent[, j])
  oX <- X_train_rent[ord, j]
  plot(oX, F_hat_rent[ord, j] - mean(F_hat_rent[, j]), type = "l", ylim = c(-250, 650), 
       xlab = colnames(train_rent)[1+j], ylab = "estimated effect on rent")
}

#pdf("bamlssrent.pdf", width = 7, height = 3.5)
plot(bamlss_rent, ask = FALSE, ylim = c(-250, 650), ylab = "estimated effect on rent")
par(mfrow = c(1, 1))
#dev.off()
# spam #########################################################################
data(spambase)
set.seed(934)
sample <- sample(c(TRUE, FALSE), nrow(spambase), replace = TRUE, prob = c(0.8,0.2))
train_spam  <- spambase[sample, ]
test_spam   <- spambase[!sample, ]

y_train_spam <- as.factor(train_spam[, 1])
X_train_spam <- train_spam[, -1]
sX_spam <- scale(X_train_spam)
strain_spam <- cbind(y_train_spam, data.frame(sX_spam))
B_train_spam <- create.OP.basis(as.matrix(sX_spam), nbasis = 9, constant = FALSE)

df_spam <- cbind(y_train_spam, as.data.frame(B_train_spam))
names(df_spam) <- c("y", paste("V", 1:dim(B_train_spam)[2],sep = ""))
B_train_spam <- df_spam[, -1]

# formula for bamlss
names(strain_spam) <- c("is.spam", paste("V", 1:dim(train_spam[, -1])[2] , sep = ""))

V <- paste0("s(V", 1:dim(train_spam[, -1])[2], ", k = 9)")
V <- stri_paste(V, collapse = "+")
formula <- as.formula(paste("is.spam", V, sep = "~"))

for (m in 1:5){
  set.seed(934)
  st <- Sys.time()
  if (m == 1){ # Fit a model using logistic horseshoe for 4,000 samples
  bayes_spam_hs <- bayesreg(as.factor(y) ~., data = df_spam, model = "logistic", 
                            prior = "horseshoe", n.samples = 12e3)}
  if (m == 2){ # lasso instead of horseshoe prior
  bayes_spam_lasso <- bayesreg(as.factor(y) ~., data = df_spam, model = "logistic",
                             prior = "lasso", n.samples = 12e3)}
  if (m == 3){ # bayes glm
  bayes_spam_hsglm <- bayesreg(as.factor(is.spam) ~ ., strain_spam, model = "logistic", 
                               prior = "horseshoe", n.samples = 12e3)}
  if(m == 4){
    bamlss_spam <- bamlss(formula, data = strain_spam, family = "binomial", 
                        n.samples = 12e3)}
  if (m == 5) {
    glm_spam <- glm(as.factor(is.spam) ~ ., data = strain_spam, family = "binomial")}
  et <- Sys.time()
  rt <- et - st
  runtime[m, 2] <- rt
  print(rt)}

#bayes_spam_hs <- readRDS("bayes_spam_hs.RData")
#bayes_spam_lasso<- readRDS("bayes_spam_lasso.RData")
#bayes_spam_ridge <- readRDS("bayes_spam_ridge.RData")
y_test_spam <- as.factor(test_spam[, 1])
stest_spam <- test_spam[, -1]
for (i in 1:(dim(stest_spam)[2])){
  mu <- attr(sX_spam, 'scaled:center')[i]
  s <- attr(sX_spam, 'scaled:scale')[i]
  stest_spam[, i] <- (stest_spam[, i] - mu)/s
}
names(stest_spam) <- c(paste("V", 1:dim(train_spam[, -1])[2] , sep = ""))

B_test_spam <- as.data.frame(get.OP.evaluation(as.matrix(sX_spam), nbasis = 9, 
                                               as.matrix(stest_spam)))
colnames(B_test_spam) <- c(paste("V", 1:dim(B_test_spam)[2],sep = ""))


# Check how well did our predictions did by generating confusion matrix
models <- c("bayesreg GAM with horseshoe", "bayesreg GAM with lasso", 
            "bayesreg GLM with horseshoe", "bamlss GAM", "bayes GAM with ridge")
for (m in 1:5){
  if (m == 1){prob_test <- predict(bayes_spam_hs, B_test_spam, type = "prob")}
  if (m == 2){prob_test <- predict(bayes_spam_lasso, B_test_spam, type = "prob")}
  if (m == 3){prob_test <- predict(bayes_spam_hsglm, stest_spam, type = "prob")}
  if (m == 4){prob_test <- predict(bamlss_spam, stest_spam, type = "parameter")}
  if (m == 5){prob_test <- predict(bayes_spam_ridge, B_test_spam, type = "prob")}
  y_hat <- prob_test
  y_hat[prob_test < 0.5] <- 0
  y_hat[prob_test >= 0.5] <- 1
  y_hat <- as.factor(y_hat)
  
  cma <- caret::confusionMatrix(y_hat, y_test_spam, positive = NULL, 
                                dnn = c("Prediction", "Observation"))
  cm <- cma[2]$table
  #print(cm)
  print(xtable(cm))
  acc <- cma[3]$overall[1]
  app_performance[m, 2] <- acc
  print(paste0("The accuracy for ", models[m], " is: ", acc))
}

# Poisson distributed Medical Care Demand ######################################
data(DebTrivedi)
medcare <- DebTrivedi[, c(1, 6:8, 11, 13, 15, 16, 18)]
set.seed(934)
sample <- sample(c(TRUE, FALSE), nrow(medcare), replace = TRUE, prob = c(0.8,0.2))
train_med  <- medcare[sample, ]
test_med   <- medcare[!sample, ]

cont <- c(2, 4, 5, 7, 8) # continous variable index
y_train_med <- train_med[, 1]
sX_med <- scale(train_med[, cont])
strain_med <- cbind(train_med[, -cont], data.frame(sX_med))
B_train_med <- create.OP.basis(as.matrix(sX_med), nbasis = 8, constant = FALSE)

df_med <- as.data.frame(cbind(train_med[, -cont], B_train_med))
colnames(df_med)[5:dim(df_med)[2]] <- paste("V", 1:dim(B_train_med)[2],sep = "")

for (m in 1:5){
  set.seed(934)
  if (m == 1){ # Fit a model using logistic horseshoe for 4,000 samples
    bayes_med_hs <- bayesreg(ofp ~., data = df_med, model = "poisson", 
                                         prior = "hs", n.samples = 12e3, burnin = 9e3)}
  if (m == 2){ # lasso instead of horseshoe prior
    bayes_med_lasso <- bayesreg(ofp ~., data = df_med, model = "poisson", 
                                prior = "lasso", n.samples = 12e3, burnin = 9e3)}
  # Comparing to bayesreg without basis expansion
  if (m == 3){ # bayes glm
    bayes_med_hsglm <- bayesreg(ofp ~., data = strain_med, model = "poisson", 
                                prior = "horseshoe", n.samples = 12e3, burnin = 9e3)}
  if(m == 4){
    bamlss_med <- bamlss(ofp ~ s(hosp, k = 8) + s(numchron, k = 8) + s(age, k = 8) + 
                           s(school, k = 8) + s(faminc, k = 8) + health + gender + 
                           privins, data = strain_med, family = "poisson", 
                         n.samples = 12e3, burnin = 9e3)}
  if (m == 5) {glm_med <- glm(ofp ~ ., data = strain_med, family = "poisson")}}

y_test_med <- test_med[, 1]
stest_med_cont <- test_med[, cont]
for (i in 1:(dim(stest_med_cont)[2])){
  mu <- attr(sX_med, 'scaled:center')[i]
  s <- attr(sX_med, 'scaled:scale')[i]
  stest_med_cont[, i] <- (stest_med_cont[, i] - mu)/s
}
stest_med <- cbind(test_med[, -cont], data.frame(stest_med_cont))

B_test_med <- as.data.frame(get.OP.evaluation(as.matrix(sX_med), nbasis = 8, as.matrix(stest_med_cont)))
XB_test_med <- as.data.frame(cbind(test_med[, -cont], B_test_med))
colnames(XB_test_med)[5:dim(df_med)[2]] <- paste("V", 1:dim(B_test_med)[2],sep = "")

models <- c("bayesreg GAM with horseshoe", "bayesreg GAM with lasso", "bayesreg GLM", 
            "bamlss", "GLM")
for (m in 1:5){
  if (m == 1){y_hat <- predict(bayes_med_hs, XB_test_med, type = "response")}
  if (m == 2){y_hat <- predict(bayes_med_lasso, XB_test_med, type = "response")}
  if (m == 3){y_hat <- predict(bayes_med_hsglm, stest_med, type = "response")}
  if (m == 4){y_hat <- predict(bamlss_med, stest_med, type = "parameter")}
  if (m == 5){y_hat <- predict(glm_med, stest_med, type = "response")}
  
  RMSE_med <- RMSE(round(y_hat), y_test_med)
  app_performance[m, 3] <- RMSE_med
  print(paste0("The RMSE for ", models[m], " is: ", RMSE_med))
}

gamma_med <- bayes_med_hs$mu.beta[-c(1:4)]
F_hat_med <- matrix(0, dim(train_med)[1], length(gamma_med))
k = 8
for(j in 1:5){
  F_hat_med[, j] <- as.matrix(B_train_med[, ((j-1)*k + 1):((j-1)*k + k)]) %*% 
    matrix(gamma_med[((j-1)*k+1):(j*k)])
  
}

par(mfrow = c(3, 2))
for(j in (1:5)){
  ord <- order(train_med[, cont][, j])
  oX <- train_med[, cont][ord, j]
  plot(oX, F_hat_med[ord, j] - mean(F_hat_med[, j]), type = "l", ylim = c(-3, 1),
       xlab = colnames(train_med[, cont])[j], ylab = "ofp") # main = paste("effect of", colnames(train_med[, cont])[j], "on ofp")
}

xtable(app_performance, digits = c(0, 4, 4, 4))

plot(bamlss_med, ask = FALSE, ylim = c(-3, 1))
par(mfrow = c(1, 1))
