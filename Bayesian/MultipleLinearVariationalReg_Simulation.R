###Variational Regression Parameters Estimation Function###
#Package dependence for the experiment
library(matrixcalc)    # Calculation of matrix
library(ggplot2)       # Better plotting
library(dplyr)         # Data transformation
library(tidyr)         # Tidying up data
library(readxl)        # Import excel files

# Illustration of not enough samples may lead to Overfitting 
post_grd = read_excel("posting_grade.xlsx")
set.seed(210101)
samp_post_id = sample(post_grd$grade, 20)
samp_post = post_grd[is.element(post_grd$grade, samp_post_id),]
outsamp_post = post_grd[!(post_grd$grade %in% samp_post_id),]

ggplot(data=samp_post, aes(x=posting, y=grade)) + 
  geom_point(color='black')
ggplot(data=samp_post, aes(x=posting, y=grade)) + 
  geom_point(data=outsamp_post, aes(x=posting, y=grade), color='grey') + 
  geom_point(color='black')


#The first part of the experiment is to build the estimation function.
#For MLE method, default function lm() from stats package is used.
#For variational estimation method, we refer to Kapourani (2018): https://rpubs.com/cakapourani/variational-bayes-lr for the function.

###Building Function to estimate regression coefficient using variational parameters###
#Input for the function is a design matrix X for linear regression model and vector of dependent variable data y
#The function below depend on a package: matrixcalc
#In this function a0 and b0 for variaitonal distribution is initiated at 0.01
#Stopping criteria is set at either maximum iteration of 10,000 or converges when lower bound < epsilon = 0.000001.
vlr = function(X, y, lambda = 1, a_0 = .01, b_0 = .01, 
               max_iter = 10000, epsilon_conv = .000001, 
               is_verbose = FALSE)
  #The function takes input of a design matrix X and vector of target variable y.
  #Initialization for variational parameters are a0 = 0.01 and b0 = 0.01.
  #For the variance of error term, precision term lambda is used as the representation
  #Lambda value may be changed from the inputation
  #The number of maximum iteration is 10000 with the tolerance value for convergence in lower bound is 0.000001
{
  X = as.matrix(X)
  c = NCOL(X)              # Number of covariates
  n = NROW(X)              # Number of observations
  XX = crossprod(X, X)     # Compute X'X
  Xy = crossprod(X, y)     # Compute X'y
  yy = c(crossprod(y))     # Compute y'y
  
  E_a = a_0 / b_0          # Expectation of gamma distribution part in variational distribution
  a = a_0 + c / 2          # Compute alpha parameter of the gamma distribution
  L = rep(-Inf, max_iter)  # Store the lower bounds
  # Iterate to find the value of optimal variational parameters
  for (i in 2:max_iter) {
    # Covariance matrix of Normal factor from variational distirbution
    S = solve(diag(E_a, c) + lambda * XX)
    # Mean of Normal factor from variational distirbution
    m = lambda * S %*% Xy
    # Expectation of E[w'w]
    E_ww = as.vector(crossprod(m) + matrix.trace(S))
    # b_N parameter of Gamma factor from variational distirbution
    b = b_0 + 0.5 * E_ww
    # Compute expectation of E[a]
    E_a = a / b
    
    # Compute lower bound for convergence criteria
    lb_py = 0.5 * n * log(lambda / (2 * pi)) - 0.5 * lambda * yy + lambda * crossprod(m, Xy) - 0.5 * lambda * matrix.trace(XX %*% (tcrossprod(m, m) + S)) 
    lb_pw = -0.5 * c * log(2*pi) + 0.5 * c * (digamma(a) - log(b)) - 0.5 * E_a * E_ww 
    lb_pa = a_0*log(b_0) + (a_0 - 1) * (digamma(a) - log(b)) - b_0*E_a - log(gamma(a_0))
    lb_qw = -0.5 * log(det(S)) - 0.5 * c * (1 + log(2*pi)) 
    lb_qa = -lgamma(a) + (a - 1)*digamma(a) + log(b) - a
    
    L[i] = lb_py + lb_pw + lb_pa - lb_qw - lb_qa # Compute the lower bound value
    
    # Gives caution for output of the function
    # Show VB difference
    if (is_verbose) { cat("It:\t",i,"\tLB:\t", L[i],"\tLB_diff:\t", 
                          L[i] - L[i - 1],"\n")}
    # Check if lower bound decreases
    if (L[i] < L[i - 1]) { message("Lower bound decreases!\n")}
    # Check for convergence
    if (L[i] - L[i - 1] < epsilon_conv) { break }
    # Check if VB converged in the given maximum iterations
    if (i == max_iter) {warning("VB did not converge!\n")}
  }
  # Construct an object as output of the function
  obj = structure(list("Estimated Regression Coefficients" = m,
                       "Covariance Matrix" = S,
                       a = a, b = b, lambda = lambda,
                       "Design Matrix" = X, "Number of Observation" = n, 
                       "Number of Covariates" = c, L = L[2:i]),
                  class = "vlr")
  return(obj) # Write the output
}


           ###Performance Comparison of Variational vs MLE in Multiple Linear Regression
#This experiment does comparison for coefficient estimates and prediction for out-of-sample data.
#Data is generated randomly using randomization function based on normal distirbution: rnorm.
#Unless stated otherwise, the normal distirbution used to randomize the data is the standard normal distirbution N(0,1)



##Multiple Linear Regression Simulation
#Consider a 3 independent variables multiple linear regression:
#y = b0 + b1*x1 + b2*x2 + b3*x3 + e
#Refer model with: b0 = 2, b1 = 3, b2 = 4, b3 = 5, lambda = 1, as model 1
#Refer model with: b0 = 0, b1 = 3, b2 = 4, b3 = 5, lambda = 1, as model 2
#Refer model with: b0 = 1, b1 = 2, b2 = 3, b3 = 4, lambda = 1/9, as model 3
#For each of the models above, comparison of estimates and MSE from the 2 methods is done for sample size n = 100, 500, and 1000
#For each of the models and sample sizes, comparison of average and variance of the parameters estimates from the 2 methods is done based on 500 different trials

D = 3                    # Number of covariates
##Model 1
coefs = c(2,3,4,5)       # Coefficient of the regression model
lamb = 1                 # Precision parameter of the model
#n = 100
set.seed(2901)                                                       # Seed for pseudo number randomization/reproducibility of the results purpose
n = 100                                                              # Sample size
X = cbind(1, replicate(D, rnorm(n)))                                 # Generate X data using randomize from normal distirbution
y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb))                       # Generate y data using model 1 with the error randomize from normal distribution
vlr_model = vlr(X = X, y = y, is_verbose = TRUE)                     # Estimate regression coefficients using variational method
mlr_model = lm(y ~ X[,-1])                                           # Estimate regression coefficients using MLE method
vlr_model$`Estimated Regression Coefficients`                        # Regression coefficient estimates from variational method
mlr_model$coefficients                                               # Regression coefficient estimates from MLE method
(n-D-1)/sum((X%*%vlr_model$`Estimated Regression Coefficients`-y)^2) # Lambda estimates from variational method
(n-D-1)/sum((X%*%mlr_model$coefficients-y)^2)                        # Lambda estimates from MLE method
#Replicated comparisons
vm_parameters_b0 = rep(NA, 500)
vm_parameters_b1 = rep(NA, 500)
vm_parameters_b2 = rep(NA, 500)
vm_parameters_b3 = rep(NA, 500)
vm_lambda = rep(NA, 500)
mle_parameters_b0 = rep(NA, 500)
mle_parameters_b1 = rep(NA, 500)
mle_parameters_b2 = rep(NA, 500)
mle_parameters_b3 = rep(NA, 500)
mle_lambda = rep(NA, 500)
for(i in 1:500){
  X = cbind(1, replicate(D, rnorm(n)))             
  y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb)) 
  
  vm_mod = vlr(X = X, y = y, is_verbose = TRUE)
  lm_mod = lm(y~X[,-1])
  vm_parameters_b0[i] = vm_mod$`Estimated Regression Coefficients`[1] 
  vm_parameters_b1[i] = vm_mod$`Estimated Regression Coefficients`[2] 
  vm_parameters_b2[i] = vm_mod$`Estimated Regression Coefficients`[3] 
  vm_parameters_b3[i] = vm_mod$`Estimated Regression Coefficients`[4] 
  vm_lambda[i] = (n-D-1)/sum((X%*%vm_mod$`Estimated Regression Coefficients`-y)^2)
  mle_parameters_b0[i] = lm_mod$coefficients[1]                      
  mle_parameters_b1[i] = lm_mod$coefficients[2]                      
  mle_parameters_b2[i] = lm_mod$coefficients[3]                      
  mle_parameters_b3[i] = lm_mod$coefficients[4]
  mle_lambda[i] = (n-D-1)/sum((X%*%lm_mod$coefficients-y)^2)
}
#Transform the data to be visualized and comparison estimates for each parameters
#b0
df = data.frame(vm_parameters_b0, mle_parameters_b0)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16) + ylim(1.6, 2.4) + geom_hline(yintercept = 2 )
#Average of squared difference of parameters
mean(vm_parameters_b0)
mean(mle_parameters_b0)
#Variance of difference of parameters
var(vm_parameters_b0)
var(mle_parameters_b0)
#b1
df = data.frame(vm_parameters_b1, mle_parameters_b1)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16) + ylim(2.6, 3.4) + geom_hline(yintercept = 3 )
#Average of squared difference of parameters
mean(vm_parameters_b1)
mean(mle_parameters_b1)
#Variance of difference of parameters
var(vm_parameters_b1)
var(mle_parameters_b1)
#b2
df = data.frame(vm_parameters_b2, mle_parameters_b2)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16) + ylim(3.7, 4.3) + geom_hline(yintercept = 4 )
#Average of squared difference of parameters
mean(vm_parameters_b2)
mean(mle_parameters_b2)
#Variance of difference of parameters
var(vm_parameters_b2)
var(mle_parameters_b2)
#b3
df = data.frame(vm_parameters_b3, mle_parameters_b3)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16) + ylim(4.7, 5.3) + geom_hline(yintercept = 5 )
#Average of squared difference of parameters
mean(vm_parameters_b3)
mean(mle_parameters_b3)
#Variance of difference of parameters
var(vm_parameters_b3)
var(mle_parameters_b3)
#lambda
df = data.frame(vm_lambda, mle_lambda)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16)  + geom_hline(yintercept = 1 ) + ylim(0.6,1.4)
#Average of squared difference of parameters
mean(vm_lambda)
mean(mle_lambda)
#Variance of difference of parameters
var(vm_lambda)
var(mle_lambda)

#n = 500
set.seed(2901)                                                       # Seed for pseudo number randomization/reproducibility of the results purpose
n = 500                                                              # Sample size
X = cbind(1, replicate(D, rnorm(n)))                                 # Generate X data using randomize from normal distirbution
y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb))                       # Generate y data using model 1 with the error randomize from normal distribution
vlr_model = vlr(X = X, y = y, is_verbose = TRUE)                     # Estimate regression coefficients using variational method
mlr_model = lm(y ~ X[,-1])                                           # Estimate regression coefficients using MLE method
vlr_model$`Estimated Regression Coefficients`                        # Regression coefficient estimates from variational method
mlr_model$coefficients                                               # Regression coefficient estimates from MLE method
(n-D-1)/sum((X%*%vlr_model$`Estimated Regression Coefficients`-y)^2) # Lambda estimates from variational method
(n-D-1)/sum((X%*%mlr_model$coefficients-y)^2)                        # Lambda estimates from MLE method
#Replicated comparisons
vm_parameters_b0 = rep(NA, 500)
vm_parameters_b1 = rep(NA, 500)
vm_parameters_b2 = rep(NA, 500)
vm_parameters_b3 = rep(NA, 500)
vm_lambda = rep(NA, 500)
mle_parameters_b0 = rep(NA, 500)
mle_parameters_b1 = rep(NA, 500)
mle_parameters_b2 = rep(NA, 500)
mle_parameters_b3 = rep(NA, 500)
mle_lambda = rep(NA, 500)
for(i in 1:500){
  X = cbind(1, replicate(D, rnorm(n)))             
  y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb)) 
  
  vm_mod = vlr(X = X, y = y, is_verbose = TRUE)
  lm_mod = lm(y~X[,-1])
  vm_parameters_b0[i] = vm_mod$`Estimated Regression Coefficients`[1] 
  vm_parameters_b1[i] = vm_mod$`Estimated Regression Coefficients`[2] 
  vm_parameters_b2[i] = vm_mod$`Estimated Regression Coefficients`[3] 
  vm_parameters_b3[i] = vm_mod$`Estimated Regression Coefficients`[4] 
  vm_lambda[i] = (n-D-1)/sum((X%*%vm_mod$`Estimated Regression Coefficients`-y)^2)
  mle_parameters_b0[i] = lm_mod$coefficients[1]                      
  mle_parameters_b1[i] = lm_mod$coefficients[2]                      
  mle_parameters_b2[i] = lm_mod$coefficients[3]                      
  mle_parameters_b3[i] = lm_mod$coefficients[4]
  mle_lambda[i] = (n-D-1)/sum((X%*%lm_mod$coefficients-y)^2)
}
#Transform the data to be visualized and comparison estimates for each parameters
#b0
df = data.frame(vm_parameters_b0, mle_parameters_b0)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16) + ylim(1.6, 2.4) + geom_hline(yintercept = 2 )
#Average of squared difference of parameters
mean(vm_parameters_b0)
mean(mle_parameters_b0)
#Variance of difference of parameters
var(vm_parameters_b0)
var(mle_parameters_b0)
#b1
df = data.frame(vm_parameters_b1, mle_parameters_b1)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16) + ylim(2.6, 3.4) + geom_hline(yintercept = 3 )
#Average of squared difference of parameters
mean(vm_parameters_b1)
mean(mle_parameters_b1)
#Variance of difference of parameters
var(vm_parameters_b1)
var(mle_parameters_b1)
#b2
df = data.frame(vm_parameters_b2, mle_parameters_b2)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16) + ylim(3.7, 4.3) + geom_hline(yintercept = 4 )
#Average of squared difference of parameters
mean(vm_parameters_b2)
mean(mle_parameters_b2)
#Variance of difference of parameters
var(vm_parameters_b2)
var(mle_parameters_b2)
#b3
df = data.frame(vm_parameters_b3, mle_parameters_b3)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16) + ylim(4.7, 5.3) + geom_hline(yintercept = 5 )
#Average of squared difference of parameters
mean(vm_parameters_b3)
mean(mle_parameters_b3)
#Variance of difference of parameters
var(vm_parameters_b3)
var(mle_parameters_b3)
#lambda
df = data.frame(vm_lambda, mle_lambda)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16)  + geom_hline(yintercept = 1 ) + ylim(0.6,1.4)
#Average of squared difference of parameters
mean(vm_lambda)
mean(mle_lambda)
#Variance of difference of parameters
var(vm_lambda)
var(mle_lambda)

#n = 1000
set.seed(2901)                                                       # Seed for pseudo number randomization/reproducibility of the results purpose
n = 1000                                                             # Sample size
X = cbind(1, replicate(D, rnorm(n)))                                 # Generate X data using randomize from normal distirbution
y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb))                       # Generate y data using model 1 with the error randomize from normal distribution
vlr_model = vlr(X = X, y = y, is_verbose = TRUE)                     # Estimate regression coefficients using variational method
mlr_model = lm(y ~ X[,-1])                                           # Estimate regression coefficients using MLE method
vlr_model$`Estimated Regression Coefficients`                        # Regression coefficient estimates from variational method
mlr_model$coefficients                                               # Regression coefficient estimates from MLE method
(n-D-1)/sum((X%*%vlr_model$`Estimated Regression Coefficients`-y)^2) # MSE of regression with coefficient estimates from variational method
(n-D-1)/sum((X%*%mlr_model$coefficients-y)^2)                        # MSE of regression with coefficient estimates from MLE method
#Replicated comparisons
vm_parameters_b0 = rep(NA, 500)
vm_parameters_b1 = rep(NA, 500)
vm_parameters_b2 = rep(NA, 500)
vm_parameters_b3 = rep(NA, 500)
vm_lambda = rep(NA, 500)
mle_parameters_b0 = rep(NA, 500)
mle_parameters_b1 = rep(NA, 500)
mle_parameters_b2 = rep(NA, 500)
mle_parameters_b3 = rep(NA, 500)
mle_lambda = rep(NA, 500)
for(i in 1:500){
  X = cbind(1, replicate(D, rnorm(n)))             
  y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb)) 
  
  vm_mod = vlr(X = X, y = y, is_verbose = TRUE)
  lm_mod = lm(y~X[,-1])
  vm_parameters_b0[i] = vm_mod$`Estimated Regression Coefficients`[1] 
  vm_parameters_b1[i] = vm_mod$`Estimated Regression Coefficients`[2] 
  vm_parameters_b2[i] = vm_mod$`Estimated Regression Coefficients`[3] 
  vm_parameters_b3[i] = vm_mod$`Estimated Regression Coefficients`[4] 
  vm_lambda[i] = (n-D-1)/sum((X%*%vm_mod$`Estimated Regression Coefficients`-y)^2)
  mle_parameters_b0[i] = lm_mod$coefficients[1]                      
  mle_parameters_b1[i] = lm_mod$coefficients[2]                      
  mle_parameters_b2[i] = lm_mod$coefficients[3]                      
  mle_parameters_b3[i] = lm_mod$coefficients[4]
  mle_lambda[i] = (n-D-1)/sum((X%*%lm_mod$coefficients-y)^2)
}
#Transform the data to be visualized and comparison estimates for each parameters
#b0
df = data.frame(vm_parameters_b0, mle_parameters_b0)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16) + ylim(1.6, 2.4) + geom_hline(yintercept = 2 )
#Average of squared difference of parameters
mean(vm_parameters_b0)
mean(mle_parameters_b0)
#Variance of difference of parameters
var(vm_parameters_b0)
var(mle_parameters_b0)
#b1
df = data.frame(vm_parameters_b1, mle_parameters_b1)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16) + ylim(2.6, 3.4) + geom_hline(yintercept = 3 )
#Average of squared difference of parameters
mean(vm_parameters_b1)
mean(mle_parameters_b1)
#Variance of difference of parameters
var(vm_parameters_b1)
var(mle_parameters_b1)
#b2
df = data.frame(vm_parameters_b2, mle_parameters_b2)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16) + ylim(3.7, 4.3) + geom_hline(yintercept = 4 )
#Average of squared difference of parameters
mean(vm_parameters_b2)
mean(mle_parameters_b2)
#Variance of difference of parameters
var(vm_parameters_b2)
var(mle_parameters_b2)
#b3
df = data.frame(vm_parameters_b3, mle_parameters_b3)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16) + ylim(4.7, 5.3) + geom_hline(yintercept = 5 )
#Average of squared difference of parameters
mean(vm_parameters_b3)
mean(mle_parameters_b3)
#Variance of difference of parameters
var(vm_parameters_b3)
var(mle_parameters_b3)
#lambda
df = data.frame(vm_lambda, mle_lambda)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2, position=position_jitter(h=0.01,w=0.01)) +
  theme_classic(base_size = 16)  + geom_hline(yintercept = 1 ) + ylim(0.6,1.4)
#Average of squared difference of parameters
mean(vm_lambda)
mean(mle_lambda)
#Variance of difference of parameters
var(vm_lambda)
var(mle_lambda)



##Model 2
coefs = c(0,3,4,5)       # Coefficient of the regression model
lamb = 1                 # Precision parameter of the model
#n = 100
set.seed(2901)                                                       # Seed for pseudo number randomization/reproducibility of the results purpose
n = 100                                                              # Sample size
X = cbind(1, replicate(D, rnorm(n, mean = 100)))                     # Generate X data using randomize from normal distirbution
y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb))                       # Generate y data using model 1 with the error randomize from normal distribution
vlr_model = vlr(X = X, y = y, is_verbose = TRUE)                     # Estimate regression coefficients using variational method
mlr_model = lm(y ~ X[,-1])                                           # Estimate regression coefficients using MLE method
vlr_model$`Estimated Regression Coefficients`                        # Regression coefficient estimates from variational method
mlr_model$coefficients                                               # Regression coefficient estimates from MLE method
(n-D-1)/sum((X%*%vlr_model$`Estimated Regression Coefficients`-y)^2) # Lambda estimates from variational method
(n-D-1)/sum((X%*%mlr_model$coefficients-y)^2)                        # Lambda estimates from MLE method
#Replicated comparisons
vm_parameters_b0 = rep(NA, 500)
vm_parameters_b1 = rep(NA, 500)
vm_parameters_b2 = rep(NA, 500)
vm_parameters_b3 = rep(NA, 500)
vm_lambda = rep(NA, 500)
mle_parameters_b0 = rep(NA, 500)
mle_parameters_b1 = rep(NA, 500)
mle_parameters_b2 = rep(NA, 500)
mle_parameters_b3 = rep(NA, 500)
mle_lambda = rep(NA, 500)
set.seed(10223)
for(i in 1:500){
  X = cbind(1, replicate(D, rnorm(n, mean = 100)))             
  y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb)) 
  
  vm_mod = vlr(X = X, y = y, is_verbose = TRUE)
  lm_mod = lm(y~X[,-1])
  vm_parameters_b0[i] = vm_mod$`Estimated Regression Coefficients`[1]
  vm_parameters_b1[i] = vm_mod$`Estimated Regression Coefficients`[2] 
  vm_parameters_b2[i] = vm_mod$`Estimated Regression Coefficients`[3] 
  vm_parameters_b3[i] = vm_mod$`Estimated Regression Coefficients`[4] 
  vm_lambda[i] = (n-D-1)/sum((X%*%vm_mod$`Estimated Regression Coefficients`-y)^2)
  mle_parameters_b0[i] = lm_mod$coefficients[1]
  mle_parameters_b1[i] = lm_mod$coefficients[2]                      
  mle_parameters_b2[i] = lm_mod$coefficients[3]                      
  mle_parameters_b3[i] = lm_mod$coefficients[4]
  mle_lambda[i] = (n-D-1)/sum((X%*%lm_mod$coefficients-y)^2)
}
#Transform the data to be visualized and comparison estimates for each parameters
#b0
df = data.frame(vm_parameters_b0, mle_parameters_b0)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(-50, 50)  + geom_hline(yintercept = 0 )
#Average of estimates parameters
mean(vm_parameters_b0)
mean(mle_parameters_b0)
#Variance of estimates parameters
var(vm_parameters_b0)
var(mle_parameters_b0)
#b1
df = data.frame(vm_parameters_b1, mle_parameters_b1)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(2.6, 3.4) + geom_hline(yintercept = 3 )
#Average of estimates parameters
mean(vm_parameters_b1)
mean(mle_parameters_b1)
#Variance of estimates parameters
var(vm_parameters_b1)
var(mle_parameters_b1)
#b2
df = data.frame(vm_parameters_b2, mle_parameters_b2)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(3.7, 4.3) + geom_hline(yintercept = 4 )
#Average of estimates parameters
mean(vm_parameters_b2)
mean(mle_parameters_b2)
#Variance of estimates parameters
var(vm_parameters_b2)
var(mle_parameters_b2)
#b3
df = data.frame(vm_parameters_b3, mle_parameters_b3)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(4.7, 5.3) + geom_hline(yintercept = 5 )
#Average of estimates parameters
mean(vm_parameters_b3)
mean(mle_parameters_b3)
#Variance of estimates parameters
var(vm_parameters_b3)
var(mle_parameters_b3)
#lambda
df = data.frame(vm_lambda, mle_lambda)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16)  + geom_hline(yintercept = 1 ) + ylim(0.7,1.3)
#Average of estimates parameters
mean(vm_lambda)
mean(mle_lambda)
#Variance of estimates parameters
var(vm_lambda)
var(mle_lambda)

#n = 500
set.seed(2901)                                                       # Seed for pseudo number randomization/reproducibility of the results purpose
n = 500                                                              # Sample size
X = cbind(1, replicate(D, rnorm(n, mean = 100)))                     # Generate X data using randomize from normal distirbution
y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb))                       # Generate y data using model 1 with the error randomize from normal distribution
vlr_model = vlr(X = X, y = y, is_verbose = TRUE)                     # Estimate regression coefficients using variational method
mlr_model = lm(y ~ X[,-1])                                           # Estimate regression coefficients using MLE method
vlr_model$`Estimated Regression Coefficients`                        # Regression coefficient estimates from variational method
mlr_model$coefficients                                               # Regression coefficient estimates from MLE method
(n-D-1)/sum((X%*%vlr_model$`Estimated Regression Coefficients`-y)^2) # Lambda estimates from variational method
(n-D-1)/sum((X%*%mlr_model$coefficients-y)^2)                        # Lambda estimates from MLE method
#Replicated comparisons
vm_parameters_b0 = rep(NA, 500)
vm_parameters_b1 = rep(NA, 500)
vm_parameters_b2 = rep(NA, 500)
vm_parameters_b3 = rep(NA, 500)
vm_lambda = rep(NA, 500)
mle_parameters_b0 = rep(NA, 500)
mle_parameters_b1 = rep(NA, 500)
mle_parameters_b2 = rep(NA, 500)
mle_parameters_b3 = rep(NA, 500)
mle_lambda = rep(NA, 500)
set.seed(10223)
for(i in 1:500){
  X = cbind(1, replicate(D, rnorm(n, mean = 100)))             
  y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb)) 
  
  vm_mod = vlr(X = X, y = y, is_verbose = TRUE)
  lm_mod = lm(y~X[,-1])
  vm_parameters_b0[i] = vm_mod$`Estimated Regression Coefficients`[1]
  vm_parameters_b1[i] = vm_mod$`Estimated Regression Coefficients`[2] 
  vm_parameters_b2[i] = vm_mod$`Estimated Regression Coefficients`[3] 
  vm_parameters_b3[i] = vm_mod$`Estimated Regression Coefficients`[4] 
  vm_lambda[i] = (n-D-1)/sum((X%*%vm_mod$`Estimated Regression Coefficients`-y)^2)
  mle_parameters_b0[i] = lm_mod$coefficients[1]
  mle_parameters_b1[i] = lm_mod$coefficients[2]                      
  mle_parameters_b2[i] = lm_mod$coefficients[3]                      
  mle_parameters_b3[i] = lm_mod$coefficients[4]
  mle_lambda[i] = (n-D-1)/sum((X%*%lm_mod$coefficients-y)^2)
}
#Transform the data to be visualized and comparison estimates for each parameters
#b0
df = data.frame(vm_parameters_b0, mle_parameters_b0)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(-50, 50)  + geom_hline(yintercept = 0 )
#Average of estimates parameters
mean(vm_parameters_b0)
mean(mle_parameters_b0)
#Variance of estimates parameters
var(vm_parameters_b0)
var(mle_parameters_b0)
#b1
df = data.frame(vm_parameters_b1, mle_parameters_b1)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(2.6, 3.4) + geom_hline(yintercept = 3 )
#Average of estimates parameters
mean(vm_parameters_b1)
mean(mle_parameters_b1)
#Variance of estimates parameters
var(vm_parameters_b1)
var(mle_parameters_b1)
#b2
df = data.frame(vm_parameters_b2, mle_parameters_b2)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(3.7, 4.3) + geom_hline(yintercept = 4 )
#Average of estimates parameters
mean(vm_parameters_b2)
mean(mle_parameters_b2)
#Variance of estimates parameters
var(vm_parameters_b2)
var(mle_parameters_b2)
#b3
df = data.frame(vm_parameters_b3, mle_parameters_b3)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(4.7, 5.3) + geom_hline(yintercept = 5 )
#Average of estimates parameters
mean(vm_parameters_b3)
mean(mle_parameters_b3)
#Variance of estimates parameters
var(vm_parameters_b3)
var(mle_parameters_b3)
#lambda
df = data.frame(vm_lambda, mle_lambda)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16)  + geom_hline(yintercept = 1 ) + ylim(0.7,1.3)
#Average of estimates parameters
mean(vm_lambda)
mean(mle_lambda)
#Variance of estimates parameters
var(vm_lambda)
var(mle_lambda)

#n = 1000
set.seed(2901)                                                       # Seed for pseudo number randomization/reproducibility of the results purpose
n = 1000                                                             # Sample size
X = cbind(1, replicate(D, rnorm(n, mean = 100)))                     # Generate X data using randomize from normal distirbution
y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb))                       # Generate y data using model 1 with the error randomize from normal distribution
vlr_model = vlr(X = X, y = y, is_verbose = TRUE)                     # Estimate regression coefficients using variational method
mlr_model = lm(y ~ X[,-1])                                           # Estimate regression coefficients using MLE method
vlr_model$`Estimated Regression Coefficients`                        # Regression coefficient estimates from variational method
mlr_model$coefficients                                               # Regression coefficient estimates from MLE method
(n-D-1)/sum((X%*%vlr_model$`Estimated Regression Coefficients`-y)^2) # Lambda estimates from variational method
(n-D-1)/sum((X%*%mlr_model$coefficients-y)^2)                        # Lambda estimates from MLE method
#Replicated comparisons
vm_parameters_b0 = rep(NA, 500)
vm_parameters_b1 = rep(NA, 500)
vm_parameters_b2 = rep(NA, 500)
vm_parameters_b3 = rep(NA, 500)
vm_lambda = rep(NA, 500)
mle_parameters_b0 = rep(NA, 500)
mle_parameters_b1 = rep(NA, 500)
mle_parameters_b2 = rep(NA, 500)
mle_parameters_b3 = rep(NA, 500)
mle_lambda = rep(NA, 500)
set.seed(10223)
for(i in 1:500){
  X = cbind(1, replicate(D, rnorm(n, mean = 100)))             
  y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb)) 
  
  vm_mod = vlr(X = X, y = y, is_verbose = TRUE)
  lm_mod = lm(y~X[,-1])
  vm_parameters_b0[i] = vm_mod$`Estimated Regression Coefficients`[1]
  vm_parameters_b1[i] = vm_mod$`Estimated Regression Coefficients`[2] 
  vm_parameters_b2[i] = vm_mod$`Estimated Regression Coefficients`[3] 
  vm_parameters_b3[i] = vm_mod$`Estimated Regression Coefficients`[4] 
  vm_lambda[i] = (n-D-1)/sum((X%*%vm_mod$`Estimated Regression Coefficients`-y)^2)
  mle_parameters_b0[i] = lm_mod$coefficients[1]
  mle_parameters_b1[i] = lm_mod$coefficients[2]                      
  mle_parameters_b2[i] = lm_mod$coefficients[3]                      
  mle_parameters_b3[i] = lm_mod$coefficients[4]
  mle_lambda[i] = (n-D-1)/sum((X%*%lm_mod$coefficients-y)^2)
}
#Transform the data to be visualized and comparison estimates for each parameters
#b0
df = data.frame(vm_parameters_b0, mle_parameters_b0)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(-50, 50)  + geom_hline(yintercept = 0 )
#Average of estimates parameters
mean(vm_parameters_b0)
mean(mle_parameters_b0)
#Variance of estimates parameters
var(vm_parameters_b0)
var(mle_parameters_b0)
#b1
df = data.frame(vm_parameters_b1, mle_parameters_b1)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(2.6, 3.4) + geom_hline(yintercept = 3 )
#Average of estimates parameters
mean(vm_parameters_b1)
mean(mle_parameters_b1)
#Variance of estimates parameters
var(vm_parameters_b1)
var(mle_parameters_b1)
#b2
df = data.frame(vm_parameters_b2, mle_parameters_b2)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(3.7, 4.3) + geom_hline(yintercept = 4 )
#Average of estimates parameters
mean(vm_parameters_b2)
mean(mle_parameters_b2)
#Variance of estimates parameters
var(vm_parameters_b2)
var(mle_parameters_b2)
#b3
df = data.frame(vm_parameters_b3, mle_parameters_b3)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(4.7, 5.3) + geom_hline(yintercept = 5 )
#Average of estimates parameters
mean(vm_parameters_b3)
mean(mle_parameters_b3)
#Variance of estimates parameters
var(vm_parameters_b3)
var(mle_parameters_b3)
#lambda
df = data.frame(vm_lambda, mle_lambda)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16)  + geom_hline(yintercept = 1 ) + ylim(0.7,1.3)
#Average of estimates parameters
mean(vm_lambda)
mean(mle_lambda)
#Variance of estimates parameters
var(vm_lambda)
var(mle_lambda)





##Model 3
coefs = c(1,2,3,4)       # Coefficient of the regression model
lamb = 1/9               # Precision parameter of the model
#n = 100
set.seed(2901)                                                       # Seed for pseudo number randomization/reproducibility of the results purpose
n = 100                                                              # Sample size
X = cbind(1, rnorm(n), rnorm(n, mean = 100, sd = 2), 
          rnorm(n, mean = 10, sd = 10))                              # Generate X data using randomize from normal distirbution
y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb))                       # Generate y data using model 1 with the error randomize from normal distribution
vlr_model = vlr(X = X, y = y, is_verbose = TRUE)                     # Estimate regression coefficients using variational method
mlr_model = lm(y ~ X[,-1])                                           # Estimate regression coefficients using MLE method
vlr_model$`Estimated Regression Coefficients`                        # Regression coefficient estimates from variational method
mlr_model$coefficients                                               # Regression coefficient estimates from MLE method
(n-D-1)/sum((X%*%vlr_model$`Estimated Regression Coefficients`-y)^2) # Lambda estimates from variational method
(n-D-1)/sum((X%*%mlr_model$coefficients-y)^2)                        # Lambda estimates from MLE method
#Replicated comparisons
vm_parameters_b0 = rep(NA, 500)
vm_parameters_b1 = rep(NA, 500)
vm_parameters_b2 = rep(NA, 500)
vm_parameters_b3 = rep(NA, 500)
vm_lambda = rep(NA, 500)
mle_parameters_b0 = rep(NA, 500)
mle_parameters_b1 = rep(NA, 500)
mle_parameters_b2 = rep(NA, 500)
mle_parameters_b3 = rep(NA, 500)
mle_lambda = rep(NA, 500)
set.seed(102223)
for(i in 1:500){
  X = cbind(1, rnorm(n), rnorm(n, mean = 100, sd = 2), 
            rnorm(n, mean = 10, sd = 10))                              
  y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb))  
  
  vm_mod = vlr(X = X, y = y, is_verbose = TRUE)
  lm_mod = lm(y~X[,-1])
  vm_parameters_b0[i] = vm_mod$`Estimated Regression Coefficients`[1] 
  vm_parameters_b1[i] = vm_mod$`Estimated Regression Coefficients`[2] 
  vm_parameters_b2[i] = vm_mod$`Estimated Regression Coefficients`[3] 
  vm_parameters_b3[i] = vm_mod$`Estimated Regression Coefficients`[4] 
  vm_lambda[i] = (n-D-1)/sum((X%*%vm_mod$`Estimated Regression Coefficients`-y)^2)
  mle_parameters_b0[i] = lm_mod$coefficients[1]                      
  mle_parameters_b1[i] = lm_mod$coefficients[2]                      
  mle_parameters_b2[i] = lm_mod$coefficients[3]                      
  mle_parameters_b3[i] = lm_mod$coefficients[4]
  mle_lambda[i] = (n-D-1)/sum((X%*%lm_mod$coefficients-y)^2)
}
#Transform the data to be visualized and comparison estimates for each parameters
#b0
df = data.frame(vm_parameters_b0, mle_parameters_b0)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 1 )
#Average of estimates parameters
mean(vm_parameters_b0)
mean(mle_parameters_b0)
#Variance of estimates parameters
var(vm_parameters_b0)
var(mle_parameters_b0)
#b1
df = data.frame(vm_parameters_b1, mle_parameters_b1)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(1, 3) + geom_hline(yintercept = 2 )
#Average of estimates parameters
mean(vm_parameters_b1)
mean(mle_parameters_b1)
#Variance of estimates parameters
var(vm_parameters_b1)
var(mle_parameters_b1)
#b2
df = data.frame(vm_parameters_b2, mle_parameters_b2)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(2, 4) + geom_hline(yintercept = 3 )
#Average of estimates parameters
mean(vm_parameters_b2)
mean(mle_parameters_b2)
#Variance of estimates parameters
var(vm_parameters_b2)
var(mle_parameters_b2)
#b3
df = data.frame(vm_parameters_b3, mle_parameters_b3)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(3.5, 4.5) + geom_hline(yintercept = 4 )
#Average of estimates parameters
mean(vm_parameters_b3)
mean(mle_parameters_b3)
#Variance of estimates parameters
var(vm_parameters_b3)
var(mle_parameters_b3)
#lambda
df = data.frame(vm_lambda, mle_lambda)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16)  + geom_hline(yintercept = 1/9 )
#Average of estimates parameters
mean(vm_lambda)
mean(mle_lambda)
#Variance of estimates parameters
var(vm_lambda)
var(mle_lambda)

#n = 500
set.seed(2901)                                                       # Seed for pseudo number randomization/reproducibility of the results purpose
n = 500                                                              # Sample size
X = cbind(1, replicate(D, rnorm(n)))                                 # Generate X data using randomize from normal distirbution
y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb))                       # Generate y data using model 1 with the error randomize from normal distribution
vlr_model = vlr(X = X, y = y, is_verbose = TRUE)                     # Estimate regression coefficients using variational method
mlr_model = lm(y ~ X[,-1])                                           # Estimate regression coefficients using MLE method
vlr_model$`Estimated Regression Coefficients`                        # Regression coefficient estimates from variational method
mlr_model$coefficients                                               # Regression coefficient estimates from MLE method
(n-D-1)/sum((X%*%vlr_model$`Estimated Regression Coefficients`-y)^2) # Lambda estimates from variational method
(n-D-1)/sum((X%*%mlr_model$coefficients-y)^2)                        # Lambda estimates from MLE method
#Replicated comparisons
vm_parameters_b0 = rep(NA, 500)
vm_parameters_b1 = rep(NA, 500)
vm_parameters_b2 = rep(NA, 500)
vm_parameters_b3 = rep(NA, 500)
vm_lambda = rep(NA, 500)
mle_parameters_b0 = rep(NA, 500)
mle_parameters_b1 = rep(NA, 500)
mle_parameters_b2 = rep(NA, 500)
mle_parameters_b3 = rep(NA, 500)
mle_lambda = rep(NA, 500)
set.seed(102223)
for(i in 1:500){
  X = cbind(1, rnorm(n), rnorm(n, mean = 100, sd = 2), 
            rnorm(n, mean = 10, sd = 10))                              
  y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb))  
  
  vm_mod = vlr(X = X, y = y, is_verbose = TRUE)
  lm_mod = lm(y~X[,-1])
  vm_parameters_b0[i] = vm_mod$`Estimated Regression Coefficients`[1] 
  vm_parameters_b1[i] = vm_mod$`Estimated Regression Coefficients`[2] 
  vm_parameters_b2[i] = vm_mod$`Estimated Regression Coefficients`[3] 
  vm_parameters_b3[i] = vm_mod$`Estimated Regression Coefficients`[4] 
  vm_lambda[i] = (n-D-1)/sum((X%*%vm_mod$`Estimated Regression Coefficients`-y)^2)
  mle_parameters_b0[i] = lm_mod$coefficients[1]                      
  mle_parameters_b1[i] = lm_mod$coefficients[2]                      
  mle_parameters_b2[i] = lm_mod$coefficients[3]                      
  mle_parameters_b3[i] = lm_mod$coefficients[4]
  mle_lambda[i] = (n-D-1)/sum((X%*%lm_mod$coefficients-y)^2)
}
#Transform the data to be visualized and comparison estimates for each parameters
#b0
df = data.frame(vm_parameters_b0, mle_parameters_b0)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 1 )
#Average of estimates parameters
mean(vm_parameters_b0)
mean(mle_parameters_b0)
#Variance of estimates parameters
var(vm_parameters_b0)
var(mle_parameters_b0)
#b1
df = data.frame(vm_parameters_b1, mle_parameters_b1)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(1, 3) + geom_hline(yintercept = 2 )
#Average of estimates parameters
mean(vm_parameters_b1)
mean(mle_parameters_b1)
#Variance of estimates parameters
var(vm_parameters_b1)
var(mle_parameters_b1)
#b2
df = data.frame(vm_parameters_b2, mle_parameters_b2)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(2, 4) + geom_hline(yintercept = 3 )
#Average of estimates parameters
mean(vm_parameters_b2)
mean(mle_parameters_b2)
#Variance of estimates parameters
var(vm_parameters_b2)
var(mle_parameters_b2)
#b3
df = data.frame(vm_parameters_b3, mle_parameters_b3)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(3.5, 4.5) + geom_hline(yintercept = 4 )
#Average of estimates parameters
mean(vm_parameters_b3)
mean(mle_parameters_b3)
#Variance of estimates parameters
var(vm_parameters_b3)
var(mle_parameters_b3)
#lambda
df = data.frame(vm_lambda, mle_lambda)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16)  + geom_hline(yintercept = 1/9 )
#Average of estimates parameters
mean(vm_lambda)
mean(mle_lambda)
#Variance of estimates parameters
var(vm_lambda)
var(mle_lambda)

#n = 1000
set.seed(2901)                                                       # Seed for pseudo number randomization/reproducibility of the results purpose
n = 1000                                                             # Sample size
X = cbind(1, replicate(D, rnorm(n)))                                 # Generate X data using randomize from normal distirbution
y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb))                       # Generate y data using model 1 with the error randomize from normal distribution
vlr_model = vlr(X = X, y = y, is_verbose = TRUE)                     # Estimate regression coefficients using variational method
mlr_model = lm(y ~ X[,-1])                                           # Estimate regression coefficients using MLE method
vlr_model$`Estimated Regression Coefficients`                        # Regression coefficient estimates from variational method
mlr_model$coefficients                                               # Regression coefficient estimates from MLE method
(n-D-1)/sum((X%*%vlr_model$`Estimated Regression Coefficients`-y)^2) # Lambda estimates from variational method
(n-D-1)/sum((X%*%mlr_model$coefficients-y)^2)                        # Lambda estimates from MLE method
#Replicated comparisons
vm_parameters_b0 = rep(NA, 500)
vm_parameters_b1 = rep(NA, 500)
vm_parameters_b2 = rep(NA, 500)
vm_parameters_b3 = rep(NA, 500)
vm_lambda = rep(NA, 500)
mle_parameters_b0 = rep(NA, 500)
mle_parameters_b1 = rep(NA, 500)
mle_parameters_b2 = rep(NA, 500)
mle_parameters_b3 = rep(NA, 500)
mle_lambda = rep(NA, 500)
set.seed(102223)
for(i in 1:500){
  X = cbind(1, rnorm(n), rnorm(n, mean = 100, sd = 2), 
            rnorm(n, mean = 10, sd = 10))                              
  y = X %*% coefs  + rnorm(n, sd = sqrt(1/lamb))  
  
  vm_mod = vlr(X = X, y = y, is_verbose = TRUE)
  lm_mod = lm(y~X[,-1])
  vm_parameters_b0[i] = vm_mod$`Estimated Regression Coefficients`[1] 
  vm_parameters_b1[i] = vm_mod$`Estimated Regression Coefficients`[2] 
  vm_parameters_b2[i] = vm_mod$`Estimated Regression Coefficients`[3] 
  vm_parameters_b3[i] = vm_mod$`Estimated Regression Coefficients`[4] 
  vm_lambda[i] = (n-D-1)/sum((X%*%vm_mod$`Estimated Regression Coefficients`-y)^2)
  mle_parameters_b0[i] = lm_mod$coefficients[1]                      
  mle_parameters_b1[i] = lm_mod$coefficients[2]                      
  mle_parameters_b2[i] = lm_mod$coefficients[3]                      
  mle_parameters_b3[i] = lm_mod$coefficients[4]
  mle_lambda[i] = (n-D-1)/sum((X%*%lm_mod$coefficients-y)^2)
}
#Transform the data to be visualized and comparison estimates for each parameters
#b0
df = data.frame(vm_parameters_b0, mle_parameters_b0)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 1 )
#Average of estimates parameters
mean(vm_parameters_b0)
mean(mle_parameters_b0)
#Variance of estimates parameters
var(vm_parameters_b0)
var(mle_parameters_b0)
#b1
df = data.frame(vm_parameters_b1, mle_parameters_b1)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(1, 3) + geom_hline(yintercept = 2 )
#Average of estimates parameters
mean(vm_parameters_b1)
mean(mle_parameters_b1)
#Variance of estimates parameters
var(vm_parameters_b1)
var(mle_parameters_b1)
#b2
df = data.frame(vm_parameters_b2, mle_parameters_b2)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(2, 4) + geom_hline(yintercept = 3 )
#Average of estimates parameters
mean(vm_parameters_b2)
mean(mle_parameters_b2)
#Variance of estimates parameters
var(vm_parameters_b2)
var(mle_parameters_b2)
#b3
df = data.frame(vm_parameters_b3, mle_parameters_b3)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + ylim(3.5, 4.5) + geom_hline(yintercept = 4 )
#Average of estimates parameters
mean(vm_parameters_b3)
mean(mle_parameters_b3)
#Variance of estimates parameters
var(vm_parameters_b3)
var(mle_parameters_b3)
#lambda
df = data.frame(vm_lambda, mle_lambda)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16)  + geom_hline(yintercept = 1/9 )
#Average of estimates parameters
mean(vm_lambda)
mean(mle_lambda)
#Variance of estimates parameters
var(vm_lambda)
var(mle_lambda)







##Polynomial Regression Simulation
#Consider a 3rd degree polynomial regression:
#y = b0 + b1*x + b2*x^2 + b3*x^3 + e
D = 3                    # Number of covariates
##Model 4
#For this model, true value for the coefficients are set as: b0 = 2, b1 = 0.3, b2 = 5, b3 = 4, lambda = 1/16
coefs_pol = c(0,0.5,5,10)       # Coefficient of the regression model
lamb = 1/16                       # Precision parameter of the model
#n = 100 
set.seed(2901)                                                       # Seed for pseudo number randomization/reproducibility of the results purpose
n = 100    
#Assuming error term is comes from normal distribution N(0,1) and independent variable X is generated from normal distribution N(25, 3^2).
#Thus, the data-generating process is:
set.seed(2707)
X_1 = rnorm(100, mean = 25, sd = 3)
X_pol = cbind(1, X_1, X_1^2, X_1^3)
y_pol = X_pol %*% coefs_pol + rnorm(100)
plot("x" = X_1, "y" = y_pol, xlab = "x", ylab = "y")  # Plot of the generated data x and y
lm_pol = lm(y_pol ~ I(X_1) + I(X_1^2) + I(X_1^3))
vm_pol = vlr(X = X_pol, y = y_pol)
lm_pol$coefficients
vm_pol$`Estimated Regression Coefficients`
#Replicated comparisons
vm_parameters_b0 = rep(NA, 500)
vm_parameters_b1 = rep(NA, 500)
vm_parameters_b2 = rep(NA, 500)
vm_parameters_b3 = rep(NA, 500)
vm_lambda = rep(NA, 500)
mle_parameters_b0 = rep(NA, 500)
mle_parameters_b1 = rep(NA, 500)
mle_parameters_b2 = rep(NA, 500)
mle_parameters_b3 = rep(NA, 500)
mle_lambda = rep(NA, 500)
set.seed(102223)
for(i in 1:500){
  X_1 = rnorm(100, mean = 25, sd = 3)
  X_pol = cbind(1, X_1, X_1^2, X_1^3)
  y_pol = X_pol %*% coefs_pol + rnorm(100)
  
  lm_pol = lm(y_pol ~ I(X_1) + I(X_1^2) + I(X_1^3))
  vm_pol = vlr(X = X_pol, y = y_pol)
  vm_parameters_b0[i] = vm_pol$`Estimated Regression Coefficients`[1] 
  vm_parameters_b1[i] = vm_pol$`Estimated Regression Coefficients`[2] 
  vm_parameters_b2[i] = vm_pol$`Estimated Regression Coefficients`[3] 
  vm_parameters_b3[i] = vm_pol$`Estimated Regression Coefficients`[4] 
  vm_lambda[i] = (n-D-1)/sum((X%*%vm_mod$`Estimated Regression Coefficients`-y)^2)
  mle_parameters_b0[i] = lm_pol$coefficients[1]                      
  mle_parameters_b1[i] = lm_pol$coefficients[2]                      
  mle_parameters_b2[i] = lm_pol$coefficients[3]                      
  mle_parameters_b3[i] = lm_pol$coefficients[4]
  mle_lambda[i] = (n-D-1)/sum((X%*%lm_mod$coefficients-y)^2)
}
#Transform the data to be visualized and comparison estimates for each parameters
#b0
df = data.frame(vm_parameters_b0, mle_parameters_b0)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 0 )
#Average of estimates parameters
mean(vm_parameters_b0)
mean(mle_parameters_b0)
#Variance of estimates parameters
var(vm_parameters_b0)
var(mle_parameters_b0)
#b1
df = data.frame(vm_parameters_b1, mle_parameters_b1)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 0.5 )
#Average of estimates parameters
mean(vm_parameters_b1)
mean(mle_parameters_b1)
#Variance of estimates parameters
var(vm_parameters_b1)
var(mle_parameters_b1)
#b2
df = data.frame(vm_parameters_b2, mle_parameters_b2)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 5 )
#Average of estimates parameters
mean(vm_parameters_b2)
mean(mle_parameters_b2)
#Variance of estimates parameters
var(vm_parameters_b2)
var(mle_parameters_b2)
#b3
df = data.frame(vm_parameters_b3, mle_parameters_b3)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 10 ) + ylim(9.990, 10.015)
#Average of estimates parameters
mean(vm_parameters_b3)
mean(mle_parameters_b3)
#Variance of estimates parameters
var(vm_parameters_b3)
var(mle_parameters_b3)
#lambda
df = data.frame(vm_lambda, mle_lambda)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16)  + geom_hline(yintercept = 1/16 )
#Average of estimates parameters
mean(vm_lambda)
mean(mle_lambda)
#Variance of estimates parameters
var(vm_lambda)
var(mle_lambda)

##Model 5
#For this model, true value for the coefficients are set as: b0 = 2, b1 = 0.3, b2 = 5, b3 = 4, lambda = 1/16
coefs_pol = c(100,10,0.5,10)     # Coefficient of the regression model
lamb = 1                         # Precision parameter of the model
#n = 100 
set.seed(2901)                                                       # Seed for pseudo number randomization/reproducibility of the results purpose
n = 100    
#Assuming error term is comes from normal distribution N(0,1) and independent variable X is generated from normal distribution N(25, 3^2).
#Thus, the data-generating process is:
set.seed(2707)
X_1 = rnorm(100, mean = 25, sd = 3)
X_pol = cbind(1, X_1, X_1^2, X_1^3)
y_pol = X_pol %*% coefs_pol + rnorm(100)
plot("x" = X_1, "y" = y_pol, xlab = "x", ylab = "y")  # Plot of the generated data x and y
lm_pol = lm(y_pol ~ I(X_1) + I(X_1^2) + I(X_1^3))
vm_pol = vlr(X = X_pol, y = y_pol)
lm_pol$coefficients
vm_pol$`Estimated Regression Coefficients`
#Replicated comparisons
vm_parameters_b0 = rep(NA, 500)
vm_parameters_b1 = rep(NA, 500)
vm_parameters_b2 = rep(NA, 500)
vm_parameters_b3 = rep(NA, 500)
vm_lambda = rep(NA, 500)
mle_parameters_b0 = rep(NA, 500)
mle_parameters_b1 = rep(NA, 500)
mle_parameters_b2 = rep(NA, 500)
mle_parameters_b3 = rep(NA, 500)
mle_lambda = rep(NA, 500)
set.seed(102223)
for(i in 1:500){
  X_1 = rnorm(100, mean = 25, sd = 3)
  X_pol = cbind(1, X_1, X_1^2, X_1^3)
  y_pol = X_pol %*% coefs_pol + rnorm(100)
  
  lm_pol = lm(y_pol ~ I(X_1) + I(X_1^2) + I(X_1^3))
  vm_pol = vlr(X = X_pol, y = y_pol)
  vm_parameters_b0[i] = vm_pol$`Estimated Regression Coefficients`[1] 
  vm_parameters_b1[i] = vm_pol$`Estimated Regression Coefficients`[2] 
  vm_parameters_b2[i] = vm_pol$`Estimated Regression Coefficients`[3] 
  vm_parameters_b3[i] = vm_pol$`Estimated Regression Coefficients`[4] 
  vm_lambda[i] = (n-D-1)/sum((X%*%vm_mod$`Estimated Regression Coefficients`-y)^2)
  mle_parameters_b0[i] = lm_pol$coefficients[1]                      
  mle_parameters_b1[i] = lm_pol$coefficients[2]                      
  mle_parameters_b2[i] = lm_pol$coefficients[3]                      
  mle_parameters_b3[i] = lm_pol$coefficients[4]
  mle_lambda[i] = (n-D-1)/sum((X%*%lm_mod$coefficients-y)^2)
}
#Transform the data to be visualized and comparison estimates for each parameters
#b0
df = data.frame(vm_parameters_b0, mle_parameters_b0)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 100 )
#Average of estimates parameters
mean(vm_parameters_b0)
mean(mle_parameters_b0)
#Variance of estimates parameters
var(vm_parameters_b0)
var(mle_parameters_b0)
#b1
df = data.frame(vm_parameters_b1, mle_parameters_b1)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 10 )
#Average of estimates parameters
mean(vm_parameters_b1)
mean(mle_parameters_b1)
#Variance of estimates parameters
var(vm_parameters_b1)
var(mle_parameters_b1)
#b2
df = data.frame(vm_parameters_b2, mle_parameters_b2)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 0.5 )
#Average of estimates parameters
mean(vm_parameters_b2)
mean(mle_parameters_b2)
#Variance of estimates parameters
var(vm_parameters_b2)
var(mle_parameters_b2)
#b3
df = data.frame(vm_parameters_b3, mle_parameters_b3)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 5 )
#Average of estimates parameters
mean(vm_parameters_b3)
mean(mle_parameters_b3)
#Variance of estimates parameters
var(vm_parameters_b3)
var(mle_parameters_b3)
#lambda
df = data.frame(vm_lambda, mle_lambda)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16)  + geom_hline(yintercept = 1 )
#Average of estimates parameters
mean(vm_lambda)
mean(mle_lambda)
#Variance of estimates parameters
var(vm_lambda)
var(mle_lambda)

##Cross-validation
set.seed(10223)
X_1 = rnorm(100, mean = 25, sd = 3)
X_pol = cbind(1, X_1, X_1^2, X_1^3)
y_pol = X_pol %*% coefs_pol + rnorm(100)
dat_pol = data.frame(y = y_pol, x0 = X_pol[,1], x = X_pol[,2], x2 = X_pol[,3], x3 = X_pol[,4])
#Partitioning data
#1st fold
pol_test1 = as.matrix(dat_pol[1:20,])
lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(1:20),])
vm_pol = vlr(X = X_pol[-c(1:20),], y = y_pol[-c(1:20)])
#Compute mean squared prediction error in test data
mean((pol_test1[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test1[,1])^2) 
mean((pol_test1[,-1]%*%lm_pol$coefficients-pol_test1[,1])^2)
#2nd fold
pol_test2 = as.matrix(dat_pol[21:40,])
lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(21:40),])
vm_pol = vlr(X = X_pol[-c(21:40),], y = y_pol[-c(21:40)])
mean((pol_test2[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test2[,1])^2) 
mean((pol_test2[,-1]%*%lm_pol$coefficients-pol_test2[,1])^2)
#3rd fold
pol_test3 = as.matrix(dat_pol[41:60,])
lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(41:60),])
vm_pol = vlr(X = X_pol[-c(41:60),], y = y_pol[-c(41:60)])
mean((pol_test3[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test3[,1])^2) 
mean((pol_test3[,-1]%*%lm_pol$coefficients-pol_test3[,1])^2)
#4th fold
pol_test4 = as.matrix(dat_pol[61:80,])
lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(61:80),])
vm_pol = vlr(X = X_pol[-c(61:80),], y = y_pol[-c(61:80)])
mean((pol_test4[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test4[,1])^2) 
mean((pol_test4[,-1]%*%lm_pol$coefficients-pol_test4[,1])^2)
#5th fold
pol_test5 = as.matrix(dat_pol[81:100,])
lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(81:100),])
vm_pol = vlr(X = X_pol[-c(81:100),], y = y_pol[-c(81:100)])
mean((pol_test5[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test5[,1])^2) 
mean((pol_test5[,-1]%*%lm_pol$coefficients-pol_test5[,1])^2)

##Replicated cross validation
var_win = rep(NA, 500)
set.seed(10223)
for(i in 1:500){
  var_win[i] = 0
  fold_win = rep(NA, 5)
  X_1 = rnorm(100, mean = 25, sd = 3)
  X_pol = cbind(1, X_1, X_1^2, X_1^3)
  y_pol = X_pol %*% coefs_pol + rnorm(100)
  dat_pol = data.frame(y = y_pol, x0 = X_pol[,1], x = X_pol[,2], x2 = X_pol[,3], x3 = X_pol[,4])
  pol_test1 = as.matrix(dat_pol[1:20,])
  pol_test2 = as.matrix(dat_pol[21:40,])
  pol_test3 = as.matrix(dat_pol[41:60,])
  pol_test4 = as.matrix(dat_pol[61:80,])
  pol_test5 = as.matrix(dat_pol[81:100,])
  lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(1:20),])
  vm_pol = vlr(X = X_pol[-c(1:20),], y = y_pol[-c(1:20)])
  fold_win[1] = ifelse(mean((pol_test1[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test1[,1])^2)<=mean((pol_test1[,-1]%*%lm_pol$coefficients-pol_test1[,1])^2),
                         1,0)
  lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(21:40),])
  vm_pol = vlr(X = X_pol[-c(21:40),], y = y_pol[-c(21:40)])
  fold_win[2] = ifelse(mean((pol_test2[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test2[,1])^2)<=mean((pol_test2[,-1]%*%lm_pol$coefficients-pol_test2[,1])^2),
                       1, 0)
  lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(41:60),])
  vm_pol = vlr(X = X_pol[-c(41:60),], y = y_pol[-c(41:60)])
  fold_win[3] = ifelse(mean((pol_test3[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test3[,1])^2)<=mean((pol_test3[,-1]%*%lm_pol$coefficients-pol_test3[,1])^2),
                       1,0)
  lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(61:80),])
  vm_pol = vlr(X = X_pol[-c(61:80),], y = y_pol[-c(61:80)])
  fold_win[4] = ifelse(mean((pol_test4[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test4[,1])^2)<mean((pol_test4[,-1]%*%lm_pol$coefficients-pol_test4[,1])^2),
                       1,0)
  lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(81:100),])
  vm_pol = vlr(X = X_pol[-c(81:100),], y = y_pol[-c(81:100)])
  fold_win[5] = ifelse(mean((pol_test5[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test5[,1])^2)<mean((pol_test5[,-1]%*%lm_pol$coefficients-pol_test5[,1])^2),
                       1,0)
  var_win[i] = sum(fold_win)
}
mean(var_win) #average of 1.922 out of 5 folds in which variational method has better or equal accuracy to MLE

#Changing lambda value in the variational estimation function
vm_parameters_b0 = rep(NA, 500)
vm_parameters_b1 = rep(NA, 500)
vm_parameters_b2 = rep(NA, 500)
vm_parameters_b3 = rep(NA, 500)
vm_lambda = rep(NA, 500)
mle_parameters_b0 = rep(NA, 500)
mle_parameters_b1 = rep(NA, 500)
mle_parameters_b2 = rep(NA, 500)
mle_parameters_b3 = rep(NA, 500)
mle_lambda = rep(NA, 500)
set.seed(102223)
for(i in 1:500){
  X_1 = rnorm(100, mean = 25, sd = 3)
  X_pol = cbind(1, X_1, X_1^2, X_1^3)
  y_pol = X_pol %*% coefs_pol + rnorm(100)
  
  lm_pol = lm(y_pol ~ I(X_1) + I(X_1^2) + I(X_1^3))
  vm_pol = vlr(X = X_pol, y = y_pol, lambda = 100)
  vm_parameters_b0[i] = vm_pol$`Estimated Regression Coefficients`[1] 
  vm_parameters_b1[i] = vm_pol$`Estimated Regression Coefficients`[2] 
  vm_parameters_b2[i] = vm_pol$`Estimated Regression Coefficients`[3] 
  vm_parameters_b3[i] = vm_pol$`Estimated Regression Coefficients`[4] 
  vm_lambda[i] = (n-D-1)/sum((X%*%vm_mod$`Estimated Regression Coefficients`-y)^2)
  mle_parameters_b0[i] = lm_pol$coefficients[1]                      
  mle_parameters_b1[i] = lm_pol$coefficients[2]                      
  mle_parameters_b2[i] = lm_pol$coefficients[3]                      
  mle_parameters_b3[i] = lm_pol$coefficients[4]
  mle_lambda[i] = (n-D-1)/sum((X%*%lm_mod$coefficients-y)^2)
}
#Transform the data to be visualized and comparison estimates for each parameters
#b0
df = data.frame(vm_parameters_b0, mle_parameters_b0)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 100 )
#Average of estimates parameters
mean(vm_parameters_b0)
mean(mle_parameters_b0)
#Variance of estimates parameters
var(vm_parameters_b0)
var(mle_parameters_b0)
#b1
df = data.frame(vm_parameters_b1, mle_parameters_b1)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 10 )
#Average of estimates parameters
mean(vm_parameters_b1)
mean(mle_parameters_b1)
#Variance of estimates parameters
var(vm_parameters_b1)
var(mle_parameters_b1)
#b2
df = data.frame(vm_parameters_b2, mle_parameters_b2)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 0.5 )
#Average of estimates parameters
mean(vm_parameters_b2)
mean(mle_parameters_b2)
#Variance of estimates parameters
var(vm_parameters_b2)
var(mle_parameters_b2)
#b3
df = data.frame(vm_parameters_b3, mle_parameters_b3)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16) + geom_hline(yintercept = 5 )
#Average of estimates parameters
mean(vm_parameters_b3)
mean(mle_parameters_b3)
#Variance of estimates parameters
var(vm_parameters_b3)
var(mle_parameters_b3)
#lambda
df = data.frame(vm_lambda, mle_lambda)
df_long = df %>% 
  mutate(rep = row_number()) %>% 
  gather(method, value, -rep)
#Plot the estimations from the 2 methods in a graphic using ggplot function
ggplot(df_long, aes(x = rep, y = value)) +
  geom_point(aes(color = method, shape = method), size = 2) +
  theme_classic(base_size = 16)  + geom_hline(yintercept = 1 )
#Average of estimates parameters
mean(vm_lambda)
mean(mle_lambda)
#Variance of estimates parameters
var(vm_lambda)
var(mle_lambda)




var_win = rep(NA, 500)
set.seed(10223)
for(i in 1:500){
  var_win[i] = 0
  fold_win = rep(NA, 5)
  X_1 = rnorm(100, mean = 25, sd = 3)
  X_pol = cbind(1, X_1, X_1^2, X_1^3)
  y_pol = X_pol %*% coefs_pol + rnorm(100)
  dat_pol = data.frame(y = y_pol, x0 = X_pol[,1], x = X_pol[,2], x2 = X_pol[,3], x3 = X_pol[,4])
  pol_test1 = as.matrix(dat_pol[1:20,])
  pol_test2 = as.matrix(dat_pol[21:40,])
  pol_test3 = as.matrix(dat_pol[41:60,])
  pol_test4 = as.matrix(dat_pol[61:80,])
  pol_test5 = as.matrix(dat_pol[81:100,])
  lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(1:20),])
  vm_pol = vlr(X = X_pol[-c(1:20),], y = y_pol[-c(1:20)], lambda = 100)
  fold_win[1] = ifelse(mean((pol_test1[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test1[,1])^2)<=mean((pol_test1[,-1]%*%lm_pol$coefficients-pol_test1[,1])^2),
                       1,0)
  lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(21:40),])
  vm_pol = vlr(X = X_pol[-c(21:40),], y = y_pol[-c(21:40)], lambda = 100)
  fold_win[2] = ifelse(mean((pol_test2[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test2[,1])^2)<=mean((pol_test2[,-1]%*%lm_pol$coefficients-pol_test2[,1])^2),
                       1, 0)
  lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(41:60),])
  vm_pol = vlr(X = X_pol[-c(41:60),], y = y_pol[-c(41:60)], lambda = 100)
  fold_win[3] = ifelse(mean((pol_test3[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test3[,1])^2)<=mean((pol_test3[,-1]%*%lm_pol$coefficients-pol_test3[,1])^2),
                       1,0)
  lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(61:80),])
  vm_pol = vlr(X = X_pol[-c(61:80),], y = y_pol[-c(61:80)], lambda = 100)
  fold_win[4] = ifelse(mean((pol_test4[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test4[,1])^2)<mean((pol_test4[,-1]%*%lm_pol$coefficients-pol_test4[,1])^2),
                       1,0)
  lm_pol = lm(y ~ x + x2 + x3, data = dat_pol[-c(81:100),])
  vm_pol = vlr(X = X_pol[-c(81:100),], y = y_pol[-c(81:100)], lambda = 100)
  fold_win[5] = ifelse(mean((pol_test5[,-1]%*%vm_pol$`Estimated Regression Coefficients`-pol_test5[,1])^2)<mean((pol_test5[,-1]%*%lm_pol$coefficients-pol_test5[,1])^2),
                       1,0)
  var_win[i] = sum(fold_win)
}
mean(var_win)   #2.632 out of 5 folds in which variational method has better or equal accuracy to MLE


##Applying to Survey data
lm_sur = lm(grade ~ posting + I(posting^2), data = samp_post)
plot(lm_sur, 1)                  # Homoskedasticity 
shapiro.test(lm_sur$residuals)   # Normal residual
anova(lm_sur)    # MSE of the model is 82
X_sur = cbind(1, samp_post$posting, samp_post$posting^2)
vm_sur = vlr(X = X_sur, y = samp_post$grade)
lm_sur$coefficients
vm_sur$`Estimated Regression Coefficients`
#Cross-validation against out-of-sample data
X_sur2 = as.matrix(cbind(1, outsamp_post$posting, outsamp_post$posting^2))

mean((X_sur2%*%vm_sur$`Estimated Regression Coefficients`-outsamp_post$grade)^2) 
var((X_sur2%*%vm_sur$`Estimated Regression Coefficients`-outsamp_post$grade)^2)
mean((X_sur2%*%lm_sur$coefficients-outsamp_post$grade)^2)
var((X_sur2%*%lm_sur$coefficients-outsamp_post$grade)^2)
outsamp.y_sur = data.frame(id = c(1:99),
                           true_value = outsamp_post$grade, 
                           MLE_predict = X_sur2%*%lm_sur$coefficients,
                           VM_predict = X_sur2%*%vm_sur$`Estimated Regression Coefficients`)
head(outsamp.y_sur, 10)

ggplot(data = outsamp.y_sur)+
  geom_point(aes(x = id, y = true_value), color = "black")+
  geom_point(aes(x = id, y = MLE_predict), color = "red")+
  geom_point(aes(x = id, y = VM_predict), color = "blue", shape = 1)