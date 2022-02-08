### Regularization ###
library(glmnet)
library(dplyr)
library(ggplot2)

# To illustrate what does it do, lets imagine we want to build a regression model of predictors X
#the target Y. We simplify on both variables being a numerical continuous variables.
#Lets create dummy data for this
set.seed(919191)
x = runif(100, 0, 1) #randomize 100 observations from a uniform distribution on [0,1]
y = sin(2*pi*x) + rnorm(100, 0, 1) #target is a combination from function of X and gaussian "noise"

#Divide into train and test dataset
idx = sample(1:length(x), 0.4*length(x))
x_train = x[idx]
x_test = x[-idx]
y_train = y[idx]
y_test = y[-idx]
train = as.data.frame(cbind(x_train, y_train))
test = as.data.frame(cbind(x_test, y_test))
colnames(train) = c("x", "y")
colnames(test) = c("x", "y")

#Fit several polynomial regression models with order represented as M
list.of.fits <- list()
for (i in 1:9) {
  fit.name = paste0("M", i)
  list.of.fits[[fit.name]] = glm(y~poly(x, i), data = train)
}

#Compare by plot
library(gridExtra)
grid.arrange(train %>% ggplot(aes(x, y), xl ) + geom_point() +
               stat_smooth(method = lm, formula = y ~ x, col = "red", se = FALSE)+ 
               ylab("Y")+ xlab ("X") + labs(title = "M = 1"),
             train %>% ggplot(aes(x, y), xl ) + geom_point() +
               stat_smooth(method = lm, formula = y ~ x + I(x^2), col = "red", se = FALSE) +
               ylab("Y")+ xlab ("X") + labs(title = "M = 2"),
             train %>% ggplot(aes(x, y), xl ) + geom_point() +
               stat_smooth(method = lm, formula = y ~ x + I(x^2) + I(x^3), col = "red", se = FALSE) +
               ylab("Y")+ xlab ("X") + labs(title = "M = 3"),
             train %>% ggplot(aes(x, y), xl ) + geom_point() +
               stat_smooth(method = lm, formula = y ~ x + I(x^2) + I(x^3) + I(x^4), col = "red", se = FALSE) +
               ylab("Y")+ xlab ("X") + labs(title = "M = 4"),
             nrow = 2)
grid.arrange(train %>% ggplot(aes(x, y), xl ) + geom_point() +
               stat_smooth(method = lm, formula = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5), col = "red", se = FALSE) +
               ylab("Y")+ xlab ("X") + labs(title = "M = 5"),
             train %>% ggplot(aes(x, y), xl ) + geom_point() +
               stat_smooth(method = lm, formula = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6), col = "red", se = FALSE) +
               ylab("Y")+ xlab ("X") + labs(title = "M = 6"),
             train %>% ggplot(aes(x, y), xl ) + geom_point() +
               stat_smooth(method = lm, formula = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7), col = "red", se = FALSE) +
               ylab("Y")+ xlab ("X") + labs(title = "M = 7"),
             train %>% ggplot(aes(x, y), xl ) + geom_point() +
               stat_smooth(method = lm, formula = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8), col = "red", se = FALSE) +
               ylab("Y")+ xlab ("X") + labs(title = "M = 8"),
             nrow = 2)
grid.arrange(train %>% ggplot(aes(x,y), xl) + geom_point() +
               stat_smooth(method = lm, formula = y ~ x, col = "red", se = FALSE) +
               ylab("Y") + xlab("X") + labs(title = "M = 1"),
             train %>% ggplot(aes(x, y), xl) + geom_point() +
               stat_smooth(method = lm, formula = y ~ x + I(x^2) + I(x^3) + I(x^4), col = "red", se = FALSE) +
               ylab("Y") + xlab("X") + labs(title = "M = 4"),
             train %>% ggplot(aes(x, y), xl ) + geom_point() +
               stat_smooth(method = lm, formula = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7), col = "red", se = FALSE) +
               ylab("Y")+ xlab ("X") + labs(title = "M = 7"),
             train %>% ggplot(aes(x, y), xl ) + geom_point() +
               stat_smooth(method = lm, formula = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9), col = "red", se = FALSE) +
               ylab("Y")+ xlab ("X") + labs(title = "M = 9"),
             nrow = 2
)

#Comparing from the previous plots, we can see that higher order polynomial seems to fit the train data better.
#Of course if we follow this, higher order polynomial is a rather plausible model for this.
#But before jumping to that conclusion, why dont we take a look of each model's performance on the "test data"
#Lets compare each of the previous model by MSE in the test data
results = data.frame()
for (i in 1:9) {
  fit.name = i
  predicted =  predict(list.of.fits[[fit.name]], newdata = test)
  mse = mean((y_test - predicted)^2)
  temp = data.frame(M = i, mse=mse, fit.name = paste0("M",fit.name))
  results = rbind(results, temp)
}
#View the results
results
#Turns out, the trend found from previous plot comparison seems to be "rejected" as mse in higher order polynomial doesnt
#decrease but rather increase. We can also compare this on the actual equation which is y = sin(2*pi*x).
#Compare M = 9 with the actual equation
grid.arrange(train %>% ggplot(aes(x, y), xl ) +
               stat_smooth(method = lm, formula = y ~ x + I(x^2) + I(x^3), 
                           col = "red", se = FALSE) +
               geom_function(fun = function(x) sin(2*pi*x), colour = "blue", inherit.aes = FALSE) +
               ylab("Y")+ xlab ("X") + labs(title = "Polynomial Regression Order 3 VS Actual Function"),
             train %>% ggplot(aes(x, y), xl ) +
               stat_smooth(method = lm, formula = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9), 
                           col = "red", se = FALSE) +
               geom_function(fun = function(x) sin(2*pi*x), colour = "blue", inherit.aes = FALSE) +
               ylab("Y")+ xlab ("X") + labs(title = "Polynomial Regression Order 9 VS Actual Function"),
             nrow = 2, ncol = 2
) 
#These actions of "comparing" performance of model in case of training and test data is known as "cross validation"
#and for this specific problem, we want to use a regularization method on the regression so that we can find a model
#which performs "okayish" in train data but not falling off greatly when fitted into a test data.
#To set such constraint in the model, regularization apply whats called "penalty term".
#This penalty works as to adjust to the previous condition : stabilizing performance on train and test data.