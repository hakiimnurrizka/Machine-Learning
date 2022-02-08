### Exploratory Data Analysis on client Data ###
library(readxl)
library(ggplot2)
library(corrplot)
library(dplyr)
library(gridExtra)

## Load the data
client_train = read_excel("~/Github/Machine-Learning/client-data.xlsx", sheet = "client_train")
client_test = read_excel("~/Github/Machine-Learning/client-data.xlsx", sheet = "client_test")
View(client_train)
client_train = client_train[,-1] #Remove id as it serves no purpose in exploration nor analysis
client_test = client_test[,-1]
str(client_train) #Data structure, variables type
client_train$default.payment.next.month = as.factor(client_train$default.payment.next.month)
#To give context, our goal is to analyze characteristics of 2 type of clients based on payment status : 
#(1) Skip a payment at least once (bad)
#(2) Never skip a payment (good)
#Then build a model to predict clients' payment status on the test data.

#Lets exclude pay_ as this is more of a time series categorical variable
client_train1 = as.data.frame(apply(client_train[,6:11], 2, as.factor))
client_test1 = as.data.frame(apply(client_test[,6:11], 2, as.factor)) #We save it for another methods to analyze

## Cleaning the empty and descriptive statistics
client_train = na.omit(client_train) #eliminate rows with NA, do this unless you want to do imputation
summary(client_train) #This is not going to help in analysis, lets instead do summary based on status payment
client_train %>% group_by(default.payment.next.month) %>%
  summarise(mean_limit_bal = mean(LIMIT_BAL), mean_age = mean(AGE),
            mean_bill1 = mean(BILL_AMT1), mean_bill2 = mean(BILL_AMT2),
            mean_bill3 = mean(BILL_AMT3), mean_bill4 = mean(BILL_AMT4),
            mean_bill5 = mean(BILL_AMT5), mean_bill6 = mean(BILL_AMT6),
            mean_pay1 = mean(PAY_AMT1), mean_pay2 = mean(PAY_AMT2),
            mean_pay3 = mean(PAY_AMT3), mean_pay4 = mean(PAY_AMT4),
            mean_pay5 = mean(PAY_AMT5), mean_pay6 = mean(PAY_AMT6))
client_train %>% group_by(SEX) %>%
  summarise(mean_limit_bal = mean(LIMIT_BAL), mean_age = mean(AGE),
            mean_bill1 = mean(BILL_AMT1), mean_bill2 = mean(BILL_AMT2),
            mean_bill3 = mean(BILL_AMT3), mean_bill4 = mean(BILL_AMT4),
            mean_bill5 = mean(BILL_AMT5), mean_bill6 = mean(BILL_AMT6),
            mean_pay1 = mean(PAY_AMT1), mean_pay2 = mean(PAY_AMT2),
            mean_pay3 = mean(PAY_AMT3), mean_pay4 = mean(PAY_AMT4),
            mean_pay5 = mean(PAY_AMT5), mean_pay6 = mean(PAY_AMT6))
client_train %>% group_by(EDUCATION) %>%
  summarise(mean_limit_bal = mean(LIMIT_BAL), mean_age = mean(AGE),
            mean_bill1 = mean(BILL_AMT1), mean_bill2 = mean(BILL_AMT2),
            mean_bill3 = mean(BILL_AMT3), mean_bill4 = mean(BILL_AMT4),
            mean_bill5 = mean(BILL_AMT5), mean_bill6 = mean(BILL_AMT6),
            mean_pay1 = mean(PAY_AMT1), mean_pay2 = mean(PAY_AMT2),
            mean_pay3 = mean(PAY_AMT3), mean_pay4 = mean(PAY_AMT4),
            mean_pay5 = mean(PAY_AMT5), mean_pay6 = mean(PAY_AMT6))
client_train %>% group_by(MARRIAGE) %>%
  summarise(mean_limit_bal = mean(LIMIT_BAL), mean_age = mean(AGE),
            mean_bill1 = mean(BILL_AMT1), mean_bill2 = mean(BILL_AMT2),
            mean_bill3 = mean(BILL_AMT3), mean_bill4 = mean(BILL_AMT4),
            mean_bill5 = mean(BILL_AMT5), mean_bill6 = mean(BILL_AMT6),
            mean_pay1 = mean(PAY_AMT1), mean_pay2 = mean(PAY_AMT2),
            mean_pay3 = mean(PAY_AMT3), mean_pay4 = mean(PAY_AMT4),
            mean_pay5 = mean(PAY_AMT5), mean_pay6 = mean(PAY_AMT6))

## Visualization
#Lets use plot for given credit and age against :
#Payment status
client_train %>% ggplot(aes(x = default.payment.next.month, y = LIMIT_BAL)) + geom_boxplot() + 
  ylab("Amount of Given Credit") + xlab("Payment Status") + ggtitle("Given Credit VS Payment Status")
client_train %>% ggplot(aes(x = default.payment.next.month, y = AGE)) + geom_boxplot() + 
  ylab("Age") + xlab("Payment Status") + ggtitle("Age VS Payment Status")
#Sex
client_train %>% ggplot(aes(x = SEX, y = LIMIT_BAL)) + geom_boxplot() + 
  ylab("Amount of Given Credit") + xlab("Sex") + ggtitle("Given Credit VS Sex")
client_train %>% ggplot(aes(x = SEX, y = AGE)) + geom_boxplot() + 
  ylab("Age") + xlab("Sex") + ggtitle("Age VS Sex")
#Education
client_train %>% ggplot(aes(x = EDUCATION, y = LIMIT_BAL)) + geom_boxplot() + 
  ylab("Amount of Given Credit") + xlab("Education") + ggtitle("Given Credit VS Education")
client_train %>% ggplot(aes(x = EDUCATION, y = AGE)) + geom_boxplot() + 
  ylab("Age") + xlab("Education") + ggtitle("Age VS Education")
#Marriage
client_train %>% ggplot(aes(x = MARRIAGE, y = LIMIT_BAL)) + geom_boxplot() + 
  ylab("Amount of Given Credit") + xlab("Education") + ggtitle("Given Credit VS Marital Status")
client_train %>% ggplot(aes(x = MARRIAGE, y = AGE)) + geom_boxplot() + 
  ylab("Age") + xlab("Education") + ggtitle("Age VS Marital Status")
#From previous plot, we have seen that some predictors such as level of education against given credit
#and age against education may have association. Thus using including these predictors together in a model
#such as viewing them in multivariate case is suggested. This is so that we can capture the covariance
#between predictors. Multicollinearity is also a concern that may happen in the case of our previous
#consideration.
#For the numeric variables, we can further study with a more compact statistical methods such as PCA
#that is to capture a denser components as to simplify our model.
#Next, lets explore for the categoric predictors

## Payment status against sex, education, and marriage
client_train %>% group_by(default.payment.next.month) %>%
  count(SEX)
chisq_sex = chisq.test(table(client_train$default.payment.next.month, client_train$SEX)) 
chisq_sex #There is a significant association between sex and payment status.
corrplot(chisq_sex$residuals, is.corr = FALSE)
#From the plot above, there is a positive association between skipping payment with male also female with never skipping
client_train %>% group_by(default.payment.next.month) %>%
  count(EDUCATION)
chisq_edu = chisq.test(table(client_train$default.payment.next.month, client_train$EDUCATION))
chisq_edu #There is a significant association between education and payment status.
corrplot(chisq_edu$residuals, is.corr = FALSE)
#From the plot above, strongest positive association with skipping payment seems to appear on high school 
#level education while other type of education has strongest postive association with never skip a payment.
client_train %>% group_by(default.payment.next.month) %>%
  count(MARRIAGE)
chisq_marry = chisq.test(table(client_train$default.payment.next.month, client_train$MARRIAGE))
chisq_marry #There doesnt seems to be any assocation between marital status and payment status
#From exploring the categoric variables, it is found that sex and marital status seems to have significant
#association with payment status.
#Lets see the association for the three predictors.
library(vcd)
threepred_client = xtabs(~MARRIAGE+SEX+EDUCATION, data = client_train)
ftable(threepred_client)
mantelhaen.test(threepred_client)
chisq.test(table(client_train$EDUCATION, client_train$SEX))
chisq.test(table(client_train$MARRIAGE, client_train$SEX))
chisq.test(table(client_train$MARRIAGE, client_train$EDUCATION))
#Turns out 3 categorical variables have significant association with each other.
#Similar to previous numerical predictors, it is recommended to look at the predictors as a unity.
#Some methods to simplify these categorical variables are available such as correspondence analysis.