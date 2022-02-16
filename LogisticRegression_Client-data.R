### Logistic Regressions on Client Data ###
library(caret)
library(FSelectorRcpp)
library(smbinning)
library(glmnet)
library(factoextra)
library(InformationValue)

## Data preparation
# Load data
client_train = read_excel("~/Github/Machine-Learning/client-data.xlsx", sheet = "client_train")
client_test = read_excel("~/Github/Machine-Learning/client-data.xlsx", sheet = "client_test")
# Cleaning the data
View(client_train)
client_train = client_train[,-1] #Remove id as it serves no purpose in exploration nor analysis
client_test = client_test[,-1]
str(client_train) #Data structure, variables type
#Lets use client_train1 as the data train with all used variables.
#To select these used variables, we'll use feature selection.
client_train1 = as.data.frame(client_train)
client_test1 = as.data.frame(client_test)
names(client_train1)[24] <- "default_payment"
#Re-scale the Pay_ variables into just 3 categories : -1, 0, 1
client_train1[, 6:11] = as.data.frame(sapply(data.frame(client_train1[,6:11]), function(x) replace(x, x>1, 1)))
client_train1[, 6:11] = as.data.frame(sapply(data.frame(client_train1[,6:11]), function(x) replace(x, x< -1, -1)))
client_test1[, 6:11] = as.data.frame(sapply(data.frame(client_test1[,6:11]), function(x) replace(x, x>1, 1)))
client_test1[, 6:11] = as.data.frame(sapply(data.frame(client_test1[,6:11]), function(x) replace(x, x< -1, -1)))
cat_client = c("SEX", "EDUCATION", "MARRIAGE", "PAY_0", "PAY_2", "PAY_3", "PAY_4", "PAY_5", "PAY_6")
cont_client = c("LIMIT_BAL", "AGE", "BILL_AMT1", "BILL_AMT2", "BILL_AMT3", "BILL_AMT4", "BILL_AMT5", 
                "BILL_AMT6", "PAY_AMT1", "PAY_AMT2", "PAY_AMT3", "PAY_AMT4", "PAY_AMT5", "PAY_AMT6")
#Segregate based on the type of predictors : categorical and continuous
client_train1[cat_client] = as.data.frame(lapply(client_train1[cat_client], as.factor)) #Change into factor/categorical variable
client_train1$default_payment = as.numeric(client_train1$default_payment) #This is done because of next feature selection
str(client_train1)
#Goal is to analyze characteristics of 2 types of clients based on payment status : 
#(1) Skip a payment at least once (bad)
#(2) Never skip a payment (good)
#Then build a model to predict clients' payment status on the test data.

## Variable selection for logistic model
# Information Value (IV)
#Shows how much information a variable contain to characterizing the target variable.
iv_table = data.frame(VARS = c(cat_client, cont_client), IV = numeric(23))

for(i in cat_client) {
  smb = smbinning.factor(client_train1, y = 'default_payment', x = i)
  if(class(smb) != 'character'){
    iv_table[iv_table$VARS == i, "IV"] = smb$iv
  }
} #Compute IV for categorical predictors

for(j in cont_client){
  smb = smbinning(client_train1, y = 'default_payment', x = j)
  if(class(smb) != 'character'){
    iv_table[iv_table$VARS == j, "IV"] = smb$iv
  }
} #Compute IV for continuous predictors
iv_table = iv_table[order(-iv_table$IV),]
iv_table
#From Siddqi(2006) IV < 0.02 is not useful for logistic model. As in these variables by itself is not
#useful to characterize the difference between good and bad clients.
#But for now we'll only eliminate SEX and MARRIAGE as we'll do extraction method for BILL_AMT and PAY_AMT
#Alternatively from above, we can also use information gain to decide which variable to eliminate.
information_gain(default_payment~., client_train1) #We'll still only eliminate SEX and MARRIAGE
client_train1 = client_train1[,-c(2,4)]
client_train1$default_payment = as.factor(client_train1$default_payment)

## Variable extraction
#PCA
#We'll proceed on using extraction method. Intuitively, we'll simplify : bill_amt and pay_amt as the continuous
#variables using PCA first.
pc_train = prcomp(client_train1[,10:21], scale= TRUE)
fviz_eig(pc_train) #We'll keep 3 components so we maintain ~70% explained variance 
client_train1 = cbind(client_train1[,-c(10:21)], pc_train$x[,1:3])

#Before we fit the model into our data, we'll split training data into the actual training data to build
#the model upon while the other is used for test. This is so we can do cross validation and improving 
#model performance later on.
table(client_train1$default_payment) #Since we have unbalanced response, we'll split on each category of client
input_1 = client_train1[which(client_train1$default_payment == 1), ]
input_0 = client_train1[which(client_train1$default_payment == 0), ]
set.seed(12345)
idx_1 = sample(1:nrow(input_1), 0.7*nrow(input_1))
idx_0 = sample(1:nrow(input_0), 0.7*nrow(input_0))
model_train = input_0[idx_0,]
model_train = rbind(model_train, input_1[idx_1,])
model_test = input_0[-idx_0,]
model_test = rbind(model_test, input_1[-idx_1,])

## Fitting the model
#After previous feature selection, we proceed with predictors : Limit_bal, Edu, Age, Pay_, Bill_Amt, Pay_Amt
log_client1 = glm(default_payment~., data = model_train, family = "binomial")
summary(log_client1)
#Warning of algorithm did not converge, we'll use regularization to penalized the fit
#We'll proceed with penalized logistic regression by applying lasso,ridge, and elastic-net regression
penalized_model = model.matrix(~.-1, model_train[,-10])
log_client2 = cv.glmnet(x = penalized_model, y = model_train$default_payment, 
                     intercept = FALSE, family = "binomial", alpha = 0, type.measure = "deviance") #Ridge
list.of.fits = list()
for(i in 0:5){
  fit.name = paste0("alpha", i/5)
  list.of.fits[[fit.name]] = cv.glmnet(x = penalized_model, y = model_train$default_payment,
                                       intercept = FALSE, family = "binomial", alpha = i/5, 
                                       type.measure = "deviance")
} #Fitting 6 models, alpha = 0 being ridge, alpha = 1 being lasso, and others being elastic-net
#Lets do cross validation on the model_test data
penalized_test = model.matrix(~.-1, model_test[,-10])
pred = predict(log_client2, newx = penalized_test, type = "response")
opt_cut0 = optimalCutoff(model_test$default_payment, pred)
misClassError(model_test$default_payment, pred, threshold = opt_cut0)
plotROC(model_test$default_payment, pred)
Concordance(model_test$default_payment, pred)
sensitivity(model_test$default_payment, pred, threshold = opt_cut0)
specificity(model_test$default_payment, pred, threshold = opt_cut0)
confusionMatrix(model_test$default_payment, pred, threshold = opt_cut0)

result_client2 = data.frame()
for(i in 0:5){
  fit.name = paste0("alpha", i/5)
  pred = predict(list.of.fits[[fit.name]], newx = penalized_test, type = "response")
  opt_cut = optimalCutoff(model_test$default_payment, pred)
  misclas = misClassError(model_test$default_payment, pred, threshold = opt_cut)
  concord = Concordance(model_test$default_payment, pred)
  Sens = sensitivity(model_test$default_payment, pred, threshold = opt_cut)
  Spec = specificity(model_test$default_payment, pred, threshold = opt_cut)
  temp = data.frame(Alpha = i/5, Misclasification.rate = misclas , Concordance = concord$Concordance, 
                    Sensitivity = Sens, Specificity = Spec)
  result_client2 = rbind(result_client2, temp)
}
result_client2
#From result above, we can get 100% classification on model_test using elastic-net and lasso regression on the logistic model
#thus, we chose either of these as the "proposed" model.
#To try resampling the training and see if the model behave consistently, below is the code
idx_1 = sample(1:nrow(input_1), 0.7*nrow(input_1))
idx_0 = sample(1:nrow(input_0), 0.7*nrow(input_0))
model_train = input_0[idx_0,]
model_train = rbind(model_train, input_1[idx_1,])
model_test = input_0[-idx_0,]
model_test = rbind(model_test, input_1[-idx_1,])
penalized_model = model.matrix(~.-1, model_train[,-10])
penalized_test = model.matrix(~.-1, model_test[,-10])
list.of.fits = list()
for(i in 0:5){
  fit.name = paste0("alpha", i/5)
  list.of.fits[[fit.name]] = cv.glmnet(x = penalized_model, y = model_train$default_payment,
                                       intercept = FALSE, family = "binomial", alpha = i/5, 
                                       type.measure = "deviance")}
result_client2 = data.frame()
for(i in 0:5){
  fit.name = paste0("alpha", i/5)
  pred = predict(list.of.fits[[fit.name]], newx = penalized_test, type = "response")
  opt_cut = optimalCutoff(model_test$default_payment, pred)
  misclas = misClassError(model_test$default_payment, pred, threshold = opt_cut)
  concord = Concordance(model_test$default_payment, pred)
  Sens = sensitivity(model_test$default_payment, pred, threshold = opt_cut)
  Spec = specificity(model_test$default_payment, pred, threshold = opt_cut)
  temp = data.frame(Alpha = i/5, Misclasification.rate = misclas , Concordance = concord$Concordance, 
                    Sensitivity = Sens, Specificity = Spec)
  result_client2 = rbind(result_client2, temp)}
result_client2