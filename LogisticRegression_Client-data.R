### Logistic Regressions on Client Data ###
library(caret)
library(randomForest)
library(FSelector)
library(FSelectorRcpp)
library(smbinning)

## Load the data
client_train = read_excel("~/Github/Machine-Learning/client-data.xlsx", sheet = "client_train")
client_test = read_excel("~/Github/Machine-Learning/client-data.xlsx", sheet = "client_test")

## Cleaning the data
View(client_train)
client_train = client_train[,-1] #Remove id as it serves no purpose in exploration nor analysis
client_test = client_test[,-1]
str(client_train) #Data structure, variables type
# Lets use client_train1 as the data train with all used variables 
# To select these used variables, we'll use feature selection
client_train1 = as.data.frame(client_train)
names(client_train1)[24] <- "default_payment"
# Rescale the Pay_ variables into just 3 categories : -1, 0, 1
client_train1[, 6:11] = as.data.frame(sapply(data.frame(client_train1[,6:11]), function(x) replace(x, x>1, 1)))
client_train1[, 6:11] = as.data.frame(sapply(data.frame(client_train1[,6:11]), function(x) replace(x, x< -1, -1)))
cat_client = c("SEX", "EDUCATION", "MARRIAGE", "PAY_0", "PAY_2", "PAY_3", "PAY_4", "PAY_5", "PAY_6")
cont_client = c("LIMIT_BAL", "AGE", "BILL_AMT1", "BILL_AMT2", "BILL_AMT3", "BILL_AMT4", "BILL_AMT5", 
                "BILL_AMT6", "PAY_AMT1", "PAY_AMT2", "PAY_AMT3", "PAY_AMT4", "PAY_AMT5", "PAY_AMT6")
#Segregate based on the type of predictors : categorical and continuous
client_train1[cat_client] = as.data.frame(lapply(client_train1[cat_client], as.factor)) #Change into factor/categorical variable
client_train1$default_payment = as.numeric(client_train1$default_payment)
client_train1$default_payment = client_train1$default_payment-1
str(client_train1)
#Goal is to analyze characteristics of 2 types of clients based on payment status : 
#(1) Skip a payment at least once (bad)
#(2) Never skip a payment (good)
#Then build a model to predict clients' payment status on the test data.

## Variable selection for logistic model
#IV
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

information_gain(default.payment.next.month~., client_train)
