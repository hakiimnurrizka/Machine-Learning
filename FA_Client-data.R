### Client Data Analysis using Factor Analysis ###
library(readxl)
library(psych)
library(corrplot)
library(ggplot2)
library(car)
library(nFactors)
library(REdaS)
library(dplyr)

##Prepare the data
client_train = read_excel("~/Github/Machine-Learning/client-data.xlsx", 
                          sheet = "client_train")
client_train1 = client_train[,-c(1:12)]

##Correlation and sampling adequation
#This acts as pre-analysis to identify whether factor analysis is a valid method to analyze the data or not
corrplot(cor(client_train1[,-13]), method = "number") #correlation matrix
KMO(cor(client_train1[,-13])) #Kaiser-Meyer-Olkin to measure factorability
KMOS(client_train1[,-13]) #Another function for KMO and MSA
cortest.bartlett(client_train1[,-13])
det(cor(client_train1[,-13]))
det(cor(client_train1[,-c(6,12,13)]))
View(client_train1[,-c(6,12,13)])

##Extraction using PCA
library(factoextra)
pc_client = prcomp(client_train1[,-c(6,12,13)], scale. = T)
View(pc_client$rotation)
fviz_eig(pc_client)

#Based on previous PCA, lets move on with the premise that there are 2 factors that make up client data
#We can also determine how many factors using following method
fa_client = fa(client_train1[,-c(6,12,13)], nfactors = ncol(client_train1[,-c(6,12,13)]), rotate = "none")
n_fa = length(fa_client$e.values)
scree_fa = data.frame(factor_n = as.factor(1:n_fa), eigenvalues = fa_client$e.values)
ggplot(scree_fa, aes(x = factor_n, y = eigenvalues, group = 1))+
  geom_point()+geom_line()+xlab("number of factors")+ylab("initial eigenvalue")+
  labs(title = "scree plot", subtitle = "(based on unreduced correlation matrix)") #either using elbow method or by setting certain percentage of total variance (calculate using eigen value)

fa.parallel(client_train1[,-c(6,12,13)]) #Paralel analysis

##Apply factor analysis
fa_client =  fa(r = client_train1[,-c(6,12,13)], nfactors = 2, covar = FALSE, SMC = TRUE,
                 fm = "pa", max.iter = 100, rotate = "varimax")
print(fa_client)

factanal(client_train1[,-c(6,12,13)], factors = 5, rotation ="varimax", scores = c("regression"))
fa.diagram(fa_client)

#The result from factor analysis then can be used to other model such as regression.
#Obviously, we are to use the factor analysis transformation to predict using logistic regression.
#Prepare the test data
client_test = read_excel("~/Github/Machine-Learning/client-data.xlsx", 
                         sheet = "client_train")
client_test.FA = predict.psych(fa_client, old.data = client_train1[,-c(6,12,13)], data = client_test[,c(13:17,19:23)])
client_test2 = as.data.frame(cbind(client_test[,-c(1,7:24)],client_test.FA))

### Logistic model using FA results
library(ROCR)
library(caret)
client_train_new1 = as.data.frame(cbind(client_train[,c(2:6,25)], fa_client$scores))
client_train_new1[,2] = as.factor(client_train_new1[,2])
client_train_new1[,3] = as.factor(client_train_new1[,3])
client_train_new1[,4] = as.factor(client_train_new1[,4])
client_train_new1[,6] = as.factor(client_train_new1[,6])
client_log1 = glm(default.payment.next.month~0++LIMIT_BAL+EDUCATION+MARRIAGE+PA1+PA2, data = client_train_new1, family = binomial)
summary(client_log1)

## Check the accuracy using specificity and sensitiviy
predict_client_stat = predict(client_log1, newdata = client_train_new1[,c(1,3,4,7,8)], type = "response")
predictions = prediction(predict_client_stat, client_train_new1$default.payment.next.month) 
sens = data.frame(x=unlist(performance(predictions, "sens")@x.values), 
                  y=unlist(performance(predictions, "sens")@y.values))
spec = data.frame(x=unlist(performance(predictions, "spec")@x.values), 
                  y=unlist(performance(predictions, "spec")@y.values))

sens %>% ggplot(aes(x,y)) + 
  geom_line() + 
  geom_line(data=spec, aes(x,y,col="red")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Specificity")) +
  labs(x='Cutoff', y="Sensitivity") +
  theme(axis.title.y.right = element_text(colour = "red"), legend.position="none")

#get the optimum threshold
sens = cbind(unlist(performance(predictions, "sens")@x.values), unlist(performance(predictions, "sens")@y.values))
spec = cbind(unlist(performance(predictions, "spec")@x.values), unlist(performance(predictions, "spec")@y.values))
cutoff_opt_client = sens[which.min(apply(sens, 1, function(x) min(colSums(abs(t(spec) - x))))), 1]

#Confusion matrix
predict_client_stat = ifelse(predict_client_stat>cutoff_opt_client,1,0)
conf_mat_client = table(predict_client_stat, client_train_new1$default.payment.next.month)
conf_mat_client
specificity(conf_mat_client)
sensitivity(conf_mat_client)

#Apply model and cutoff on test data
predict_client_test2 = ifelse(predict(client_log1, newdata = client_test2[,c(1,3,4,6,7)], type = "response")>cutoff_opt_client,1,0)
client_test2$predict.paymentFA = predict_client_test2
View(head(client_test1[,-c(2, 5:11)], 100))