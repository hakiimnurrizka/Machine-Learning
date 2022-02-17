### Categorical data analysis ###
## Based on vcd package
library(vcd)
library(dplyr)
library(lessR)
library(plotrix)
library(corrplot)
library(ca)

## Data preparation
# Load data and picking the categorical variables
client_train = read_excel("~/Github/Machine-Learning/client-data.xlsx", sheet = "client_train")
client_train1 = client_train[,c(2:4, 6:11, 24)]
names(client_train1)[10] <- "default_payment"
client_train1[, 4:9] = as.data.frame(sapply(data.frame(client_train1[,4:9]), function(x) replace(x, x>1, 1)))
client_train1[, 4:9] = as.data.frame(sapply(data.frame(client_train1[,4:9]), function(x) replace(x, x< -1, -1)))
client_train1[,4:10] = as.data.frame(lapply(client_train1[,4:10], as.factor))
str(client_train1)

## Chart and plot
#For categorical variables, we can use pie chart to see the proportion of each category
PieChart(SEX, data = client_train1, hole = 0.1, fill = c("lightblue", "pink"), main = "Sex Proportion")
pie3D(table(client_train1$EDUCATION),  labels = c("Graduate school", "High school", "Others", "University"),
      labelcex = 1, col = c("green", "blue", "red", "purple"))
barplot(table(client_train1$MARRIAGE), main="Count of each marital status", xlab="Marriage status",  
        ylab="Count", border="black", col="lightblue")
barplot(table(client_train1$SEX, client_train1$default_payment),
        main = "Default status for each sex", xlab = "Default Payment Status", col = c("red","green")) + 
  legend("topright", c("Female","Male"), fill = c("red","green"))
## Frequency and table form
#We'll create the frequency for for each non-PAY_ predictors and the default_payment status
table(client_train1$SEX, client_train1$default_payment) #The cross-table
freq_sex = data.frame(expand.grid(sex=c("female", "male"), default_status=c("0", "1")), 
                      count=c(9778,6167,4723,3332)) #manual
mat_edu = table(client_train1$EDUCATION, client_train1$default_payment)
dimnames(mat_edu) = list(Education = c("Graduate school", "High school", "Others", "University"), default_status = c("0", "1"))
View(mat_edu) #From the cross-table
freq_edu = data.frame(expand.grid(education=c("High school", "University", "Graduate school", "Others"), 
                                  default_status=c("0", "1")), 
                      count=c(2430, 7309, 5896, 310, 1466, 3929, 2595, 65))
mat_marry = table(client_train1$MARRIAGE, client_train1$default_payment)
dimnames(mat_marry) = list(marital_status = c("Married", "Others", "Single"), default_status = c("0", "1"))
mat_marry
#3-way table
rm(sex.edu)
sex.marry = xtabs(~SEX+default_payment+MARRIAGE, data = client_train1)
ftable(sex.marry)

## Independence test of payment status against sex, education, and marriage
#2way
chisq_sex = chisq.test(table(client_train$default.payment.next.month, client_train$SEX)) 
chisq_sex #There is a significant association between sex and payment status.
corrplot(chisq_sex$residuals, is.corr = FALSE)
#From the plot above, there is a positive association between skipping payment with male also female with never skipping
chisq_edu = chisq.test(table(client_train$default.payment.next.month, client_train$EDUCATION))
chisq_edu #There is a significant association between education and payment status.
corrplot(chisq_edu$residuals, is.corr = FALSE)
#From the plot above, strongest positive association with skipping payment seems to appear on high school 
#level education while other type of education has strongest postive association with never skip a payment.
chisq_marry = chisq.test(table(client_train$default.payment.next.month, client_train$MARRIAGE))
chisq_marry
#There doesnt seems to be any assocation between marital status and payment status
#3way
mantelhaen.test(sex.marry, conf.level = )
sex.marry
oddsratio(sex.marry, log = FALSE) #comparison of odds ratio for each marriage status
#Correspondence
#Another way to explore the relationship between row and collumn
ca_sex = ca(margin.table(sex.marry, c(1,2)))
ca(margin.table(sex.marry, c(1,2)))
plot(ca_sex)
