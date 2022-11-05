###Some Basic Hypothesis Tests###
##Data
aquifer = read.csv("aquifer.csv")
#This data is about characteristics of aquifer based on the type of the limestone
#The characteristics are defined as : density of the layer, porosity of the layer (ratio of the volume of pores to the volume of bulk rock),
#insoluble residual contained in the limestone, and grain axis length 
View(aquifer)


##Descriptive Statistics
library(dplyr)
group_by(aquifer, type) %>%
  summarise(
    count = n(),#Count number of observation for each type of limestone
    mean_density = mean(density, na.rm = TRUE),#Compute mean of density for each type of limestone
    sd_density = sd(density, na.rm = TRUE),#Compute SD of density for each type of limestone
    mean_porosity = mean(porosity, na.rm = TRUE),#Compute mean of porosity for each type of limestone
    sd_porosity = sd(porosity, na.rm = TRUE),#Compute SD of porosity for each type of limestone
    mean_residue = mean(residue, na.rm = TRUE),#Compute mean of residue for each type of limestone
    sd_residue = sd(residue, na.rm = TRUE),#Compute SD of residue for each type of limestone
    mean_glength = mean(glength, na.rm = TRUE),#Compute mean of grain length for each type of limestone
    sd_glength = sd(glength, na.rm = TRUE),#Compute SD of grain length for each type of limestone
  )


##Visualization
library(ggplot2)
library(ggpubr)
#Density with box plot
ggboxplot(aquifer, x = "type", y = "density", 
          color = "type",
          order = c("Axeann", "Bellefonte", "Nittany", "Stonehenge"),
          ylab = "Density", xlab = "Limestone Type")
#Porosity with box plot
ggboxplot(aquifer, x = "type", y = "porosity", 
          color = "type",
          order = c("Axeann", "Bellefonte", "Nittany", "Stonehenge"),
          ylab = "Porosity", xlab = "Limestone Type")
#Residue with mean plot
ggline(aquifer, x = "type", y = "residue", 
       add = c("mean_se"), 
       order = c("Axeann", "Bellefonte", "Nittany", "Stonehenge"),
       ylab = "Unsoluble Residue", xlab = "Limestone Type")
#Grain length with mean plot
ggline(aquifer, x = "type", y = "glength", 
       add = c("mean_se", "jitter"), 
       order = c("Axeann", "Bellefonte", "Nittany", "Stonehenge"),
       ylab = "Grain Axis Length", xlab = "Limestone Type")


##Hypothesis Testing##

###One variable

##One population
#Lets focus only on one of the variable that is density and one population that is limestone type Axeann
Axeann = aquifer[aquifer$type == "Axeann",] #Filter data by type of limestone
View(Axeann)
den.ax = Axeann[,2]#Choose second collumn only
#Suppose we have a claim that says density of aquifer land which limestone type is Axeann has an average of 2.7g/cm3
#To test whether this claim is true or not, we can do 2 type of parametric test.
#First is is test of mean difference using Z statistics, secondly using one population t test.
#For Z statistics, we can only this method IF we know the true value of population's variance.
#Lets say the value of population varince for bulk density of Axeann aquifer is 0.003g/cm3
xbar = mean(den.ax)#Sample mean from density of axeann aquifer
mu = 2.7#Hypothesized population
n = length(den.ax)#Number of observations in sample
z.stat = sqrt(n)*(xbar - mu)/0.003#Z statistics
z.stat
p.val.z = 2*pnorm(abs(z.stat), lower.tail = FALSE)#P value for 2 way hypothesis test for the Z statistics
p.val.z#Atomic number, very close to zero
library(BSDA)
z.test(den.ax, mu = 2.7, sigma.x = 0.003)
#From the test above, null hypothesis of mu = 2.7 at the significance level alpha = 0.05 is rejected. 
#Thus, we can conclude that the average bulk density of axeann aquifer is not equal to 2.7g/cm3.
#The second test can be used in the case of unknown true value of population's variance.
t.stat = sqrt(n)*(xbar - mu)/sd(den.ax)#T statistics
t.stat
2*pt(abs(t.stat), df= 8,lower.tail = FALSE)
t.test(den.ax, mu = 2.7)
#From the test above, null hypothesis of mu = 2.7 at the significance level alpha = 0.05 is not rejected. 
#Thus, we can conclude that the average bulk density of axeann aquifer is equal to 2.7g/cm3.

##Two Population
#Suppose now we want compare between 2 type of aquifer limestone, we can do the test for two population
#such as t test for two population.
#Lets consider bulk density for two types of limestone : Axeann and Bellefonte.
#Lets consider equal variance of both sample
den.bel = aquifer[aquifer$type == "Axeann", 2]
ybar = mean(den.bel)
sp = (8*var(den.ax)+8*var(den.bel))/16
t.stat2 = (xbar-ybar)/sqrt(sp*2/9)
2*pt(abs(t.stat2), df= 16,lower.tail = FALSE)
t.test(den.ax, den.bel, var.equal = T)
#From the test above, null hypothesis of mu1 = mu2 at the significance level alpha = 0.05 is not rejected. 
#Thus, we can conclude that the average bulk density of axeann aquifer is equal to average bulk density of bellefonte aquifer.
#Then,consider unequal variance for both sample
t.stat2.2 = (xbar-ybar)/sqrt(var(den.ax)/9+var(den.bel)/9)
v = (var(den.ax)/9+var(den.bel)/9)^2/(((var(den.ax)/9)^2+(var(den.bel)/9)^2)/8)
2*pt(abs(t.stat2), df= v,lower.tail = FALSE)
t.test(den.ax, den.bel, var.equal = F)
#From the test above, null hypothesis of mu1 = mu2 at the significance level alpha = 0.05 is not rejected. 
#Thus, we can conclude that the average bulk density of axeann aquifer is equal to average bulk density of bellefonte aquifer.

##More Than Two Population
#A more general case, that is lets consider we have p population and a varaible to specify some numeric measurement.
#That is for example we have our complete data of aquifer
View(aquifer)
#And we want to compare the bulk density between 4 types of limestone
ggboxplot(aquifer, x = "type", y = "density", 
          color = "type",
          order = c("Axeann", "Bellefonte", "Nittany", "Stonehenge"),
          ylab = "Density", xlab = "Limestone Type")
#We can use one way analysis of variance
anov1 = aov(density ~ type, data = aquifer)
summary(anov1)#ANOVA Table
#With the null hypothesis H0: mu1 = mu2 = mu3 = mu4 at the significance level alpha = 0.05 is rejected.
#Thus, we can conclude that there is at least two types of limestone that have different bulk density mean. 

###Two Variables
##Correlation
#One of the simplest form of analysis between 2 numerical variables is by learning how one variable change in
#accordance to the other variable change.
#Lets take an example with the data aquifer for variable porosity and residue
View(aquifer[,3:4])
plot(aquifer[,3:4])
#There is faint pattern that when the value of porosity increases, the value of residuea also increase
#We can measure how "strong" this pattern/relationship using correlation coeffiecient
x = aquifer[, 3]
y = aquifer[, 4]
r = sum((x-mean(x))*(y-mean(y)))/sqrt(sum((x-mean(x))^2)*sum((y-mean(y))^2))
r
cor(x,y)    
#Correlation coeffiecient values are ranged within [-1,1], that is in this case the correlation between
#porosity and residue fall within "positive side". That is in another words, the correlation is positive
#meaning that the previous identification of both variable change together in the same direction is somewhat
#true.
#Obviously we can do something to attest how significant this claim is using hypothesis test.
#Lets say we claim that the correlation between porosity and residue is significantly positive, that is
#in the hypothesis testing we want to test null H0 : r <= 0 vs alternative H1 : r > 0 (a one way hypothesis test)
cor.test(x,y, alternative = "greater")#We specify our alternative hypothesis
#The t statistics for the coefficient correlation r is 2.54 which correponds to p value of one way test 0.008.
#That is we can reject null at the significance level 0.05 and conclude that the claim of "positive association"
#between porosity and residue of aquifer is true.

##Simple Linear Regression
#A dependence statistical analysis, that is in this method we are not only interested in the pattern of change
#but we are also interested in predicting the value of the variable.
#In this analysis, we first assign "role" to variables, that is : 
#(1) outcome/target/dependent variable which we will use notation y
#(2) explanatory/predictor/independent variable which we will use notation x
#Lets see the scatterplot between porosity and residue once again
plot(aquifer[,3:4])
#Now imagine we want to "fit" a straight line to represent this scatterplot
abline(lm(residue~porosity, data = aquifer))#This line is called regression line
#To create this line, we need to make sure that we accomodate the variation in the data.
#A simpler intuition is that we want to create a line that overall has smallest "total distance" with the points.
#a more technical approach, regression line can be written in the form of line equation with "slope" and "intercept"
#y = b0 + b1 * x, here we call b0 the intercept (where it crosses y-axis) and b1 the slope.
#The problem is that, we build this equation not just simply by connecting the dots, rather we want to 
#find a line that has total minimum distance to each point. 
#Skipping the more complicated technical details, to get this regression line we simply do estimation
#for b0 and b1 using our data. The estimation method is called least square.
b1 = sum((x-mean(x))*(y-mean(y)))/sum((x-mean(x))^2)
b0 = mean(y) - b1*mean(x)
b1
b0
#Thus we get our regression line for predicting residue (y) using porosity (x) is
#y = 6.709685 + 1.663072 * x
#we can use the regression line above to predict residue of aquifer with certain value of porosity.
#Say, we want to predict residue of aquifer with porosity 5
pred = 6.709685 + 1.663072 * 5
pred
#So we predict by using the regression line above that the residue of aquifer with porosity 5 is equal to 15.025%
mod1 = lm(y ~  x)
summary(mod1)
#By using lm function we can also get a more specific details about our regression line, included is quality of
#the regression and significance of the regression. To do prediction using this lm function can be done via
#predict function.
predict(mod1, newdata = data.frame(x=5), interval = "confidence")
#To see the fitted and confidence interval lines
ggplot(aquifer, aes(x=porosity, y=residue)) + geom_point() + stat_smooth(formula = y ~ x, method = "lm")