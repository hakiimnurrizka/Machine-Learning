### Bayesian approach ###
##in regard to describe the difference between frequentist and bayesian approach
##suppose we are to analyze penalty events. this illustration will also provide 
##links that is the essentials to bayesian methods : prior, likelihood, and posterior

#say we have 10 shots from observation on the penalty kicks
shots = c("goal", "goal", "goal", "miss", "miss",
          "goal", "goal", "miss", "miss", "goal")
shotsNum = as.numeric(shots=='goal')
N = length(shots)                      # sample size
nGoal = sum(shots=='goal')             # number of shots made
nMiss = sum(shots=='miss')             # number of those miss

#the goal is to estimate in bayesian approach the parameter theta that is the 
#probability of scoring a goal.

#when doing this in frequentist approach, one can simply use the binomial distribution
#then try fitting the appropriate parameter so that the mean results close to the
#observed.

#next we proceed on the bayesian estimation, as it stated before, the bayesian approach 
#does not "fit" certain value to parameter, instead we assume certain distribution
#that we belief represent the "initial information".
#to keep it limited and easier to understand, we'll only consider 10 values for our 
#theta estimation which these value falls between 0 to 1.
theta = seq(from=1/(N+1), to=N/(N+1), length=10)

#For prior distribution here we use several options
# triangular as in Kruschke text example
pTheta = pmin(theta, 1-theta)

# uniform
pTheta = dunif(theta)

# beta prior with mean = .5
pTheta = dbeta(theta, 10, 10)

# Normalize so that values sum to 1
pTheta = pTheta/sum(pTheta)

#For now, we are proceeding with the first prior that  is triangular distribution.
#After applying the prior, we'll calculate the likelihood of our data given the prior
#To calculat this we are using the frequentist approach belief of binomial distribution
pDataGivenTheta = choose(N, nGoal) * theta^nGoal * (1-theta)^nMiss

#The last one is to compute the posterior based on previous prior and likelihood
pData = sum(pDataGivenTheta*pTheta)  # marginal probability of the data
pThetaGivenData = pDataGivenTheta*pTheta / pData  # Bayes theorem

#Examine what we've got
options(digits = 4)
bayes1 = as.data.frame(cbind(theta, pTheta, pDataGivenTheta, pThetaGivenData))
colnames(bayes1) = c("Theta", "Prior", "Likelihood", "Posterior")
View(bayes1)
#Select the meanas the parameter estimator
posteriorMean = sum(pThetaGivenData*theta)
posteriorMean #using this approach, given the data, itis found that the bayes estimation for theta yield 0.56 result
par(mfrow = c(2,2))
#the next one will proceed with bigger sample size 
#for prior, we'll use beta(10,10) distribution
theta2 = seq(from=1/(N+1), to=N/(N+1), length=100)
pTheta2 = dbeta(theta2, 10,10)/10
plot(theta2, pTheta2, type = "l", xlab = "Theta", ylab = "Density",
     main = "Prior")
polygon(theta2, pTheta2, col = "blue", border = "blue")

pDataGivenTheta2 = choose(N, nGoal) * theta2^nGoal * (1-theta2)^nMiss
plot(theta2, pDataGivenTheta2, type = "l", xlab = "Theta", ylab = "Density",
     main = "Likelihood")
polygon(theta2, pDataGivenTheta2, col = "red", border = "red")

pData2 = sum(pDataGivenTheta2*pTheta2)  
pThetaGivenData2 = pDataGivenTheta2*pTheta2*10 / pData2  
plot(theta2, pThetaGivenData2, type = "l", xlab = "Theta", ylab = "Density",
     main = "Posterior")
polygon(theta2, pThetaGivenData2, col = "green", border = "green")

par(mfrow = c(1,1))
plot(theta2, pTheta2, type = "l", xlab = "Theta", ylab = "Density",
     main = "Comparison prior-likelihood-posterior")
lines(theta2, pDataGivenTheta2, col = "red")
lines(theta2, pThetaGivenData2, col = "blue")
legend("topleft", legend = c("Prior", "Likelihood", "Posterior"),
       text.col = c(1,2,4), bty = "n")
sum(pThetaGivenData2*theta2/10) #2nd simulation yield 0.533 as the estimation for theta given the data

#Next we'll proceed into a simulation to show a case where ones' prior belief can be
#broadly useful for statistical inference

###Prevalence of certain disease's infection (taken from Hoff, PD (2009))
#Suppose we are interested in the prevalence of certain disease in a city.
#Higher prevalence should be handled with more health precaution.
#Random sample of 20 individuals is taken,checked for the infection
#From this illustration, we are interested in theta as the proportion of the infected.
#If theta is known, it's reasonable enough to use the binomial distribution such that
#Y | theta ~ bin(20, theta)
#Lets simulate this distribution given several values for theta

##Simulating sampling model
n = 20
theta = c(0.05,0.1, 0.2)
J = length(theta)
prob = array(dim=c(n,J))

for (j in 1:J){
  for (i in 1:n){
    prob[i,j]<-dbinom(i,n,theta[j])
    i2<-i+.1
    i3<-i2+.1}}

plot(prob[,1], ylab="Probability of getting infected", xlab="Number of infected",
     main = "Comparison of binomial sampling model")
segments(1:n, 0, 1:n, prob[,1])
#Adding the prob for theta2 = .1
points(prob[,2], col="red")
segments(1:n, 0, 1:n, prob[,2], col="red")
#Adding the prob for theta3 = .2
points(prob[,3], col="blue")
segments(1:n, 0, 1:n, prob[,3], col="blue")
colour = c("black", "red", "blue")
labels = c("0.05", "0.1","0.2")
legend("topright", inset=.05, title="Theta",labels, lwd=2, lty=c(1, 1, 1), col=colour, bty="n")

#Next, suppose there is a study saying the infection rate for the disease is around
#0.05 to 0.2 with an average 0.1 between observed cities.
#This information suggests that the prior distribution that is to be used
#should have substantial amount of probability for the interval (0.05 , 0.2) and an
#expected value of theta is around 0.1.
#Such distribution os of course not as unique and in fact there are wide range of 
#option for distributions which satisfy the previous constraints.
#Thus, we can choose the prior based on the mathematical convenience.
#In this case, we will proceed with beta as prior distribution, specifically
#beta(2,20) in which the expectation and most probable value is 0.1 and 0.05
#generating a random number from a prior beta(a,b) distribution
nrand = 100000
a = 2
b = 20
dt = rbeta(nrand,a,b)
plot(density(dt))

#calculating the probability of .05<theta<.2
pbeta(.2,a,b)-pbeta(0.05,a,b)
#calculating the probability of theta<.1
pbeta(.1,a,b)

#generating a random number from a posterior 
#beta(a1,b1) distribution
#for sample of size 20(n) with 0(y) success
y = 0; n = 20
a1 = a+y 
b1 = b+n+y
dt.post = rbeta(nrand,a1,b1)
mean(dt.post) #posterior mean of theta

#Graph overly prior and posterior density
plot(density(dt), ylim=c(0,20)) #prior density
lines(density(dt.post), col="red")  #posterior density
#Add legend
legend("topright", legend=c("Prior", "Posterior"),
       col=c("black", "red"), lty=1, cex=0.8, bty="n")