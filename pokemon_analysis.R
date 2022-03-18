### Pokemon Data Analysis with clustering and discriminant analysis###
library(readxl)
library(factoextra)
library(dplyr)
library(corrplot)
library(psych)
library(ggplot2)

##Data
#Preprocessing done by eliminating individual id/name and sparse variables 
pokemon = read_excel("~/Github/Multivariate-Statistics-R/pokemon.xlsx")
View(pokemon)
str(pokemon)
# Lets try to cluster the pokemon based on their capability to fight against different type of monsters

##Variable extraction
pc_poke = prcomp(pokemon[,1:18], scale. = TRUE)
fviz_eig(pc_poke)
summary(pc_poke) #it looks like the information is quite spreaded among several first variables
fviz_pca_var(pc_poke, repel = TRUE, col.var = "cos2")
#We'll try to use factor analysis and assume the capalibity to fight against each type of monster is influenced
#by certain factor(s)
#First lets take a look at paralel analysis
fa.parallel(pokemon[,1:18]) #Based on parallel analysis the recommended number of components for PCA is 6
#while for factor analysis it is suggested that 7 factors are needed to extract the information from the data.
#Judging by this, PCA is more preferred since we ended up in a simpler data.
#But for the sake of studying, we'll proceed on applying factor analysis on data.
#We'll only use exploratory factor analysis.
corrplot(cor(pokemon[,1:18]))
KMO(pokemon[,1:18]) 
cortest.bartlett(pokemon[,c(1:18)])
det(cor(pokemon[,c(1,3,5,7,12,15,16,18)]))
#Judging by above result, it doesnt look like factor analysis is a good approach to extract the information
#from the data. This is mostly due to the sparse or a lot of elements in the correlation is close to zero.
#We propose to eliminate variables by MSA.
KMO(pokemon[,c(3,5,7,12,15,16,18)])
cortest.bartlett(pokemon[,c(3,5,7,12,15,16,18)])
det(cor(pokemon[,c(3,5,7,12,15,16,18)]))
#We find a decent KMO value close to 0.6
fa_poke = fa(pokemon[,c(3,5,7,12,15,16,18)], nfactors = ncol(pokemon[,c(3,5,7,12,15,16,18)]), rotate = 'none')
n_fa_poke = length(fa_poke$values)
scree_fa = data.frame(
  n_fa = as.factor(1:n_fa_poke),
  eigenval = fa_poke$values
)
scree_fa %>% ggplot(aes(n_fa, eigenval, group = 1)) + geom_point() + geom_line() + xlab("number of factors") +
  ylab("initial eigen values") + labs(title = "Scree Plot", subtitle = "Based on unreduced correlation matrix")
fa.parallel(pokemon[,c(3,5,7,12,15,16,18)])
#After analyzing the data using factor analysis, we still find that PCA is more preferred to simplify the data.
#Next, we save the transformed data by using first 6 components of PCA (from previous result on PCA)
pokemon1 = as.data.frame(pc_poke$x[,1:6])
pokemon1 = na.omit(pokemon1)

##Clustering
#First, lets take a look at distance matrix and how conventional hierarchical clustering will treat the data.
hier1_poke = hclust(dist(pokemon1))
fviz_nbclust(pokemon1, FUNcluster = hcut, method = "silhouette")
pokemon1$hier1_clust = cutree(hier1_poke, 2)-1
mean(pokemon$is_legendary == pokemon1$hier1_clust)#good accuracy to predict legendary/not.
#But there is a catch, legendary pokemon is quite rare which means the classification accuracy above might just
#come from accuracy to predict non-legendary.
table(pokemon$is_legendary, pokemon1$hier1_clust)
#And yes it can be seen that the classification rule using basic hierarchical clustering on transformed 
#"against_" variables only has 20% correctness on classifying legendary pokemon.
#Of course we can improve it by adding more information from other variables we didnt include.
#But lets do a bit more exploring on clustering model for the transformed
#Lets take a look using kmeans clustering
km1_poke = kmeans(pokemon1[,1:6], centers = 2)
pokemon1$km1_clust = (km1_poke$cluster-2)*-1
mean(pokemon$is_legendary == pokemon1$km1_clust)
table(pokemon$is_legendary, pokemon1$km1_clust)
#Classification of legendary using kmeans cluster seems to perform worse than hierarchical clustering
#Lets check how it actually handle the clustering
fviz_nbclust(pokemon1[,1:6], FUNcluster = kmeans, method = "silhouette", nboot = 1000)
#Optimal cluster is 8, which means kmeans has different result on clustering the transformed data than
#the previous hierarchical clustering.