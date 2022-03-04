### Clustering of car data sales ###
library(haven)
library(ggplot2)
library(dplyr)
library(reshape2)
library(factoextra)
library(cluster)
library(mclust)
library(dendextend)
library(tidyverse)
library(fpc)

## Preparing data
model_car_sales = read_sas("~/Github/Machine-Learning/model_car_sales.sas7bdat", NULL)
model_car_sales = na.exclude(model_car_sales)
model_car_sales[,c(1,4:7)]
car_sales_plot = melt(model_car_sales[,c(3:7)])
ggplot(car_sales_plot, aes(x = DEALER_CODE, value, color = variable)) + geom_point()
#Since UTE type has very low average sales across dealers, we'll ignore it from clustering feature
#This is becausewe are trying to find segments with lowest number with each cluster characterize by which type
#of car(s) is their highest sales. Thus, we'll cluster by three variables : hatch, wagon, adn sedan
#Lets see the data in its principal components
pc_car = prcomp(model_car_sales[,5:7])
fviz_eig(pc_car) #around 80% of total variance can be explained by first 2 components
data.frame(pc_car$x) %>% ggplot(aes(x = PC1, y = PC2)) + geom_point() #doesnt show clear pattern of segmentation
#we'll try clustering without PC

## Clustering 
#Kmeans
k3m = kmeans(model_car_sales[,5:7], centers = 3)
fviz_cluster(k3m, data = model_car_sales[,5:7]) 
#we'll try to use some methods to find optimum k in kmeans
#Elbow (WSS)
set.seed(222)
kmax = 10
wss = sapply(1:kmax, function(k){
               kmeans(model_car_sales[,5:7], k, nstart = 50, iter.max = 100)$tot.withinss
             })
wss
plot(1:kmax, wss, type = "b", pch =19, xlab = "Number of clusters", ylab = "Total within-clusters sum squares")
#4 seems to be good choice based on plot above
k4m = kmeans(model_car_sales[,5:7], centers = 4, iter.max = 100)
fviz_cluster(k4m, model_car_sales[,5:7]) #we'll keep this result as k4m column
model_car_sales1 = cbind(model_car_sales[,c(3,5:7)], Km_cluster = k4m$cluster)

#Hierarchical
h_clus = hclust(dist(model_car_sales[,5:7]), method = "complete")
#we'll asses which method is the best for hierarchical agglomerative
meth = c("average", "single", "complete", "ward")
names(meth) = c("average", "single", "complete", "ward")
ac = function(x){
  agnes(model_car_sales[,5:7], method =x)$ac
}
map_dbl(meth, ac) #ward is the mest method for agglomerative hierarchical clustering
h_clus2 = hclust(dist(model_car_sales[,5:7]), method = "ward.D2")
plot(h_clus2, hang = -1) #4 clusters seems to be good option
rect.hclust(h_clus2, 4)
hier_car = cutree(h_clus2, k = 4) 
fviz_cluster(list(data = model_car_sales[,5:7], cluster = hier_car))
#we'll try to look for optimal number of cluster using several methods
#Elbow
fviz_nbclust(model_car_sales[,5:7], FUN = hcut, method = "wss")
#Silhouette
fviz_nbclust(model_car_sales[,5:7], FUN = hcut, method = "silhouette")
#Gap statistic
gap_hier_car = clusGap(model_car_sales[,5:7], FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_hier_car)
#Based on the previous 3 methods, 4 cluster is the best number of cluster for hierarchicalclustering

#Model-based clustering with normal mixture modelling
#We choose number of cluster k based on bayesian inference criterion (BIC)
k_clus = Mclust(model_car_sales[,5:7], G = 1:10, modelNames = mclust.options("emModelNames"))
summary(k_clus)
k_clus$BIC
plot(k_clus) #optimum k is 9 based on BIC
fviz_cluster(k_clus)

#Density-Based Spatial Clustering with Noise (DBSCAN)
