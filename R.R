install.packages("cluster")    
install.packages("factoextra")
install.packages("ggplot2")
install.packages("BBmisc")
install.packages("plsgenomics")
install.packages("dendextend")
install.packages("circlize")
install.packages("d3heatmap")
install.packages("plsgenomics")
install.packages("cluster")
install.packages("fclus") 
install.packages("FactoMineR")
install.packages("sjPlot")
install.packages("dplyr")
install.packages("stats")
install.packages("readxl")
install.packages("e1071")
install.packages("ggplot2")
install.packages("ggrepel")






#############################################################################
#########################Hierarchical clustering in R (iris dataset)#########
#############################################################################


library(cluster)    
library(factoextra)

?hclust  # we need the distance matrix

?dist # distances between the rows
dim(iris)

head(iris)

tail(iris)

data=iris

data[, 3:4]


#############################################################################
#########################Euclidean distance + complete linkage###############
#############################################################################

dd=dist(data[, 3:4], method = "euclidean")
round(dd,2)

clusters = hclust(dd, method = "complete")
clusters
plot(clusters)


?cutree

cut = cutree(clusters, 3)

cut

table(cut, data$Species)


rect.hclust(clusters , k = 3, border = 2:6)


#############################################################################
#########################Plot the dendrogram#################################
#############################################################################

library(ggplot2)
names(data)
names(data)[3]
names(data)[4]

attach(data)

p = ggplot(data, aes(Petal.Length, Petal.Width))

p + geom_point(aes(colour = factor(Species)), size = 4) + 
  ggtitle("Real Iris Categories")

p = ggplot(data, aes(Petal.Length, Petal.Width))

p + geom_point(aes(colour = factor(cut)), size = 4) + ggtitle("Clustering Results")


#############ADDITIONAL PLOT
fviz_cluster(list(data = data[, 3:4], cluster = cut))




#############################################################################
######################### Minkowski distance + average linkage ##############
#############################################################################

dd=dist(data[, 3:4], method = "minkowski",  p = 3)
round(dd,2)

clusters = hclust(dd,method = "average")
clusters
plot(clusters)

cut = cutree(clusters, 3)
cut
table(cut, data$Species)
rect.hclust(clusters , k = 3, border = 2:6)

#############################################################################
#########################Plot ###############################################
#############################################################################


p = ggplot(data, aes(Petal.Length, Petal.Width))

p + geom_point(aes(colour = factor(Species)), size = 4) + ggtitle("Real Iris Categories")

p = ggplot(data, aes(Petal.Length, Petal.Width))

p + geom_point(aes(colour = factor(cut)), size = 4) + ggtitle("Clustering Results")


#############################################################################
######################### Manhattan distance + simple linkage ##############
#############################################################################

dd=dist(data[, 3:4], method = "manhattan")
dd

clusters = hclust(dd,method = "single")
clusters
plot(clusters)

cut = cutree(clusters, 3)
cut
table(cut, data$Species)
rect.hclust(clusters , k = 3, border = 2:6)

p = ggplot(data, aes(Petal.Length, Petal.Width))

p + geom_point(aes(colour = factor(Species)), size = 4) + ggtitle("Real Iris Categories")

p = ggplot(data, aes(Petal.Length, Petal.Width))

p + geom_point(aes(colour = factor(cut)), size = 4) + ggtitle("Clustering Results")

##########################NORMALIZATION############

#There are two things we must do before starting to cluster datasets with many different variables:
# 1) Scaling: each observation’s feature values are represented as coordinates in n-dimensional space (n is the number of features) and then the distances between these coordinates are calculated. If these coordinates are not normalized, then it may lead to false results.
#Alternatives:
###############################  Min-max normalization   
standardize <- function(x){(x-min(x))/(max(x)-min(x))} 

#or the function 

#normalize()
###############################  Standardisation
#x(s)=x(i)-mean(x)/sd(x) 

#or the function 

#scale()
library(BBmisc)
View(data)
nor=data[,-5]
normalize(nor)


#2) Missing Value imputation: 
#if there are missing values, we have two alternatives:

#Missing values imputation:
  ?impute

#########Delete observations with missing values:
  ?na.omit # look at the slides of the Statistical Learning course
any(is.na(data)) #check if any


##########################################################################
##########################Clustering patients affected by leukemia#######
##########################################################################

install.packages("plsgenomics")
library(plsgenomics)

data(leukemia)    # patient x genes

attach(leukemia)
leukemia$gene.names

class(leukemia)  # not a dataframe
dim(leukemia$X)
head(leukemia$X)
leukemia$X[,1]

dd=dist(leukemia$X, method = "manhattan")
dd

clusters = hclust(dd,method = "average")
clusters
plot(clusters)

cut = cutree(clusters, 2)
cut
table(cut, leukemia$Y)
rect.hclust(clusters, k = 2, border = 2:3)
###################

leukemia$Y==cut
sum(leukemia$Y==cut)/length(leukemia$Y)


#NOTE:
#This is trivial check.
#We are assuming that we do not know the number of groups.
#So, we are imagining we do not know leukemia$Y.
#If we knew, we could do supervised classification!
  
#NOTE n.2
#We do not plot the graph because there are too many genes.
#Possible solution: PCA and plot of the first and second PCs.


##################################################################################
#####################Plot groups using the first 2 dimensions#####################
##################################################################################

fviz_cluster(list(data = leukemia$X, cluster = cut))


#################################################################################
############################Divisive hierarchical clustering in R################
#################################################################################


??diana

data=iris
data[, 3:4]

dd=dist(data[, 3:4], method = "euclidean")
round(dd,2)

clusters <- diana(data)

print(clusters)

plot(clusters)

###########################################################################
########################Selecting the number of groups#####################
###########################################################################

###################Elbow method########################
fviz_nbclust 
data=data[, 3:4]   # 2 variables from iris dataset
p1 <- fviz_nbclust(data, FUN = hcut, method = "wss", 
                   k.max = 10) +
  ggtitle("(A) Elbow method")

p1
###################Silhouette method########################

p2 <- fviz_nbclust(data, FUN = hcut, method = "silhouette", 
                   k.max = 10) +
  ggtitle("(B) Silhouette method")

gridExtra::grid.arrange(p1, p2, nrow = 1)

p2

###################################VISUALIZATION
data("leukemia")
##########################THE HEATMAP
?heatmap
heatmap(leukemia$X,scale = "column",
        xlab = "Genes", ylab =  "Patients",
        main = "Leukemia heatmap")


###########################CIRCLE DENDROGRAM
data=leukemia$X
library(dendextend)
library(circlize)
?as.dendrogram
dd=dist(data[, 3:4], method = "euclidean")
clusters = hclust(dd,method = "ave")
dend <- as.dendrogram(clusters)
dend <- color_branches(dend, k=3) #, groupLabels=iris_species)
circlize_dendrogram(dend)




###########################################################################
###################################K-means-TWO PREDICTORS##################
#########################################################################

data=iris[,3:4]
?scale
data=scale(data)
data=as.data.frame(data)
head(iris)
head(data)

?kmeans
kk=kmeans(data, centers = 3, iter.max = 10, nstart = 1)
kk
kk$cluster
kk$iter
kk$tot.withinss

p = ggplot(data, aes(Petal.Length, Petal.Width))+
  geom_point(aes(colour = factor(Species)), size = 4) + 
  ggtitle("Real Iris Categories")
#dev.off()
p = ggplot(data, aes(Petal.Length, Petal.Width))

p + geom_point(aes(colour = factor(kk$cluster)), size = 4) + ggtitle("K-means Results")

?fviz_cluster

fviz_cluster(kk, data,ellipse.type = "norm")



###########################################################################
###################################K-means-TWO PREDICTORS##################
#########################################################################
library(plsgenomics)
data(leukemia)    
scaled_leukemia=scale(leukemia$X)
scaled_leukemia=as.data.frame(scaled_leukemia)
attach(scaled_leukemia)
leukemia_kk=kmeans(leukemia$X, 
                   centers = 5, iter.max = 50, nstart = 1)
leukemia_kk
leukemia_kk$cluster
leukemia_kk$iter
leukemia_kk$tot.withinss
fviz_cluster(leukemia_kk, leukemia$X)

 ###############Check the number of groups

p1 <- fviz_nbclust(scaled_leukemia, kmeans, method = "wss", 
                   k.max = 10) +
  ggtitle("(A) Elbow method")
p2 <- fviz_nbclust(scaled_leukemia, kmeans, method = "silhouette", 
                   k.max = 10) +
  ggtitle("(B) Silhouette method")

gridExtra::grid.arrange(p1, p2, nrow = 1)


###########################################################
################Clustering Leukemia dataset using 2 groups



leukemia_kk=kmeans(leukemia$X, centers = 2, iter.max = 50, nstart = 1)

leukemia_kk
leukemia_kk$cluster
leukemia_kk$iter
leukemia_kk$tot.withinss

fviz_cluster(leukemia_kk, leukemia$X)

##################################BES

library(cluster)
library(fclus) 
library(factoextra)
library(FactoMineR)

library(sjPlot)
library(dplyr)
library(ggplot2)
library(stats)
library(e1071)


library(readxl)
BES_2015 <- read_excel("BES_2015.xlsx")
View(BES_2015)

data=data.frame(BES_2015)   
attach(data)

dim(data)
summary(data)
data
rownames(data)=data[,1]
data=data[,-1]
data2015=dplyr::select(data,  contains("2015"))
data2015
sjp.corr(data2015,corr.method="pearson")

# Note that r -> 1 gives increasingly crisper clusterings  whereas 
#r -> Inf leads to complete fuzzyness. 

?fanny

fannyx=fanny(data2015,2)
fannyx

plot(fannyx)
summary(fannyx)


fannyx$silinfo
fannyx$membership
fannyx$coeff
fannyx$clustering
fannyx$diss
fannyx$data
fannyx$convergence
fannyx$call
fannyx$k.crisp


##############################################################################
################«Crisp» Cluster Plot of the Italian Regions###################
##############################################################################

library(e1071)
library(factoextra)
library(cluster)
?fviz_cluster
?fanny
head(data2015)

fanni=fanny(data2015, 2, diss = FALSE, memb.exp = 2, metric = "euclidean", 
            stand = FALSE, maxit = 500)
plot(fanni)

fviz_cluster(fanni,repel=TRUE)    # repel avoid overlapping
summary(fanni)
fanni$membership

#################################################################
##############«Fuzzy» Cluster Plot of the Italian Regions########
#################################################################

library(stats)
library(FactoMineR)
?PCA
pc=PCA(data2015,ncp=2,scale.unit = TRUE)
pc$eig
pc$var







coo=pc$ind$coord
dataf<-data.frame(coo[,1],coo[,2])
x1=coo[,1]
x2=coo[,2]
colnames(dataf)<-c("x1","x2")
Name=rownames(dataf)

library(ggplot2)
library(ggrepel)

?dist
distanze=dist(data2015,method = "euclidean")  # can be avoided if we use 							the data

library(ggrepel)

fuzzyplot<-ggplot(data2015, aes(x=x1, y=x2, 
      colour=fanni$membership[,1]))+
  geom_point(size=1)  +
  #geom_text(aes(label=Name),hjust=0, vjust=0, check_overlap = TRUE)+
  scale_colour_gradient("Membership degree for group 1",low="green",high="red") +
  ggtitle("Fuzzy Clustering") +
  geom_text(aes(label = rownames(data2015)), size = 4)
  

fuzzyplot




data(diana)
