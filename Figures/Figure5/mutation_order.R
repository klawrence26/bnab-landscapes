library(dplyr)
library(stringr)
library(ggplot2)
library(ggpmisc)
library(reshape2)
library(nlmrt)
library(dplyr)
library(ggpmisc)
library(matrixStats)
library(rethinking)
library(rstan)
library(robustbase)
library(RColorBrewer)

########6261#######
matrix_pathways <- read.csv('data/6261_probability_random_median_moderate.csv',
                            quote="",header=TRUE)
matrix_pathways.melted <- melt(matrix_pathways,id=c("Mutation"))
matrix_pathways.melted$Mutation <- factor(matrix_pathways.melted$Mutation,
                                          levels = c("Mut 1", "Mut 2", "Mut 3", "Mut 4", "Mut 5", "Mut 6",
                                                     "Mut 7", "Mut 8", "Mut 9", "Mut 10", "Mut 11"))
levels(matrix_pathways.melted$Mutation) <- c("29","35","65","66","69","82","83","84","85","87","112.1")
levels(matrix_pathways.melted$variable) <- c('#1','#2','#3','#4','#5',"#6","#7","#8","#9","#10","#11")
mycolors <- c("#66c2bd",
              "#1b9d9b", 
              "#7570B3",
              "#E78AC3",
              "#E7298A",
              "#A6D854",
              "#66A61E",
              "#FFD92F",
              "#E6AB02",
              "#B3B3B3",
              "#666666")

ggplot(data = matrix_pathways.melted,aes(x=variable,y=value,fill=Mutation))+
  geom_bar(position="fill", stat="identity")+
  xlab("Mutation order")+
  ylab("Relative Likelihood")+
  scale_fill_manual(values=mycolors)+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

########9114#######
matrix_pathways <- read.csv('data/9114_probability_random_median_moderate.csv',quote="",header=TRUE)


matrix_pathways.melted <- melt(matrix_pathways,id=c("Mutation"))
matrix_pathways.melted$Mutation <- factor(matrix_pathways.melted$Mutation,
                                          levels = c("Mut 1","Mut 2","Mut 3","Mut 4","Mut 5","Mut 6",
                                                     "Mut 7","Mut 8","Mut 9","Mut 10","Mut 11",
                                                     "Mut 12","Mut 13","Mut 14","Mut 15","Mut 16"))

levels(matrix_pathways.melted$variable) <- c('1','2','3','4','5',"6","7","8",
                                             "9","10","11","12","13","14","15","16")
levels(matrix_pathways.melted$Mutation) <- c("30","35","36","57","64","65","66","79","82",
                                              "83","84","85","92","95","103","113")
library(RColorBrewer)

df<-data.frame()
for(x in 1:16){
  for(y in 1:16){
    newrow<-c(x,y,sample(1:1000,1))
    df<-rbind(df,newrow)
  }
}
colnames(df)<-c('X','Y','Val')


brewer.pal(n=8,"Set2")
brewer.pal(n=8,"Pastel2")
brewer.pal(n=8,"Dark2")
mycolors <- c("#66c2bd", 
              "#1b9d9b", 
              "#fc7864", 
              "#d94602", 
              "#8DA0CB",
              "#7570B3",
              "#E78AC3",
              "#E7298A",
              "#A6D854",
              "#66A61E",
              "#FFD92F",
              "#E6AB02",
              "#E5C494",
              "#A6761D",
              "#B3B3B3",
              "#666666")

ggplot(data = matrix_pathways.melted,aes(x=variable,y=value,fill=Mutation))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=mycolors)+
  xlab("Mutation order")+
  ylab("Relative Likelihood")+
  theme_classic()+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=7))

