
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RandomWalkRestartMH")

#install.packages("tidyverse")
#install.packages(RandomWalkRestartMH)
library(RandomWalkRestartMH)
#install.packages("tidyverse")
library(tidyverse)
library(igraph)

#### TCGA-LUNG-train data
sg_data <- read.csv('score_sg_all_tcga_brca_all_cosmic.csv', header=TRUE,na.strings="")
string_data <- read.csv('Gene_interaction_list1.csv', header=TRUE,na.strings="")
colnames(string_data)[3] <- 'weight'
delta = 0.3
string_graph = graph_from_data_frame(string_data, directed=F, vertices=NULL)
string_MultiplexObject = create.multiplex(list(STRING=string_graph))
AdjMatrix_string <- compute.adjacency.matrix(string_MultiplexObject, delta=delta)
AdjMatrixNorm_string <- normalize.multiplex.adjacency(AdjMatrix_string)
#seeds= c("RNF43","NOTCH3","POLR1F")
string_node_list = get.data.frame(string_graph, what='vertices')[,1]
seed_data <- read.csv('seed_tcga_brca_all_cosmic.csv', header=FALSE,na.strings="")
#probs <- list()
#probs <- vector(mode='list', length=7)
l = list()
probs= NULL
p=NULL
l1 =list()
l2 = list()
l3 = list()
l4 = list()

#unique_gene_1<-unique(string_data$Gene1)

#unique_gene_2<-unique(string_data$Gene2)

Col_gene_list<-read.csv("Cosmic.csv", header=TRUE,na.strings="")

len_genes_list<-length(Col_gene_list)

final_gene_list<-unique(Col_gene_list$Genes)

score_data = data.frame(matrix(nrow = 1009, ncol = 739))

# assign column names
colnames(score_data) = final_gene_list



for(i in 1:nrow(seed_data)) { 
  g = seed_data[i,]
  l[[i]] = g[which(g != "")]
  #seed_data[i, ] <- seed_data[i, ] 
  rwr_result <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_string, string_MultiplexObject,l[[i]])
  probs[i] <- rwr_result$RWRM_Results["Score"]
  p[i] <- rwr_result$RWRM_Results["NodeNames"]
  #colnames(probs)[i] <- paste0("Col_", i) 
  df<- as.data.frame(rbind(probs[[i]],p[[i]]))
  ### here data should be prepared
  l2[[i]]= df
  df1<- as.data.frame(t(df))
  l3[[i]]= df1
  colnames(l3[[i]]) <- c("Score","NodeName")
}

for(i in 1:1009)
  
{
  
  l4[[i]] =l3[[i]][order(l3[[i]]$Score,decreasing=TRUE),]
}

#top_200= list()
#for(i in 1:length(l4))
#{
  #top_200[[i]] <- l4[[i]][1:200,]
#}
merge_data = list()
sum_data =list()
len_of_list =list()

for(i in 1:length(l4))
{
      
  merge_data[[i]] <- merge(l4[i], sg_data, by.x = "NodeName", by.y = "NodeName")
  len_of_list[[i]] <- length(merge_data[[i]]$NodeName)
}

## accessing elements
#score_data$
#nrow(score_data)
#ncol(score_data)
#length(merge_data)
##length(merge_data[[1]])
#nrow(merge_data)
#nrow(merge_data[[1]])

list_of_head=list()
list_of_head<-colnames(score_data)
list_of_head


len_head<-length(list_of_head)

#list_of_head[[1]]

#for(x in 1:len_head){
  #print(x)
#}

#n<-"CBL"
# find index value of "Swift" using which()
#pos<-which(list_of_head ==n) # 2


pos<-0

for(i in 1:nrow(score_data)){
  a<-merge_data[[i]]
  #nrow(a)
  for(j in 1:nrow(a)){
    n<-a[j,1]
    s<-a[j,2]
    ## here number of columns are getting increased , we need to check here
    if (n %in% list_of_head) {
      pos<-which(list_of_head ==n)
      print(pos)
      score_data[i,pos]<-s
      #print("Item is present in the List.")
    }
    else
    {
      print("not present")
      #score_data[i,pos]<-0
    }
    
  }
}

write.csv(score_data,"tcga_BRCA_matabric_cosmic_all.csv")
  

