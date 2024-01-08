import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
def run_r_script(input_csv_sg, input_csv_seed, output_csv):
    r_script = f"""
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RandomWalkRestartMH")

library(RandomWalkRestartMH)
library(tidyverse)
library(igraph)

sg_data <- read.csv('{input_csv_sg}', header=TRUE, na.strings="")
string_data <- read.csv('Gene_interaction_list1.csv', header=TRUE, na.strings="")
# print(string_data)
colnames(string_data)[3] <- 'weight'
delta = 0.3
string_graph = graph_from_data_frame(string_data, directed=F, vertices=NULL)
string_MultiplexObject = create.multiplex(list(STRING=string_graph))
AdjMatrix_string <- compute.adjacency.matrix(string_MultiplexObject, delta=delta)
AdjMatrixNorm_string <- normalize.multiplex.adjacency(AdjMatrix_string)

seed_data <- read.csv('{input_csv_seed}', header=FALSE, na.strings="")
l = list()
probs= NULL
p=NULL
l1 =list()
l2 = list()
l3 = list()
l4 = list()

Col_gene_list <- read.csv("Cosmic.csv", header=TRUE, na.strings="")
len_genes_list <- length(Col_gene_list)
final_gene_list <- unique(Col_gene_list$Genes)

score_data = data.frame(matrix(nrow = 1009, ncol = 739))
colnames(score_data) = final_gene_list

for(i in 1:nrow(seed_data)) {{
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
}}

for(i in 1:1009)

{{

  l4[[i]] =l3[[i]][order(l3[[i]]$Score,decreasing=TRUE),]
}}

merge_data = list()
sum_data =list()
len_of_list =list()

for(i in 1:length(l4))
{{

  merge_data[[i]] <- merge(l4[i], sg_data, by.x = "NodeName", by.y = "NodeName")
  len_of_list[[i]] <- length(merge_data[[i]]$NodeName)
}}

list_of_head=list()
list_of_head<-colnames(score_data)
list_of_head
len_head<-length(list_of_head)


pos<-0

for(i in 1:nrow(score_data)){{
  a<-merge_data[[i]]
  #nrow(a)
  for(j in 1:nrow(a)){{
    n<-a[j,1]
    s<-a[j,2]
    ## here number of columns are getting increased , we need to check here
    if (n %in% list_of_head) {{
      pos<-which(list_of_head ==n)
      print(pos)
      score_data[i,pos]<-s
      #print("Item is present in the List.")
    }}
    else
    {{
      print("not present")
      #score_data[i,pos]<-0
    }}

  }}
}}

result_path <- '{output_csv}'
write.csv(score_data, result_path)
result_path
"""

    return str(robjects.r(r_script))

input_csv_sg = "score_sg_all_tcga_brca_all.csv"
input_csv_seed = "seed_tcga_brca_all.csv"
output_csv = "output_tcga_brca_cosmic_all.csv"
result_path = run_r_script(input_csv_sg, input_csv_seed, output_csv)
