# This script creates and plots networks of TCRs which share sequence similarity as measured 
# by the string kernel, looking for similar triplets.    
# required libraries are : 
library(kernlab)
library(igraph)
library(stringdist)

################# CDR3 CLUSTERING - NP16 CLONES ALPHA ###################

cd8_np16 <- read.csv(("/t1-data/user/lfelce/TCR_analysis/cd8_np16_tcr_clones.csv"))

# this script requires two sets of input; a set of  expanded CDR3 sequences (each CDR3 should only occur once) in the set of TCRs to be analysed, 
# in a variable called CDR_exp; and all the CDR3s (each CDR3 should only occur once) from the same set, in a variable called CDR_all.    
CDR_exp <- cd8_np16 %>% group_by(CDR3_aa.x) %>% filter(n() >= 2)
CDR_exp <- CDR_exp %>% distinct(CDR3_aa.x, .keep_all = TRUE)
CDR_exp <- CDR_exp[,"CDR3_aa.x"]
colnames(CDR_exp) <- "CDR3_alpha_exp"

CDR_all <- cd8_np16 %>% distinct(CDR3_aa.x, .keep_all = TRUE)
CDR_all <- as.data.frame(CDR_all[,"CDR3_aa.x"])
colnames(CDR_all) <- "CDR3_alpha_all"

#nc<-c()
#network_sum<-c()

# initialize the string kernel, which measures the pairwise distance between string vectors, and outputs a kernel matrix 
# containing the string distance between each CDR3 in set 1 and each CDR3 in set 2.  
sk <- stringdot(type="spectrum", length=3, normalized=TRUE)

# calculate string distance between expanded clones and rest
km_CDR <- kernelMatrix(sk,CDR_exp,CDR_all)

#remove all the CDR3s which are remotely connected to any of the inputs
colmax <- apply(km_CDR,2,max)
# optionally plot a histogram of distances
#hist(colmax)

# chose top 1500 TCRS and discard the rest. This is an arbitrary cutoff, which retains all the similar TCRs
# but does not produce too large a distance matrix 
l_colmax<-min(1500,length(colmax))
rem<-order(colmax,decreasing = TRUE)[1:l_colmax]

# now work out full kernel  matrix for the remaining 1500 CDRs
km_CDR_all<-kernelMatrix(sk,CDR_all[rem],CDR_all[rem])

#replace those > 0.98 with 0; this removes any identical pairs of CDR3s. 
ones <- which(km_CDR_all>0.98)
km_CDR_all[ones]<-0

# optionally plot a histogram of distances again to check 
# hist(km_CDR_all)

# convert the distance matrix into an adjacentcy  matrix; any two TCRs which are closer than 
# an arbitrary threshold th, are given a distance value of 1.  Any which are less than th, are given a distance of 0.
# the size of the threshold is arbitrary, and needs to be set so a suitable control set of CDR3s  give few clusters.

thc<-0.85
connect<-(km_CDR_all>=thc)
km_CDR_adj<-km_CDR_all
km_CDR_adj[which(connect)]<-1
km_CDR_adj[which(!connect)]<-0

# nc is the number of non-zero entries in the set; i.e. number of edges     
nc <- length(which(km_CDR_adj>0))/2

# hist(km_CDR_adj)
# convert the adjacency  matrix into a graph using igraph
km_graph<-graph.adjacency(km_CDR_adj, mode="undirected",weighted = TRUE,add.colnames=NULL,diag=FALSE)
summary(km_graph)

# names of vertices are names of CDR3s
V(km_graph)$name<-CDR_all[rem]

# remove isolated nodes; new smaller graph is called km_graph_s
isolates <- function(g){return(which(igraph::degree(g)==0))}
km_graph_s<-delete.vertices(km_graph,isolates(km_graph))
summary(km_graph_s)

#define which vertices are also exp clones
inter <- match(CDR_exp,V(km_graph_s)$name)
exp_v <- inter[which(!is.na(inter))]

names_g<-V(km_graph_s)$name
# collect summary data from the graph
# number of vertices
n_nodes1 <- length(V(km_graph_s))

# size of each cluster
cluster_size <- clusters(km_graph_s)$csize

# number of clusters
noclusters <- clusters(km_graph_s)$no

# degree of each vertix 
degree <- igraph::degree(km_graph_s,V(km_graph_s))
degree_dist <- degree.distribution(km_graph_s)

# which node(vertex)  is in which cluster 
membership <- clusters(km_graph_s)$membership

# capture netowrk layout so can plot reproducibly

coords.exp <- layout.kamada.kawai(km_graph_s)

# Create the plot
plot(km_graph_s, vertex.label=NA,layout=coords.exp,edge.width=1,vertex.color="black",edge.color="grey",vertex.size=4)
title(main="expanded CDR3s", cex.main=1, font.main=1)

#can add text to plot using mtext
#mtext(paste0(nam," (thr = ",thc," c.thr = ",c.thr,")"), outer = TRUE, cex = 1)
