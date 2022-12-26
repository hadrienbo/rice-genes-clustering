library(WGCNA)
minModuleSize = 50
hclust.method='average'
merge.hclust.method='average'
merge.threshold = 0.2
pamStage=T
corOptions = "use = 'p'"
cor.analysis=F
RsquaredCut=0.85
deepSplit = 3


data2wgcna <- read.csv('data2wgcna.csv',row.names = 1)
## part of data2wgcna (z-score)
#                     AT1G01010  AT1G01060   AT1G01120
# T1_biorep1_techrep1 0.9410511 -0.8421927 -0.52460297
# T1_biorep1_techrep2 1.5308907 -0.8424077 -0.50793792
# T1_biorep1_techrep3 1.8280669 -0.8429868 -0.50471229
# T1_biorep2_techrep1 1.7955671 -0.6875331 -0.74422657
# T1_biorep2_techrep2 1.2606236 -0.6891777 -0.70320075
# T1_biorep2_techrep3 1.6421986 -0.6876681 -0.65506606
# T1_biorep3_techrep1 1.5436608 -0.6772809 -0.62068630
# T1_biorep3_techrep2 1.9468675 -0.6768155 -0.62072068
# T1_biorep3_techrep3 2.0969749 -0.6761502 -0.64546329
# T2_biorep1_techrep1 0.8201197 -0.8418698 -0.08430271
###set powers

powers = c(c(1:10), seq(from = 11, to=50, by=3))
sft = pickSoftThreshold(data2wgcna, powerVector = powers, verbose = 5,
                        networkType = 'signed',RsquaredCut=0.85)

##soft threshold plot
cex1<-sft$powerEstimate


# ###make a plot
par(mfrow = c(1,2));
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels= sft$fitIndices$Power,cex=0.7,col="red");
# this line corresponds to using an R^2 cut-off of h
# abline(h=sft$fitIndices$SFT.R.sq[sft$fitIndices$Power==sft$powerEstimate],col="red")
h <- sft$fitIndices[sft$fitIndices$Power==sft$powerEstimate,2]
abline(h=h,col="red")
mtext(text = round(h,2),at=h, side=2,col="red",cex=0.7)
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels= sft$fitIndices$Power, cex=0.7,col="red")
h <- sft$fitIndices[sft$fitIndices$Power==sft$powerEstimate,5]
abline(h=h,col="red")
mtext(text = round(h,2),at=h, side=2,col="red",cex=0.7)


##########################################################################################################
#Calculate the TOM
softPower<-sft$powerEstimate

##adjacency matrix
adjacency = adjacency(data2wgcna, power = softPower,type='signed',corOptions="use = 'p'")

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency,TOMType = 'signed');
rownames(TOM) <- colnames(TOM) <- rownames(adjacency)

labels <- rownames(adjacency)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = hclust.method);

##########################################################################################################
#####Clustering results

###dynamic cluster
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,pamStage=pamStage,
                            deepSplit = deepSplit, pamRespectsDendro = T,
                            minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
cluster.results <- data.frame(labels=labels,color=dynamicColors,cluster=dynamicMods,row.names = labels)
head(cluster.results) ## the cluster column is the cluster number of different genes
#               labels     color cluster
# AT1G01010 AT1G01010       tan      12
# AT1G01060 AT1G01060       red       6
# AT1G01120 AT1G01120    salmon      13
# AT1G01140 AT1G01140       red       6
# AT1G01170 AT1G01170 lightcyan      16
# AT1G01180 AT1G01180    yellow       4

##########################################################################################################
###---> merge clusters
### merge to reduce the cluster numbers, if not necessary, skip the steps below

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


##eigen before merge
eigen.beforemerge = moduleEigengenes(data2wgcna, colors = dynamicColors)


MEs = eigen.beforemerge$eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree.beforemerge = hclust(as.dist(MEDiss), method = hclust.method);
# Plot the result

plot(METree.beforemerge, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=merge.threshold, col = "red")


##do the merge
merge = mergeCloseModules(data2wgcna, dynamicColors, cutHeight = merge.threshold, verbose = 3)


# The merged module colors
mergedColors = merge$colors;

cat(paste0(' Module number before merge: ',length(table(dynamicColors))))
cat(paste0(' Module number after merge: ',length(table(mergedColors))))

##make new plots

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

###----add new cluster results
cluster.results$new_color <- mergedColors
new_cluster.idx <- by(cluster.results,INDICES = cluster.results$new_color,function(x){
  data.frame(x,new_cluster=min(x$cluster))
})
cluster.results <- do.call(rbind,new_cluster.idx)

##re rank the clustrs
new_cluster_rerank<- dplyr::dense_rank(cluster.results$new_cluster)
if(any(cluster.results$new_cluster==0))
  new_cluster_rerank <- new_cluster_rerank-1

cluster.results$new_cluster_rerank <- new_cluster_rerank
rownames(cluster.results) <- NULL

head(cluster.results)

# labels color cluster new_color new_cluster new_cluster_rerank
# AT1G01940  blue       2     black           2                  2
# AT1G02370  blue       2     black           2                  2
# AT1G02630  blue       2     black           2                  2
# AT1G02670 black       7     black           2                  2
# AT1G02870  blue       2     black           2                  2
# AT1G03110  blue       2     black           2                  2
