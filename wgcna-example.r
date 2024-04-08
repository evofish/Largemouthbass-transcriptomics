#this is the WGCNA for Large mouth bass North 

#installing WGCNA again
#install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
#install.packages("BiocManager")
#BiocManager::install("WGCNA")
#BiocManager::install("DESeq2")

setwd('~/Desktop/large-mouth-bass/wgcna/')
library('DESeq2')
library('pheatmap')
library('RColorBrewer')
library('ggplot2')
library('WGCNA')
library('tidyverse')
library('flashClust')

options(stringsAsFactors = FALSE)

x<-read.csv('North-allcounts.txt', row.names="Geneid", header=T, sep='\t')
head(x)
nrow(x)#114287
x$filtr<-apply(x, 1, function(k) mean(k > 10)) > 0.5
x<-x[!(x$filtr=="FALSE"),]
nrow(x) #17048

#This is the code to only analyze the 50% with the highest variation
#counts$variance = apply(counts, 1, var)
#counts2 = counts[counts$variance >= quantile(counts$variance, c(.50)), ] #50% most variable genes
#counts2$variance <- NULL
#counts=counts2

x$filtr=NULL
#outlier removal if necessary

metaData <- read.csv('northern-traits.csv', header = TRUE, sep = ",", row.names = "sample") 
head(metaData)
nrow(metaData)

totalCounts=colSums(x) 
min(totalCounts) #361818
max(totalCounts) #3263825
mean(totalCounts) #1666593

#countdata<-apply(counts,1,function(x) {all(x>0)})
#conditions for LRT

dds <- DESeqDataSetFromMatrix(countData=x,
                              colData=metaData,
                              design=~temperature+age+tank )
dds<-dds[rowSums(counts(dds))>0,]
dds <- DESeq(dds, test="LRT", reduced=~tank)  
summary(dds) #12946
res <-results(dds)
table(res$padj<0.01)
#FALSE  TRUE 
# 9912  7136 #genes differentiated by age and by temp

#to obtain variance stabilized data folow: 
head(res)
vsd=getVarianceStabilizedData(dds)
head(vsd)
write.csv(vsd, file="North-WGCNA-variancestabilized-Jan23.csv", quote=FALSE)

# Check that the data has the correct format for many functions operating on multiple sets:
#WGCNA THE DATA NEEDS TO BE TRANSPOSED to make dendrogram

gsg = goodSamplesGenes(vsd, verbose = 5);
gsg$allOK

dim(vsd) #17048    48
xx=as.data.frame(t(vsd)) #The stabilized data has to be transposed! 
dim(xx)
rownames(xx) #make sure the rownames are actually 

sampleTree = hclust(dist(xx), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

traitData = read.csv("northern-traits-wgcna.csv",header=TRUE);
dim(traitData)
names(traitData)

flosamples = rownames(xx);
traitRows = match(flosamples, traitData$sample);
datTraits = traitData[traitRows, -1];
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage();

#store the data in an R function 

save(xx, datTraits, file = "NOR-wgcna-step1.RData")

##To Find blocks of modules 8, use the data that is not transposed 
# Choose a set of soft-thresholding powers
#This was not done... not sure how to set nSets. 
lnames= load(file="NOR-wgcna-step1.RData")
lnames
checkSets=TRUE
nSets = checkSets(multiExpr)$nSets

powers = c(c(1:10), seq(from = 15, to=40, by=2))

powerTables = vector(mode = "list", length=nSets);
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();

# Call the network topology analysis function
sft = pickSoftThreshold(xx, powerVector = powers, verbose = 10)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#the curve is initially stabilized at 24!!!! - with xx dataframe
net = blockwiseModules(xx, power = 24, corType="bicor", networkType = "signed",
                       TOMType = "signed", minModuleSize = 50,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE,
                       deepSplit = 2,
                       pamRespectsDendro = F,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "FLOR-test",
                       verbose = 5, maxBlockSize = 5000)
names(net)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
consMEs = net$MEs;
consTree = net$dendrograms[[1]]

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

save(consMEs, moduleLabels, moduleColors, consTree, file = "Nor-NetworkConstruction-auto.RData")
######## 
#Clusters to All Samples (not very clear pattern)
module_df <- data.frame(gene_id = names(net$colors),
                        colors = labels2colors(net$colors))
module_df
module_df[1:5,]
write.table(module_df, file ="FLO-gene_modules.txt",sep = "\t", quote=FALSE)

MEs0 <- moduleEigengenes(xx, mergedColors)$eigengenes
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
MEs0$treatment = row.names(MEs0)

mME = MEs0 %>% pivot_longer(-treatment) %>% mutate(
  name = gsub("ME", "", name),
  name = factor(name, levels = module_order))

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

########
#Cluster Dendrogram Dissimilarity and Merged Dynamic

softPower = 24;
adjacency = adjacency(xx, power = softPower, type='signed');
TOM = TOMsimilarity(adjacency, TOMType='signed');
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
#pdf(file=Dendrogram_signed_BM10.pdf, width=20, height=20)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
#dev.off()
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
#dynamicMods
#0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#1995 5079 2543 1854 1206  979  895  543  447  318  293  226  221  142  123  110   74

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#DynamicColors
#black         blue        brown         cyan        green  greenyellow 
#543         2543         1854          123          979          226 
#grey    lightcyan      magenta midnightblue         pink       purple 
#1995           74          318          110          447          293 
#red       salmon          tan    turquoise       yellow 
#895          142          221         5079         1206 
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
#pdf(file=Dendrogram_signed_BM10_colors.pdf, width=20, height=20)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(xx, colors = dynamicColors)
MEList$eigengenes #gives you the eigenes by sample 
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
#pdf(file=ClusteringEigengenes.pdf, width=20, height=20)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#dev.off()

MEDissThres = 0.2 #decided to choose this one instead of 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(xx, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#sizeGrWindow(12, 9)
#pdf(file = "DendroAndColors_sft6_bm10.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# MERGING: Rename to moduleColors that are similar
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "Florida-RNA_WGCNA_networkConstruct_signed-Mar2023.RData")
# Define numbers of genes and samples
nGenes = ncol(xx);
nSamples = nrow(xx);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(xx, moduleColors)$eigengenes 
MEs = orderMEs(MEs0)

######################
#correlations of traits with eigengenes

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitCor 

## Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
par(mar=c(0.8,0.8,0.8,0.8))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

######
##EXTRACTING DATA FROM MODULES FOR GO
annot=read.table(file="nor-GO-final.txt", header=T, sep='\t')
head(annot)
probes = names(xx)
probes2annot = match(probes,annot$V1)
probes2annot

datGS.Traits=data.frame(cor(xx,datTraits,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(xx,moduleColors)$eigengenes
datKME=signedKME(xx, datME, outputColumnName="MM.")

datOutput=data.frame(ProbeID=names(xx), annot[probes2annot,],moduleColors,datKME,datGS.Traits)
dim(datOutput)

head(datOutput)
write.table(datOutput, "Nor-GOAnnotated-Modules-2023.txt", sep='\t') 
#this table was edited in excel to make the final input for the analysis of enrichment using the Fisher's test

datOutput[complete.cases(datOutput),]
nrow(datOutput) #17048

ww=datOutput[c(2,4)]
ww2=ww[complete.cases(ww),] 
nrow(ww2) #11264
head(ww2)
write.table(ww2, file='nor-isogroup-to-module.txt', sep='\t', quote=F)
#this contains all the isogroups with the corresponding Module and their module. 

### for GO-MWU 
## for this we now need to see the correspondance between Modules and the complete table of GO-cats for Florida

head(ww2)
goannot<-read.table("nor-GO-final.txt", header = T, sep='\t')
head(goannot) 
go_iso_module<-merge(ww2, goannot, by='V1', all=T)
head(go_iso_module)
nrow(goannot)
nrow(go_iso_module)
write.table(go_iso_module, file='nor-all-GO-module.txt', sep='\t')
