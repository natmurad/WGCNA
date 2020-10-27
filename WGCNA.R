
############################ Preparing the Data

#Display the current directory
getwd();
#If necessary, change the pathway.
workingDir = ".";
setwd(workingDir);
setwd("YOUR DIRECTORY HERE")
dir()

#Loading WGCNA
library(WGCNA);

options(stringsAsFactors = FALSE);
#Enable multithread
enableWGCNAThreads();

#Read the data
#Rows are the sample and columns the genes  
expressiondata = read.csv("expressiondataEx.csv", sep = ";");

#remove gene name column
expression = as.data.frame(expressiondata[, -c(1)]) 
expression = t(expression)

#Column 1 -  gene names
colnames(expression) = expressiondata$EST
rownames(expression) = names(expressiondata)[-c(1)];

####Filtering missing values
gsg = goodSamplesGenes(expression, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optional: print the gene names and samples removed
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(expression)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(expression)[!gsg$goodSamples], collapse = ", ")));
  # Removing samples and genes unwanted
  expression = expression[gsg$goodSamples, gsg$goodGenes]
}

#dim(expressiondata)
#dim(expression)

#selecting data for test
#expression <- expression[,1:5000]

#Group data in a dendogram to check outliers
sampleTree = hclust(dist(expression), method = "average");

#Plot a sample tree: Open the output in a 12:9 inchs size window
dev.off()
sizeGrWindow(9,6)

pdf(file = "YOUR DIRECTORY/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#Plot a line showing the cut-off
abline(h = 47000, col = "red");

# Determine clusters below the line
#help("cutreeStatic")
clust = cutreeStatic(sampleTree, cutHeight = 47000, minSize = 10)
table(clust)
#clust

#Cluster 1 contains the samples we want to keep.
keepSamples = (clust==1)
expression0 = expression[keepSamples, ]
dim(expression0)
#dim(expression0)
nGenes = ncol(expression0)
nSamples = nrow(expression0)

#Read phenotypic data
traitData = read.csv("phenoTraits.csv", sep = ";");
dim(traitData)
names(traitData)
Samples = rownames(expression0);

traitRows = match(Samples, traitData$EST);
datTraits = traitData[traitRows, -1];
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage();

#Regrouping samples
sampleTree2 = hclust(dist(expression0), method = "average")
#Converting phenotypic characters in a color representation: white means low value, red means high value
#and gray missing value
traitColors = numbers2colors(datTraits, signed = FALSE);

#Plot a sample dendogram with the colors below
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
save(expression0, datTraits, file = "SpodopteradataInput.RData")


#######BUILDING THE NETWORK AND DETECTING MODULES

#expression0 = expression0[, 1:10000]
#lnames = load(file = "SpodopteradataInput.RData");
#The variable lnames contains the names of loaded variables.
#lnames

#Defining normalization power factor through the graphs

#Chosing a threshold set for the power
powers = c(c(1:10), seq(from = 12, to=20, by=2))

#Calling the function of network topology analysis
sft = pickSoftThreshold(expression0, powerVector = powers, verbose = 5)

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

#Index the scale free topology adjust as a function of the power soft thresholding.
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

#This line corresponds to use a cut-off R² of h
abline(h=0.90,col="red")

#Connectivity mean as a function of soft power thresholding
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#This line corresponds to use a cut-off R² of h
abline(h=0.90,col="red")

#Constructing the network and detecting modules
############ ESSA PARTE ESTOURA A MEMORIA

softPower = 7; #Chosen in the graphs before
#dim(expression0)
#expression0 = expression0[, 1:10000]
adjacency = adjacency(expression0, power = softPower); #Calculating the adjacency matrix
#adjacency[1:20, 1:20]

#Transforming the adjacency matrix in a topological overlap
TOM = TOMsimilarity(adjacency); #Calculating the topological overlap matrix
dissTOM = 1-TOM ##Calculating the dissimilarity

#Calling the hierarchical grouping function
geneTree = hclust(as.dist(dissTOM), method = "average");

#Plotting the clustering tree resulting (Dendogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

#Big modules are wanted then a high minimum module size must be chosen
minModuleSize = 100;

#Identification of module using dynamic tree cut-off
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

#Converting numeric names in colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

#Plot a dendogram and colors below
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#Clustering the modules

#Calculating eigengenes
MEList = moduleEigengenes(expression0, colors = dynamicColors)
MEs = MEList$eigengenes

#Calculating the module dissimilarity eigengenes
MEDiss = 1-cor(MEs);

#Clustering the eigengenes modules
METree = hclust(as.dist(MEDiss), method = "average");
#Plotting the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#Grouping the clusters from a cut-off
MEDissThres = 0.25
#Plotting a cut-off line
abline(h=MEDissThres, col = "red")

#Calling an automatic function of grouping (merge)
merge = mergeCloseModules(expression0, dynamicColors, cutHeight = MEDissThres, verbose = 3)

#Grouping module colors
mergedColors = merge$colors;

#Eigengenes of new grouped modules
mergedMEs = merge$newMEs;
#getwd()
sizeGrWindow(12, 9)
pdf(file = "YOUR DIRECTORY/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


#Renaming the module colors
moduleColors = mergedColors

#Building numeric labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

#Saving the module colors and labels to ensuing use
save(MEs, moduleLabels, moduleColors, geneTree, file = "Spodoptera-02-networkConstruction-stepByStep.RData")


#Dealing with a big data set: Constructing a block-wise network and detecting modules
bwnet = blockwiseModules(expression0, maxBlockSize = 5000,
                         power = 10, TOMType = "unsigned", minModuleSize = 100,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "SpodopteraTOM-blockwise",
                         verbose = 3)



#Loding the results pf single-block analysis
#load(file = "Spodoptera-02-networkConstruction-auto.RData");
# Re-labeling blockwise modules
bwLabels = matchLabels(bwnet$colors, moduleLabels);

#Converting the labels in colors to plot
bwModuleColors = labels2colors(bwLabels)


#Opening graphic window
sizeGrWindow(6,6)

#Plot the dendogram and the modules colors below to the block 1
#(Replacing the number of the block you want to plot in bwnet$dendrograms[[X]]
#e bwModuleColors[bwnet$blockGenes[[X]]])
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#Plotting more than one together
sizeGrWindow(12,9)
plotDendroAndColors(geneTree,
                    cbind(moduleColors, bwModuleColors),
                    c("Single block", "2 blocks"),
                    main = "Single block gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

singleBlockMEs = moduleEigengenes(expression0, moduleColors)$eigengenes;
blockwiseMEs = moduleEigengenes(expression0, bwModuleColors)$eigengenes;

single2blockwise = match(names(singleBlockMEs), names(blockwiseMEs))
signif(diag(cor(blockwiseMEs[, single2blockwise], singleBlockMEs)), 3)
dev.off();

#Relating modules to characteristics and identifying important genes
#Defining the number of genes and samples
nGenes = ncol(expression0);
nSamples = nrow(expression0);

#Recalculating MEs with label colors
MEs0 = moduleEigengenes(expression0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(8,4)

#Displaying correlations and its p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

#Displaying the correlation values in a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#Defining the variable Peso10dias containing the column Peso10dias of datTrait
Peso10dias = as.data.frame(datTraits$Peso10dias);
names(Peso10dias) = "Peso10d"

#names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(expression0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(expression0, Peso10dias, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(Peso10dias), sep="");
names(GSPvalue) = paste("p.GS.", names(Peso10dias), sep="");

module = "greenyellow" #########################putting the color below the plot
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Peso 10 dias",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#colnames(expression0)
#Display the gene names inside the module
colnames(expression0)[moduleColors=="greenyellow"] 

#Identifying most important genes for one determined characteristic inside of the cluster
geneInfo0 = data.frame(EST = colnames(expression0),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)


modOrder = order(-abs(cor(MEs, Peso10dias, use = "p")));
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Peso10d));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")

#Exporting the network to a cytoscape format
#Recalculating topological overlap, if necessary
#TOM = TOMsimilarityFromExpr(expression0, power = 10);
#Select the modules
#modules = c("brown", "red"); #chose modules that u want to export
#Select the gene modules
genes = colnames(expression0)

#if you want export specific colors, substitute the second modulecolors by above modules
inModule = is.finite(match(moduleColors, moduleColors));
modGenes = genes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];

#Select the corresponding topologic overlap 
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modGenes, modGenes)
#########modTOMSignificantes = which(modTOM>0.1)
#####warnings()

#Organize the genes by importance inside the module
genes = colnames(expression0)
sum(is.na(genes))
#It must return 0.


#Create the dataframe since the beginning
geneInfo0 = data.frame(ESTs = genes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

#Order the modules by the significance by a character Ex: peso10days
modOrder = order(-abs(cor(MEs, Peso10dias, use = "p")));

#Add information of the members of the module in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

#Order the genes of geneinfo variable first by the color of the module, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Peso10d));
geneInfo = geneInfo0[geneOrder, ]

#write the file with the ordered values
write.csv(geneInfo, file = "geneInfo.csv")

#Export the network in list files os n edges that cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = "CytoscapeEdgeFile.txt",
                               nodeFile = "CytoscapeNodeFile.txt",
                               weighted = TRUE,
                               threshold = 0.4,
                               nodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])

