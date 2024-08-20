# ========================================================
# Example script for how snRNA-seq data (10x Genomics, 
# filtered feature barcode matrix output) was processed 
# through initial QC, PCA, UMAP, and clustering. 
# ======================================================== 

source("~/00_scOnline_sourceCode/sconline_code.R")
source("~/00_scOnline_sourceCode/Human_mouse_alignment_code.R")
source("~/00_scOnline_sourceCode/integrated_analysis_functions.R")

libList <- gsub('.h5','',list.files('~/insert_directory_here/')) # filtered_feature_bc_matrix data is available for download from GEO (GSE271421)

dataList <- list()
seuratList <- list()

# Generate counts matrix and Seurat object for each library
for (i in 1:length(libList)) {
  dataList[[i]] <- Read10X_h5(paste0("~/insert_directory_here/",libList[[i]],".h5"))
  seuratList[[i]] <- CreateSeuratObject(counts = dataList[[i]], project = "thryoid")
}

for (i in 1:length(libList)) {
  seuratList[[i]]$libName <- libList[[i]]
}

#Add mitochondrial QC metrics to meta data
for (i in 1:length(libList)) {
  seuratList[[i]]$percent.mt <- PercentageFeatureSet(seuratList[[i]], pattern = "^mt-")
}

# Basic QC on both seurat object lists
seuratList <- lapply(seuratList, function(x) {
  x <- subset(x, subset = nFeature_RNA > 200 & percent.mt < 5) # 
  return(x)
})

# add orig.ident to distinguish the replicates from the same donor prior to merging
# add donor_id to facilitate merge
# add condition and treatment 
seuratList <- lapply(1:length(libList), function(i) {
  obj <- seuratList[[i]]
  obj@meta.data$orig.ident <- libList[[i]]
  obj@meta.data$donor_id <- ... # annotate with mouse ID (ie merge sequencing replicates) according to metadata
  obj@meta.data$THRB <- ... # annotate as "DN-THRB", "WT-THRB", or "N/A" according to metadata
  obj@meta.data$treatment <- ... # annotate as "T3" or "C" according to metadata
  return(obj)
})

# Use Seurat merge() function to merge Seurat objects for sequencing replicates (ie Seurat objects with same donor_id but different orig.ident)
completeMergedSeuratList <- qread("~/insert_directory_here/completeMergedSeuratListFinal.qs")

# Create list of Single Cell Experiment objects from the list of Seurat objects
sceList <- list() 
sceList <- lapply(completeMergedSeuratList,function(x) {
  x <- .myExpSetCreatorFn(inputExpData=x@assays$RNA@counts,organism="mouse",minExpCells=0,
                          inputPdata=as.data.frame(x@meta.data),inputFdata=NULL,addExtraAnno=T,
                          server=T,redownload_files=F,ncores=3)
  return(x)
})

data_thyroidAnalysis=sceList
# set ds_batch and anno_batch meta.data
data_thyroidAnalysis <- lapply(data_thyroidAnalysis, function(x) {
  x$ds_batch <- x$donor_id
  x$anno_batch <- x$donor_id
  return(x)
})

# annotate iCre and DNThrb - they will have NA in rowData(data_thyroidAnalysis[[i]])$ensembl_gene_id
customAnnot <- lapply(data_thyroidAnalysis,function(x){
  i1 = which(rownames(x) == "iCre")
  i2 = which(rownames(x) == "DNThrb")
  if(length(i1) == 0) i1 = "NA"
  if(length(i2) == 0) i2 = "NA"
  return(c(i1,i2))
})

#we first remove genes with no ensembl_gene_id
data_thyroidAnalysis <- lapply(1:length(data_thyroidAnalysis), function(i) {
  obj <- data_thyroidAnalysis[[i]]
  rowData(obj)$ensembl_gene_id[which(rownames(obj) == "iCre")] <- "iCre" # replace NA for iCre so it doesn't get removed
  rowData(obj)$ensembl_gene_id[which(rownames(obj) == "DNThrb")] <- "DNThrb" # replace NA for DNThrb so it doesn't get removed
  print(rownames(obj)[which(is.na(rowData(obj)$ensembl_gene_id))])
  obj <- obj[!is.na(rowData(obj)$ensembl_gene_id),] # remove remaining NA values
  return(obj)
})

# second, we remove genes with duplicate ensemble_gene_ids
data_thyroidAnalysis <- lapply(1:length(data_thyroidAnalysis), function (i) {
  print(paste("Removing duplicated gene names:",sum(duplicated(rowData(data_thyroidAnalysis[[i]])$ensembl_gene_id)),"genes"))
  obj <- data_thyroidAnalysis[[i]]
  obj <- obj[!duplicated(rowData(obj)$ensembl_gene_id),]
  row.names(obj) <- rowData(obj)$ensembl_gene_id
  return(obj)
})

input_data = list(data=data_thyroidAnalysis,data_m=NULL)

.ArgList=.myArgCreatorFn(HVG_list=3,prefix="TH_completeMerged",nPCs=25,indScaling=F,includeHuman=F,includeMouse=T,saveDir="~/test",ncores=6)
names(.ArgList)
.ArgList$saveDir

# Select highly variable genes
data=.myHighVarGeneSlFn(input_data,dataorganism="Mouse",argList = .ArgList)

# Set HVG threshold based on the number of desired variable genes
hvg <- .sconline.select_HVG_thr(data=data,nGene=c(2000,4000))
.ArgList$HVG_list=hvg[,1]

# Example of running PCA, UMAP, and Harmony (change parameters as specified in Methods) 
library(irlba)
.myPCAfn(data, argList = .ArgList, UMI_cor_thr = .ArgList$UMI_cor_thr)
dir(.ArgList$saveDir) # PCA results are saved in saveDir directory from .ArgList:
.myUMAPgeneratorFn(argList=.ArgList,pca_UMAP=T,umap.method='uwot') # calculate UMAP representations from PCA results 
.myQCfn(argList=.ArgList,pca_UMAP=T,colname_list=c("anno_batch")) # generate some basic QC plots
.myClusteringFn(argList=.ArgList,pca_level = T,group.singletons=F) # Seurat clustering of the PCA space

# Run harmony to account for batch effects (ie from "anno_batch")
res=parallel::mclapply(.ArgList$HVG_list,.myHarmonyFn_1var,batch_column="anno_batch",argList=.ArgList,mc.cores = 7)

# Generate UMAP, clustering and QC plots for harmony results
.myUMAPgeneratorFn(argList=.ArgList,umap.method='uwot',pca_UMAP=F)
.myQCfn(argList=.ArgList,pca_UMAP=F,colname_list=c("anno_batch"))
.myClusteringFn(argList=.ArgList,pca_level = F,group.singletons=F)

# Convert PCA or harmony solution back to Seurat to perform marker analysis
pcaSlnToSeurat = .mySeuratMakerFn(inputExpData = input_data, argList = .ArgList, select_pca = T, select_HVG = .ArgList$HVG_list) 
harmonySlnToSeurat = .mySeuratMakerFn(inputExpData = input_data, argList = .ArgList, select_pca = F, select_HVG = .ArgList$HVG_list, select_harmony_index = ) #incorporates UMAP and pca coordinates based on the chosen solution

# Remove clusters with high MT% and OXPHOS%
Idents(pcaSlnToSeurat)=pcaSlnToSeurat$anno_cluster_res.0.6
artifactIDs = as.character(c(0,23,30,32))
filtered_pcaSlnToSeurat <- subset(pcaSlnToSeurat,idents = artifactIDs,invert = TRUE)

# Reran PCA/clustering steps above on cleaned up dataset with 100 nPCs (again using 2000 highly variable genes) to generate final cleaned dataset for major cell class annotation
