# ========================================================
# Example script for transferring subtype annotations from 
# a motor cortex reference ("source") dataset for non-
# neuronal cell classes. Harmony was performed to align 
# source and target datasets and account for batch effects
# prior to performing the label transfer algorithm. Final 
# QC (i.e. doublet removal) was performed based on final 
# embeddings and subtype annotation. 
# ======================================================== 

source("~/00_scOnline_sourceCode/sconline_code.R")
source("~/00_scOnline_sourceCode/Human_mouse_alignment_code.R")
source("~/00_scOnline_sourceCode/integrated_analysis_functions.R")
source("~/02_labelTransfer/labelTransfer_plots.R")
library(Matrix)

# Load "source" dataset (in Single Cell Experiment format) with annotations from published transcriptomic atlas of mouse motor cortex
astro_SCE_source <- ... # load reference data with annotations to transfer
data_source=astro_SCE_source

# Load "target" dataset: subset of our study dataset for cell type of interest (in this example, astrocytes)
astro_Seurat_target <- ... # load target data with original embeddings 
data_target=astro_Seurat_target

#remove anno_cluster_res columns so that new clustering is performed
data_target@meta.data=data_target@meta.data[,!grepl("anno_cluster_res.",colnames(data_target@meta.data))]
#convert to SCE
data_target=.myExpSetCreatorFn(inputExpData=data_target@assays$RNA@counts,organism="mouse",minExpCells=0,
                               inputPdata=as.data.frame(data_target@meta.data),inputFdata=NULL,addExtraAnno=F,
                               server=T,redownload_files=F,ncores=10)

# Make sure datasets are annotated for any needed harmony batch correction
data_source$ds_batch="yao"
data_source$donor_id = data_source$anno_batch
data_target$ds_batch="thyroid"
# data_target$anno_batch=data_target$anno_batch #already defined in the dataset

# First, we remove genes with no ensembl_gene_id (have to edit this to not lose custom annotated transgenes)
sum(is.na(rowData(data_source)$ensembl_gene_id))
sum(is.na(rowData(data_target)$ensembl_gene_id))
rowData(data_target[which(rownames(data_target) == "iCre")])$ensembl_gene_id <- "iCre"
rowData(data_target[which(rownames(data_target) == "DNThrb")])$ensembl_gene_id <- "DNThrb"
data_target=data_target[!is.na(rowData(data_target)$ensembl_gene_id),]

# Second, we remove the genes w duplicate ensemble_gene_ids
print(paste("Removing duplicated gene names:",sum(duplicated(rowData(data_source)$ensembl_gene_id)),"genes from the source and",sum(duplicated(rowData(data_target)$ensembl_gene_id)),"from the target dataset."))
data_source=data_source[!duplicated(rowData(data_source)$ensembl_gene_id),]
data_target=data_target[!duplicated(rowData(data_target)$ensembl_gene_id),]

# Third, we convert the gene names to ensembl gene ids
row.names(data_source)=rowData(data_source)$ensembl_gene_id
row.names(data_target)=rowData(data_target)$ensembl_gene_id

print(paste("Number of genes in common between the two datasets",length(intersect(row.names(data_source),row.names(data_target)))))

# Lastly, we only include genes that are in common between the two datasets (okay to exclude transgenes for purposes of cell subtype annotation)
data_source=data_source[which(row.names(data_source) %in% row.names(data_target)),]
data_target=data_target[which(row.names(data_target) %in% row.names(data_source)),]

# Split each dataset based on donor_id:
data_source_list=.mySplitObject(data_source,"anno_batch")
names(data_source_list)
data_target_list=.mySplitObject(data_target,"anno_batch")
names(data_target_list)

# Join source and target datasets as input to integrative analysis
input_data=list(data=c(data_source_list,data_target_list),data_m=NULL)
.ArgList=.myArgCreatorFn(HVG_list=3,prefix="astroAlignedLT_harmony",nPCs=100,indScaling=F,includeHuman=F,includeMouse=T,saveDir="~/test",ncores=10)
data=.myHighVarGeneSlFn(input_data,dataorganism="Mouse",argList = .ArgList)
hvg=.sconline.select_HVG_thr(data=data,nGene=c(1500,2000))
.ArgList$HVG_list=hvg[,1]

# First perform PCA, which for non-neuronal subtypes was insufficient for aligning source and target datasets
.myPCAfn(data, argList = argList, UMI_cor_thr = argList$UMI_cor_thr)
.myUMAPgeneratorFn(argList=argList,pca_UMAP=T,umap.method='uwot') # calculate UMAP representations from PCA results
.myQCfn(argList=argList,pca_UMAP=T,colname_list=c("anno_batch")) # generate some basic QC plots
.myClusteringFn(argList=argList,pca_level = T,group.singletons=F) #Seurat clustering of the PCA space

# Harmony was utilized for non-neuronal cell classes to account for batch effects and obtain improved cell subtype annotation 
res=parallel::mclapply(argList$HVG_list,.myHarmonyFn_2var,batch_column=c("anno_batch","ds_batch"),argList=argList,mc.cores = 7)
.myUMAPgeneratorFn(argList=argList,umap.method='uwot',pca_UMAP=F) # non
.myQCfn(argList=argList,pca_UMAP=F,colname_list=c("anno_batch")) #QC analysis of the harmony results
.myClusteringFn(argList=argList,pca_level = F,group.singletons=F,cluster_resolution=c(0.6,1)) #Seurat clustering of the Harmony-corrected PC space

# Generate Seurat object from optimal solution
seuratObj = .mySeuratMakerFn(inputExpData = input_data, argList = .ArgList, select_pca = F, select_HVG = .ArgList$HVG_list[1],select_harmony_index=2)

## Once source and target datasets are properly aligned, we can run label transfer algorithm

# Split Seurat into source and target to extract embeddings for each
source_AlignSeurat=subset(seuratSubset,subset=ds_batch=="yao")
# source_rowNames=as.list(unique(source_AlignSeurat$anno_batch))
target_AlignSeurat=subset(seuratSubset,subset=ds_batch=="thyroid")
my_pca_source=source_AlignSeurat@reductions$pca@cell.embeddings
my_pca_target=target_AlignSeurat@reductions$pca@cell.embeddings

# Perform label transfer
res_astro_aligned_subset=.myLabelTransfer_aligned(pca_source=my_pca_source,
                                                  meta_source=my_meta_source,
                                                  pca_target=my_pca_target,
                                                  meta_target=my_meta_target,
                                                  source_label_col="cluster_label",
                                                  target_label_col="anno_cluster_res.1",
                                                  source_data=my_data_source,target_data=my_data_target,
                                                  nPCs=NULL) # default is 30 nPCs

# Generate plots to help evaluate quality of cell subtype annotation transfer from reference dataset
astro_genes=c("ENSMUSG00000024411","ENSMUSG00000020932","ENSMUSG00000030495","ENSMUSG00000020914")
astro_gene_names=c("Aqp4","Gfap","Slc7a10","Top2a")
plots=labelTransfer_plots(res=res_astro_aligned_subset,
                          clusteringRes='anno_cluster_res.1',
                          confidenceLevel=0.8,
                          markerGenes=c(astro_genes),
                          markerGeneNames=c(astro_gene_names))
cowplot::plot_grid(plotlist=plots)

# Perform final QC fine-tuning by removing doublet clusters (or, in the case of immune cells, reactive microglia)
seuratObj=res_astro_aligned_subset$seurat_obj
Idents(seuratObj)=seuratObj$anno_cluster_res.0.6
data_astro=subset(seuratObj,subset=anno_cluster_res.0.6=="6",invert=TRUE)
