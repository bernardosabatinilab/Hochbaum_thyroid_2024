# ========================================================
# Example script for transferring subtype annotations from 
# a motor cortex reference ("source") dataset for glutama-
# tergic and GABAergic neurons. Embeddings from original 
# clustering was sufficient for a high-quality annotation.
# ======================================================== 

source("~/00_scOnline_sourceCode/sconline_code.R")
source("~/00_scOnline_sourceCode/integrated_analysis_functions.R")
source("~/02_labelTransfer/labelTransfer_plots.R")

# Load "source" dataset (in Single Cell Experiment format) with annotations from published transcriptomic atlas of mouse motor cortex
exN_SCE_source <- ... # load reference data with annotations to transfer

# Load "target" dataset: subset of our study dataset for cell type of interest (in this example, glutamatergic neurons)
exN_Seurat_target <- ... # load target data with original embeddings 
# Convert to Single Cell Experiment format 
exN_SCE_target=.myExpSetCreatorFn(inputExpData=exN_Seurat_target@assays$RNA@counts,organism="mouse",minExpCells=0,
                                  inputPdata=as.data.frame(exN_Seurat_target@meta.data),inputFdata=NULL,addExtraAnno=T,
                                  server=T,redownload_files=F,ncores=10)

data_source <- exN_SCE_source
data_target <- exN_SCE_target

# If needed, remove genes from source and target datasets with no ensembl_gene_id (have to edit this to not lose custom annotated transgenes)
data_source=data_source[!is.na(rowData(data_source)$ensembl_gene_id),]
rowData(data_target[which(rownames(data_target) == "iCre")])$ensembl_gene_id <- "iCre"
rowData(data_target[which(rownames(data_target) == "DNThrb")])$ensembl_gene_id <- "DNThrb"
data_target=data_target[!is.na(rowData(data_target)$ensembl_gene_id),]

# Second, we remove the genes with duplicate ensemble_gene_ids
print(paste("Removing duplicated gene names:",sum(duplicated(rowData(data_source)$ensembl_gene_id)),"genes from the source and",sum(duplicated(rowData(data_target)$ensembl_gene_id)),"from the target dataset."))
data_source=data_source[!duplicated(rowData(data_source)$ensembl_gene_id),]
data_target=data_target[!duplicated(rowData(data_target)$ensembl_gene_id),]

# Third, we convert the gene names to ensembl gene ids
row.names(data_source)=rowData(data_source)$ensembl_gene_id
row.names(data_target)=rowData(data_target)$ensembl_gene_id

print(paste("Number of genes in common between the two datasets",length(intersect(row.names(data_source),row.names(data_target)))))

#Lastly, we only include genes that are in common between the two datasets (okay to exclude transgenes for purposes of cell subtype annotation)
data_source=data_source[which(row.names(data_source) %in% row.names(data_target)),]
data_target=data_target[which(row.names(data_target) %in% row.names(data_source)),]

# matching the row names between the two datasets
data_source=data_source[match(row.names(data_target),row.names(data_source)),]

# running scOnline label transfer function
res_prop_harmony=.myLabelTransfer_harmony(dataset_source = data_source,dataset_target = data_target,indScaling=T,source_label_col = "subclass_label",target_label_col = 'anno_cluster_res.1',return_seurat_obj=T,covariates = NULL,calculate_depth_per_gene = F)

# Generate plots to help evaluate quality of cell subtype annotation transfer from reference dataset
exNMarkerGenes=c("ENSMUSG00000029151","ENSMUSG00000021743","ENSMUSG00000016918","ENSMUSG00000028871")
exNMarkerGeneNames=c("Slc30a3","Fezf2","Sulf1","Rspo1")
plots=labelTransfer_plots(res=res_prop_harmony,clusteringRes='anno_cluster_res.1',confidenceLevel=0.8,markerGenes=exNMarkerGenes,markerGeneNames=exNMarkerGeneNames)
cowplot::plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]])