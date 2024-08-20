# ========================================================
# Example script for how the data were annotated for gene
# markers of major cell classes were annotated from 
# initial QC and clustering results. 
# ======================================================== 

source("~/00_scOnline_sourceCode/sconline_code.R")
source("~/00_scOnline_sourceCode/Human_mouse_alignment_code.R")
source("~/00_scOnline_sourceCode/integrated_analysis_functions.R")

pcaSlnToSeurat <- qread("~/insert_directory_here/pcaSlnToSeurat.qs") # load in final, QC-ed Seurat object from 01_dataPreprocess.R

Idents(pcaSlnToSeurat)=pcaSlnToSeurat$anno_cluster_res.0.6

p1 <- DimPlot(pcaSlnToSeurat, reduction = "umap", group.by = "anno_batch")
p2 <- DimPlot(pcaSlnToSeurat, reduction = "umap", label = TRUE)
p1+p2

# Visualize QC metrics for removal of low quality cells 
qc_plot <- FeaturePlot(pcaSlnToSeurat, features = c("QC_MT.pct","QC_OXPHOS.pct"),cols = c("lightgrey","blue"),pt.size = 0.5,combine=FALSE)
cowplot::plot_grid(plotlist = qc_plot,ncol=3)

# Marker genes explored for major cell class annotation
# initialMarkerGenes <- c("ENSMUSG00000024411", "ENSMUSG00000029648","ENSMUSG00000036192", "ENSMUSG00000042589","ENSMUSG00000021743", "ENSMUSG00000070570","ENSMUSG00000019935",
#                         "ENSMUSG00000070880", "ENSMUSG00000026787", "ENSMUSG00000054675", "ENSMUSG00000046160",
#                         "ENSMUSG00000039830")
initialMarkerGenesNames <- c("Aqp4","Flt1","Rorb","Cux2","Fezf2","Slc17a7","Slc17a8","Gad1","Gad2","Tmem119","Olig1","Olig2")

# exNMarkerGenes <- c("ENSMUSG00000036192", "ENSMUSG00000042589","ENSMUSG00000030500", "ENSMUSG00000030500")
exNMarkerGenesNames <- c("Rorb","Cux2","Slc17a6","Slc17a7")

# micro_genes=c("ENSMUSG00000026395","ENSMUSG00000054675", "ENSMUSG00000024610","ENSMUSG00000026712","ENSMUSG00000070570")
micro_gene_names=c("Ptprc","Tmem119","Cd74","Mrc1","Slc17a7")

endo_gene_names=c("Flt1","Dcn","Kcnj8","Acta2","Ctla2a","Cytl1","Adamts19","Spp1","Col15a1","Osr1","Cd74")
# endo_genes=c("ENSMUSG00000029648","ENSMUSG00000019929","ENSMUSG00000030247","ENSMUSG00000035783","ENSMUSG00000044258","ENSMUSG00000062329",
#              "ENSMUSG00000053441","ENSMUSG00000029304","ENSMUSG00000028339","ENSMUSG00000048387","ENSMUSG00000024610")

p3 <- FeaturePlot(pcaSlnToSeurat, features = c(initialMarkerGenesNames),cols = c("lightgrey","blue"), max.cutoff = "q98", pt.size = 0.5, order=TRUE,combine=FALSE,slot = "data") 
for (i in 1:(length(c(initialMarkerGenesNames)))) {
  p3[[i]] <- p3[[i]] + labs(title=c(initialMarkerGenesNames)[[i]])
}
cowplot::plot_grid(plotlist = p3,ncol=3)

# Example of subsetting to major cell classes, excluding clusters with high %MT and high %OXPHOS
Idents(pcaSlnToSeurat)=pcaSlnToSeurat$anno_cluster_res.1
ExN_clusters <- as.character(c(1,2,6,9,12,28,
                               4,11,16,
                               17,3,
                               8,55))
ExN_subset <- subset(pcaSlnToSeurat,idents = ExN_clusters)