# ========================================================
# Example script for performing pseudocell differential 
# expression within a specific cell class to compare effect
# of Cre+ vs Cre- within T3+ DN-THR samples
# ======================================================== 

source("~/00_scOnline_sourceCode/sconline_code.R")
source("~/03_differentialExpression/robustnessCalc.R")

data_seurat <- ... # load in Seurat object with cell subtype annotations and optimal PCA/harmony solution
# Annotate Cre+ and Cre-
cre_ExN=as.data.frame(data_seurat@assays$RNA@counts["iCre",])
colnames(cre_ExN)="cre_counts"
cre_ExN$creStatus=ifelse(cre_ExN$cre_counts >= 1, "cre_pos", ifelse(cre_ExN$cre_counts == 0, "cre_neg", cre_ExN$cre_counts))
data_seurat$creStatus=cre_ExN$creStatus
# Convert to SCE
data_sce=.myExpSetCreatorFn(inputExpData=data_seurat@assays$RNA@counts,organism="Mouse",minExpCells=0,inputPdata=as.data.frame(data_ExN_addPCA@meta.data[,!grepl("QC",colnames(data_ExN_addPCA@meta.data))]),inputFdata=NULL,addExtraAnno=T,server=T,redownload_files=F,ncores=6)

pca_embeddings=data_seurat@reductions$pca@cell.embeddings
ps=suppressWarnings(.sconline.PseudobulkGeneration(argList=NULL, parsing.col.names = c("donor_id","anno","creStatus"),
                                                     pseudocell.size=30,inputExpData=data_sce,min_size_limit=15,
                                                     inputPhenoData=as.data.frame(colData(data_sce)),
                                                     inputEmbedding=pca_embeddings,nPCs=30,ncores=9,
                                                     rand_pseudobulk_mod=T,organism="Mouse"))

ps$QC_Gene_total_log=log2(ps$QC_Gene_total_count)
ps$QC_Gene_unique_log=log2(ps$QC_Gene_unique_count)

ps_list=.mySplitObject(ps,'anno') #splitting the objects based on the cluster annotations


###### T3+ DN-THR Cre+ vs Cre- comparison ######
#remove all cells with irrelevant donor_ids 
ps_list=lapply(ps_list,function(x) {
  x=x[,x$THRB %in% c("DN-THRB")] 
  return(x)
})
print(unique(ps_list$inferred_L4.5.IT$donor_id)) # Check for correct donor_ids based on metadata
#remove all cells without +T3 treatment
ps_list=lapply(ps_list,function(x) {
  x=x[,x$treatment=="T3"]
  return(x)
})
print(unique(ps_list$inferred_L4.5.IT$donor_id)) # Check for correct donor_ids based on metadata

#Next, we perform DE analysis using LimmaTrend method -- compare Cre+ to Cre- in T3+ DN-THRB mice
covList=c('creStatus',"QC_Gene_total_log",'QC_MT.pct') # no sex needed as covariate because all male
bkgGene_pct_thr=0.01 # or make 0.01?
rand_var="donor_id"
quantile_norm=F

DE_res=lapply(names(ps_list),function(x_name){
  x=ps_list[[x_name]]
  res_arranged=NULL
  #requiring existence of at least 6 pseudocells for the DE analysis
  if(length(unique(x$creStatus))>1&ncol(x)>5){
    print(paste0("Calculating DE genes for ",x_name))
    
    if(!is.null(bkgGene_pct_thr)){
      tmp_bkg_genes=counts(data_sce)[,which(data_sce$anno==x_name)]
      tmp_bkg_genes=rowSums(tmp_bkg_genes>0)/sum(data_sce$anno==x_name)
      tmp_bkg_genes=row.names(data_sce)[tmp_bkg_genes>=bkgGene_pct_thr]
      
    } else {
      tmp_bkg_genes=NULL
    }
    print(paste0("Number of bkg genes used in analysis is ",length(tmp_bkg_genes)))
    
    res=.sconline.fitLimmaFn(inputExpData=x,covariates=covList,randomEffect=rand_var,bkg_genes=tmp_bkg_genes,quantile.norm=quantile_norm,prior.count=1)
    print(paste0("Consensus correlation is ",res$dc$consensus.correlation))

    if(sum(grepl("creStatus",colnames(res$model)))>0){
      contr_thyroid <- makeContrasts(creStatuscre_pos - creStatuscre_neg, levels = res$model) # specify comparison of interest
      fit2_thyroid=contrasts.fit(res$fit,contrasts=contr_thyroid)
      fit2_thyroid=eBayes(fit2_thyroid,robust = T,trend=T)
      #running topTable for each contrast
      res_thyroid_TH=topTable(fit2_thyroid,number=dim(fit2_thyroid)[1], adjust.method = "BH", coef="creStatuscre_pos - creStatuscre_neg");
      res_thyroid_TH$gene=row.names(res_thyroid_TH)
      res_thyroid_TH$blocked_analysis=res$blocked_analysis
      res_thyroid_TH$block.cor=res$dc$consensus.correlation 
      res_thyroid_TH$cluster=unique(as.character(x$anno))
      res_thyroid_TH = res_thyroid_TH %>% relocate("gene")
      res_arranged=res_thyroid_TH
    }
  }
  return(res_arranged)
})