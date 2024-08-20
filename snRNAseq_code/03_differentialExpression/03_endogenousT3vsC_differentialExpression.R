# ========================================================
# Example script for performing pseudocell differential 
# expression within a specific cell class to compare effect
# of T3 vs. C treatment on endogenous samples
# ======================================================== 

source("~/00_scOnline_sourceCode/sconline_code.R")
source("~/00_scOnline_sourceCode/Human_mouse_alignment_code.R")
source("~/00_scOnline_sourceCode/integrated_analysis_functions.R")
source("~/03_differentialExpression/robustnessCalc.R")

data_seurat <- ... # load in Seurat object with cell subtype annotations and optimal PCA/harmony solution
pca_embeddings=data_seurat@reductions$pca@cell.embeddings
data_sce <- ... # load in SCE object with cell subtype annotations and optimal PCA/harmony solution

ps=suppressWarnings(.sconline.PseudobulkGeneration(argList=NULL, parsing.col.names = c("donor_id","anno"),
                                                     pseudocell.size=30,inputExpData=data_sce,min_size_limit=15,
                                                     inputPhenoData=as.data.frame(colData(data_sce)),
                                                     inputEmbedding=pca_embeddings,nPCs=30,ncores=9,
                                                     rand_pseudobulk_mod=T,organism="Mouse"))

ps$QC_Gene_total_log=log2(ps$QC_Gene_total_count)
ps$QC_Gene_unique_log=log2(ps$QC_Gene_unique_count)

ps_list=.mySplitObject(ps,'anno') #splitting the objects based on the cluster annotations


###### Endogenous T3 vs. C comparison ######
#remove all cells with irrelevant donor_ids for the endogenous T3 vs C comparison
ps_list=lapply(ps_list,function(x) {
  x=x[,x$THRB %in% c("N/A")] #THRB = N/A
  return(x)
})
print(unique(ps_list$inferred_L4.5.IT$donor_id)) # Check correct donor_ids [1] "13" "14" "15" "17" "18" "20" "21" "22"

#Next, we perform DE analysis using LimmaTrend method -- compare T3 to C in N/A mice
covList=c('treatment',"QC_Gene_total_log",'QC_MT.pct') # no sex needed as covariate because all male
bkgGene_pct_thr=0.01 
rand_var="donor_id"
quantile_norm=F

DE_res=lapply(names(ps_list),function(x_name){
  x=ps_list[[x_name]]
  res_arranged=NULL
  #requiring existence of at least 6 pseudocells for the DE analysis
  if(length(unique(x$treatment))>1&ncol(x)>5){
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

    if(sum(grepl("treatment",colnames(res$model)))>0){
      contr_thyroid <- makeContrasts(treatmentT3 - treatmentC, levels = res$model) # specify comparison of interest
      fit2_thyroid=contrasts.fit(res$fit,contrasts=contr_thyroid)
      fit2_thyroid=eBayes(fit2_thyroid,robust = T,trend=T)
      #running topTable for each contrast
      res_thyroid_TH=topTable(fit2_thyroid,number=dim(fit2_thyroid)[1], adjust.method = "BH", coef="treatmentT3 - treatmentC");
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