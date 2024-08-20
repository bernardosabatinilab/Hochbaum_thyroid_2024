# Performs orthogonal robustness analysis on SCE input object to assess confidence in DE genes
robustnessCalc=function(input_data) {
  input_exp_data=input_data
  res_robustness=list()
  # input_exp_data=input_exp_data[,!is.na(input_exp_data$anno)]
  print(paste0("#cells in input sce data object: ", dim(input_exp_data)[2]))
  input_exp_data=input_exp_data[,input_exp_data$THRB %in% c("N/A")]
  print(paste0("#cells in input sce data object once subsetting to only ",unique(input_exp_data$THRB)," cells is: ",dim(input_exp_data)[2]))
  
  print(paste0("Calculating robustness for ",unique(input_exp_data$amandaCellClass)))
  
  tmp_data=input_exp_data # NO NEED TO SUBSET BY ANNOTATION
  tmp_data_seurat=.extraExport2SeuratFn(inputData=input_exp_data,project="scRNA")
  tmp_data_seurat=Seurat::NormalizeData(tmp_data_seurat,verbose=F)
  
  if(sum(unique(input_exp_data$treatment) %in% c("T3","C"))==2){
    robustness_T3=.sconline.RobustFC(inputData=tmp_data,batch_col="donor_id",contrast_col="treatment",contrast_1="T3",contrast_2="C",sex_col=NULL,ncores=5,groupLevel=F)
    robustness_T3=do.call("rbind",robustness_T3)
    robustness_cell.count.1=aggregate(cell.count.1~gene,data=robustness_T3,mean)
    robustness_FCscore=aggregate(score_logFC~gene,data=robustness_T3,mean)
    robustness_PCTscore=aggregate(score_pct~gene,data=robustness_T3,mean)
    robustness_meanRefCount=aggregate(ref_count~gene,data=robustness_T3,mean)
    robustness_T3=merge(robustness_FCscore,robustness_PCTscore,by="gene")
    robustness_T3=merge(robustness_T3,robustness_meanRefCount,by="gene")
    robustness_T3=merge(robustness_T3,robustness_cell.count.1,by="gene")
    robustness_T3=robustness_T3[order(robustness_T3$score_pct,decreasing = T),]
    robustness_T3$treatment="T3"
    
    tmp_res_seurat=.myEvalMarkers(object=tmp_data_seurat, cells.1=colnames(tmp_data_seurat)[which(tmp_data_seurat@meta.data$treatment=="T3")], cells.2=colnames(tmp_data_seurat)[which(tmp_data_seurat@meta.data$treatment=="C")], slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cells.1.weight.col=NULL,cluster_name=NULL)
    colnames(tmp_res_seurat)=paste0("seurat_",colnames(tmp_res_seurat))
    tmp_res_seurat$gene=row.names(tmp_res_seurat)
    robustness_T3=merge(robustness_T3,tmp_res_seurat,by="gene",all.x=T)
    res_robustness=robustness_T3
  }
  return(res_robustness)
}
