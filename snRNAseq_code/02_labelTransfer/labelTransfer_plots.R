labelTransfer_plots=function(res,clusteringRes,confidenceLevel,markerGenes,markerGeneNames) {
  seuratObj=res$seurat_obj
  print(class(seuratObj))
  
  # Generate cell assignment heatmap
  p1 <- .mycellAssignHeatmap_binary(input_labelTransfer_object=res,confidenceLevel=confidenceLevel,target_cluster_col=clusteringRes)
  
  # UMAP split by batch label
  # Idents(seuratObj)=seuratObj$anno_cluster_res.0.6
  p2 = DimPlot(seuratObj,group.by=clusteringRes,split.by ="batch_label",label=TRUE)
  
  # Generate inferred labels that include NA for assignments under confidence threshold
  inferredLbls=seuratObj@meta.data[,grepl("inferred_",colnames(seuratObj@meta.data))]
  inferredLblsIndx=apply(inferredLbls,1,function(x) { # apply the function over rows in lbls
    if(sum(!is.na(x))>0){ # if 
      if(max(x,na.rm = T)>confidenceLevel){
        which(x==max(x))[1]
      }else{NA}
    }else{NA}
  })
  amandaLbls=colnames(inferredLbls)[inferredLblsIndx]
  seuratObj@meta.data$amanda_inferred=amandaLbls
  
  # Plot umap with confidence threshold
  data_obj=seuratObj
  Idents(data_obj)=data_obj$anno
  p3=DimPlot(data_obj,split.by ="batch_label")
  
  # Marker plot to check quality of label transfer
  seuratObj_target=subset(seuratObj,subset=batch_label=="target")
  Idents(seuratObj_target)=seuratObj_target$anno
  markerPlot=FeaturePlot(seuratObj_target, features = c(markerGenes),cols = c("lightgrey","blue"), max.cutoff = "q90", pt.size = 0.5, order=TRUE,combine=FALSE) 
  for (i in 1:length(markerGenes)) {
    markerPlot[[i]] <- markerPlot[[i]] + labs(title=markerGeneNames[[i]])
  }
  p4=cowplot::plot_grid(plotlist = markerPlot)
  
  # cowplot::plot_grid(p1,p2,p3,p4)
  return(list(p1,p2,p3,p4))
}
