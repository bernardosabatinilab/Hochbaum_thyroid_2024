.myArgCreatorFn=function(HVG_list,prefix,nPCs=60,indScaling=T,includeHuman,includeMouse,saveDir,min_ds_size=50,min_ds_genes=5000,ncores=6,newRun=F,input_highly_var_genes=NULL,covariates=NULL){
  #CB var genes: "~/myBucket/varGenes.rda"
  runIndx=26 
  FinnishSbjBased=F
  DE_supportingFractionThr=0.1
  DE_n.adaptiveKernel=20
  pval_thr=0.001
  varGenesDir=NULL
  if(indScaling){
    saveDir=file.path(saveDir,paste0(prefix,"_indScaling"))
  } else {
    saveDir=file.path(saveDir,prefix)
  }
  
  saveDir=paste0(saveDir,"_nPCs",nPCs)
  
  
  #####################
  #could be important!
  conservativeMapping=F
  oldMapping=F
  MGIsymbol=T
  ######################
  
  #prefix="MGonly_"
  #prefix="PCprojection2"
  #prefix="organismConserved3"
  #runIndx=runIndx;exNonMicCells=F;conservativeMapping = conservativeMapping;oldMapping=oldMapping;MGIsymbol=MGIsymbol;ncores=ncores;sensitiveSearch=sensitiveSearch;includeHuman=includeHuman;includeMouse=includeMouse;FinnishSbjBased=FinnishSbjBased;uniformZscore=uniformZscore;DE_supportingFractionThr=DE_supportingFractionThr;DE_n.adaptiveKernel=DE_n.adaptiveKernel;DE_nPropIter=DE_nPropIter;dist_zscore_gamma=dist_zscore_gamma;dist_zscore_norm=dist_zscore_norm;dist_zscore_nbinom=dist_zscore_nbinom;regularize=regularize;geoMean=geoMean;prefix=prefix;newRun = newRun;inputDf=NULL;pseudocell_count=pseudocell_count;external_DE_path=external_DE_path;external_DE_name=external_DE_name;saveDir=saveDir
  if(!is.null(varGenesDir)){
    load(varGenesDir)
    input_highly_var_genes=res$varGenes
  }
  
  .ArgList=.myArgFn(runIndx=runIndx,nPCs=nPCs,conservativeMapping = conservativeMapping,oldMapping=oldMapping,MGIsymbol=MGIsymbol,ncores=ncores,includeHuman=includeHuman,includeMouse=includeMouse,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,prefix=prefix,newRun = newRun,saveDir=saveDir)
  
  .ArgList$input_highly_var_genes=input_highly_var_genes
  .ArgList$indScaling=indScaling
  .ArgList$saveDirGlobal=paste0(saveDir,"_Global")
  
  
  
  .ArgList$min_ds_size=min_ds_size
  .ArgList$min_ds_genes=min_ds_genes
  .ArgList$HVG_list=HVG_list
  
  .ArgList$allGenesFraction=0.5
  
  ##############
  .ArgList$covariates=covariates
  ###############
  
  
  if(!dir.exists(.ArgList$saveDirGlobal)){
    dir.create(.ArgList$saveDirGlobal,recursive = T)
  }
  return(.ArgList)
}

.myUMAPgeneratorFn=function(argList,pca_UMAP=T,indxList=c(1:4),umap.method='uwot',newRun=F){
  
  if(pca_UMAP){
    for(iHVG in argList$HVG_list){
      argList$HVG_count=iHVG
      if(!file.exists(.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))){
        load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
      } else {
        load(.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))
      }
      
      if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2&(!newRun)){
        stop("UMAP_1 and UMAP_2 columns already exists in the data")
      }
      
      tst=.reductionUMAPFn(pca_res[,1:argList$nPCs],umap.method='uwot')
      resUMAP=tst$embedding
      pd$UMAP_1=resUMAP[,1]
      pd$UMAP_2=resUMAP[,2]
      save(pd,file=.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))
      
    }
  } else {
    for(indx in indxList){
      for(iHVG in argList$HVG_list){
        argList$HVG_count=iHVG
        if(!file.exists(.myFilePathMakerFn(paste0("harmony_umap_",indx),argList=argList,pseudoImportant = F))){
          load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
        } else {
          load(.myFilePathMakerFn(paste0("harmony_umap_",indx),argList=argList,pseudoImportant = F))
        }
        
        if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2&(!newRun)){
          stop("UMAP_1 and UMAP_2 columns already exists in the data")
        }
        
        pca_res=qread(.myFilePathMakerFn(paste0("harmony_embedding_",indx),argList=argList,pseudoImportant = F,qsFormat=T))
        
        tst=.reductionUMAPFn(pca_res[,1:argList$nPCs],umap.method='uwot')
        resUMAP=tst$embedding
        pd$UMAP_1=resUMAP[,1]
        pd$UMAP_2=resUMAP[,2]
        save(pd,file=.myFilePathMakerFn(paste0("harmony_umap_",indx),argList=argList,pseudoImportant = F))
        
      }
    }
  }
  
  return("Done!")
}

.myQCfn=function(argList,pca_UMAP=T,indxList=c(1:4),colname_list=c("batch_merging"),complementrary_pd=NULL){
  
  myQCplotMakerFn=function(inputPd,colname,argList,pca_UMAP=T,indx=NULL){
    inputPd$facet_col=inputPd[,colname]
    p=ggplot(inputPd,aes(UMAP_1,UMAP_2))+geom_point(size=0.1)+facet_wrap(~facet_col)+theme_classic()
    if(pca_UMAP){
      ggsave(plot=p,.myFilePathMakerFn(paste0("QC_PCA_",colname),argList=argList,pseudoImportant = F,pdf=T),width = 10,height = 10)
    } else {
      ggsave(plot=p,.myFilePathMakerFn(paste0("QC_Harmony_",colname,"_",indx),argList=argList,pseudoImportant = F,pdf=T),width = 10,height = 10)
    }
    
    p=.my2dPlot2(inputPCA = inputPd[!is.na(inputPd[,colname]),],batchCol = colname,reductionCols = c('UMAP_1','UMAP_2'))
    if(pca_UMAP){
      ggsave(plot=p,.myFilePathMakerFn(paste0("QC_PCA2d_",colname),argList=argList,pseudoImportant = F,pdf=T),width = 10,height = 10)
    } else {
      ggsave(plot=p,.myFilePathMakerFn(paste0("QC_Harmony2d_",colname,"_",indx),argList=argList,pseudoImportant = F,pdf=T),width = 10,height = 10)
    }
    
    return("Done")
    
  }
  
  if(!is.null(complementrary_pd)){
    colname_list=colnames(complementrary_pd)
  }
  
  for(iHVG in unique(argList$HVG_list)){
    argList$HVG_count=iHVG
    for(indx in indxList){
      for(icol in colname_list){
        if(pca_UMAP){
          if(!file.exists(.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))){
            stop("UMAP coordinates need to be generated first")
          }
          
          #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
          load(.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))
          if(!is.null(complementrary_pd)){
            pd=pd[match(row.names(complementrary_pd),row.names(pd)),]
            complementrary_pd$UMAP_1=pd$UMAP_1
            complementrary_pd$UMAP_2=pd$UMAP_2
            pd=complementrary_pd
          }
          tmp_res=myQCplotMakerFn(inputPd=pd,colname=icol,argList=argList,pca_UMAP=T,indx=NULL)
        } else {
          if(!file.exists(.myFilePathMakerFn(paste0("harmony_umap_",indx),argList=argList,pseudoImportant = F))){
            stop("UMAP coordinates need to be generated first")
          }
          #pca_res=qread(.myFilePathMakerFn(paste0("harmony_embedding_",indx),argList=argList,pseudoImportant = F,qsFormat=T))
          load(.myFilePathMakerFn(paste0("harmony_umap_",indx),argList=argList,pseudoImportant = F))
          if(!is.null(complementrary_pd)){
            pd=pd[match(row.names(complementrary_pd),row.names(pd)),]
            complementrary_pd$UMAP_1=pd$UMAP_1
            complementrary_pd$UMAP_2=pd$UMAP_2
            pd=complementrary_pd
          }
          tmp_res=myQCplotMakerFn(inputPd=pd,colname=icol,argList=argList,pca_UMAP=F,indx=indx)
        }
        
      }
    }
  }
  
  
  return("Done")
}

.myClusteringFn=function(argList,pca_level=T,indxList=c(1:4),newRun=F,cluster_resolution=c(0.6,1),group.singletons=T){
  require(cowplot)
  myInternalFn=function(inputPd,inputPCA,pca_level,argList,indx,newRun=F,group.singletons=T,cluster_resolution=cluster_resolution){
    if(sum(grepl("anno_cluster_res",colnames(inputPd)))==0|newRun){
      
      harmony_net=Seurat:::FindNeighbors.default(object=inputPCA[row.names(inputPd),1:argList$nPCs])
      clusters=.netFindClusters(inputGraph=harmony_net[["snn"]],resolution=cluster_resolution, algorithm = 1,group.singletons = group.singletons,modularity.fxn = 1)
      colnames(clusters)=paste0("anno_cluster_",colnames(clusters))
      inputPd=cbind(inputPd,clusters)
    } else {
      warning("anno_cluster_res already exists, skipping the seurat clustering")
    }
    
    for(icol in colnames(inputPd)[grepl("anno_cluster_res",colnames(inputPd))]){
      p=.myDimPlotFn(object=inputPd, dimCols = c("UMAP_1", "UMAP_2"), attCol = icol, 
                     pt.size = NULL, label = TRUE, label.size = 4, 
                     repel = FALSE, cells.highlight = NULL, cols.highlight = "#DE2D26", 
                     sizes.highlight = 1, na.value = "grey50", ncol = NULL, combine = TRUE)
      
      if(pca_level){
        ggsave(plot=p,.myFilePathMakerFn(paste0("QC_PCA_seurat",gsub("anno_cluster_res","cluster",icol)),argList=argList,pseudoImportant = F,pdf=T),width = 10,height = 10)
      } else {
        ggsave(plot=p,.myFilePathMakerFn(paste0("QC_Harmony_seurat",gsub("anno_cluster_res","cluster",icol),"_",indx),argList=argList,pseudoImportant = F,pdf=T),width = 10,height = 10)
      }
      
      tmp_pd=inputPd
      tmp_pd$anno_cluster_res=inputPd[,icol]
      p1=ggplot(tmp_pd,aes(anno_cluster_res,QC_MT.pct,fill=anno_cluster_res))+geom_violin()+theme_classic()+theme(legend.position = "none")
      p2=ggplot(tmp_pd,aes(anno_cluster_res,QC_IEG.pct,fill=anno_cluster_res))+geom_violin()+theme_classic()+theme(legend.position = "none")
      p3=ggplot(tmp_pd,aes(anno_cluster_res,QC_top50_pct,fill=anno_cluster_res))+geom_violin()+theme_classic()+theme(legend.position = "none")
      p=p1+p2+p3+ plot_layout(nrow = 3, byrow = FALSE)
      if(pca_level){
        ggsave(plot=p,.myFilePathMakerFn(paste0("QC_PCA_seurat",gsub("anno_cluster_res","cluster",icol),"_MTpct"),argList=argList,pseudoImportant = F,pdf=T),width = 10,height = 10)
      } else {
        ggsave(plot=p,.myFilePathMakerFn(paste0("QC_Harmony_seurat",gsub("anno_cluster_res","cluster",icol),"_MTpct_",indx),argList=argList,pseudoImportant = F,pdf=T),width = 10,height = 10)
      }
      
      
    }
    
    
    
    if(sum(!is.na(inputPd$anno_cellState))>0){
      if(sum(is.na(inputPd$anno_cellState))>0){
        inputPd$anno_cellState[is.na(inputPd$anno_cellState)]="N/A"
      }
      inputPd=inputPd[!is.na(inputPd$anno_cellState),]
      p=.myDimPlotFn(object=inputPd, dimCols = c("UMAP_1", "UMAP_2"), attCol = 'anno_cellState', 
                     pt.size = NULL, label = TRUE, label.size = 4, 
                     repel = FALSE, cells.highlight = NULL, cols.highlight = "#DE2D26", 
                     sizes.highlight = 1, na.value = "grey50", ncol = NULL, combine = TRUE)+theme(legend.position="none")
      if(pca_level){
        ggsave(plot=p,.myFilePathMakerFn(paste0("QC_PCA_cellState"),argList=argList,pseudoImportant = F,pdf=T),width = 10,height = 10)
      } else {
        ggsave(plot=p,.myFilePathMakerFn(paste0("QC_Harmony_cellState_",indx),argList=argList,pseudoImportant = F,pdf=T),width = 10,height = 10)
      }
      
    }
    
    return(inputPd)
    
  }
  
  for(iHVG in argList$HVG_list){
    argList$HVG_count=iHVG
    load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    
    if(pca_level){
      if(file.exists(.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))){
        load(.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))
      }
      
      pd=myInternalFn(inputPd = pd,pca_level = T,inputPCA = pca_res,argList = argList,newRun = newRun,cluster_resolution=cluster_resolution,group.singletons=group.singletons)
      save(pd,file=.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))
    } else {
      for(indx in indxList){
        pca_res=qread(.myFilePathMakerFn(paste0("harmony_embedding_",indx),argList=argList,pseudoImportant = F,qsFormat=T))
        if(file.exists(.myFilePathMakerFn(paste0("harmony_umap_",indx),argList=argList,pseudoImportant = F))){
          load(.myFilePathMakerFn(paste0("harmony_umap_",indx),argList=argList,pseudoImportant = F))
        }
        
        pd=myInternalFn(inputPd = pd,indx = indx,inputPCA = pca_res,pca_level = F,argList = argList,newRun = newRun,cluster_resolution=cluster_resolution,group.singletons=group.singletons)
        save(pd,file=.myFilePathMakerFn(paste0("harmony_umap_",indx),argList=argList,pseudoImportant = F))
      }
    }
    
  }
  
  
  return("Done")
}

.mySeuratMakerFn=function(inputExpData,argList,select_pca,select_HVG,select_harmony_index){
  argList$HVG_count=select_HVG
  
  
  harmony_res=NULL
  arranged_data=.mycBindFn(inputExpData$data)
  
  arranged_data=.extraExport2SeuratFn(inputData=arranged_data,project="scRNA")
  
  load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
  
  if(select_pca){
    if(file.exists(.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))){
      load(.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))
    }
    
  } else {
    pca_res=qread(.myFilePathMakerFn(paste0("harmony_embedding_",select_harmony_index),argList=argList,pseudoImportant = F,qsFormat=T))
    if(file.exists(.myFilePathMakerFn(paste0("harmony_umap_",select_harmony_index),argList=argList,pseudoImportant = F))){
      load(.myFilePathMakerFn(paste0("harmony_umap_",select_harmony_index),argList=argList,pseudoImportant = F))
    }
  }
  
  pd=pd[match(colnames(arranged_data),row.names(pd)),]
  pca_res=pca_res[match(colnames(arranged_data),row.names(pca_res)),]
  if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2){
    umapdata=as.matrix(pd[,c("UMAP_1","UMAP_2")])
    umap.reduction <- CreateDimReducObject(embeddings = umapdata, 
                                           key = "UMAP_", assay = "RNA", global = TRUE)
    arranged_data[["umap"]]=umap.reduction
  }
  
  if(sum(!grepl("^PC_",colnames(pca_res)))>0){
    colnames(pca_res)=gsub("PC_","",colnames(pca_res))
  }
  
  reduction.data <- CreateDimReducObject(embeddings = pca_res, assay = "RNA", key = "PC_")
  arranged_data[["pca"]] <- reduction.data
  
  if(!is.null(harmony_res)&F){
    harmony_res=pca_res[match(colnames(arranged_data),row.names(pca_res)),]
    
    if(sum(!grepl("^PC_",colnames(harmony_res)))>0){
      colnames(harmony_res)=gsub("PC_","",colnames(harmony_res))
    }
    
    reduction.data <- CreateDimReducObject(embeddings = harmony_res, assay = "RNA", key = "PC_")
    arranged_data[["harmony_res"]] <- reduction.data
  }
  
  pd=pd[,!colnames(pd) %in% c("UMAP_1","UMAP_2")]
  #if(length(setdiff(colnames(pd),colnames(arranged_data@meta.data)))>0){
  #arranged_data@meta.data=cbind(arranged_data@meta.data,pd)
  #}
  
  arranged_data@meta.data=pd
  
  return(arranged_data)
  
}

.myHarmonyFn_1var=function(iHVG,argList,batch_column){
  argList$HVG_count=iHVG
  load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
  #pd=pd[,c('batch_merging',"ds_batch")]
  gc()
  harmony_gridSearch=function(indx,inputDf,pca_res,pd,argList){
    res=harmony::HarmonyMatrix(pca_res, pd, c('ds_batch'),theta=c(inputDf$theta1[indx]),lambda=c(inputDf$lambda1[indx]), do_pca = FALSE, verbose=FALSE)
    qsave(res,file=.myFilePathMakerFn(paste0("harmony_embedding_",indx),argList=argList,pseudoImportant = F,qsFormat=T))
    return("Done")
  }
  
  dfSettings=NULL
  dfSettings=rbind(dfSettings,data.frame(theta1=2,lambda1=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=4,lambda1=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=6,lambda1=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=8,lambda1=1))
  
  #indx=1;inputDf=dfSettings;pca_res=pca_res[,1:30]
  res=parallel::mclapply(1:nrow(dfSettings),harmony_gridSearch,inputDf=dfSettings,pca_res=pca_res[,1:argList$nPCs],argList=argList,pd=pd,mc.cores = 7)
  return(res)
}

.myHarmonyFn_2var=function(iHVG,argList,batch_columns){
  argList$HVG_count=iHVG
  load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
  #pd=pd[,c('batch_merging',"ds_batch")]
  gc()
  harmony_gridSearch=function(indx,inputDf,pca_res,pd,argList){
    res=harmony::HarmonyMatrix(pca_res, pd, batch_columns,theta=c(inputDf$theta1[indx],inputDf$theta2[indx]),lambda=c(inputDf$lambda1[indx],inputDf$lambda2[indx]), do_pca = FALSE, verbose=FALSE)
    qsave(res,file=.myFilePathMakerFn(paste0("harmony_embedding_",indx),argList=argList,pseudoImportant = F,qsFormat=T))
    return("Done")
  }
  
  dfSettings=NULL
  dfSettings=rbind(dfSettings,data.frame(theta1=2,theta2=2,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=4,theta2=2,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=6,theta2=2,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=8,theta2=2,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=2,theta2=4,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=2,theta2=6,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=2,theta2=8,lambda1=1,lambda2=1))
  
  dfSettings=dfSettings[c(1,2,3,5),]
  
  #indx=1;inputDf=dfSettings;pca_res=pca_res[,1:30]
  res=parallel::mclapply(1:nrow(dfSettings),harmony_gridSearch,inputDf=dfSettings,pca_res=pca_res[,1:argList$nPCs],argList=argList,pd=pd,mc.cores = 7)
  return(res)
}

.myMarkerPlot=function(argList,input_data,marker_genes,pca_level=T,indxList=c(1:4),newRun=F,fig_width=49,fig_height=40){
  
  for(iHVG in argList$HVG_list){
    argList$HVG_count=iHVG
    load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    
    if(pca_level){
      if(file.exists(.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))){
        load(.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))
      }
      
      tmp=input_data$data
      names(tmp)=unlist(lapply(tmp,function(x) as.character(x$anno_batch[1])))
      pd2=pd
      
      x=unname(unlist(lapply(tmp,function(x) colnames(x))))
      if(length(setdiff(pd2$sample,x))==0){
        row.names(pd2)=pd2$sample
      }
      
      if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename=paste0("PCA_markers"),pdf=T),argList)))){
        p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(marker_genes,function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6,ncores=argList$ncores)
        #ggsave(plot=p,file="~/myBucket/torm1.png",width = 49,height = 40)
        ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("PCA_markers"),pdf=T),argList)),width = fig_width,height = fig_height)
        
      }
      
    } else {
      for(indx in indxList){
        pca_res=qread(.myFilePathMakerFn(paste0("harmony_embedding_",indx),argList=argList,pseudoImportant = F,qsFormat=T))
        if(file.exists(.myFilePathMakerFn(paste0("harmony_umap_",indx),argList=argList,pseudoImportant = F))){
          load(.myFilePathMakerFn(paste0("harmony_umap_",indx),argList=argList,pseudoImportant = F))
        }
        
        tmp=input_data$data
        names(tmp)=unlist(lapply(tmp,function(x) as.character(x$anno_batch[1])))
        pd2=pd
        
        x=unname(unlist(lapply(tmp,function(x) colnames(x))))
        if(length(setdiff(pd2$sample,x))==0){
          row.names(pd2)=pd2$sample
        }
        
        if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename=paste0("harmony_markers_",indx),pdf=T),argList)))){
          p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(marker_genes,function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6,ncores=argList$ncores)
          #ggsave(plot=p,file="~/myBucket/torm1.png",width = 49,height = 40)
          ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("harmony_markers_",indx),pdf=T),argList)),width = fig_width,height = fig_height)
          
        }
        
      }
    }
    
    
    
    
  }
  
  
  
  return("Done")
}

.myDotPlot_normData=function (norm_data,cluster_anno,cluster_col,features, cols = c("blue","red"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                              scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA,re_order_markers=T) {
  
  require(cowplot)
  scale.func <- switch(EXPR = scale.by, size = scale_size, radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  data.features=as.data.frame(t(norm_data))
  cluster_anno=as.data.frame(cluster_anno)
  cluster_anno=cluster_anno[match(row.names(data.features),row.names(cluster_anno)),]
  if(sum(is.na(cluster_anno[,cluster_col]))>0){
    stop("Clustering information was not found for some of the cells")
  } else {
    data.features$id=factor(cluster_anno[,cluster_col])
  }
  
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (re_order_markers) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = cor(t(mat))))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  color.by <- "avg.exp.scaled"
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                     color = color.by)) + scale.func(range = c(0, dot.scale), 
                                                                                                                                     limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                               axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = "Identity") + theme_cowplot()
  
  if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  
  return(plot)
}


.myDotPlot_slMarkers_listData=function(argList,input_data,markers,pca_level=T,new_run=T,gene_name_col="gene_short_name",indxList=c(1:4),clustering_resolution=0.6,
                                       cols = c("blue","red"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                                       scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA,fig_width=10,fig_height=17){
  
  .local_cbindFn=function(data1,data2){
    data1=as.matrix(data1)
    data2=as.matrix(data2)
    if(length(setdiff(row.names(data1),row.names(data2)))>0){
      tmp_data2=matrix(0,nrow=length(setdiff(row.names(data1),row.names(data2))),ncol=ncol(data2))
      row.names(tmp_data2)=setdiff(row.names(data1),row.names(data2))
      colnames(tmp_data2)=colnames(data2)
      data2=rbind(data2,tmp_data2)
      rm(tmp_data2)
    }
    
    if(length(setdiff(row.names(data2),row.names(data1)))>0){
      tmp_data1=matrix(0,nrow=length(setdiff(row.names(data2),row.names(data1))),ncol=ncol(data1))
      row.names(tmp_data1)=setdiff(row.names(data2),row.names(data1))
      colnames(tmp_data1)=colnames(data1)
      data2=rbind(data1,tmp_data1)
      rm(tmp_data1)
    }
    
    data2=data2[match(row.names(data1),row.names(data2)),]
    
    res=cbind(data1,data2)
    res=as(res,Class = "dgCMatrix")
    return(res)
  }
  
  tmp=input_data$data
  
  fd=NULL
  if(!is.null(gene_name_col)){
    fd=lapply(input_data$data,function(x){
      x=x@assays$RNA@meta.features
      x$rwname=row.names(x)
      x=x[,c('rwname',gene_name_col)]
      x=x[which(toupper(x[,gene_name_col]) %in% toupper(markers)),]
      return(x)
    })
    fd=do.call("rbind",fd)
    fd=fd[!duplicated(fd$rwname),]
    fd=fd[!duplicated(fd[,gene_name_col]),]
    markers=markers[markers %in% fd[,gene_name_col]]
    fd=fd[match(markers,fd[,gene_name_col]),]
    markers=fd$rwname
  }
  
  if(length(markers)<1){
    stop("None of the markers were identified!")
  }
  
  norm_data=tmp[[1]]@assays$RNA@data[markers,,drop=F]
  not_found_counter=0
  if(length(tmp)>1){
    for(i in 2:length(tmp)){
      if(sum(row.names(tmp[[i]]@assays$RNA@data) %in% markers)==0){
        not_found_counter=not_found_counter+1
      }
      norm_data=.local_cbindFn(norm_data,tmp[[i]]@assays$RNA@data[markers,,drop=F])
    }
  }
  
  markers=markers[which(toupper(markers) %in% toupper(row.names(norm_data)))]
  norm_data=norm_data[match(toupper(markers),toupper(row.names(norm_data))),]
  if(!is.null(fd)){
    fd=fd[match(toupper(row.names(norm_data)),toupper(fd$rwname)),]
    row.names(norm_data)=fd[,gene_name_col]
    markers=row.names(norm_data)
  }
  
  if(not_found_counter>0){
    warning(paste("Marker genes were not found in",not_found_counter,"batches!"))
  }
  
  for(iHVG in argList$HVG_list){
    argList$HVG_count=iHVG
    load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    
    if(pca_level){
      if(file.exists(.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))){
        load(.myFilePathMakerFn("pca_umap",argList=argList,pseudoImportant = F))
      }
      
      
      
      if(sum(colnames(pd)==paste0("anno_cluster_res.",clustering_resolution))==0){
        stop("specified resolution value could not be found!")
      }
      
      
      if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename=paste0("pca_dotPlot"),pdf=T),argList)))|new_run){
        clust_count=as.data.frame(table(pd[,paste0("anno_cluster_res.",clustering_resolution)]))
        clust_count=clust_count[which(clust_count[,2]>10),]
        pd=pd[pd[,paste0("anno_cluster_res.",clustering_resolution)] %in% as.character(clust_count[,1]),]
        norm_data=norm_data[,which(colnames(norm_data) %in% row.names(pd))]
        p=.myDotPlot_normData(norm_data=norm_data,cluster_anno=pd,cluster_col=paste0("anno_cluster_res.",clustering_resolution),features=markers, cols = cols, col.min = col.min, col.max = col.max, dot.min = dot.min, dot.scale = dot.scale, 
                              scale = scale, scale.by = scale.by, scale.min = scale.min, scale.max = scale.max)
        #ggsave(plot=p,file="~/myBucket/torm1.png",width = 49,height = 40)
        p1=p+theme(axis.text.x=element_text(angle=90,vjust=0.6,hjust=1))
        ggsave(plot=p1,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("pca_dotPlot"),pdf=T),argList)),width = fig_width,height = fig_height)
        
      }
      
    } else {
      for(indx in indxList){
        pca_res=qread(.myFilePathMakerFn(paste0("harmony_embedding_",indx),argList=argList,pseudoImportant = F,qsFormat=T))
        if(file.exists(.myFilePathMakerFn(paste0("harmony_umap_",indx),argList=argList,pseudoImportant = F))){
          load(.myFilePathMakerFn(paste0("harmony_umap_",indx),argList=argList,pseudoImportant = F))
        }
        
        if(sum(colnames(pd)==paste0("anno_cluster_res.",clustering_resolution))==0){
          stop("specified resolution value could not be found!")
        }
        
        
        
        if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename=paste0("harmony_dotPlot_",indx),pdf=T),argList)))|new_run){
          clust_count=as.data.frame(table(pd[,paste0("anno_cluster_res.",clustering_resolution)]))
          clust_count=clust_count[which(clust_count[,2]>10),]
          pd=pd[pd[,paste0("anno_cluster_res.",clustering_resolution)] %in% as.character(clust_count[,1]),]
          norm_data=norm_data[,which(colnames(norm_data) %in% row.names(pd))]
          p=.myDotPlot_normData(norm_data=norm_data,cluster_anno=pd,cluster_col=paste0("anno_cluster_res.",clustering_resolution),features=markers, cols = cols, col.min = col.min, col.max = col.max, dot.min = dot.min, dot.scale = dot.scale, 
                                scale = scale, scale.by = scale.by, scale.min = scale.min, scale.max = scale.max)
          #ggsave(plot=p,file="~/myBucket/torm1.png",width = 49,height = 40)
          p1=p+theme(axis.text.x=element_text(angle=90,vjust=0.6,hjust=1))
          
          #ggsave(plot=p,file="~/myBucket/torm1.png",width = 49,height = 40)
          ggsave(plot=p1,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("harmony_dotPlot_",indx),pdf=T),argList)),width = fig_width,height = fig_height)
        }
      }
    }
  }
  return("Done")
}