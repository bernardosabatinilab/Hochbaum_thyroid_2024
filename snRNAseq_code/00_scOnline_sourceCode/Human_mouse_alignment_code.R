.myArgFn=function(runIndx,nPCs=30,conservativeMapping=F,oldMapping=F,MGIsymbol=F,ncores=7,includeHuman,includeMouse,FinnishSbjBased,DE_supportingFractionThr,DE_n.adaptiveKernel,Leng200=F,Jakel200=F,prefix=NULL,newRun=F,saveDir=NULL){
  
  .dfSetting=.myRunSettingsComplete_mouse()
  
  .commonExpressed=.dfSetting$commonExpressed[runIndx]
  .includeHuman=includeHuman
  .includeMouse=includeMouse
  .nPCs=nPCs
  .HVG_count=.dfSetting$HVG_count[runIndx] #4
  .HVG_method="vst" #"vst" and "sctransform"
  
  .covariates=.dfSetting$covariates[runIndx]
  .removeHighExp=.dfSetting$removeHighExp[runIndx]
  .UMI_cor_thr=.dfSetting$UMI_cor_thr[runIndx]
  
  if(!is.null(.covariates)){
    if(all(.covariates=='3')){
      .covariates=c('QC_Gene_total_count','QC_top50_pct','QC_Gene_unique_count')
    } else if(all(.covariates=='2')) {
      .covariates=c('QC_top50_pct','QC_Gene_unique_count')
    } else if(all(.covariates=='1')){
      .covariates=c('QC_Gene_unique_count')
    } else if(all(.covariates=="")){
      .covariates=NULL
    }
  }
  
  
  
  #"depth_per_gene"#'QC_Gene_unique_count'#c('QC_Gene_total_count','QC_top50_pct',)#, "organism","depth_per_gene" 'QC_Gene_total_count','QC_top50_pct','QC_Gene_unique_count'
  .indScaling=.dfSetting$indScaling[runIndx]
  .indScaling=T
  .dfSetting$indScaling[runIndx]=T
  if(is.null(prefix)){
    prefix=""
  }
  
  if(sum(c(.includeHuman,.includeMouse))==2){
    prefix=paste0(prefix,"human-mouse")
  } else if(.includeMouse){
    prefix=paste0(prefix,"mouse")
  } else if(.includeHuman){
    prefix=paste0(prefix,"human")
  }
  
  if(.dfSetting$do.liger[runIndx]){
    prefix=paste0(prefix,"_liger")
  }
  
  if(Leng200){
    prefix=paste0(prefix,"_Leng200")
  }
  if(Jakel200){
    prefix=paste0(prefix,"_Jakel200")
  }
  
  if(is.null(saveDir)){
    .saveDir=.mySaveDirMaker(paste0("tmpBucket/results/",prefix),nPCs = .nPCs,cov=.covariates)
  } else {
    .saveDir=saveDir
  }
  
  .saveDirGlobal=paste0(.saveDir,"-Global/")
  
  argList=list(commonExpressed=.dfSetting$commonExpressed[runIndx],conservativeMapping=conservativeMapping,oldMapping=oldMapping,MGIsymbol=MGIsymbol,Leng200=Leng200,Jakel200=Jakel200,includeHuman=.includeHuman,includeMouse=.includeMouse,nPCs=nPCs,HVG_count=.dfSetting$HVG_count[runIndx],HVG_method="vst",covariates=.covariates,indScaling=.indScaling,saveDir=.saveDir,saveDirGlobal=.saveDirGlobal,ncores=ncores,excludeHighExp=.removeHighExp,UMI_cor_thr=.UMI_cor_thr,onlyRankBased=.dfSetting$onlyRankBased[runIndx],varScore.thr=.dfSetting$varScore.thr[runIndx],do.liger=.dfSetting$do.liger[runIndx],allGenesFraction=.dfSetting$allGenesFraction[runIndx],DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,prefix=prefix)
  
  .extraCreateDirFn(argList)
  
  .extraNewRun(argList=argList,newRun=newRun)
  argList$newRun=newRun
  
  return(argList)
  
}

.myLigerScaleNotCenter=function (object, remove.missing = TRUE, chunk = 1000, verbose = TRUE) {
  require(rliger)
  if (class(object@raw.data[[1]])[1] == "H5File") {
    hdf5_files = names(object@raw.data)
    vargenes = object@var.genes
    for (i in 1:length(hdf5_files)) {
      if (verbose) {
        message(hdf5_files[i])
      }
      chunk_size = chunk
      if (object@h5file.info[[i]][["format.type"]] == "AnnData") {
        genes = object@raw.data[[i]][["raw.var"]][]$index
      }
      else {
        genes = object@h5file.info[[i]][["genes"]][]
      }
      num_cells = object@h5file.info[[i]][["barcodes"]]$dims
      num_genes = length(genes)
      num_entries = object@h5file.info[[i]][["data"]]$dims
      prev_end_col = 1
      prev_end_data = 1
      prev_end_ind = 0
      gene_vars = rep(0, num_genes)
      gene_means = object@raw.data[[i]][["gene_means"]][1:num_genes]
      gene_sum_sq = object@raw.data[[i]][["gene_sum_sq"]][1:num_genes]
      gene_inds = which(genes %in% vargenes)
      gene_root_mean_sum_sq = sqrt(gene_sum_sq/(num_cells - 
                                                  1))
      safe_h5_create(object = object, idx = i, dataset_name = "scale.data", 
                     dims = c(length(vargenes), num_cells), mode = h5types$double, 
                     chunk_size = c(length(vargenes), chunk_size))
      num_chunks = ceiling(num_cells/chunk_size)
      if (verbose) {
        pb = txtProgressBar(0, num_chunks, style = 3)
      }
      ind = 0
      while (prev_end_col < num_cells) {
        ind = ind + 1
        if (num_cells - prev_end_col < chunk_size) {
          chunk_size = num_cells - prev_end_col + 1
        }
        start_inds = object@h5file.info[[i]][["indptr"]][prev_end_col:(prev_end_col + 
                                                                         chunk_size)]
        row_inds = object@h5file.info[[i]][["indices"]][(prev_end_ind + 
                                                           1):(tail(start_inds, 1))]
        counts = object@norm.data[[i]][(prev_end_ind + 
                                          1):(tail(start_inds, 1))]
        scaled = sparseMatrix(i = row_inds[1:length(counts)] + 
                                1, p = start_inds[1:(chunk_size + 1)] - prev_end_ind, 
                              x = counts, dims = c(num_genes, chunk_size))
        scaled = scaled[gene_inds, ]
        scaled = as.matrix(scaled)
        root_mean_sum_sq = gene_root_mean_sum_sq[gene_inds]
        scaled = sweep(scaled, 1, root_mean_sum_sq, "/")
        rownames(scaled) = genes[gene_inds]
        scaled = scaled[vargenes, ]
        scaled[is.na(scaled)] = 0
        scaled[scaled == Inf] = 0
        object@raw.data[[i]][["scale.data"]][1:length(vargenes), 
                                             prev_end_col:(prev_end_col + chunk_size - 1)] = scaled
        num_read = length(counts)
        prev_end_col = prev_end_col + chunk_size
        prev_end_data = prev_end_data + num_read
        prev_end_ind = tail(start_inds, 1)
        if (verbose) {
          setTxtProgressBar(pb, ind)
        }
      }
      object@scale.data[[i]] = object@raw.data[[i]][["scale.data"]]
      if (verbose) {
        setTxtProgressBar(pb, num_chunks)
        cat("\n")
      }
    }
    names(object@scale.data) <- names(object@raw.data)
  }
  else {
    
    if(F){
      for(i in 1:length(object@norm.data)){
        x=as.matrix(object@norm.data[[i]])[match(object@var.genes,row.names(object@norm.data[[i]])),]
        x[is.na(x)]=0
        x=as(x, "dgCMatrix")
        tst=rliger:::scaleNotCenterFast(t(x))
      }
    }
    
    object@scale.data <- lapply(1:length(object@norm.data), 
                                function(i) {
                                  x=as.matrix(object@norm.data[[i]])[match(object@var.genes,row.names(object@norm.data[[i]])),]
                                  x[is.na(x)]=0
                                  x=as(x, "dgCMatrix")
                                  rliger:::scaleNotCenterFast(t(x))
                                })
    object@scale.data <- lapply(object@scale.data, function(x) {
      as.matrix(x)
    })
    names(object@scale.data) <- names(object@norm.data)
    for (i in 1:length(object@scale.data)) {
      object@scale.data[[i]][is.na(object@scale.data[[i]])] <- 0
      rownames(object@scale.data[[i]]) <- colnames(object@raw.data[[i]])
      colnames(object@scale.data[[i]]) <- object@var.genes
    }
    if (remove.missing) {
      object <- removeMissingObs(object, slot.use = "scale.data", 
                                 use.cols = FALSE, verbose = verbose)
    }
  }
  return(object)
}

.myPCAfn_archive=function(data,argList,UMI_cor_thr){
  library(future)
  plan("multiprocess", workers = 5)
  plan()
  options(future.globals.maxSize = 1000 * 1024^4)
  
  myPrepDR=function (scaledData, features, verbose = TRUE) {
    
    data.use <- scaledData
    if (nrow(x = data.use) == 0) {
      stop("Data has not been scaled. Please run ScaleData and retry")
    }
    features.keep <- unique(x = features[features %in% rownames(x = data.use)])
    if (length(x = features.keep) < length(x = features)) {
      features.exclude <- setdiff(x = features, y = features.keep)
      if (verbose) {
        warning(paste0("The following ", length(x = features.exclude), 
                       " features requested have not been scaled (running reduction without them): ", 
                       paste0(features.exclude, collapse = ", ")))
      }
    }
    features <- features.keep
    features.var <- apply(X = data.use[features, ], MARGIN = 1, 
                          FUN = var)
    features.keep <- features[features.var > 0]
    if (length(x = features.keep) < length(x = features)) {
      features.exclude <- setdiff(x = features, y = features.keep)
      if (verbose) {
        warning(paste0("The following ", length(x = features.exclude), 
                       " features requested have zero variance (running reduction without them): ", 
                       paste0(features.exclude, collapse = ", ")))
      }
    }
    features <- features.keep
    features <- features[!is.na(x = features)]
    data.use <- data.use[features, ]
    return(data.use)
  }
  
  #reRunCheck=tryCatch({load(.myFilePathMakerFn("pca_anno",argList = argList));F}, error=function(e) {return(T)})
  
  if(T){
    mySeuratFn2=function(inputData,varFeatures){
      inputData@assays$RNA@var.features=varFeatures
      inputData = ScaleData(inputData, verbose = FALSE,features=varFeatures)
      return(inputData)
    }
    
    tmpInd=c()
    if(sum(!is.null(argList$covariates))>0&!is.null(data$data_m)){
      for(i in argList$covariates[!is.null(argList$covariates)]){
        tmpInd=c(tmpInd,which(is.na(data$data_m@meta.data[,i])))
      }
      if(length(tmpInd)>0){
        data$data_m=data$data_m[,-unique(tmpInd)]
      }
    }
    
    if(argList$indScaling){
      if(!is.null(data$data_m)){
        dataList=SplitObject(data$data_m, split.by = "anno_batch")
      } else {
        dataList=data$data
      }
      varFeatures=c()
      if(!is.null(argList$HVG_list)){
        for(iHVG in argList$HVG_list){
          varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))])
        }
      } else {
        varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(argList$HVG_count-1))])
      }
      
      varFeatures=unique(varFeatures)
      
      dataList=parallel::mclapply(dataList,mySeuratFn2,varFeatures=varFeatures,mc.cores = argList$ncores)
      gc()
      .mycBindFillFn=function(mat1,mat2){
        
        mat1c=setdiff(row.names(mat2),row.names(mat1))
        mat2c=setdiff(row.names(mat1),row.names(mat2))
        if(length(mat1c)>0){
          mat1cc=matrix(0,nrow=length(mat1c),ncol=ncol(mat1))
          row.names(mat1cc)=mat1c
          mat1=rbind(mat1,mat1cc)
        }
        if(length(mat2c)>0){
          mat2cc=matrix(0,nrow=length(mat2c),ncol=ncol(mat2))
          row.names(mat2cc)=mat2c
          mat2=rbind(mat2,mat2cc)
        }
        mat2=mat2[match(row.names(mat1),row.names(mat2)),]
        mat=cbind(mat1,mat2)
        return(mat)
      }
      
      resScaled=dataList[[1]]@assays$RNA@scale.data
      for(i in 2:length(dataList)){
        resScaled=.mycBindFillFn(mat1=resScaled,mat2=dataList[[i]]@assays$RNA@scale.data)
      }
      gc()
      
      if(!is.null(data$data_m)){
        resScaled=resScaled[,match(colnames(data$data_m),colnames(resScaled))]
        if(!all(colnames(data$data_m)==colnames(resScaled),na.rm = F)){
          stop("Error in the matching!")
        }
        if(!all(row.names(data$data_m)[row.names(data$data_m) %in% data$varFeatures]==row.names(resScaled))){
          stop("Error in the matching!")
        }
        data$data_m@assays$RNA@scale.data=resScaled
        if(!is.null(argList$HVG_list)){
          for(iHVG in argList$HVG_list){
            argList$HVG_count=iHVG
            tmp_varFeatures=data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))]
            data$data_m = RunPCA(data$data_m,features=tmp_varFeatures)
            pca_res=data$data_m@reductions$pca@cell.embeddings
            pd=data$data_m@meta.data
            save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
          }
        } else {
          tmp_varFeatures=data$varFeatures$Gene[which(data$varFeatures$Freq>(argList$HVG_count-1))]
          data$data_m = RunPCA(data$data_m,features=tmp_varFeatures)
          pca_res=data$data_m@reductions$pca@cell.embeddings
          pd=data$data_m@meta.data
          save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
        }
        
      } else {
        pd=dataList[[1]]@meta.data
        pd$sample=colnames(dataList[[1]])
        for(i in 2:length(dataList)){
          tmp=dataList[[i]]@meta.data
          tmp$sample=colnames(dataList[[i]])
          pd=rbind.fill(pd,tmp)
        }
        
        organism=rep("human",nrow(pd))
        organism[grepl("^mouse",tolower(pd$anno_batch))]="mouse"
        
        if(length(which(grepl("organism",colnames(pd))))>0){
          if(sum(is.na(pd$organism))>0){
            pd=pd[,-which(grepl("organism",colnames(pd)))]
            pd$organism=organism
          }
        } else {
          pd$organism=organism
        }
        
        
        pd_all=pd
        gc()
        if(!is.null(argList$HVG_list)){
          for(iHVG in argList$HVG_list){
            argList$HVG_count=iHVG
            tmp_varFeatures=data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))]
            
            #pca_res=.extraPCAfn(object=resScaled[which(row.names(resScaled) %in% tmp_varFeatures),], npcs = 30,findElbowPoint=F,seed.use = 42,reduction.key = "PC_",weight.by.var=T, approx = TRUE)
            #pca_res=pca_res$embeddings
            pca_res=Seurat:::RunPCA.default(object=myPrepDR(scaledData=resScaled, features=tmp_varFeatures, verbose = TRUE))
            pca_res=pca_res@cell.embeddings
            
            pd=pd_all[pd_all$sample %in% row.names(pca_res),]
            row.names(pd)=pd$sample
            pca_res=pca_res[match(row.names(pd),row.names(pca_res)),]
            save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
            gc()
          }
        } else {
          tmp_varFeatures=data$varFeatures$Gene[which(data$varFeatures$Freq>(argList$HVG_count-1))]
          
          #pca_res=.extraPCAfn(object=resScaled[which(row.names(resScaled) %in% tmp_varFeatures),], npcs = 30,findElbowPoint=F,seed.use = 42,reduction.key = "PC_",weight.by.var=T, approx = TRUE)
          #pca_res=pca_res$embeddings
          pca_res=Seurat:::RunPCA.default(object=myPrepDR(scaledData=resScaled, features=tmp_varFeatures, verbose = TRUE))
          pca_res=pca_res@cell.embeddings
          pd=dataList[[1]]@meta.data
          
          pd=pd_all[pd_all$sample %in% row.names(pca_res),]
          row.names(pd)=pd$sample
          pca_res=pca_res[match(row.names(pd),row.names(pca_res)),]
          save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
        }
      }
    }
    
    if(!argList$indScaling){
      
      varFeatures=c()
      if(!is.null(argList$HVG_list)){
        for(iHVG in argList$HVG_list){
          varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))])
        }
      } else {
        varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(argList$HVG_count-1))])
      }
      
      data$data_m <- ScaleData(data$data_m, verbose = FALSE,features=varFeatures, vars.to.regress =argList$covariates)
      
      if(!is.null(argList$HVG_list)){
        for(iHVG in argList$HVG_list){
          argList$HVG_count=iHVG
          tmp_varFeatures=data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))]
          data$data_m = RunPCA(data$data_m,features=tmp_varFeatures)
          pca_res=data$data_m@reductions$pca@cell.embeddings
          pd=data$data_m@meta.data
          save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
        }
      }
        
      
    }
  }
  
  return("done")
}

.myHarmonyFn=function(iHVG,res=NULL,runIndx=c(1,2,3,5),cluster_prior =NULL){
  .ArgList$HVG_count=iHVG
  load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F))
  
  organism=rep("Human",nrow(pd))
  organism[grepl("mouse",pd$anno_batch)]="Mouse"
  pd$organism=organism
  
  pd=pd[,c('anno_batch',"organism","ds_batch")]
  gc()
  harmony_gridSearch=function(indx,res=NULL,inputDf,pca_res,pd,cluster_prior =cluster_prior){
    if(is.null(res)){
      res=harmony::HarmonyMatrix(pca_res, pd, c('anno_batch',"organism"),theta=c(inputDf$theta1[indx],inputDf$theta2[indx]),lambda=c(inputDf$lambda1[indx],inputDf$lambda2[indx]), do_pca = FALSE, verbose=FALSE)
    }
    return(res)
  }
  
  dfSettings=NULL
  dfSettings=rbind(dfSettings,data.frame(theta1=2,theta2=2,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=4,theta2=2,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=6,theta2=2,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=8,theta2=2,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=2,theta2=4,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=2,theta2=6,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=2,theta2=8,lambda1=1,lambda2=1))
  
  #indx=1;inputDf=dfSettings;pca_res=pca_res[,1:30]
  res=parallel::mclapply(runIndx,res=res,harmony_gridSearch,inputDf=dfSettings,pca_res=pca_res[,1:.ArgList$nPCs],pd=pd,cluster_prior=cluster_prior,mc.cores = 7)
  return(res)
}

.myHarmonyFn_modified=function(iHVG,res=NULL,runIndx=c(1,2,3,5)){
  .ArgList$HVG_count=iHVG
  load(gsub("_pseudocell200","",.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F)))
  refined_anno=qread("~/myBucket/torm_pd_ExN_refined.qs")
  pd=pd[pd$sample %in% refined_anno$sample,]
  pca_res=pca_res[match(pd$sample,row.names(pca_res)),]
  
  if(sum(colnames(pd)=="organism")==0|sum(is.na(pd$organism))>0){
    organism=rep("Human",nrow(pd))
    organism[grepl("mouse",pd$anno_batch)]="Mouse"
    pd$organism=organism
  }
  
  pd=pd[,c('anno_batch',"organism","ds_batch")]
  gc()
  harmony_gridSearch=function(indx,res=NULL,inputDf,pca_res,pd){
    if(is.null(res)){
      res=harmony::HarmonyMatrix(pca_res, pd, c('ds_batch',"organism"),theta=c(inputDf$theta1[indx],inputDf$theta2[indx]),lambda=c(inputDf$lambda1[indx],inputDf$lambda2[indx]), do_pca = FALSE, verbose=FALSE)
    }
    return(res)
  }
  
  dfSettings=NULL
  dfSettings=rbind(dfSettings,data.frame(theta1=2,theta2=2,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=4,theta2=2,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=6,theta2=2,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=8,theta2=2,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=2,theta2=4,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=2,theta2=6,lambda1=1,lambda2=1))
  dfSettings=rbind(dfSettings,data.frame(theta1=2,theta2=8,lambda1=1,lambda2=1))
  
  #indx=1;inputDf=dfSettings;pca_res=pca_res[,1:30]
  res=parallel::mclapply(runIndx,res=res,harmony_gridSearch,inputDf=dfSettings,pca_res=pca_res[,1:.ArgList$nPCs],pd=pd,mc.cores = 7)
  return(res)
}


.myLigerRunFn=function(iHVG=0,pd,argList,varFeatures,inputLiger){
  if(iHVG!=0){
    argList$HVG_count=iHVG
    tmp_varFeatures=varFeatures$Gene[which(varFeatures$Freq>(iHVG-1))]
    inputLiger@var.genes=unique(as.character(tmp_varFeatures))
  } else {
    argList$HVG_count=0
  }
  
  if(!file.exists(.myFilePathMakerFn("pca_anno",argList = argList,pseudoImportant = F))){
    
    
    #object=inputLiger;remove.missing=T; chunk = 1000; verbose = TRUE
    inputLiger <- .myLigerScaleNotCenter(inputLiger,remove.missing=T, chunk = 1000, verbose = TRUE)
    inputLiger <- optimizeALS(inputLiger,k = argList$nPCs)
    inputLiger <- quantile_norm(inputLiger)
    
    harmony_embeddings = inputLiger@H.norm
    .mySaveFn(harmony_embeddings,file=.myFilePathMakerFn("harmony-embeddings",argList=argList))
    
    pca_res=do.call("rbind",inputLiger@H)
    save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList = argList,pseudoImportant = F))
    
  }
  
}

.AnnoFn_graph=function(i,res,.ArgList,newRun,cluster_resolution=0.6,offsetCount,changeSaveDir=T,Liger=F,group.singletons=T,dsLevel=F,save_dir=NULL){
  internalFn=function(j,res,argList,pd,newRun=F,offsetCount=0,changeSaveDir=T,Liger=F,group.singletons=T,dsLevel=F,cluster_resolution,save_dir){
    require(future)
    plan("multiprocess", workers = 4)
    plan()
    options(future.globals.maxSize = 1000 * 1024^4)
    
    if(sum(names(argList)=="Leng200")==0){
      argList$Leng200=F
    }
    if(sum(names(argList)=="Jakel200")==0){
      argList$Jakel200=F
    }
    if(sum(names(argList)=="projection")==0){
      argList$projection=F
    }
    
    if(!is.null(save_dir)){
      argList$saveDir=save_dir
      
      
      if(!dir.exists(argList$saveDir)){
        dir.create(argList$saveDir,recursive = T)
      }
    } else if(changeSaveDir){
      .saveDir=paste0("~/data/",argList$cellClass)
      if(.ArgList$includeHuman&.ArgList$includeMouse){
        .saveDir=paste0(.saveDir,"_HumanMouse")
      } else if(.ArgList$includeHuman) {
        .saveDir=paste0(.saveDir,"_Human")
      } else if(.ArgList$includeMouse){
        .saveDir=paste0(.saveDir,"_Mouse")
      } else {
        .saveDir=paste0(.saveDir,"_unknown")
      }
      if(!argList$indScaling){
        .saveDir=paste0(.saveDir,"_notIndScaling")
      } else {
        .saveDir=paste0(.saveDir,"_indScaling")
      }
      if(Liger){
        argList$saveDir=paste0(.saveDir,"_liger")
      } else {
        argList$saveDir=.saveDir
      }
      
      if(dsLevel){
        .saveDir=paste0(.saveDir,"_dsLevel")
      }
      if(argList$projection){
        .saveDir=paste0(.saveDir,"_projection")
      }
      
      if(argList$nPCs!=30){
        .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
      }
      
      
      if(argList$Leng200){
        .saveDir=paste0(.saveDir,"_leng200")
      }
      if(argList$Jakel200){
        .saveDir=paste0(.saveDir,"_Jakel200")
      }
      
      argList$saveDir=.saveDir
      
      
      if(!dir.exists(argList$saveDir)){
        dir.create(argList$saveDir,recursive = T)
      }
    }
    
    #print(argList$saveDir)
    
    reRun=T
    if(!newRun){
      if(file.exists(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))){
        reRunCheck=tryCatch({load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F));F}, error=function(e) {return(T)})
        if(!reRunCheck){
          if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2){
            if(sum(is.na(pd$UMAP_1))==0&sum(is.na(pd$UMAP_2))==0){
              reRun=F
            }
          }
        }
        
      }
    }
    
    
    if(reRun|newRun){
      
      harmony_embeddings=res[[j]]
      if(sum(colnames(pd)=="sample")>0){
        if(sum(is.na(pd$sample))==0){
          pd=pd[pd$sample %in% row.names(harmony_embeddings),]
          harmony_embeddings=harmony_embeddings[match(pd$sample,row.names(harmony_embeddings)),]
        } else {
          pd=pd[row.names(pd) %in% row.names(harmony_embeddings),]
          harmony_embeddings=harmony_embeddings[match(row.names(pd),row.names(harmony_embeddings)),]
        }
        
      } else {
        pd=pd[row.names(pd) %in% row.names(harmony_embeddings),]
        harmony_embeddings=harmony_embeddings[match(row.names(pd),row.names(harmony_embeddings)),]
      }
      
      
      #harmony_net=Seurat:::FindNeighbors(object=harmony_embeddings[,-5])
      harmony_net=Seurat:::FindNeighbors(object=harmony_embeddings)
      clusters=.netFindClusters(inputGraph=harmony_net[["snn"]], algorithm = 1,resolution = cluster_resolution,group.singletons = group.singletons,modularity.fxn = 1)
      
      
      tst=.reductionUMAPFn(harmony_embeddings,umap.method='uwot')
      resUMAP=tst$embedding
      
      pd$UMAP_1=resUMAP[,1]
      pd$UMAP_2=resUMAP[,2]
      if(sum(is.na(pd$organism))>0){
        organism=rep("human",nrow(pd))
        organism[grepl("^mouse_",pd$anno_batch)]="mouse"
        pd$organism=organism
      }
      
      colnames(clusters)=paste0("anno_cluster_",colnames(clusters))
      pd=cbind(pd,clusters)
      save(pd,file=.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))
      
      p=.myDimPlotFn(object=pd, dimCols = c("UMAP_1", "UMAP_2"), attCol = 'anno_cluster_res', 
                     pt.size = NULL, label = TRUE, label.size = 4, 
                     repel = FALSE, cells.highlight = NULL, cols.highlight = "#DE2D26", 
                     sizes.highlight = 1, na.value = "grey50", ncol = NULL, combine = TRUE)
      ggsave(plot = p,.myFilePathMakerFn(paste0("hm_harmony",j+offsetCount,"_clusters"),argList=argList,pdf=T),width = 12,height = 8)
      
      p1=ggplot(pd,aes(anno_cluster_res,QC_MT.pct,fill=anno_cluster_res))+geom_violin()+theme_classic()+theme(legend.position = "none")
      p2=ggplot(pd,aes(anno_cluster_res,QC_IEG.pct,fill=anno_cluster_res))+geom_violin()+theme_classic()+theme(legend.position = "none")
      p3=ggplot(pd,aes(anno_cluster_res,QC_top50_pct,fill=anno_cluster_res))+geom_violin()+theme_classic()+theme(legend.position = "none")
      p=p1+p2+p3+ plot_layout(nrow = 3, byrow = FALSE)
      ggsave(plot = p,.myFilePathMakerFn(paste0("hm_harmony",j+offsetCount,"_clusterQC"),argList=argList,pdf=T),width = 10,height = 8)
      
      
    }
    
    return(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))
  }
  
  .ArgList$HVG_count=as.numeric(gsub("HVG_","",names(res)[i]))
  load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F))
  
  pd$ds_batch[which(grepl("human_Tushar_iNPH",pd$ds_batch))]="human_Tushar_iNPH"
  
  res2=mclapply(1:length(res[[i]]),internalFn,res=res[[i]],argList=.ArgList,pd=pd,newRun=newRun,offsetCount=offsetCount,cluster_resolution=cluster_resolution,changeSaveDir=changeSaveDir,Liger=Liger,dsLevel=dsLevel,group.singletons=group.singletons,mc.cores = length(res[[i]]),save_dir=save_dir)
  return(res2)
}

#offsetCount=0;doStriatum=F;changeSaveDir=T;Liger=F
.AnnoFn=function(j,res,data_m,argList,ncores,offsetCount=0,doStriatum=F,changeSaveDir=T,markers_path=NULL,Liger=F,dsLevel=F){
  require(ggplot2)
  require(Seurat)
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  if(sum(names(argList)=="projection")==0){
    argList$projection=F
  }
  
  if(changeSaveDir){
    .saveDir=paste0("~/data/",argList$cellClass)
    if(.ArgList$includeHuman&.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_HumanMouse")
    } else if(.ArgList$includeHuman) {
      .saveDir=paste0(.saveDir,"_Human")
    } else if(.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_Mouse")
    } else {
      .saveDir=paste0(.saveDir,"_unknown")
    }
    if(!argList$indScaling){
      .saveDir=paste0(.saveDir,"_notIndScaling")
    } else {
      .saveDir=paste0(.saveDir,"_indScaling")
    }
    if(Liger){
      argList$saveDir=paste0(.saveDir,"_liger")
    } else {
      argList$saveDir=.saveDir
    }
    
    if(dsLevel){
      .saveDir=paste0(.saveDir,"_dsLevel")
    }
    if(argList$projection){
      .saveDir=paste0(.saveDir,"_projection")
    }
    
    if(argList$nPCs!=30){
      .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
    }
    
    
    if(argList$Leng200){
      .saveDir=paste0(.saveDir,"_leng200")
    }
    if(argList$Jakel200){
      .saveDir=paste0(.saveDir,"_Jakel200")
    }
    argList$saveDir=.saveDir
    
    
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  
  
  reRun=T
  if(T){
    if(file.exists(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))){
      reRunCheck=tryCatch({load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F));F}, error=function(e) {return(T)})
      if(!reRunCheck){
        if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2){
          if(sum(is.na(pd$UMAP_1))==0&sum(is.na(pd$UMAP_2))==0){
            reRun=F
          }
        }
      }
      
    }
  }
  
  
  if(reRun){
    stop("First step needs to be run")
  } else {
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))
  }
  
  if(sum(colnames(pd)=="ds_batch")==0){
    pd$ds_batch=pd$anno_batch
  }
  
  
  #row.names(pd)=pd$sample
  
  if(class(data_m)==class(list())){
    tmp_snames=unlist(lapply(data_m,colnames))
  } else {
    tmp_snames=colnames(data_m)
  }
  
  if(sum(tmp_snames %in% pd$sample)>sum(tmp_snames %in% row.names(pd))){
    row.names(pd)=pd$sample
  }
  
  if(length(setdiff(row.names(pd),tmp_snames))>0|length(setdiff(tmp_snames,row.names(pd)))>0){
    warning("Inconsistency between pd and expression data")
    pd=pd[row.names(pd) %in% tmp_snames,]
    if(class(data_m)==class(list())){
      tormlist=c()
      for(i in 1:length(data_m)){
        if(sum(colnames(data_m[[i]]) %in% row.names(pd))>0){
          data_m[[i]]=data_m[[i]][,colnames(data_m[[i]]) %in% row.names(pd)]
        } else {
          tormlist=c(tormlist,i)
        }
        
      }
      if(length(tormlist)>0){
        data_m=data_m[-tormlist]
      }
    } else {
      data_m=data_m[,colnames(data_m) %in% row.names(pd)]
    }
  }
  
 
  
  cat("     Initial QC figures...\n")
  if(sum(colnames(pd)=="organism")>0){
    if(length(unique(pd$organism))>1){
      p=.myDimPlotFn(object=pd, dimCols = c("UMAP_1", "UMAP_2"),attCol = "organism")
      ggsave(plot=p,.myFilePathMakerFn(paste0("hm_harmony",j+offsetCount),argList=argList,pdf=T),width = 15,height = 10)
    }
    
  }
  
  
  
  pd_Tushar=pd[which(pd$ds_batch=="human_Tushar_iNPH"),]
  if(nrow(pd_Tushar)>0){
    if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony_Tushar",j+offsetCount,"2"),pdf=T),argList)))){
      p=.my2dPlot2(inputPCA = pd_Tushar,batchCol = "anno_cellState",reductionCols = c('UMAP_1','UMAP_2'))
      ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony_Tushar",j+offsetCount,"2"),pdf=T),argList)),width = 12,height = 8)
      
    }
    
  }
  
  
  p=.my2dPlot2(inputPCA = pd,batchCol = "ds_batch",reductionCols = c('UMAP_1','UMAP_2'))
  ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"2"),pdf=T),argList)),width = 12,height = 8)
  
  if(sum(colnames(pd)=="organism")>0){
    p=.my2dPlot2(inputPCA = pd,batchCol = "organism",reductionCols = c('UMAP_1','UMAP_2'))
    ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"6"),pdf=T),argList)),width = 12,height = 8)
  }
  
  
  #p=.my2dPlot2(inputPCA = tmp,batchCol = "anno_cellState",reductionCols = c('UMAP_1','UMAP_2'))
  #ggsave(plot=p,.myFilePathMakerFn(paste0("hm_harmony",j,"3"),argList=argList,pdf=T),width = 12,height = 8)
  
  #p=.my2dPlot2(inputPCA = tmp[tmp$organism=="mouse",],batchCol = "anno_cellState",reductionCols = c('UMAP_1','UMAP_2'))
  #ggsave(plot=p,.myFilePathMakerFn(paste0("hm_harmony",j,"3mouse"),argList=argList,pdf=T),width = 12,height = 8)
  
  #p=.my2dPlot2(inputPCA = tmp[tmp$organism=="human",],batchCol = "anno_cellState",reductionCols = c('UMAP_1','UMAP_2'))
  #ggsave(plot=p,.myFilePathMakerFn(paste0("hm_harmony",j,"3human"),argList=argList,pdf=T),width = 12,height = 8)
  
  if(is.null(markers_path)){
    markerGenes=read.table("~/myBucket/markergenes/all_markers.txt",sep="\t",header = T,stringsAsFactors = F)
    
    
    markerGenes=rbind(markerGenes,data.frame(Gene=c("CD200R1","BRIP1","TOP2A","FTH1","LYVE1","GPR126","HCG22","PAQR5","PLAT","MGAM","CECR2","CD300E","TNFAIP3","IFIT2"),pattern=NA,source="Tushar",MainMarker="Yes",page="page2",stringsAsFactors = F))
    
  } else {
    markerGenes=read.table(markers_path,sep="\t",header = T,stringsAsFactors = F)
    
  }
  #markerGenes=data.frame(Gene=toupper(c("EP2","EP4","EP3","EP1","PTGS1","PBR","GYS1","egr3","hbegf","nr4a1","slc2a3","b4galt1","map3k8","zfp36")),pattern=NA,source="Tushar",MainMarker="Yes",page="page2",stringsAsFactors = F)
  
  markerGenes=markerGenes[!duplicated(markerGenes$Gene),]
  markerGenes_all=markerGenes
  
  if(F){
    all_gene_list=markerGenes$Gene
    if(nrow(pd_Tushar)>0){
      Tushar_markers=unique(pd_Tushar$anno_cellState)
      Tushar_markers=unique(unlist(strsplit(Tushar_markers,"_")))
      all_gene_list=c(all_gene_list,Tushar_markers)
    }
    
    data_m=list()
    for(ik in 1:length(data)){
      tmp=data[[ik]]
      tmp=tmp[which(toupper(rowData(tmp)$gene_short_name) %in% toupper(all_gene_list)),]
      tmp_count=apply(as.matrix(counts(tmp)),2,sum)
      tmp=tmp[,which(tmp_count>0)]
      data_m=c(data_m,tmp)
      names(data_m)[length(data_m)]=names(data)[ik]
    }
    data_m=.mycBindFn(data_m)
    data_m=.extraExport2SeuratFn(data_m)
  }
  
  gc()
  cat("     Marker figures...\n")
  for(ik in unique(markerGenes_all$page)){
    markerGenes=markerGenes_all[which(markerGenes_all$page==ik),]
    tmp=data_m
    names(tmp)=unlist(lapply(tmp,function(x) as.character(x$anno_batch[1])))
    #tmp=tmp[,colnames(tmp) %in% row.names(pd)]
    pd2=pd#[match(colnames(tmp),row.names(pd)),]
    #p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = markerGenes$Gene)
    #ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j,"4",ik),pdf=T),argList)),width = 49,height = 40)
    
    
    #inputPCA=pd2;batch_values="anno_batch";reductionCols=c("UMAP_1","UMAP_2");geneNameList=lapply(markerGenes$Gene,function(x) x);geneNameCol="gene_short_name";expData=tmp;ncolumns=6;exponetial_pwr=8;combine_figs=F;ncores=1
    
    
    #p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply("C1QA",function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6,ncores=ncores)
    #inputPCA=pd2;batch_values="anno_batch";reductionCols=c("UMAP_1","UMAP_2");geneNameList=lapply(markerGenes$Gene,function(x) x);geneNameCol="gene_short_name";expData=tmp;ncolumns=6
    
    x=unname(unlist(lapply(tmp,function(x) colnames(x))))
    if(length(setdiff(pd2$sample,x))==0){
      row.names(pd2)=pd2$sample
    }
    
    if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5",ik),pdf=T),argList)))){
      p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(markerGenes$Gene,function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6,ncores=ncores)
      #ggsave(plot=p,file="~/myBucket/torm1.png",width = 49,height = 40)
      ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5",ik),pdf=T),argList)),width = 49,height = 40)
      
    }
    
    pd_mouse=pd2[which(tolower(pd2$organism)=="mouse"),]
    if(nrow(pd_mouse)>0){
      tmp_mouse=list()
      for(ij in 1:length(tmp)){
        tmp2=tmp[[ij]]
        if(sum(colnames(tmp2) %in% row.names(pd_mouse))>0){
          tmp2=tmp2[,colnames(tmp2) %in% row.names(pd_mouse)]
          tmp_mouse=c(tmp_mouse,tmp2)
        }
      }
      names(tmp_mouse)=unlist(lapply(tmp_mouse,function(x) as.character(x$anno_batch[1])))
      
      if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5mouse",ik),pdf=T),argList)))){
        p=.my2dPlot_counts(inputPCA=pd_mouse,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(markerGenes$Gene,function(x) x),geneNameCol="gene_short_name",expData=tmp_mouse,ncolumns=6,ncores=ncores)
        ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5mouse",ik),pdf=T),argList)),width = 49,height = 40)
        
      }
      
    }
    
    
    pd_human=pd2[which(tolower(pd2$organism)=="human"),]
    if(nrow(pd_human)>0){
      tmp_human=list()
      for(ij in 1:length(tmp)){
        tmp2=tmp[[ij]]
        if(sum(colnames(tmp2) %in% row.names(pd_human))>0){
          tmp2=tmp2[,colnames(tmp2) %in% row.names(pd_human)]
          tmp_human=c(tmp_human,tmp2)
        }
        
      }
      names(tmp_human)=unlist(lapply(tmp_human,function(x) as.character(x$anno_batch[1])))
      
      if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5human",ik),pdf=T),argList)))){
        p=.my2dPlot_counts(inputPCA=pd_human,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(markerGenes$Gene,function(x) x),geneNameCol="gene_short_name",expData=tmp_human,ncolumns=6,ncores=ncores)
        ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5human",ik),pdf=T),argList)),width = 49,height = 40)
        
      }
      
    }
    
  }
  cat("     Tushar's marker figures...\n")
  if(nrow(pd_Tushar)>0){
    Tushar_markers=unique(pd_Tushar$anno_cellState)
    Tushar_markers=unique(unlist(strsplit(Tushar_markers,"_")))
    
    tmp=data_m
    #tmp=tmp[,colnames(tmp) %in% row.names(pd)]
    pd2=pd#[match(colnames(tmp),row.names(pd)),]
    #p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = Tushar_markers)
    #ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony_",j,"4_Tushar"),pdf=T),argList)),width = 49,height = 40)
    
    
    if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5","_Tushar"),pdf=T),argList)))){
      p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(Tushar_markers,function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6,ncores=ncores)
      ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5","_Tushar"),pdf=T),argList)),width = 6.7*min(length(Tushar_markers),6),height = 40/7*min(ceiling(length(Tushar_markers)/6),7))
      
    }
    
    
    pd_mouse=pd2[which(tolower(pd2$organism)=="mouse"),]
    if(nrow(pd_mouse)>0){
      tmp_mouse=list()
      for(ij in 1:length(tmp)){
        tmp2=tmp[[ij]]
        if(sum(colnames(tmp2) %in% row.names(pd_mouse))>0){
          tmp2=tmp2[,colnames(tmp2) %in% row.names(pd_mouse)]
          tmp_mouse=c(tmp_mouse,tmp2)
        }
      }
      names(tmp_mouse)=unlist(lapply(tmp_mouse,function(x) as.character(x$anno_batch[1])))
      
      if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5mouse","_Tushar"),pdf=T),argList)))){
        p=.my2dPlot_counts(inputPCA=pd_mouse,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(Tushar_markers,function(x) x),geneNameCol="gene_short_name",expData=tmp_mouse,ncolumns=6,ncores=ncores)
        ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5mouse","_Tushar"),pdf=T),argList)),width = 7*min(length(Tushar_markers),6),height = 40/7*min(ceiling(length(Tushar_markers)/6),7))
        
      }
      
    }
    pd_human=pd2[which(tolower(pd2$organism)=="human"),]
    if(nrow(pd_human)>0){
      tmp_human=list()
      for(ij in 1:length(tmp)){
        tmp2=tmp[[ij]]
        if(sum(colnames(tmp2) %in% row.names(pd_human))>0){
          tmp2=tmp2[,colnames(tmp2) %in% row.names(pd_human)]
          tmp_human=c(tmp_human,tmp2)
        }
      }
      names(tmp_human)=unlist(lapply(tmp_human,function(x) as.character(x$anno_batch[1])))
      
      if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5human","_Tushar"),pdf=T),argList)))){
        p=.my2dPlot_counts(inputPCA=pd_human,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(Tushar_markers,function(x) x),geneNameCol="gene_short_name",expData=tmp_human,ncolumns=6,ncores=ncores)
        ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5human","_Tushar"),pdf=T),argList)),width = 7*min(length(Tushar_markers),6),height = 40/7*min(ceiling(length(Tushar_markers)/6),7))
      }
      
    }
    
  }
  
  
  
  if(sum(pd$ds_batch=="OHDA")>1000|doStriatum){
    cat("     Striatum marker figures...\n")
    #markerGenes=read.table("~/myBucket/sl_hm_InN",sep="\t",header = F,stringsAsFactors = F)
    markerGenes=read.table("~/myBucket/markergenes/SPN_markers_set2",sep="\t",header = F,stringsAsFactors = F)
    markerGenes=markerGenes[,1]
    markerGenes=markerGenes[!duplicated(markerGenes)]
    
    #markerGenes=toupper("Ptgs2")
    
    tmp=data_m
    #tmp=tmp[,colnames(tmp) %in% row.names(pd)]
    pd2=pd#[which(as.character(pd$anno_cluster_res) !="23"),]#[match(colnames(tmp),row.names(pd)),]
    #p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = Tushar_markers)
    #ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony_",j,"4_Tushar"),pdf=T),argList)),width = 49,height = 40)
    
    p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(markerGenes,function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6,ncores=ncores)
    ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5","_Striatum"),pdf=T),argList)),width = 6.7*min(length(markerGenes),6),height = 40/7*max(ceiling(length(markerGenes)/6),7))
    
    
    pd_mouse=pd2[which(tolower(pd2$organism)=="mouse"),]
    if(nrow(pd_mouse)>0){
      tmp_mouse=list()
      for(ij in 1:length(tmp)){
        tmp2=tmp[[ij]]
        if(sum(colnames(tmp2) %in% row.names(pd_mouse))>0){
          tmp2=tmp2[,colnames(tmp2) %in% row.names(pd_mouse)]
          tmp_mouse=c(tmp_mouse,tmp2)
        }
      }
      names(tmp_mouse)=unlist(lapply(tmp_mouse,function(x) as.character(x$anno_batch[1])))
      
      p=.my2dPlot_counts(inputPCA=pd_mouse,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(markerGenes,function(x) x),geneNameCol="gene_short_name",expData=tmp_mouse,ncolumns=6,ncores=ncores)
      ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5mouse","_Striatum"),pdf=T),argList)),width = 7*min(length(markerGenes),6),height = 40/7*max(ceiling(length(markerGenes)/6),7))
      
    }
    
    pd_human=pd2[which(tolower(pd2$organism)=="human"),]
    if(nrow(pd_human)>0){
      tmp_human=list()
      for(ij in 1:length(tmp)){
        tmp2=tmp[[ij]]
        if(sum(colnames(tmp2) %in% row.names(pd_human))>0){
          tmp2=tmp2[,colnames(tmp2) %in% row.names(pd_human)]
          tmp_human=c(tmp_human,tmp2)
        }
      }
      names(tmp_human)=unlist(lapply(tmp_human,function(x) as.character(x$anno_batch[1])))
      
      p=.my2dPlot_counts(inputPCA=pd_human,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(markerGenes,function(x) x),geneNameCol="gene_short_name",expData=tmp_human,ncolumns=6,ncores=ncores)
      ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j+offsetCount,"5human","_Striatum"),pdf=T),argList)),width = 7*min(length(markerGenes),6),height = 40/7*max(ceiling(length(markerGenes)/6),7))
    }
    
  }
  
  return("Done!")
}

.AnnoFn_Mouse=function(j,res,argList,ncores,offsetCount=0,changeSaveDir=T,Liger=F,dsLevel=F){
  require(ggplot2)
  require(Seurat)
  
  if(argList$includeHuman){
    stop("This function is only for the mouse datasets!")
  }
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  
  
  if(sum(names(argList)=="projection")==0){
    argList$projection=F
  }
  
  if(changeSaveDir){
    .saveDir=paste0("~/data/",argList$cellClass)
    if(.ArgList$includeHuman&.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_HumanMouse")
    } else if(.ArgList$includeHuman) {
      .saveDir=paste0(.saveDir,"_Human")
    } else if(.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_Mouse")
    } else {
      .saveDir=paste0(.saveDir,"_unknown")
    }
    if(!argList$indScaling){
      .saveDir=paste0(.saveDir,"_notIndScaling")
    } else {
      .saveDir=paste0(.saveDir,"_indScaling")
    }
    if(Liger){
      argList$saveDir=paste0(.saveDir,"_liger")
    } else {
      argList$saveDir=.saveDir
    }
    
    if(dsLevel){
      .saveDir=paste0(.saveDir,"_dsLevel")
    }
    if(argList$projection){
      .saveDir=paste0(.saveDir,"_projection")
    }
    
    if(argList$nPCs!=30){
      .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
    }
    
    
    if(argList$Leng200){
      .saveDir=paste0(.saveDir,"_leng200")
    }
    if(argList$Jakel200){
      .saveDir=paste0(.saveDir,"_Jakel200")
    }
    argList$saveDir=.saveDir
    
    
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  
  
  
  reRun=T
  if(T){
    if(file.exists(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))){
      reRunCheck=tryCatch({load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F));F}, error=function(e) {return(T)})
      if(!reRunCheck){
        if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2){
          if(sum(is.na(pd$UMAP_1))==0&sum(is.na(pd$UMAP_2))==0){
            reRun=F
          }
        }
      }
      
    }
  }
  
  
  if(reRun){
    stop("First step needs to be run")
  } else {
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))
  }
  
  
  p=.my2dPlot2(inputPCA = pd,batchCol = "anno_cellState",reductionCols = c('UMAP_1','UMAP_2'))
  ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony_mouseAnno",j+offsetCount,"2"),pdf=T),argList)),width = 12,height = 8)
  
  
  
  return("Done!")
}

.AnnoFn_Alec=function(j,res,data_m,pd,argList,ncores,offsetCount=0,changeSaveDir=T,markers_path=NULL,Liger=F,dsLevel=F){
  require(ggplot2)
  require(Seurat)
  
  if(changeSaveDir){
    if(Liger){
      argList$saveDir=paste0("~/data/",argList$cellClass,"_HumanMouse_indScaling_liger")
    } else {
      argList$saveDir=paste0("~/data/",argList$cellClass,"_HumanMouse_indScaling")
    }
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  
  reRun=T
  if(T){
    if(file.exists(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))){
      reRunCheck=tryCatch({load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F));F}, error=function(e) {return(T)})
      if(!reRunCheck){
        if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2){
          if(sum(is.na(pd$UMAP_1))==0&sum(is.na(pd$UMAP_2))==0){
            reRun=F
          }
        }
      }
      
    }
  }
  
  
  if(reRun){
    stop("First step needs to be run")
  } else {
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))
  }
  
  if(sum(colnames(pd)=="ds_batch")==0){
    pd$ds_batch=pd$anno_batch
  }
  
  cat("     timepoint QC figures...\n")
  timePoints=pd$ds_batch
  if(length(timePoints)>0){
    timePoints=unlist(lapply(strsplit(timePoints,"_"),function(x) x[length(x)]))
    pd$timePoint=timePoints
    p=.myDimPlotFn(object=pd, dimCols = c("UMAP_1", "UMAP_2"),attCol = "timePoint")
    ggsave(plot=p,.myFilePathMakerFn(paste0("hm_harmony_timePoint",j+offsetCount),argList=argList,pdf=T),width = 15,height = 10)
    
    timePoints=unlist(lapply(strsplit(timePoints,"_"),function(x) x[length(x)]))
    pd$timePoint=timePoints
    p=.myDimPlotFn(object=pd, dimCols = c("UMAP_1", "UMAP_2"),attCol = "timePoint")+facet_wrap(~timePoint)
    ggsave(plot=p,.myFilePathMakerFn(paste0("hm_harmony_timePoint",j+offsetCount,"_v2"),argList=argList,pdf=T),width = 15,height = 10)
  }
  
  cat("     Alec immune clusters figures...\n")
  pd=pd[!is.na(pd$seurat_clusters),]
  immune_data=qread("~/myBucket/Alec_BCell_Unfiltered/P14_arranged.qs")
  non_immune_data=qread("~/myBucket/Alec_BCell_Unfiltered/P14_CB_all_CD45neg_arranged.qs")
  pd$seurat_clusters=as.character(pd$seurat_clusters)
  pd$seurat_clusters[pd$sample %in% colnames(immune_data)]=paste0("Immune_",pd$seurat_clusters)[pd$sample %in% colnames(immune_data)]
  pd$seurat_clusters[pd$sample %in% colnames(non_immune_data)]=paste0("CD45neg_",pd$seurat_clusters)[pd$sample %in% colnames(non_immune_data)]
  p=.myDimPlotFn(object=pd, dimCols = c("UMAP_1", "UMAP_2"), attCol = "seurat_clusters")+theme(legend.position = "none")
  ggsave(plot=p,.myFilePathMakerFn(paste0("hm_harmony_alec",j+offsetCount),argList=argList,pdf=T),width = 15,height = 10)
  
  
  return("Done!")
}

.AnnoFn2=function(j,data_m,argList,ncores,offsetCount=0,updateBICCN=T,writeMarkers=T,changeSaveDir=T,Liger=F,dsLevel=F,...){
  
  j=j+offsetCount
  reRun=T
  
  updatelbl=""
  if(updateBICCN){
    updatelbl="updatedBICCNlbl"
  }
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  
  if(sum(names(argList)=="projection")==0){
    argList$projection=F
  }
  
  if(changeSaveDir){
    .saveDir=paste0("~/data/",argList$cellClass)
    if(.ArgList$includeHuman&.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_HumanMouse")
    } else if(.ArgList$includeHuman) {
      .saveDir=paste0(.saveDir,"_Human")
    } else if(.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_Mouse")
    } else {
      .saveDir=paste0(.saveDir,"_unknown")
    }
    if(!argList$indScaling){
      .saveDir=paste0(.saveDir,"_notIndScaling")
    } else {
      .saveDir=paste0(.saveDir,"_indScaling")
    }
    if(Liger){
      argList$saveDir=paste0(.saveDir,"_liger")
    } else {
      argList$saveDir=.saveDir
    }
    
    if(dsLevel){
      .saveDir=paste0(.saveDir,"_dsLevel")
    }
    if(argList$projection){
      .saveDir=paste0(.saveDir,"_projection")
    }
    
    if(argList$nPCs!=30){
      .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
    }
    
    
    
    if(argList$Leng200){
      .saveDir=paste0(.saveDir,"_leng200")
    }
    if(argList$Jakel200){
      .saveDir=paste0(.saveDir,"_Jakel200")
    }
    argList$saveDir=.saveDir
    
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  
  
  
  
  if(T){
    if(file.exists(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F))){
      reRunCheck=tryCatch({load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F));F}, error=function(e) {return(T)})
      if(!reRunCheck){
        if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2){
          if(sum(is.na(pd$UMAP_1))==0&sum(is.na(pd$UMAP_2))==0){
            reRun=F
          }
        }
      }
      
    }
  }
  
  
  if(reRun){
    stop("First steps needs to be run")
  } else {
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F))
    
    if(class(data_m)==class(list())){
      tmp_snames=unlist(lapply(data_m,colnames))
    } else {
      tmp_snames=colnames(data_m)
    }
    
    if(sum(tmp_snames %in% pd$sample)>sum(tmp_snames %in% row.names(pd))){
      row.names(pd)=pd$sample
    }
    
    if(length(setdiff(row.names(pd),tmp_snames))>0|length(setdiff(tmp_snames,row.names(pd)))>0){
      warning("Inconsistency between pd and expression data")
      pd=pd[row.names(pd) %in% tmp_snames,]
      if(class(data_m)==class(list())){
        for(i in 1:length(data_m)){
          data_m[[i]]=data_m[[i]][,colnames(data_m[[i]]) %in% row.names(pd)]
        }
      } else {
        data_m=data_m[,colnames(data_m) %in% row.names(pd)]
      }
    }
    
    
  }
  
  if(updateBICCN){
    anno_biccn=qread("~/myBucket/SC_data/Mouse/brain/snRNA/BICCN/BICCN.integrative.lbls.arranged.qs")
    if(sum(pd$ds_batch=="mouse_snBICCN")>0){
      anno_biccn$dataset[which(anno_biccn$dataset=="10X_nuclei_v3_Broad")]="mouse_snBICCN"
      anno_biccn$dataset[which(anno_biccn$dataset=="10X_cells_v2_AIBS")]="mouse_scBICCN_v2"
      anno_biccn$dataset[which(anno_biccn$dataset=="10X_cells_v3_AIBS")]="mouse_scBICCN_v3"
      
    } else {
      anno_biccn$dataset[which(anno_biccn$dataset=="10X_nuclei_v3_Broad")]="mouse_BICCN"
      anno_biccn$dataset[which(anno_biccn$dataset=="10X_cells_v2_AIBS")]="mouse_BICCN_v2"
      anno_biccn$dataset[which(anno_biccn$dataset=="10X_cells_v3_AIBS")]="mouse_BICCN_v3"
      
    }
    
    anno_biccn$id=paste0(anno_biccn$dataset,"_",anno_biccn$cell_lbl)
    anno_biccn=anno_biccn[match(row.names(pd),anno_biccn$id),]
    pd$cluster_label[!is.na(anno_biccn$cluster_label)]=anno_biccn$cluster_label[!is.na(anno_biccn$cluster_label)]
  }
  
  pd_mouse_sn=pd[which(pd$ds_batch %in% c("mouse_BICCN","mouse_snBICCN")),]
  if(nrow(pd_mouse_sn)>0){
    p=.my2dPlot2(inputPCA = pd_mouse_sn,batchCol = "cluster_label",reductionCols = c('UMAP_1','UMAP_2'))
    ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony_Mouse_sn",j,"2",updatelbl),pdf=T),argList)),width = 12,height = 8)
    
  }
  
  pd_mouse_sc=pd[which(pd$ds_batch %in% c("mouse_BICCN_v2","mouse_scBICCN_v2")),]
  if(nrow(pd_mouse_sc)>0){
    p=.my2dPlot2(inputPCA = pd_mouse_sc,batchCol = "cluster_label",reductionCols = c('UMAP_1','UMAP_2'))
    ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony_Mouse_scv2",j,"2",updatelbl),pdf=T),argList)),width = 12,height = 8)
    
  }
  
  pd_mouse_sc=pd[which(pd$ds_batch %in% c("mouse_BICCN_v3","mouse_scBICCN_v3")),]
  if(nrow(pd_mouse_sc)>0){
    p=.my2dPlot2(inputPCA = pd_mouse_sc,batchCol = "cluster_label",reductionCols = c('UMAP_1','UMAP_2'))
    ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony_Mouse_scv3",j,"2",updatelbl),pdf=T),argList)),width = 12,height = 8)
    
  }
  
  
  cat("     Mouse Broad_snRNAseq marker figures...\n")
  if(nrow(pd_mouse_sn)>0&writeMarkers){
    mouse_markers=unique(pd_mouse_sn$cluster_label)
    mouse_markers=unique(unlist(lapply(strsplit(mouse_markers,"_"),function(x)x[1])))
    mouse_markers=toupper(unique(unlist(strsplit(mouse_markers," "))))
    
    tmp=data_m
    #tmp=tmp[,colnames(tmp) %in% row.names(pd)]
    pd2=pd#[match(colnames(tmp),row.names(pd)),]
    #p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = Tushar_markers)
    #ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony_",j,"4_Tushar"),pdf=T),argList)),width = 49,height = 40)
    
    p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(mouse_markers,function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6,ncores=ncores)
    ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j,"5","_MouseSn"),pdf=T),argList)),width = 7*min(length(mouse_markers),6),height = 40/7*max(ceiling(length(mouse_markers)/6),7))
    
    
    if(length(which(pd2$organism=="mouse"))>0){
      pd_mouse=pd2[which(pd2$organism=="mouse"),]
      tmp_mouse=list()
      for(ij in 1:length(tmp)){
        tmp2=tmp[[ij]]
        if(sum(colnames(tmp2) %in% row.names(pd_mouse))>0){
          tmp2=tmp2[,colnames(tmp2) %in% row.names(pd_mouse)]
          tmp_mouse=c(tmp_mouse,tmp2)
        }
      }
      if(length(names(tmp))>0){
        names(tmp_mouse)=names(tmp)
      }
      p=.my2dPlot_counts(inputPCA=pd_mouse,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(mouse_markers,function(x) x),geneNameCol="gene_short_name",expData=tmp_mouse,ncolumns=6,ncores=ncores)
      ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j,"5mouse","_MouseSn"),pdf=T),argList)),width = 7*min(length(mouse_markers),6),height = 40/7*max(ceiling(length(mouse_markers)/6),7))
      
    }
    
    if(length(which(pd2$organism=="human"))>0){
      pd_human=pd2[which(pd2$organism=="human"),]
      tmp_human=list()
      for(ij in 1:length(tmp)){
        tmp2=tmp[[ij]]
        if(sum(colnames(tmp2) %in% row.names(pd_human))>0){
          tmp2=tmp2[,colnames(tmp2) %in% row.names(pd_human)]
          tmp_human=c(tmp_human,tmp2)
        }
      }
      if(length(names(tmp))>0){
        names(tmp_human)=names(tmp)
      }
      p=.my2dPlot_counts(inputPCA=pd_human,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(mouse_markers,function(x) x),geneNameCol="gene_short_name",expData=tmp_human,ncolumns=6,ncores=ncores)
      ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony",j,"5human","_MouseSn"),pdf=T),argList)),width = 7*min(length(mouse_markers),6),height = 40/7*max(ceiling(length(mouse_markers)/6),7))
    }
    
  }
  
  return("Done!")
}

.Anno_dsFn=function(j,data_m,argList,ncores,offsetCount=0,changeSaveDir=T,Liger=F,FeaturePlotOnly=F,dsLevel=F){
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  
  if(sum(names(argList)=="projection")==0){
    argList$projection=F
  }
  
  if(changeSaveDir){
    .saveDir=paste0("~/data/",argList$cellClass)
    if(.ArgList$includeHuman&.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_HumanMouse")
    } else if(.ArgList$includeHuman) {
      .saveDir=paste0(.saveDir,"_Human")
    } else if(.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_Mouse")
    } else {
      .saveDir=paste0(.saveDir,"_unknown")
    }
    if(!argList$indScaling){
      .saveDir=paste0(.saveDir,"_notIndScaling")
    } else {
      .saveDir=paste0(.saveDir,"_indScaling")
    }
    if(Liger){
      argList$saveDir=paste0(.saveDir,"_liger")
    } else {
      argList$saveDir=.saveDir
    }
    
    if(dsLevel){
      .saveDir=paste0(.saveDir,"_dsLevel")
    }
    if(argList$projection){
      .saveDir=paste0(.saveDir,"_projection")
    }
    
    if(argList$nPCs!=30){
      .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
    }
    
    
    
    if(argList$Leng200){
      .saveDir=paste0(.saveDir,"_leng200")
    }
    if(argList$Jakel200){
      .saveDir=paste0(.saveDir,"_Jakel200")
    }
    argList$saveDir=.saveDir
    
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  
  
  
  if(sum(argList=="pseudocell_size")==0){
    argList$pseudocell_size=200
  }
  
  reRun=T
  if(file.exists(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))){
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))
    if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2){
      if(sum(is.na(pd$UMAP_1))==0&sum(is.na(pd$UMAP_2))==0){
        reRun=F
      }
    }
  }
  
  if(reRun){
    stop("First step needs to be run!")
  } else {
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))
  }
  
  argList$saveDir=paste0("~/data/markerPlots/",argList$cellClass,"_HumanMouse_indScaling")
  if(argList$Leng200){
    argList$saveDir=paste0(argList$saveDir,"_Leng200")
  }
  if(argList$Jakel200){
    argList$saveDir=paste0(argList$saveDir,"_Jakel200")
  }
  
  
  if(FeaturePlotOnly){
    argList$saveDir=paste0(argList$saveDir,"_featurePlot")
  }
  
  if(!dir.exists(.ArgList$saveDir)){
    dir.create(argList$saveDir,recursive = T)
  }
  
  pd_Tushar=pd[which(pd$ds_batch=="human_Tushar_iNPH"),]
  if(nrow(pd_Tushar)>0){
    
    Tushar_markers=unique(pd_Tushar$anno_cellState)
    Tushar_markers=unique(unlist(strsplit(Tushar_markers,"_")))
    
    tmp=data_m
    #tmp=tmp[,colnames(tmp) %in% row.names(pd)]
    pd2=pd#[match(colnames(tmp),row.names(pd)),]
    #p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = Tushar_markers)
    #ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony_",j,"4_Tushar"),pdf=T),.ArgList)),width = 49,height = 40)
    
    ds_list=unlist(lapply(data_m,function(x) x$ds_batch[1]))
    
    res=list()
    for(ids in unique(ds_list)){
      tmp=data_m[which(ds_list==ids)]
      if(sum(pd$ds_batch==ids)>0){
        if(length(tmp)==1){
          tmp=tmp[[1]]
          tmp=tmp[,colnames(tmp) %in% row.names(pd)]
          pd2=pd[match(colnames(tmp),row.names(pd)),]
          #inputSeurat = tmp;inputDimData = pd2;inputGenes = Tushar_markers;combine_figs = F;order=F
          p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = Tushar_markers,combine_figs = F,order=F)
        } else {
          pd2=pd[pd$ds_batch==ids,]
          
          #inputPCA=pd2;batch_values="anno_batch";reductionCols=c("UMAP_1","UMAP_2");geneNameList=lapply(Tushar_markers,function(x) x);geneNameCol="gene_short_name";expData=tmp;ncolumns=6;ncores=ncores;combine_figs = F;exponetial_pwr=8
          p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(Tushar_markers,function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6,ncores=ncores,combine_figs = F)
        }
        
        res=c(res,list(p))
        names(res)[length(res)]=ids
      }
      
    }
    
    geneNameList=unlist(lapply(res,names))
    geneNameList=unique(geneNameList)
    
    for(geneName in geneNameList){
      
      for(ik in 1:length(res)){
        p=ggplot(data=NULL)
        if(sum(names(res[[ik]])==geneName)>0){
          p=res[[ik]][[geneName]]
        }
        p=p+ggtitle(names(res)[ik])
        if(ik==1){
          p_combined=p
        } else {
          p_combined=p_combined+p
        }
      }
      
      #ggsave(plot=p_combined,file=paste0("~/myBucket/torm.png"),width = 29,height = 29)
      ggsave(plot=p_combined,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("dsPattern_",j+offsetCount,"_",geneName),pdf=T),argList)),width = 29,height = 29)
      
    }
    
  }
  
  return("Done!")
}
.Anno_dsFn2=function(j,data_m,argList,ncores,offsetCount=0,changeSaveDir=T,Liger=F,FeaturePlotOnly=F,dsLevel=F){
  
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  
  if(sum(names(argList)=="projection")==0){
    argList$projection=F
  }
  
  if(changeSaveDir){
    .saveDir=paste0("~/data/",argList$cellClass)
    if(.ArgList$includeHuman&.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_HumanMouse")
    } else if(.ArgList$includeHuman) {
      .saveDir=paste0(.saveDir,"_Human")
    } else if(.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_Mouse")
    } else {
      .saveDir=paste0(.saveDir,"_unknown")
    }
    if(!argList$indScaling){
      .saveDir=paste0(.saveDir,"_notIndScaling")
    } else {
      .saveDir=paste0(.saveDir,"_indScaling")
    }
    if(Liger){
      argList$saveDir=paste0(.saveDir,"_liger")
    } else {
      argList$saveDir=.saveDir
    }
    
    if(dsLevel){
      .saveDir=paste0(.saveDir,"_dsLevel")
    }
    if(argList$projection){
      .saveDir=paste0(.saveDir,"_projection")
    }
    
    if(argList$nPCs!=30){
      .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
    }
    
    
    
    if(argList$Leng200){
      .saveDir=paste0(.saveDir,"_leng200")
    }
    if(argList$Jakel200){
      .saveDir=paste0(.saveDir,"_Jakel200")
    }
    argList$saveDir=.saveDir
    
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  
  
  
  reRun=T
  if(file.exists(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))){
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))
    if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2){
      if(sum(is.na(pd$UMAP_1))==0&sum(is.na(pd$UMAP_2))==0){
        reRun=F
      }
    }
  }
  
  if(reRun){
    stop("First step needs to be run!")
  } else {
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))
  }
  if(changeSaveDir){
    argList$saveDir=paste0("~/data/markerPlots/",argList$cellClass,"_HumanMouse_indScaling")
    if(argList$Leng200){
      argList$saveDir=paste0(argList$saveDir,"_Leng200")
    }
    if(argList$Jakel200){
      argList$saveDir=paste0(argList$saveDir,"_Jakel200")
    }
    
    
    if(FeaturePlotOnly){
      argList$saveDir=paste0(argList$saveDir,"_featurePlot")
    }
  }
  
  if(!dir.exists(argList$saveDir)){
    dir.create(argList$saveDir,recursive = T)
  }
  
  pd_Tushar=pd[which(pd$ds_batch=="human_Tushar_iNPH"),]
  if(nrow(pd_Tushar)>0){
    
    markerGenes=read.table("myBucket/markergenes/all_markers.txt",sep="\t",header = T,stringsAsFactors = F)
    markerGenes=rbind(markerGenes,data.frame(Gene=c("CD200R1","BRIP1","TOP2A","FTH1","LYVE1","GPR126","HCG22","PAQR5","PLAT","MGAM","CECR2","CD300E","TNFAIP3","IFIT2"),pattern=NA,source="Tushar",MainMarker="Yes",page="page2",stringsAsFactors = F))
    markerGenes=markerGenes[!duplicated(markerGenes$Gene),]
    markerGenes=markerGenes[!grepl(",",markerGenes$Gene),]
    
    tmp=data_m
    #tmp=tmp[,colnames(tmp) %in% row.names(pd)]
    pd2=pd#[match(colnames(tmp),row.names(pd)),]
    #p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = Tushar_markers)
    #ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("hm_harmony_",j,"4_Tushar"),pdf=T),argList)),width = 49,height = 40)
    
    ds_list=unlist(lapply(data_m,function(x) x$ds_batch[1]))
    
    res=list()
    for(ids in unique(ds_list)){
      tmp=data_m[which(ds_list==ids)]
      if(length(tmp)==1){
        tmp=tmp[[1]]
        pd2=pd[match(colnames(tmp),row.names(pd)),]
        p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = markerGenes$Gene,combine_figs = F,order=F)
      } else {
        pd2=pd[pd$ds_batch==ids,]
        
        #inputPCA=pd2;batch_values="anno_batch";reductionCols=c("UMAP_1","UMAP_2");geneNameList=lapply(Tushar_markers,function(x) x);geneNameCol="gene_short_name";expData=tmp;ncolumns=6;ncores=ncores;combine_figs = F;exponetial_pwr=8
        p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(markerGenes$Gene,function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6,ncores=ncores,combine_figs = F)
      }
      
      res=c(res,list(p))
      names(res)[length(res)]=ids
    }
    
    geneNameList=unlist(lapply(res,names))
    geneNameList=unique(geneNameList)
    
    for(geneName in geneNameList){
      
      for(ik in 1:length(res)){
        p=ggplot(data=NULL)
        if(sum(names(res[[ik]])==geneName)>0){
          p=res[[ik]][[geneName]]
        }
        p=p+ggtitle(names(res)[ik])
        if(ik==1){
          p_combined=p
        } else {
          p_combined=p_combined+p
        }
      }
      
      ggsave(plot=p_combined,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0("dsPattern_",j+offsetCount,"_",geneName),pdf=T),argList)),width = 25,height = 25)
      
    }
    
  }
  
  return("Done!")
}

.lblTransfer=function(j,res,source_ds_name,target_ds_name,source_label_col,target_label_col=NULL,argList,offsetCount=0,n.adaptiveKernel=5,nPropIter=3,updateBICCN=T,changeSaveDir=F,Liger=F,dsLevel=F){
  harmony_embeddings=res[[j]]
  j=j+offsetCount
  
  if(is.null(target_label_col)){
    target_label_col=source_label_col
  }
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  if(sum(names(argList)=="projection")==0){
    argList$projection=F
  }
  
  if(changeSaveDir){
    .saveDir=paste0("~/data/",argList$cellClass)
    if(.ArgList$includeHuman&.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_HumanMouse")
    } else if(.ArgList$includeHuman) {
      .saveDir=paste0(.saveDir,"_Human")
    } else if(.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_Mouse")
    } else {
      .saveDir=paste0(.saveDir,"_unknown")
    }
    if(!argList$indScaling){
      .saveDir=paste0(.saveDir,"_notIndScaling")
    } else {
      .saveDir=paste0(.saveDir,"_indScaling")
    }
    if(Liger){
      argList$saveDir=paste0(.saveDir,"_liger")
    } else {
      argList$saveDir=.saveDir
    }
    
    if(dsLevel){
      .saveDir=paste0(.saveDir,"_dsLevel")
    }
    if(argList$projection){
      .saveDir=paste0(.saveDir,"_projection")
    }
    
    if(argList$nPCs!=30){
      .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
    }
    
    
    argList$saveDir=.saveDir
    if(argList$Leng200){
      argList$saveDir=paste0(argList$saveDir,"_leng200")
    }
    if(argList$Jakel200){
      argList$saveDir=paste0(argList$saveDir,"_Jakel200")
    }
    
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  
  
  reRun=T
  if(T){
    if(file.exists(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F))){
      reRunCheck=tryCatch({load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F));F}, error=function(e) {return(T)})
      if(!reRunCheck){
        if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2){
          if(sum(is.na(pd$UMAP_1))==0&sum(is.na(pd$UMAP_2))==0){
            reRun=F
          }
        }
      }
      
    }
  }
  
  
  if(reRun){
    stop("First steps needs to be run")
  } else {
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F))
  }
  
  if(updateBICCN){
    anno_biccn=qread("~/myBucket/SC_data/Mouse/brain/snRNA/BICCN/BICCN.integrative.lbls.arranged.qs")
    if(sum(pd$ds_batch=="mouse_snBICCN")>0){
      anno_biccn$dataset[which(anno_biccn$dataset=="10X_nuclei_v3_Broad")]="mouse_snBICCN"
      anno_biccn$dataset[which(anno_biccn$dataset=="10X_cells_v2_AIBS")]="mouse_scBICCN_v2"
      anno_biccn$dataset[which(anno_biccn$dataset=="10X_cells_v3_AIBS")]="mouse_scBICCN_v3"
      
    } else {
      anno_biccn$dataset[which(anno_biccn$dataset=="10X_nuclei_v3_Broad")]="mouse_BICCN"
      anno_biccn$dataset[which(anno_biccn$dataset=="10X_cells_v2_AIBS")]="mouse_BICCN_v2"
      anno_biccn$dataset[which(anno_biccn$dataset=="10X_cells_v3_AIBS")]="mouse_BICCN_v3"
      
    }
    
    anno_biccn$id=paste0(anno_biccn$dataset,"_",anno_biccn$cell_lbl)
    anno_biccn=anno_biccn[match(row.names(pd),anno_biccn$id),]
    pd$cluster_label[!is.na(anno_biccn$cluster_label)]=anno_biccn$cluster_label[!is.na(anno_biccn$cluster_label)]
  }
  
  pd=pd[which(pd$ds_batch %in% c(source_ds_name,target_ds_name)),]
  
  
  harmony_embeddings=harmony_embeddings[match(row.names(pd),row.names(harmony_embeddings)),]
  
  
  label_col=source_label_col
  training_idx=which(pd$ds_batch==source_ds_name)
  
  res=.myKnnLabelTransferFn(inputPCAembeddings=harmony_embeddings,meta_data=pd,training_idx=training_idx,label_col=source_label_col,n.adaptiveKernel=n.adaptiveKernel,nPropIter=nPropIter)
  
  resTest=res$test_labels
  
  slCols=colnames(resTest)[which(grepl("inferred_",colnames(resTest)))]
  
  plotTestDf=NULL
  for(i in slCols){
    tmp=data.frame(label=resTest[,target_label_col],inferred=gsub("inferred_","",i),weight=resTest[,i],stringsAsFactors = F)
    plotTestDf=rbind(plotTestDf,tmp)
  }
  
  plotTestDf2=aggregate(weight~label+inferred,data=plotTestDf,sum)
  
  plotTestDf2=plotTestDf2[order(plotTestDf2$weight,decreasing = T),]
  
  plotTrainingDf=NULL
  for(i in slCols){
    tmp=data.frame(label=res$training_labels[,source_label_col],inferred=gsub("inferred_","",i),weight=res$training_labels[,i],stringsAsFactors = F)
    plotTrainingDf=rbind(plotTrainingDf,tmp)
  }
  
  plotTrainingDf2=aggregate(weight~label+inferred,data=plotTrainingDf,sum)
  
  return(c(res,list(dSource=plotTrainingDf2,dTarget=plotTestDf2,source_label_col=source_label_col,target_label_col=target_label_col,propagation_res=res,harmony_embeddings=harmony_embeddings)))
  
}



.Anno_iNPH_lblTransfer=function(j,argList,data_m,res,Liger=F,changeSaveDir=T,offsetCount=0,newRun=T,dsLevel=F,plotMarkers=T){
  
  harmony_embeddings=res[[j]]
  j=j+offsetCount
  
  if(sum(names(argList)=="pseudocell_size")==0){
    argList$pseudocell_size=200
  }
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  if(sum(names(argList)=="projection")==0){
    argList$projection=F
  }
  
  if(changeSaveDir){
    .saveDir=paste0("~/data/",argList$cellClass)
    if(.ArgList$includeHuman&.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_HumanMouse")
    } else if(.ArgList$includeHuman) {
      .saveDir=paste0(.saveDir,"_Human")
    } else if(.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_Mouse")
    } else {
      .saveDir=paste0(.saveDir,"_unknown")
    }
    if(!argList$indScaling){
      .saveDir=paste0(.saveDir,"_notIndScaling")
    } else {
      .saveDir=paste0(.saveDir,"_indScaling")
    }
    if(Liger){
      argList$saveDir=paste0(.saveDir,"_liger")
    } else {
      argList$saveDir=.saveDir
    }
    
    if(dsLevel){
      .saveDir=paste0(.saveDir,"_dsLevel")
    }
    if(argList$projection){
      .saveDir=paste0(.saveDir,"_projection")
    }
    
    if(argList$nPCs!=30){
      .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
    }
    
    
    if(argList$Leng200){
      .saveDir=paste0(.saveDir,"_leng200")
    }
    if(argList$Jakel200){
      .saveDir=paste0(.saveDir,"_Jakel200")
    }
    argList$saveDir=.saveDir
    
    
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  
  
  
  reRun=T
  if(T){
    if(file.exists(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F))){
      reRunCheck=tryCatch({load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F));F}, error=function(e) {return(T)})
      if(!reRunCheck){
        if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2){
          if(sum(is.na(pd$UMAP_1))==0&sum(is.na(pd$UMAP_2))==0){
            reRun=F
          }
        }
      }
      
    }
  }
  
  
  if(reRun){
    stop("First steps needs to be run")
  } else {
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F))
  }
  
  .lblTransfer=function(target_ds_name,harmony_embeddings,pd,source_ds_name,source_label_col,target_label_col=NULL,n.adaptiveKernel=10,nPropIter=3){
    
    
    if(is.null(target_label_col)){
      target_label_col=source_label_col
    }
    
    
    
    
    
    pd=pd[which(pd$ds_batch %in% c(source_ds_name,target_ds_name)),]
    
    
    harmony_embeddings=harmony_embeddings[match(row.names(pd),row.names(harmony_embeddings)),]
    
    
    label_col=source_label_col
    training_idx=which(pd$ds_batch==source_ds_name)
    
    #inputPCAembeddings=harmony_embeddings;meta_data=pd;training_idx=training_idx;label_col=source_label_col
    res=.myKnnLabelTransferFn(inputPCAembeddings=harmony_embeddings,meta_data=pd,training_idx=training_idx,label_col=source_label_col,n.adaptiveKernel=n.adaptiveKernel,nPropIter=nPropIter)
    
    return(res)
    
  }
  
  if(!file.exists(.myFilePathMakerFn(paste0("iNPH_lblTransfer",j+offsetCount),argList=argList,pdf=F))|newRun){
    #source_ds_name="human_Tushar_iNPH";source_label_col="anno_cellState";target_label_col=NULL;n.adaptiveKernel=5;nPropIter=3
    res=parallel::mclapply(setdiff(unique(pd$ds_batch),"human_Tushar_iNPH"),.lblTransfer,harmony_embeddings=harmony_embeddings,pd=pd,source_ds_name="human_Tushar_iNPH",source_label_col="anno_cellState",target_label_col=NULL,n.adaptiveKernel=5,nPropIter=3,mc.cores = 3)
    res2=lapply(res,function(tst) tst$test_labels)
    res2=do.call("rbind",res2)
    save(res2,file=.myFilePathMakerFn(paste0("iNPH_lblTransfer",j+offsetCount),argList=argList,pdf=F))
  } else {
    load(.myFilePathMakerFn(paste0("iNPH_lblTransfer",j+offsetCount),argList=argList,pdf=F))
  }
  
  
  if(plotMarkers){
    if(!dir.exists(paste0(argList$saveDir,"/plots"))){
      dir.create(paste0(argList$saveDir,"/plots"))
    }
    
    argList$saveDir=paste0(argList$saveDir,"/plots")
    
    geneList=unique(pd$anno_cellState[which(pd$ds_batch=="human_Tushar_iNPH")])
    geneList=unique(unlist(strsplit(geneList,"_")))
    
    for(igene in geneList){
      scores=res2[,grepl("^inferred_",colnames(res2))]
      scores=scores[,grepl(igene,colnames(scores))]
      if(class(scores)!="numeric"){
        scores=rowSums(scores)
      } else {
        names(scores)=row.names(res2)
      }
      
      cluster_name=rep("Low",nrow(pd))
      cluster_name[row.names(pd) %in% names(scores)[scores>0.35]]='Intermediate_low'
      cluster_name[row.names(pd) %in% names(scores)[scores>0.5]]='Intermediate'
      cluster_name[row.names(pd) %in% names(scores)[scores>0.85]]="Cluster-high"
      tmp_count=as.data.frame(table(cluster_name))
      tmp_count[,2]=paste0(tmp_count[,2]," (",round(tmp_count[,2]/sum(tmp_count[,2]),3),")")
      tmp_count=tmp_count[match(cluster_name,as.character(tmp_count[,1])),]
      pd$cluster_name=paste0(cluster_name,"\n",tmp_count[,2])
      
      #inputDataList=data_m;single_gene=igene;inputPd=pd;cluster_col="cluster_name";normalizeExp=T;gene_name_col="gene_short_name";dataset_name_col="ds_batch"
      tst=.my2d_dotPlot(inputDataList=data_m,single_gene=igene,inputPd=pd,cluster_col="cluster_name",normalizeExp=T,gene_name_col="gene_short_name",dataset_name_col="ds_batch")+ggtitle(igene)
      ggsave(plot=tst,file=.myFilePathMakerFn(paste0(igene,"_",j+offsetCount),argList=argList,pdf=T))
      
    }
  }
  
  
  
  return(paste0("Done ","iNPH_lblTransfer",j+offsetCount,"!"))
}


.Anno_iNPH_lblTransfer2=function(j,argList,data_m,res,Liger=F,changeSaveDir=T,offsetCount=0,dsLevel=F){
  
  harmony_embeddings=res[[j]]
  j=j+offsetCount
  
  if(sum(names(argList)=="pseudocell_size")==0){
    argList$pseudocell_size=200
  }
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  if(sum(names(argList)=="projection")==0){
    argList$projection=F
  }
  
  if(changeSaveDir){
    .saveDir=paste0("~/data/",argList$cellClass)
    if(.ArgList$includeHuman&.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_HumanMouse")
    } else if(.ArgList$includeHuman) {
      .saveDir=paste0(.saveDir,"_Human")
    } else if(.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_Mouse")
    } else {
      .saveDir=paste0(.saveDir,"_unknown")
    }
    if(!argList$indScaling){
      .saveDir=paste0(.saveDir,"_notIndScaling")
    } else {
      .saveDir=paste0(.saveDir,"_indScaling")
    }
    if(Liger){
      argList$saveDir=paste0(.saveDir,"_liger")
    } else {
      argList$saveDir=.saveDir
    }
    
    if(dsLevel){
      .saveDir=paste0(.saveDir,"_dsLevel")
    }
    if(argList$projection){
      .saveDir=paste0(.saveDir,"_projection")
    }
    
    if(argList$nPCs!=30){
      .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
    }
    
    
    if(argList$Leng200){
      .saveDir=paste0(.saveDir,"_leng200")
    }
    if(argList$Jakel200){
      .saveDir=paste0(.saveDir,"_Jakel200")
    }
    argList$saveDir=.saveDir
    
    
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  
  
  
  reRun=T
  if(T){
    if(file.exists(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F))){
      reRunCheck=tryCatch({load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F));F}, error=function(e) {return(T)})
      if(!reRunCheck){
        if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2){
          if(sum(is.na(pd$UMAP_1))==0&sum(is.na(pd$UMAP_2))==0){
            reRun=F
          }
        }
      }
      
    }
  }
  
  
  if(reRun){
    stop("First steps needs to be run")
  } else {
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F))
  }
  
  if(!file.exists(.myFilePathMakerFn(paste0("iNPH_lblTransfer",j+offsetCount),argList=argList,pdf=F))){
    stop("First propagation step needs to be run!")
  } else {
    load(.myFilePathMakerFn(paste0("iNPH_lblTransfer",j+offsetCount),argList=argList,pdf=F))
  }
  
  geneList=unique(pd$anno_cellState[which(pd$ds_batch=="human_Tushar_iNPH")])
  geneList=unique(unlist(strsplit(geneList,"_")))
  
  
  total_gene_list=unlist(lapply(data_m,function(x){
    x=x@assays$RNA@meta.features$gene_short_name
    x=toupper(x)
    x=x[!is.na(x)]
    x=unique(x)
    x
  }))
  geneList=geneList[which(toupper(geneList) %in% total_gene_list)]
  
  argList$saveDir=paste0(argList$saveDir,"/plots")
  for(igene in geneList){
    
    tmp_pd=res2
    sl_cols=colnames(res2)[grepl("^inferred_",colnames(res2))&grepl(igene,toupper(colnames(res2)))]
    if(length(sl_cols)==1){
      tmp_pd$inference_score=tmp_pd[,sl_cols]
    } else if(length(sl_cols)>1){
      tmp_pd$inference_score=rowSums(tmp_pd[,sl_cols])
    } else {
      stop("Error!")
    }
      
    
    tst=.myProp_geneDensityPlot(inputDataList=data_m,single_gene=igene,inputPd=tmp_pd,lbl_inf_col="inference_score",gene_name_col="gene_short_name",dataset_name_col="ds_batch")
    #inputDataList=data_m;single_gene=igene;inputPd=pd;cluster_col="cluster_name";normalizeExp=T;gene_name_col="gene_short_name";dataset_name_col="ds_batch"
    p=tst$p_all
    ggsave(plot=p,file=.myFilePathMakerFn(paste0(igene,"_",j+offsetCount,"_all"),argList=argList,pdf=T))
    
    p=tst$p_ind
    ggsave(plot=p,file=.myFilePathMakerFn(paste0(igene,"_",j+offsetCount,"_Ind"),argList=argList,pdf=T))
    
  }
  
  return(paste0("Done ","iNPH_lblTransfer",j+offsetCount,"!"))
}


#j=1;argList=.ArgList;data_m=.data2$data;res=res[[1]];Liger=F;changeSaveDir=T;offsetCount=0
.Anno_iNPH_guided_clustering=function(j,argList,data_m,res,Liger=F,changeSaveDir=T,offsetCount=0,newRun=F,dsLevel=F){
  
  harmony_embeddings=res[[j]]
  j=j+offsetCount
  
  if(sum(names(argList)=="pseudocell_size")==0){
    argList$pseudocell_size=200
  }
  if(sum(names(argList)=="projection")==0){
    argList$projection=F
  }
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  
  if(changeSaveDir){
    .saveDir=paste0("~/data/",argList$cellClass)
    if(.ArgList$includeHuman&.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_HumanMouse")
    } else if(.ArgList$includeHuman) {
      .saveDir=paste0(.saveDir,"_Human")
    } else if(.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_Mouse")
    } else {
      .saveDir=paste0(.saveDir,"_unknown")
    }
    if(!argList$indScaling){
      .saveDir=paste0(.saveDir,"_notIndScaling")
    } else {
      .saveDir=paste0(.saveDir,"_indScaling")
    }
    if(Liger){
      argList$saveDir=paste0(.saveDir,"_liger")
    } else {
      argList$saveDir=.saveDir
    }
    
    if(dsLevel){
      .saveDir=paste0(.saveDir,"_dsLevel")
    }
    if(argList$projection){
      .saveDir=paste0(.saveDir,"_projection")
    }
    
    if(argList$nPCs!=30){
      .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
    }
    
    
    if(argList$Leng200){
      .saveDir=paste0(.saveDir,"_leng200")
    }
    if(argList$Jakel200){
      .saveDir=paste0(.saveDir,"_Jakel200")
    }
    argList$saveDir=.saveDir
    
    
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  reRun=T
  if(T){
    if(file.exists(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F))){
      reRunCheck=tryCatch({load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F));F}, error=function(e) {return(T)})
      if(!reRunCheck){
        if(sum(colnames(pd) %in% c("UMAP_1","UMAP_2"))==2){
          if(sum(is.na(pd$UMAP_1))==0&sum(is.na(pd$UMAP_2))==0){
            reRun=F
          }
        }
      }
      
    }
  }
  
  
  if(reRun){
    stop("First steps needs to be run")
  } else {
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j),argList=argList,pdf=F))
  }
  
  if(!file.exists(.myFilePathMakerFn(paste0("iNPH_lblTransfer",j+offsetCount),argList=argList,pdf=F))){
    stop("First propagation step needs to be run!")
  } else {
    load(.myFilePathMakerFn(paste0("iNPH_lblTransfer",j+offsetCount),argList=argList,pdf=F))
  }
  
  
  
  if(sum(colnames(res2)=="inferred_prop_iNPH_lbl")==0|newRun|!file.exists(.myFilePathMakerFn(paste0("iNPH_lblTransfer_total",j+offsetCount),argList=argList,pdf=F))){
    prop_annotations=res2[,grepl("^inferred",colnames(res2))]
    prop_annotations=prop_annotations[,which(colnames(prop_annotations)!="inferred_prop_iNPH_lbl")]
    prop_annotations=apply(prop_annotations,1,function(x) {
      x=x/(max(x,na.rm = T)+0.00001)
      x=which(x>0.9)
      if(length(x)==1){
        x=colnames(prop_annotations)[x]
      } else {
        x=NA
      }
      x
    })
    
    #tst=pd[which(pd$ds_batch=="human_Tushar_iNPH"),]
    #tst=setdiff(unique(tst$anno_cellState),gsub("inferred_","",prop_annotations))
    #tst
    res2$inferred_prop_iNPH_lbl=prop_annotations
    
    
    harmony_net=Seurat:::FindNeighbors.default(object=harmony_embeddings[,1:argList$nPCs])
    
    
    tst=pd[which(pd$ds_batch=="human_Tushar_iNPH"),]
    df_anno=data.frame(sample=row.names(res2),anno=res2$inferred_prop_iNPH_lbl,stringsAsFactors = F)
    df_anno=rbind(df_anno,data.frame(sample=row.names(tst),anno=paste0("inferred_",tst$anno_cellState),stringsAsFactors = F))
    
    #df_anno=data.frame(sample=row.names(tst),anno=paste0("inferred_",tst$anno_cellState),stringsAsFactors = F)
    snn=harmony_net[["snn"]]
    snn=snn[df_anno$sample,df_anno$sample[!is.na(df_anno$anno)]]
    df_anno$anno=gsub("-",".",df_anno$anno)
    snn=snn %*% as.matrix(.myOneHotFn(df_anno$anno[!is.na(df_anno$anno)]))
    
    snn=apply(as.matrix(snn),1,function(x) {
      x=x/(max(x,na.rm = T)+0.00001)
      x=which(x>0.99)
      if(length(x)==1){
        x=colnames(snn)[x]
      } else {
        x=NA
      }
      x
    })
    
    snn_inph=snn[row.names(tst)]
    
    snn=snn[row.names(res2)]
    if(!all(names(snn)==row.names(res2))){
      stop("Error in matching!")
    }
    
    inph_concordance=table(snn_inph,tst$anno_cellState)
    inph_names=colnames(inph_concordance)
    inph_concordance=sweep(inph_concordance,1,rowSums(inph_concordance),"/")
    inph_concordance2=as.data.frame(inph_concordance)
    inph_concordance=diag(inph_concordance)
    names(inph_concordance)=inph_names
    rm(inph_names)
    
    inph_concordance2$snn_inph=paste0("inferred_",as.character(inph_concordance2$snn_inph))
    p=ggplot(inph_concordance2,aes(snn_inph,Var2,fill=Freq))+geom_tile(color="black")+scale_fill_gradient2(low="white",mid="orange",high="red",midpoint = 0.5)+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1))
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("concordancePlot_v2","_",j+offsetCount),argList=argList,pdf=T),width = 20,height = 20)
    
    
    other_concordance=table(snn,res2$inferred_prop_iNPH_lbl)
    other_names=colnames(other_concordance)
    other_concordance=diag(sweep(other_concordance,1,rowSums(other_concordance),"/"))
    names(other_concordance)=gsub("inferred_","",other_names)
    rm(other_names)
    df=merge(data.frame(name=names(inph_concordance),inph=inph_concordance,stringsAsFactors = F),data.frame(name=names(other_concordance),inferred=other_concordance,stringsAsFactors = F),by="name")
    df=reshape::melt(df)
    p=ggplot(df,aes(name,value,fill=variable))+geom_bar(stat="identity",position = "dodge",alpha=0.75,color="black")+theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.8,hjust = 1,color="black"))+scale_y_continuous(breaks = seq(0,1,0.1))+viridis::scale_fill_viridis(discrete = T)+xlab("Initial label")+ylab("% retained")
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("concordancePlot","_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 6)
    
    
    res_markers=NULL
    diff_samples=which(snn!=res2$inferred_prop_iNPH_lbl)
    for(ids in 1:length(data_m)){
      tmp=data_m[[ids]]
      if(sum(colnames(tmp) %in% names(diff_samples))>0){
        tmp=tmp[,intersect(colnames(tmp),names(diff_samples))]
        for(i in 1:ncol(tmp)){
          m1=unlist(strsplit(snn[colnames(tmp)[i]],"_"))
          m2=unlist(strsplit(res2$inferred_prop_iNPH_lbl[row.names(res2)==colnames(tmp)[i]],"_"))
          m0=intersect(m1,m2)
          m1=setdiff(m1,m0)
          m2=setdiff(m2,m0)
          rm(m0)
          if(length(m1)>0&length(m2)>0){
            tmp2=tmp@assays$RNA@counts[toupper(tmp@assays$RNA@meta.features$gene_short_name) %in% toupper(c(m1)),i]
            m1=sum(tmp2>0)/length(m1)
            tmp2=tmp@assays$RNA@counts[toupper(tmp@assays$RNA@meta.features$gene_short_name) %in% toupper(c(m2)),i]
            m2=sum(tmp2>0)/length(m2)
            if((m1>0|m2>0)&m1!=m2){
              res_markers=rbind(res_markers,data.frame(batch=tmp$anno_batch[1],ds=tmp$ds_batch[1],initial_score=m1,refined_score=m2,initial_cluster=snn[colnames(tmp)[i]],refined_cluster=res2$inferred_prop_iNPH_lbl[row.names(res2)==colnames(tmp)[i]],stringsAsFactors = F))
            }
          }
          
        }
      }
    }
    
    
    tst$anno_cellState=paste0("inferred_",tst$anno_cellState)
    diff_samples=which(snn_inph!=tst$anno_cellState)
    for(ids in 1:length(data_m)){
      tmp=data_m[[ids]]
      if(sum(colnames(tmp) %in% names(diff_samples))>0){
        tmp=tmp[,intersect(colnames(tmp),names(diff_samples))]
        for(i in 1:ncol(tmp)){
          m1=unlist(strsplit(snn_inph[colnames(tmp)[i]],"_"))
          m2=unlist(strsplit(tst$anno_cellState[row.names(tst)==colnames(tmp)[i]],"_"))
          m0=intersect(m1,m2)
          m1=setdiff(m1,m0)
          m2=setdiff(m2,m0)
          rm(m0)
          
          if(length(m1)>0&length(m2)>0){
            tmp2=tmp@assays$RNA@counts[toupper(tmp@assays$RNA@meta.features$gene_short_name) %in% toupper(c(m1)),i]
            m1=sum(tmp2>0)/length(m1)
            tmp2=tmp@assays$RNA@counts[toupper(tmp@assays$RNA@meta.features$gene_short_name) %in% toupper(c(m2)),i]
            m2=sum(tmp2>0)/length(m2)
            if((m1>0|m2>0)&m1!=m2){
              res_markers=rbind(res_markers,data.frame(batch=tmp$anno_batch[1],ds=tmp$ds_batch[1],initial_score=m1,refined_score=m2,initial_cluster=snn_inph[colnames(tmp)[i]],refined_cluster=tst$anno_cellState[row.names(tst)==colnames(tmp)[i]],stringsAsFactors = F))
            }
          }
          
        }
      }
    }
    
    p=res_markers[,c("ds","initial_score","refined_score")]
    p=reshape::melt(p)
    p=ggplot(data=p,aes(ds,value,fill=variable))+geom_violin()+ 
      stat_summary(fun = "mean",geom = "point",color = "black")+theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.8,hjust = 1,color="black"))
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("clustering_byDS_newAnnotations","_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 6)
    
    p=NULL
    for(i in unique(res_markers$ds)){
      tmp=res_markers[res_markers$ds==i,]
      p=rbind(p,data.frame(ds=i,refinement=sum(tmp$initial_score<tmp$refined_score)/nrow(tmp),stringsAsFactors = F))
    }
    p=ggplot(data=p,aes(ds,refinement))+geom_bar(stat="identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.8,color="black"))+geom_hline(yintercept = 0.5,color="purple")
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("clustering_byDS_refinements","_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 6)
    
    
    p=res_markers[res_markers$ds=="human_Tushar_iNPH",c("initial_cluster","initial_score","refined_score")]
    p=reshape::melt(p)
    p=ggplot(data=p,aes(initial_cluster,value,fill=variable))+geom_violin()+ 
      stat_summary(fun = "mean",geom = "point",color = "black")+theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.8,hjust = 1,color="black"))
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("clustering_iNPH_newAnnotations","_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 6)
    
    p=res_markers[res_markers$ds!="human_Tushar_iNPH",c("initial_cluster","initial_score","refined_score")]
    p=reshape::melt(p)
    p=ggplot(data=p,aes(initial_cluster,value,fill=variable))+geom_violin()+ 
      stat_summary(fun = "mean",geom = "point",color = "black")+theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.8,hjust = 1,color="black"))
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("clustering_initialCluster_others_newAnnotations","_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 6)
    
    pdata=res_markers[res_markers$ds=="human_Tushar_iNPH",c("initial_cluster","initial_score","refined_score")]
    p=NULL
    for(i in unique(pdata$initial_cluster)){
      tmp=pdata[pdata$initial_cluster==i,]
      p=rbind(p,data.frame(initial_cluster=i,refinement=sum(tmp$initial_score<tmp$refined_score)/nrow(tmp),stringsAsFactors = F))
    }
    p=ggplot(data=p,aes(initial_cluster,refinement))+geom_bar(stat="identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.8,color="black"))+geom_hline(yintercept = 0.5,color="purple")
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("clustering_initialCluster_iNPH_refinements","_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 6)
    
    pdata=res_markers[res_markers$ds!="human_Tushar_iNPH",c("initial_cluster","initial_score","refined_score")]
    p=NULL
    for(i in unique(pdata$initial_cluster)){
      tmp=pdata[pdata$initial_cluster==i,]
      p=rbind(p,data.frame(initial_cluster=i,refinement=sum(tmp$initial_score<tmp$refined_score)/nrow(tmp),stringsAsFactors = F))
    }
    p=ggplot(data=p,aes(initial_cluster,refinement))+geom_bar(stat="identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.8,color="black"))+geom_hline(yintercept = 0.5,color="purple")
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("clustering_initialCluster_iNPH_refinements","_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 6)
    
    
    pdata=res_markers[,c("initial_cluster","initial_score","refined_score")]
    p=NULL
    for(i in unique(pdata$initial_cluster)){
      tmp=pdata[pdata$initial_cluster==i,]
      p=rbind(p,data.frame(initial_cluster=i,refinement=sum(tmp$initial_score<tmp$refined_score)/nrow(tmp),stringsAsFactors = F))
    }
    improved_clusters=p$initial_cluster[p$refinement>0.5]
    #x2=res2$inferred_prop_iNPH_lbl[(res2$inferred_prop_iNPH_lbl %in% improved_clusters)]
    #x1=snn[(res2$inferred_prop_iNPH_lbl %in% improved_clusters)]
    if(length(improved_clusters)>0){
      res2$inferred_prop_iNPH_lbl[(res2$inferred_prop_iNPH_lbl %in% improved_clusters)]=snn[(res2$inferred_prop_iNPH_lbl %in% improved_clusters)]
    }
    
    res2$inferred_prop_iNPH_lbl[is.na(res2$inferred_prop_iNPH_lbl)]=snn[is.na(res2$inferred_prop_iNPH_lbl)]
    
    res2_inhp=tst
    res2_inhp$inferred_prop_iNPH_lbl=res2_inhp$anno_cellState
    res2_inhp$inferred_prop_iNPH_lbl[(res2_inhp$inferred_prop_iNPH_lbl %in% improved_clusters)]=snn_inph[(res2_inhp$inferred_prop_iNPH_lbl %in% improved_clusters)]
    res2_t=plyr::rbind.fill(res2,res2_inhp)
    
    tmp_inph=as.data.frame(table(tst$anno_cellState))
    tmp_inph[,1]=as.character(tmp_inph[,1])
    tmp_inph[,2]=tmp_inph[,2]/sum(tmp_inph[,2])
    p=NULL
    for(i in unique(res2$ds_batch)){
      tmp=res2[which(res2$ds_batch==i),]
      tmp=as.data.frame(table(tmp$inferred_prop_iNPH_lbl))
      tmp[,2]=tmp[,2]/sum(tmp[,2])
      tmp[,1]=as.character(tmp[,1])
      tmp=merge(tmp,tmp_inph,by="Var1")
      tmp$fraction=(tmp$Freq.x-tmp$Freq.y)/(tmp$Freq.y+0.02) *100
      tmp$ds=i
      p=rbind(p,tmp[,c("ds","Var1","fraction")])
    }
    p=ggplot(p,aes(ds,Var1,fill=fraction))+geom_tile(color="black")+scale_fill_gradient2(low="darkblue",mid="white",high="darkred",midpoint = 0)+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1))+ggtitle("% difference vs. iNPH","offset of 0.02")
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("clustering_pctDiff","_",j+offsetCount),argList=argList,pdf=T),width = 20,height = 20)
    
    
    res2_inhp$anno_cellState=gsub("inferred_","",as.character(res2_inhp$anno_cellState))
    inph_concordance=table(res2_inhp$inferred_prop_iNPH_lbl,res2_inhp$anno_cellState)
    inph_names=colnames(inph_concordance)
    inph_concordance=sweep(inph_concordance,1,rowSums(inph_concordance),"/")
    inph_concordance2=as.data.frame(inph_concordance)
    names(inph_concordance)=inph_names
    rm(inph_names)
    
    p=ggplot(inph_concordance2,aes(Var2,Var1,fill=Freq))+geom_tile(color="black")+scale_fill_gradient2(low="white",mid="orange",high="red",midpoint = 0.5)+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1))
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("concordancePlot_Final","_",j+offsetCount),argList=argList,pdf=T),width = 20,height = 20)
    
    
    p=.my2dPlot2(inputPCA = res2_t,batchCol = "inferred_prop_iNPH_lbl",reductionCols = c('UMAP_1','UMAP_2'))
    ggsave(plot = p,file=.myFilePathMakerFn(paste0("clustering_umap_Final_v1","_",j+offsetCount),argList=argList,pdf=T),width = 30,height = 30)
    
    p=.myDimPlotFn(object=res2_t, dimCols = c("UMAP_1", "UMAP_2"),attCol = "inferred_prop_iNPH_lbl")
    ggsave(plot = p,file=.myFilePathMakerFn(paste0("clustering_umap_Final_v2","_",j+offsetCount),argList=argList,pdf=T),width = 20,height = 12)
    
    tmp_res2=res2_t[which(res2_t$ds_batch=="human_Tushar_iNPH"),]
    tmp_res2$anno_cellState=gsub("inferred_","",as.character(tmp_res2$anno_cellState))
    p=.myDimPlotFn(object=tmp_res2, dimCols = c("UMAP_1", "UMAP_2"),attCol = "anno_cellState")
    ggsave(plot = p,file=.myFilePathMakerFn(paste0("clustering_umap_inph","_",j+offsetCount),argList=argList,pdf=T),width = 20,height = 12)
    
    
    #data_m
    save(res2,file=.myFilePathMakerFn(paste0("iNPH_lblTransfer",j+offsetCount),argList=argList,pdf=F))
    
    save(res2_t,file=.myFilePathMakerFn(paste0("iNPH_lblTransfer_total",j+offsetCount),argList=argList,pdf=F))
    
  }
  
  load(.myFilePathMakerFn(paste0("iNPH_lblTransfer_total",j+offsetCount),argList=argList,pdf=F))
  geneList=unique(pd$anno_cellState[which(pd$ds_batch=="human_Tushar_iNPH")])
  geneList=unique(unlist(strsplit(geneList,"_")))
  
  
  total_gene_list=unlist(lapply(data_m,function(x){
    x=x@assays$RNA@meta.features$gene_short_name
    x=toupper(x)
    x=x[!is.na(x)]
    x=unique(x)
    x
  }))
  geneList=geneList[which(toupper(geneList) %in% total_gene_list)]
  
  argList$saveDir=paste0(argList$saveDir,"/plots")
  
  tmpSeuratFn=function(inputData){
    inputData=Seurat::NormalizeData(inputData,verbose=F)
    return(inputData)
  }
  data_m=parallel::mclapply(data_m,tmpSeuratFn,mc.cores = 5)
  
  
  
  for(igene in geneList){
    
    #if(!file.exists(.myFilePathMakerFn(paste0("clustering_",igene,"_",j+offsetCount),argList=argList,pdf=T))|newRun)
    if(T){
      tmp2=res2_t
      tmp2$cluster_name=tmp2$inferred_prop_iNPH_lbl
      tmp2=tmp2[!is.na(tmp2$cluster_name),]
      tmp2$cluster_name=gsub("inferred_","",as.character(tmp2$cluster_name))
      
      row.names(tmp2)=tmp2$sample
      #inputDataList=data_m;single_gene=igene;inputPd=tmp2;cluster_col="cluster_name";normalizeExp=F;gene_name_col="gene_short_name";dataset_name_col="ds_batch";return_plot=F;ignore_missing_anno=T
      tmp=.my2d_dotPlot(inputDataList=data_m,single_gene=igene,inputPd=tmp2,cluster_col="cluster_name",normalizeExp=F,gene_name_col="gene_short_name",dataset_name_col="ds_batch",ignore_missing_anno=T,return_plot=F,ncores=10)#+ggtitle(igene)
      
      inputPd=pd[which(pd$ds_batch=="human_Tushar_iNPH"),]
      tmp2=.my2d_dotPlot(inputDataList=data_m,single_gene=igene,inputPd=inputPd,cluster_col="anno_cellState",normalizeExp=F,gene_name_col="gene_short_name",dataset_name_col="ds_batch",ignore_missing_anno=T,return_plot=F,ncores=10)
      
      if(!is.null(tmp2$p_data)){
        if(nrow(tmp2$p_data)>0){
          tmp2$p_data$dataset=paste0("org_",tmp2$p_data$dataset)
          p_data=rbind(tmp$p_data,tmp2$p_data)
          p=tmp2$plot
          
          p_data$features.plot=igene
          p_data$avg.exp.scaled[which(is.nan(p_data$avg.exp.scaled))]=0
          tmp_avg.exp=aggregate(avg.exp~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
          tmp_pct.exp=aggregate(pct.exp~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
          tmp_avg.exp.scaled=aggregate(avg.exp.scaled~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
          
          tmp_data=merge(tmp_avg.exp,tmp_avg.exp.scaled,by=c("id","features.plot","dataset"),all=T)
          tmp_data=merge(tmp_data,tmp_pct.exp,by=c("id","features.plot","dataset"),all=T)
          
          
          tmp_feature.plot=tmp_data$features.plot
          tmp_data$features.plot=tmp_data$dataset
          tmp_data$gene=tmp_feature.plot
          p$data=tmp_data
          p=p+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5))+xlab("Dataset")+ggtitle(igene)
          
          ggsave(plot=p,file=.myFilePathMakerFn(paste0("clustering_",igene,"_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 15)
          
        }
      }
      
      
    }
    
  }
  
  
  
  return(paste0("Done ","iNPH_lblTransfer",j+offsetCount,"!"))
}


.Anno_otherAD=function(j,argList,Liger=F,changeSaveDir=T,offsetCount=0,newRun=F,dsLevel=F){
  
  j=j+offsetCount
  
  if(sum(names(argList)=="pseudocell_size")==0){
    argList$pseudocell_size=200
  }
  if(sum(names(argList)=="projection")==0){
    argList$projection=F
  }
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  
  if(changeSaveDir){
    .saveDir=paste0("~/data/",argList$cellClass)
    if(.ArgList$includeHuman&.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_HumanMouse")
    } else if(.ArgList$includeHuman) {
      .saveDir=paste0(.saveDir,"_Human")
    } else if(.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_Mouse")
    } else {
      .saveDir=paste0(.saveDir,"_unknown")
    }
    if(!argList$indScaling){
      .saveDir=paste0(.saveDir,"_notIndScaling")
    } else {
      .saveDir=paste0(.saveDir,"_indScaling")
    }
    if(Liger){
      argList$saveDir=paste0(.saveDir,"_liger")
    } else {
      argList$saveDir=.saveDir
    }
    
    if(dsLevel){
      .saveDir=paste0(.saveDir,"_dsLevel")
    }
    if(argList$projection){
      .saveDir=paste0(.saveDir,"_projection")
    }
    
    if(argList$nPCs!=30){
      .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
    }
    
    
    if(argList$Leng200){
      .saveDir=paste0(.saveDir,"_leng200")
    }
    if(argList$Jakel200){
      .saveDir=paste0(.saveDir,"_Jakel200")
    }
    argList$saveDir=.saveDir
    
    
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  
  
  reRun=T
  if(T){
    if(file.exists(.myFilePathMakerFn(paste0("iNPH_lblTransfer_total",j+offsetCount),argList=argList,pdf=F))){
      reRunCheck=tryCatch({load(.myFilePathMakerFn(paste0("iNPH_lblTransfer_total",j+offsetCount),argList=argList,pdf=F));F}, error=function(e) {return(T)})
      if(!reRunCheck){
        if(sum(colnames(res2_t) %in% c("UMAP_1","UMAP_2"))==2){
          if(sum(is.na(res2_t$UMAP_1))==0&sum(is.na(res2_t$UMAP_2))==0){
            reRun=F
          }
        }
      }
      
    }
  }
  
  
  if(reRun){
    stop("First steps needs to be run")
  } else {
    load(.myFilePathMakerFn(paste0("iNPH_lblTransfer_total",j+offsetCount),argList=argList,pdf=F))
  }
  pd=res2_t
  
  if(sum(colnames(pd)=="inferred_prop_iNPH_lbl")==0){
    stop("cluster assignments couldn't be found")
  }
  
  
  #mathys et al
  
  tmp_pd=pd[which(pd$ds_batch=="human_Mathys_31042697"),]
  tmp_pd$anno_cellState=tmp_pd$Subcluster
  if(length(unique(tmp_pd$anno_cellState))>1){
    tmp_tbl=table(tmp_pd$inferred_prop_iNPH_lbl,tmp_pd$anno_cellState)
    tmp_tbl=sweep(tmp_tbl,1,rowSums(tmp_tbl),"/")
    tmp_tbl=as.data.frame(tmp_tbl)
    p=ggplot(tmp_tbl,aes(Var2,Var1,fill=Freq))+geom_tile(color="black")+scale_fill_gradient(low="white",high="red")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1))
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("clustering_mathys_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 15)
    
  }
  
  #Leng et al
  leng_anno=qread("~/myBucket/SC_data/Human/brain/snRNA/Leng_GSE147528/Anno_Leng_processed.qs")
  tmp_pd=pd[which(pd$ds_batch=="human_Leng_GSE147528"),]
  leng_anno$orig=paste0("human_Leng_GSE147528_",as.character(leng_anno$orig))
  leng_anno=leng_anno[match(tmp_pd$sample,leng_anno$orig),]
  leng_anno$clusterCellType=as.character(leng_anno$clusterCellType)
  leng_anno$clusterCellType[leng_anno$clusterCellType=="Exc"]="ExN"
  leng_anno$clusterCellType[leng_anno$clusterCellType=="Inh"]="InN"
  leng_anno$clusterCellType[leng_anno$clusterCellType=="Micro"]="MG"
  tmp_pd$anno_cellState=leng_anno$clusterAssignment
  tmp_pd$anno_cellState[which(leng_anno$clusterCellType!=tmp_pd$anno_cellClass)]=NA
  tmp_pd=tmp_pd[!is.na(tmp_pd$anno_cellState),]
  tmp_pd$anno_cellState=as.character(tmp_pd$anno_cellState)
  
  
  if(length(unique(tmp_pd$anno_cellState))>1){
    tmp_tbl=table(tmp_pd$inferred_prop_iNPH_lbl,tmp_pd$anno_cellState)
    tmp_tbl=sweep(tmp_tbl,1,rowSums(tmp_tbl),"/")
    tmp_tbl=as.data.frame(tmp_tbl)
    p=ggplot(tmp_tbl,aes(Var2,Var1,fill=Freq))+geom_tile(color="black")+scale_fill_gradient(low="white",high="red")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1))
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("clustering_leng","_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 15)
    
    p=.my2dPlot2(inputPCA = tmp_pd,batchCol = "anno_cellState",reductionCols = c('UMAP_1','UMAP_2'))
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("clustering_leng_umap1","_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 15)
    
    p=.myDimPlotFn(object=tmp_pd, dimCols = c("UMAP_1", "UMAP_2"),attCol = "anno_cellState")
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("clustering_leng_umap2","_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 15)
    
  }
  
  return(paste0("Done ","iNPH_lblTransfer",j+offsetCount,"!"))
}

.Anno_MGanno=function(j,argList,Liger=F,changeSaveDir=T,offsetCount=0,newRun=F,dsLevel=F){
  
  j=j+offsetCount
  
  if(sum(names(argList)=="pseudocell_size")==0){
    argList$pseudocell_size=200
  }
  if(sum(names(argList)=="projection")==0){
    argList$projection=F
  }
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  
  if(changeSaveDir){
    .saveDir=paste0("~/data/",argList$cellClass)
    if(.ArgList$includeHuman&.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_HumanMouse")
    } else if(.ArgList$includeHuman) {
      .saveDir=paste0(.saveDir,"_Human")
    } else if(.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_Mouse")
    } else {
      .saveDir=paste0(.saveDir,"_unknown")
    }
    if(!argList$indScaling){
      .saveDir=paste0(.saveDir,"_notIndScaling")
    } else {
      .saveDir=paste0(.saveDir,"_indScaling")
    }
    if(Liger){
      argList$saveDir=paste0(.saveDir,"_liger")
    } else {
      argList$saveDir=.saveDir
    }
    
    if(dsLevel){
      .saveDir=paste0(.saveDir,"_dsLevel")
    }
    if(argList$projection){
      .saveDir=paste0(.saveDir,"_projection")
    }
    
    if(argList$nPCs!=30){
      .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
    }
    
    
    if(argList$Leng200){
      .saveDir=paste0(.saveDir,"_leng200")
    }
    if(argList$Jakel200){
      .saveDir=paste0(.saveDir,"_Jakel200")
    }
    argList$saveDir=.saveDir
    
    
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  
  
  reRun=T
  if(T){
    if(file.exists(.myFilePathMakerFn(paste0("iNPH_lblTransfer_total",j+offsetCount),argList=argList,pdf=F))){
      reRunCheck=tryCatch({load(.myFilePathMakerFn(paste0("iNPH_lblTransfer_total",j+offsetCount),argList=argList,pdf=F));F}, error=function(e) {return(T)})
      if(!reRunCheck){
        if(sum(colnames(res2_t) %in% c("UMAP_1","UMAP_2"))==2){
          if(sum(is.na(res2_t$UMAP_1))==0&sum(is.na(res2_t$UMAP_2))==0){
            reRun=F
          }
        }
      }
      
    }
  }
  
  
  if(reRun){
    stop("First steps needs to be run")
  } else {
    load(.myFilePathMakerFn(paste0("iNPH_lblTransfer_total",j+offsetCount),argList=argList,pdf=F))
  }
  pd=res2_t
  
  if(sum(colnames(pd)=="inferred_prop_iNPH_lbl")==0){
    stop("cluster assignments couldn't be found")
  }
  
  
  tmp=pd
  tmp$anno_cellState2=tmp$anno_cellState
  tmp$anno_cellState2[grepl("Microglia_",tmp$anno_cellState)]="MG"
  tmp$anno_cellState2[grepl("mp_",tmp$anno_cellState)]="Macro"
  tmp$anno_cellState2[grepl("Macro",tmp$anno_cellState)]="Macro"
  tmp$anno_cellState2[grepl("Macro",tmp$anno_cellState)]="Macro"
  tmp$anno_cellState2[grepl("mgl_",tmp$anno_cellState)]="MG"
  tmp$anno_cellState2[grepl("^Mic",tmp$anno_cellState)]="MG"
  tmp$anno_cellState2[grepl("^homeo",tmp$anno_cellState)]="MG"
  tmp$anno_cellState2[grepl("top2a",tolower(tmp$anno_cellState))|grepl("brip1",tolower(tmp$anno_cellState))|grepl("prolif",tolower(tmp$anno_cellState))]="prolif"
  tmp$anno_cellState2[which(tmp$anno_cellState %in% c("inferred_Microglia_GPNMB_LPL","ARM"))]="DAM"
  tmp$anno_cellState2[which(tmp$anno_cellState %in% c("mg","H1/2","H1/2M","H1M","H2M","ARM","TRM","MHC.high"))]="MG"
  #mathys et al
  tmp=tmp[!is.na(tmp$anno_cellState),]
  p=.myDimPlotFn(object=tmp, dimCols = c("UMAP_1", "UMAP_2"),attCol = "anno_cellState2")
  ggsave(plot=p,file=.myFilePathMakerFn(paste0("Major_cellState","_",j+offsetCount),argList=argList,pdf=T),width = 15,height = 15)
  
  #ggsave(plot=p,file="~/myBucket/torm.pdf",width = 15,height = 15)
  
  return(paste0("Done ","iNPH_lblTransfer",j+offsetCount,"!"))
}


#iHVG=4;argList=.ArgList;Liger=F;changeSaveDir=T;offsetCount=0;newRun=F;dsLevel=F
.myiNPH_QCfn=function(iHVG,argList,Liger=F,changeSaveDir=T,offsetCount=0,newRun=F,dsLevel=F){
  require(ggplot2)
  require(Seurat)
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  if(sum(names(argList)=="projection")==0){
    argList$projection=F
  }
  
  if(changeSaveDir){
    .saveDir=paste0("~/data/",argList$cellClass)
    if(.ArgList$includeHuman&.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_HumanMouse")
    } else if(.ArgList$includeHuman) {
      .saveDir=paste0(.saveDir,"_Human")
    } else if(.ArgList$includeMouse){
      .saveDir=paste0(.saveDir,"_Mouse")
    } else {
      .saveDir=paste0(.saveDir,"_unknown")
    }
    if(!argList$indScaling){
      .saveDir=paste0(.saveDir,"_notIndScaling")
    } else {
      .saveDir=paste0(.saveDir,"_indScaling")
    }
    if(Liger){
      argList$saveDir=paste0(.saveDir,"_liger")
    } else {
      argList$saveDir=.saveDir
    }
    
    if(dsLevel){
      .saveDir=paste0(.saveDir,"_dsLevel")
    }
    if(argList$projection){
      .saveDir=paste0(.saveDir,"_projection")
    }
    
    if(argList$nPCs!=30){
      .saveDir=paste0(.saveDir,"_nPC",argList$nPCs)
    }
    
    
    if(argList$Leng200){
      .saveDir=paste0(.saveDir,"_leng200")
    }
    if(argList$Jakel200){
      .saveDir=paste0(.saveDir,"_Jakel200")
    }
    argList$saveDir=.saveDir
    
    
    if(!dir.exists(argList$saveDir)){
      dir.create(argList$saveDir,recursive = T)
    }
  }
  
  
  
  argList$HVG_count=iHVG
  for(j in 1:4){
    load(.myFilePathMakerFn(paste0("UMAP_anno_hm_harmony",j+offsetCount),argList=argList,pdf=F))
    
    load(.myFilePathMakerFn(paste0("iNPH_lblTransfer_total",j+offsetCount),argList=argList,pdf=F))
    
    tst=table(pd$anno_cellState,pd$anno_cluster_res)
    tst=sweep(tst,1,rowSums(tst),"/")
    
    
    
    tst=table(res2_t$inferred_prop_iNPH_lbl,res2_t$anno_cluster_res)
    tst=sweep(tst,1,rowSums(tst),"/")
    
    pd=res2_t[res2_t$ds_batch=="human_Tushar_iNPH",]
    table(pd$anno_cellState,pd$anno_cluster_res)
    
    table(pd$anno_cellState,pd$inferred_prop_iNPH_lbl=="inferred_Astrocyte_HGF_CHI3L1"|pd$anno_cluster_res=="6")
    
    
    if(sum(colnames(pd)=="ds_batch")==0){
      pd$ds_batch=pd$anno_batch
    }
    
    pd$organism="human"
    pd$organism[grepl("mouse",pd$ds_batch)]="mouse"
    p=.myDimPlotFn(object=pd, dimCols = c("UMAP_1", "UMAP_2"), attCol = 'organism')
    ggsave(plot = p,.myFilePathMakerFn(paste0("hm_harmony_organism",j+offsetCount,"_clusters"),argList=argList,pdf=T),width = 12,height = 8)
    
    pd=pd[pd$ds_batch=="human_Tushar_iNPH",]
    
    tst=table(pd$anno_cellState,pd$anno_cluster_res)
    tst=sweep(tst,1,rowSums(tst),"/")
    tst=as.data.frame(tst)
    colnames(tst)=c("x","y","fraction")
    p=ggplot(tst,aes(x=y,y=x,fill=fraction))+geom_tile(color="black")+scale_fill_gradient2(low="white",high="red")+theme_classic()
    ggsave(plot = p,.myFilePathMakerFn(paste0("hm_harmony_confusionMat",j+offsetCount,"_clusters"),argList=argList,pdf=T),width = 12,height = 8)
    
    
    p=.myDimPlotFn(object=pd, dimCols = c("UMAP_1", "UMAP_2"), attCol = 'inferred_prop_iNPH_lbl')
    ggsave(plot = p,"~/myBucket/torm.pdf",width = 12,height = 8)
    
    p=.myDimPlotFn(object=pd[pd$anno!="other"&!grepl("REMO",pd$anno),], dimCols = c("UMAP_1", "UMAP_2"), attCol = 'anno')
    pd$anno[pd$anno=="other"]=pd$inferred_prop_iNPH_lbl[pd$anno=="other"]
    pd=pd[pd$anno_cluster_res!="8",]
    pd=pd[!grepl("REMOVE",pd$anno),]
    tst=as.data.frame(table(pd$anno_cluster_res))
    tst=tst[tst[,2]>10,]
    pd=pd[pd$anno_cluster_res %in% tst[,1],]
    
    
    
  }
  
  
}


