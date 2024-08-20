# Saved files -> variable name

## | file name prefix     | variable name                           |
## |----------------------+-----------------------------------------|
## | "UMAP_anno"          | pd                                      |
## | "pca_centroids"      | pca_centroid                            |
## | "harmony-embeddings" | harmony_embeddings                      |
## | "exp_merged"         | expData                                 |
## | "umap_pseudocells"      | umap_pseudocells                           |
## | "pca_anno"           | pca_res                                 |
## | "UMAP_anno"          | pd=as.matrix(pd[,c("UMAP_1","UMAP_2")]) |
## | "res_DE_wZscore"     | res_arranged                            |
## | "res_dataset_array"  | data                                    |
## | "res_meta"           | data                                    |
## |----------------------+-----------------------------------------|
## | global/varGenes      | tmp                                     |
## | pca/pca_Anno         | pca_res, pd                             |
## | pca_centroids        | pca_centroid                            |
## | kmeans_res_clustes   | res_clusters                            |
## | UMAP_res             | resUMAP                                 |
## | centroid_scores      | resDistances                            |
## | res_DE               | res_arranged                            |
## |----------------------+-----------------------------------------|

#set the main_source and utilities source code addresses
#change the harmony path in .sconline.embeddingFn(), if desired

#library(googleCloudStorageR)
#gcs_source("vgazesta/code/mySC.R")
#gcs_source("vgazesta/code/my_SC_metaAnalysis.R")

source("~/code/mySC.R")
source("~/code/my_SC_metaAnalysis.R")


#library(rliger)
library(qs)

#Accessory functions
.extra_sconline.scatterPlot_summary2d=function(object,reductionCols,n=300){
  
  if(sum(colnames(object) %in% reductionCols)!=2){
    stop("Specified reduction columns couldn't be found ")
  }
  
  
  yval=object[,reductionCols]
  xval=yval[,1]
  yval=yval[,2]
  
  rescoord=expand.grid(seq(min(xval),max(xval),length.out = n),seq(min(yval),max(yval),length.out = n))
  colnames(rescoord)=reductionCols
  
  knet=RANN::nn2(data=rescoord,query = object[,reductionCols],k=1,eps=0)
  counts=as.data.frame(table(knet$nn.idx[,1]))
  counts$Var1=as.character(counts$Var1)
  rescoord$id=1:nrow(rescoord)
  rescoord=merge(rescoord,counts,by.x="id",by.y="Var1")
  rescoord=rescoord[,-which(colnames(rescoord)=="id")]
  
  return(rescoord)
}

#assay = NULL; rev.pca = FALSE; weight.by.var = TRUE;
#verbose = TRUE; ndims.print = 1:5; nfeatures.print = 30; 
#reduction.key = "PC_"; seed.use = 42; approx = TRUE
.mySeuratRunPCA=function (object,projection_data=NULL, assay = NULL, npcs = 50, rev.pca = FALSE, weight.by.var = TRUE, 
          verbose = TRUE, ndims.print = 1:5, nfeatures.print = 30, 
          reduction.key = "PC_", seed.use = 42, approx = TRUE, ...) {
  require(Seurat)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  projection_pca=NULL
  if (rev.pca) {
    npcs <- min(npcs, ncol(x = object) - 1)
    pca.results <- irlba(A = object, nv = npcs, ...)
    total.variance <- sum(Seurat:::RowVar(x = t(x = object)))
    sdev <- pca.results$d/sqrt(max(1, nrow(x = object) - 
                                     1))
    if (weight.by.var) {
      feature.loadings <- pca.results$u %*% diag(pca.results$d)
    }
    else {
      feature.loadings <- pca.results$u
    }
    cell.embeddings <- pca.results$v
  } else {
    total.variance <- sum(Seurat:::RowVar(x = object))
    if (approx) {
      npcs <- min(npcs, nrow(x = object) - 1)
      pca.results <- irlba(A = t(x = object), nv = npcs, ...)
      feature.loadings <- pca.results$v
      
      if(!is.null(projection_data)){
        projection_pca=t(projection_data) %*% feature.loadings
        projection_pca=sweep(projection_pca,2,pca.results$d,"/")
      }
      
      
      sdev <- pca.results$d/sqrt(max(1, ncol(object) -1))
      if (weight.by.var) {
        cell.embeddings <- pca.results$u %*% diag(pca.results$d)
        if(!is.null(projection_data)){
          projection_pca = projection_pca %*% diag(pca.results$d)
        }
        
      }
      else {
        cell.embeddings <- pca.results$u
      }
    } else {
      npcs <- min(npcs, nrow(x = object))
      pca.results <- prcomp(x = t(object), rank. = npcs, ...)
      feature.loadings <- pca.results$rotation
      sdev <- pca.results$sdev
      if (weight.by.var) {
        cell.embeddings <- pca.results$x
      }
      else {
        cell.embeddings <- pca.results$x/(pca.results$sdev[1:npcs] * 
                                            sqrt(x = ncol(x = object) - 1))
      }
    }
  }
  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:npcs)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- Seurat::CreateDimReducObject(embeddings = cell.embeddings, 
                                         loadings = feature.loadings, assay = assay, stdev = sdev, 
                                         key = reduction.key, misc = list(total.variance = total.variance))
  if (verbose) {
    msg <- capture.output(print(x = reduction.data, dims = ndims.print, 
                                nfeatures = nfeatures.print))
    message(paste(msg, collapse = "\n"))
  }
  
  return(list(reduction.data=reduction.data,projection_pca=projection_pca))
}


.myPCAfn=function(data, argList,projection_data=NULL,saveFiles=T,...){
  
  library(future)
  require(purrr)
  require(irlba)
  plan("multicore", workers = 5)
  plan()
  options(future.globals.maxSize = 1000 * 1024^4)
  
  UMI_cor_thr=argList$UMI_cor_thr
  
  myPrepDR=function (scaledData, features, verbose = TRUE,projection_state=F) {
    
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
    # TODO jonah parallize buyt make sure chunked
    features.var <- apply(X = data.use[features, ], MARGIN = 1,
                          FUN = var)
    if(!projection_state){
      features.keep <- features[features.var > 0]
    } else {
      features.keep <- features
    }
    
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
  
  mySeuratFn2_org=function(inputData,varFeatures){
    plan("sequential")
    inputData@assays$RNA@var.features=varFeatures
    inputData = ScaleData(inputData, verbose = FALSE,features=varFeatures)
    return(inputData)
  }
  
  mySeuratFn2_archive=function(inputData,varFeatures,scaleData=T){
    plan("sequential")
    inputData@assays$RNA@var.features=varFeatures
    if(scaleData){
      inputData = Seurat:::ScaleData.default(inputData@assays$RNA@data, verbose = FALSE,features=varFeatures)
    } else {
      inputData=inputData@assays$RNA@data[row.names(inputData) %in% varFeatures,]
    }
    
    return(inputData)
  }
  
  internal_scaleFn=function(inputData,varFeatures,scaleData=T){
    plan("sequential")
    
    if(scaleData){
      inputData = Seurat:::ScaleData.default(counts(inputData), verbose = FALSE,features=varFeatures)
    } else {
      inputData=counts(inputData)[row.names(inputData) %in% varFeatures,]
    }
    
    return(inputData)
  }
  
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
  
  if(is.null(argList$HVG_list)){
    argList$HVG_list=argList$HVG_count
  }
  #argList$HVG_list=unique(c(argList$HVG_list,argList$HVG_count))
  
  reRunCheck=F
  if(sum(names(argList)=="newRun")>0){
    if(argList$newRun){
      reRunCheck=T
    }
  }
  if(!saveFiles){
    reRunCheck=T
    if(length(argList$HVG_list)>1){
      stop("Only one HVG count threshold can be specified")
    }
  }
  
  if(!reRunCheck){
    
    for(iHVG in argList$HVG_list){
      argList$HVG_count=iHVG
      tmpCheck= tryCatch({load(.myFilePathMakerFn("pca_anno",argList = argList,pseudoImportant = F));F}, error=function(e) {return(T)})
      if(!tmpCheck){
        argList$HVG_list=setdiff(argList$HVG_list,iHVG)
      }
    }
    
  }
  
  pca_final_res="Done"
  
  if(reRunCheck|length(argList$HVG_list)>0){
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
      dataList=data$data
      
      varFeatures=c()
      if(!is.null(argList$HVG_list)){
        for(iHVG in argList$HVG_list){
          varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))])
        }
      } else {
        varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(argList$HVG_count-1))])
      }
      
      varFeatures=unique(varFeatures)
      
      # browser()
      dataList2=parallel::mclapply(dataList,internal_scaleFn,varFeatures=varFeatures,mc.cores = argList$ncores)
      sl_ind=which(!unlist(lapply(dataList2,function(x) class(x)[1])) %in% c("array","matrix"))
      while(length(sl_ind)>0){
        gc()
        datalist3=parallel::mclapply(dataList[sl_ind],internal_scaleFn,varFeatures=varFeatures,mc.cores = 2)
        dataList2[sl_ind]=datalist3
        sl_ind=which(!unlist(lapply(dataList2,function(x) class(x)[1])) %in% c("array","matrix"))
        rm(datalist3)
      }
      dataList=dataList2
      rm(dataList2)
      
      if(!is.null(projection_data)){
        projectionDataList=parallel::mclapply(projection_data$data,internal_scaleFn,varFeatures=varFeatures,mc.cores = argList$ncores)
      }
      gc()
      
      resScaled=list()
      all_genes=unique(unlist(lapply(dataList,function(x) row.names(x))))
      for(i in 1:length(dataList)){
        tmp=dataList[[i]]
        tmp=tmp[match(all_genes,row.names(tmp)),,drop=F]
        if(sum(is.na(tmp))>0){
          tmp[is.na(tmp)]=0
        }
        
        row.names(tmp)=all_genes
        resScaled=c(resScaled,list(tmp))
        #dataList[[i]]=tmp
      }
      resScaled=do.call("cbind",resScaled)
      
      
      if(!is.null(projection_data)){
        projectionScaled=list()
        for(i in 1:length(projectionDataList)){
          tmp=projectionDataList[[i]]
          tmp=tmp[match(all_genes,row.names(tmp)),]
          tmp[is.na(tmp)]=0
          row.names(tmp)=all_genes
          projectionScaled=c(projectionScaled,list(tmp))
          
        }
        projectionScaled=do.call("cbind",projectionScaled)
        row.names(projectionScaled)=all_genes
      }
      
      pd=lapply(1:length(data$data),function(i){
        tmp=as.data.frame(colData(data$data[[i]]))
        tmp$sample=colnames(dataList[[i]])
        tmp
      })
      pd=do.call(eval(parse(text='plyr::rbind.fill')), pd)
      
      if(!is.null(projection_data)){
        pd_projection=lapply(1:length(projection_data$data),function(i){
          tmp=as.data.frame(colData(projection_data$data[[i]]))
          tmp$sample=colnames(projection_data$data[[i]])
          tmp
        })
        pd_projection=do.call(eval(parse(text='plyr::rbind.fill')), pd_projection)
        pd_all_projection=pd_projection
      }
      
      pd_all=pd
      gc()
      if(is.null(argList$HVG_list)){
        argList$HVG_list=argList$HVG_count
      }
      
      if(!saveFiles){
        if(length(argList$HVG_count)==1){
          argList$HVG_list=argList$HVG_count
        }
        
      }
      
      if(is.null(argList$input_highly_var_genes)){
        for(iHVG in argList$HVG_list){
          
          argList$HVG_count=iHVG
          tmp_varFeatures=data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))]
          
          #pca_res=.extraPCAfn(object=resScaled[which(row.names(resScaled) %in% tmp_varFeatures),], npcs = 30,findElbowPoint=F,seed.use = 42,reduction.key = "PC_",weight.by.var=T, approx = TRUE)
          #pca_res=pca_res$embeddings
          
          if(F){
            #sklearn: too slow
            library(reticulate)
            pd <- import("pandas", delay_load = TRUE)
            np <- import("numpy", delay_load = TRUE)
            skl_lr <- import("sklearn.decomposition", delay_load = TRUE)
            
            
            
            ipc=skl_lr$IncrementalPCA(n_components=50L, batch_size=10L)
            for(i in 1:length(dataList)){
              ipc = ipc$partial_fit(t(dataList[[i]]))
            }
          }
          
          if(F){
            #rsvd
            tst=rsvd::rsvd(A=t(myPrepDR(scaledData=resScaled,
                                        features=tmp_varFeatures, verbose = TRUE)),k=50)
          }
          
          object=myPrepDR(scaledData=resScaled,
                          features=tmp_varFeatures, verbose = TRUE)
          
          if(!is.null(projection_data)){
            projection_data=myPrepDR(scaledData=projectionScaled[row.names(object),],
                                     features=tmp_varFeatures, verbose = TRUE,projection_state = T)
          }
          
          projection_pca=NULL
          if(F){
            tst=Seurat:::RunPCA.default(object=myPrepDR(scaledData=resScaled,
                                                            features=tmp_varFeatures, verbose = TRUE),
                                            npcs = max(50,argList$nPCs+10))
            tst=tst@cell.embeddings
          } else {
            pca_res=.mySeuratRunPCA(object=object, npcs = max(50,argList$nPCs+10),projection_data=projection_data)
            projection_pca=pca_res$projection_pca
            pca_res=pca_res$reduction.data
          }
          
          pca_res=pca_res@cell.embeddings
          
          if(!is.null(projection_pca)){
            pca_res=rbind(pca_res,projection_pca)
            pd=plyr::rbind.fill(pd_all,pd_all_projection)
          } else {
            pd=pd_all
          }
          
          pd=pd[pd$sample %in% row.names(pca_res),]
          row.names(pd)=pd$sample
          pca_res=pca_res[match(row.names(pd),row.names(pca_res)),]
          if(saveFiles){
            save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
          } else {
            pca_final_res=list(pd=pd,pca_res=pca_res)
          }
          
          gc()
        }
      } else {
        print("Highly variable genes are provided in argList, using it!")
        {
          
          tmp_varFeatures=as.character(argList$input_highly_var_genes)
          
          #pca_res=.extraPCAfn(object=resScaled[which(row.names(resScaled) %in% tmp_varFeatures),], npcs = 30,findElbowPoint=F,seed.use = 42,reduction.key = "PC_",weight.by.var=T, approx = TRUE)
          #pca_res=pca_res$embeddings
          pca_res=Seurat:::RunPCA.default(object=myPrepDR(scaledData=resScaled,
                                                          features=tmp_varFeatures, verbose = TRUE,
                                                          ...),
                                          npcs = max(50,argList$nPCs+10))
          pca_res=pca_res@cell.embeddings
          
          pd=pd_all[pd_all$sample %in% row.names(pca_res),]
          row.names(pd)=pd$sample
          pca_res=pca_res[match(row.names(pd),row.names(pca_res)),]
          if(saveFiles){
            save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
          } else {
            pca_final_res=list(pd=pd,pca_res=pca_res)
          }
          
          gc()
        }
      }
      
      
    }
    
    if(!argList$indScaling){
      
      if(!is.null(projection_data)){
        stop("Projection is not yet implemented!")
      }
      
      dataList=data$data
      varFeatures=c()
      {
        if(is.null(argList$input_highly_var_genes)){
          if(!is.null(argList$HVG_list)){
            for(iHVG in argList$HVG_list){
              varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))])
            }
          } else {
            varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(argList$HVG_count-1))])
          }
          varFeatures=unique(varFeatures)
        } else {
          varFeatures=argList$input_highly_var_genes
        }
        
        
        resNorm=list()
        dataList=parallel::mclapply(dataList,internal_scaleFn,varFeatures=varFeatures,scaleData=F,mc.cores = argList$ncores)
        all_genes=unique(unlist(lapply(dataList,function(x) row.names(x))))
        for(i in 1:length(dataList)){
          tmp=as.matrix(dataList[[i]])
          tmp=tmp[match(all_genes,row.names(tmp)),,drop=F]
          if(sum(is.na(tmp))>0){
            tmp[is.na(tmp),,drop=F]=0
          }
          
          row.names(tmp)=all_genes
          resNorm=c(resNorm,list(tmp))
          #dataList[[i]]=tmp
        }
        resNorm=do.call("cbind",resNorm)
        resScaled=Seurat:::ScaleData.default(as(resNorm,"dgCMatrix"))
        
        
        
        pd=lapply(1:length(dataList),function(i){
          tmp=as.data.frame(colData(data$data[[i]]))
          tmp$sample=colnames(dataList[[i]])
          tmp
        })
        pd=do.call(eval(parse(text='plyr::rbind.fill')), pd)
        
        
        pd_all=pd
        gc()
        if(is.null(argList$HVG_list)){
          argList$HVG_list=argList$HVG_count
        }
        
        if(!saveFiles){
          argList$HVG_list=argList$HVG_count
        }
        
        if(is.null(argList$HVG_list)){
          argList$HVG_list=HVG_count
        }
        if(is.null(argList$input_highly_var_genes)){
          for(iHVG in argList$HVG_list){
            
            argList$HVG_count=iHVG
            tmp_varFeatures=data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))]
            
            #pca_res=.extraPCAfn(object=resScaled[which(row.names(resScaled) %in% tmp_varFeatures),], npcs = 30,findElbowPoint=F,seed.use = 42,reduction.key = "PC_",weight.by.var=T, approx = TRUE)
            #pca_res=pca_res$embeddings
            
            if(F){
              #sklearn: too slow
              library(reticulate)
              pd <- import("pandas", delay_load = TRUE)
              np <- import("numpy", delay_load = TRUE)
              skl_lr <- import("sklearn.decomposition", delay_load = TRUE)
              
              
              
              ipc=skl_lr$IncrementalPCA(n_components=50L, batch_size=10L)
              for(i in 1:length(dataList)){
                ipc = ipc$partial_fit(t(dataList[[i]]))
              }
            }
            
            if(F){
              #rsvd
              tst=rsvd::rsvd(A=t(myPrepDR(scaledData=resScaled,
                                          features=tmp_varFeatures, verbose = TRUE)),k=50)
            }
            
            object=myPrepDR(scaledData=resScaled,
                            features=tmp_varFeatures, verbose = TRUE)
            
            if(!is.null(projection_data)){
              projection_data=myPrepDR(scaledData=projectionScaled[row.names(object),],
                                       features=tmp_varFeatures, verbose = TRUE,projection_state = T)
            }
            
            projection_pca=NULL
            if(F){
              tst=Seurat:::RunPCA.default(object=myPrepDR(scaledData=resScaled,
                                                          features=tmp_varFeatures, verbose = TRUE),
                                          npcs = max(50,argList$nPCs+10))
              tst=tst@cell.embeddings
            } else {
              pca_res=.mySeuratRunPCA(object=object, npcs = max(50,argList$nPCs+10),projection_data=projection_data)
              projection_pca=pca_res$projection_pca
              pca_res=pca_res$reduction.data
            }
            
            pca_res=pca_res@cell.embeddings
            
            if(!is.null(projection_pca)){
              pca_res=rbind(pca_res,projection_pca)
              pd=plyr::rbind.fill(pd_all,pd_all_projection)
            } else {
              pd=pd_all
            }
            
            pd=pd[pd$sample %in% row.names(pca_res),]
            row.names(pd)=pd$sample
            pca_res=pca_res[match(row.names(pd),row.names(pca_res)),]
            if(saveFiles){
              save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
            } else {
              pca_final_res=list(pd=pd,pca_res=pca_res)
            }
            
            gc()
          }
        } else {
          print("Highly variable genes are provided in argList, using it!")
          {
            
            tmp_varFeatures=as.character(argList$input_highly_var_genes)
            
            #pca_res=.extraPCAfn(object=resScaled[which(row.names(resScaled) %in% tmp_varFeatures),], npcs = 30,findElbowPoint=F,seed.use = 42,reduction.key = "PC_",weight.by.var=T, approx = TRUE)
            #pca_res=pca_res$embeddings
            pca_res=Seurat:::RunPCA.default(object=myPrepDR(scaledData=resScaled,
                                                            features=tmp_varFeatures, verbose = TRUE,
                                                            ...),
                                            npcs = max(50,argList$nPCs+10))
            pca_res=pca_res@cell.embeddings
            
            pd=pd_all[pd_all$sample %in% row.names(pca_res),]
            row.names(pd)=pd$sample
            pca_res=pca_res[match(row.names(pd),row.names(pca_res)),]
            if(saveFiles){
              save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
            } else {
              pca_final_res=list(pd=pd,pca_res=pca_res)
            }
            
            gc()
          }
        }
      }
      
      
      
    }
  }
  
  
  return(pca_final_res)
}

.extra_sconline_LigerToExpSet=function(inputData,organism,limit=NULL){
  require(rliger)
  require(scater)
  require(scran)
  
  expList=list()
  for(i in 1:length(inputData@raw.data)){
    tmp=inputData@raw.data[[i]]
    if(!is.null(limit)){
      if(ncol(tmp)>limit){
        tmp=tmp[,sample(ncol(tmp),min(limit,ncol(tmp)))]
      }
    }
    exp=SingleCellExperiment(assays = list(counts = tmp),colData = data.frame(sampleName=colnames(tmp),stringsAsFactors = F),rowData=data.frame(gene=row.names(tmp),stringsAsFactors = F))
    expList=c(expList,list(exp))
    names(expList)[length(expList)]=names(inputData@raw.data)[i]
  }
  
  expList=.extra_sconline_cBindFn(expList,batchNames = names(expList))
  
  pd=as.data.frame(inputData@cell.data)
  
  if(length(inputData@clusters)>0){
    if(!all(row.names(pd)==names(inputData@clusters))){
      stop("Error!")
    }
    pd$clusters=as.character(inputData@clusters)
  }
  
  pd=pd[match(expList$sampleName,row.names(pd)),]
  if(!all(row.names(pd)==expList$sample)){
    stop("Error!")
  }
  
  
  if(sum(colnames(pd)=="anno_batch")>0){
    pd$anno_batch_prv=pd$anno_batch
  }
  
  pd$anno_batch=expList$anno_batch
  
  
  row.names(pd)=as.character(expList$sampleName)
  
  data=.myExpSetCreatorFn(inputExpData=counts(expList),organism=organism,minExpCells=0,inputPdata=pd,inputFdata=as.data.frame(rowData(expList)),addExtraAnno=F,server=T)
  return(data)
}

.extra_sconline_SplitObject=function(object,colName){
  if(class(object)=="Seurat"){
    res=Seurat::SplitObject(object,split.by = colName)
  } else{
    pd=colData(object)[,colName]
    res=list()
    for(i in unique(pd)){
      res=c(res,list(object[,which(pd==i)]))
      names(res)[length(res)]=i
    }
  }
  return(res)
}

.extra_sconline_cBindFn=function(inputList,batchNames=NULL,verbose=F){
  #cbinds multiple singleCellExpression datasets with differring number of rows.
  #inputList: the list of datasets to be merged
  #batchNames: the batch name to be assigned to each dataset. length(batchNames)==length(inputList)
  res_m=""
  if(!is.null(batchNames)){
    if(length(inputList)==1){
      res_m=inputList[[1]]
    } else if(length(inputList)==2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = batchNames[1:2])
    } else if(length(inputList)>2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = batchNames[1:2])
      for(i in 3:length(inputList)){
        res_m=.extracBindDetailFn(x1=res_m,x2=inputList[[i]],batchNames=c("",batchNames[i]))
      }
    }
  } else {
    if(length(inputList)==1){
      res_m=inputList[[1]]
    } else if(length(inputList)==2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = c("",""))
    } else if(length(inputList)>2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = c("",""))
      for(i in 3:length(inputList)){
        
        #x1=res_m;x2=inputList[[i]];batchNames=c("","")
        res_m=.extracBindDetailFn(x1=res_m,x2=inputList[[i]],batchNames=c("",""))
        if(verbose){
          print(paste("dataset:",i,"; nrow:",nrow(res_m),"; ncol",ncol(res_m)))
        }
      }
    }
  }
  
  print("batch information is in the anno_batch variable")
  return(res_m)
  
}

.extra_sconline.exp_creatorFn=function(argList,inputExpData,batch_variable,organism="Mouse",addAnno=F){
  
  if(class(inputExpData)!=class(list())){
    if(class(inputExpData)=="Seurat"){
      inputExpData=.myExpSetCreatorFn(inputExpData=inputExpData@assays$RNA@counts,
                                      organism=organism,
                                      minExpCells=0,
                                      inputPdata=as.data.frame(inputExpData@meta.data),
                                      inputFdata=as.data.frame(inputExpData@assays$RNA@meta.features),
                                      addExtraAnno=addAnno,server=T)
    } else if(class(inputExpData)=="liger"){
      inputExpData=.extra_sconline_LigerToExpSet(inputData=data,organism=organism,limit=NULL)
    } else if(class(inputExpData)!="SingleCellExperiment"){
      stop("Unrecognized input expression data format!")
    }
    
    if(sum(duplicated(row.names(inputExpData)))>0){
      print(paste(sum(duplicated(row.names(inputExpData))),"Duplicate gene ids were found! duplicates were randomly removed from the data"))
      inputExpData=inputExpData[!duplicated(row.names(inputExpData)),]
    }
    
    inputExpData$anno_batch=as.character(colData(inputExpData)[,batch_variable])
    #inputExpData$anno_batch=gsub("-","_",inputExpData$anno_batch)
    if(class(inputExpData)!=class(list())){
      inputExpData=.myReadData_spliterFn(inputData=inputExpData,removeHighExp=argList$excludeHighExp)
    }
  } else {
    if(length(inputExpData)>2&sum(names(inputExpData)=="data_m")==0){
      inputExpData=list(data=inputExpData,data_m=NULL)
    }
  }
  
  
  
  return(inputExpData)
}

.extra_sconline.AffinityFn=function(sim_mat){
  
  thr=apply(sim_mat,1,function(x){
    x=x[!is.na(x)]
    y=NA
    if(length(x)>2){
      x=x[order(x,decreasing = T)]
      y=x[3]
    }
    y
  })
  thr=median(thr,na.rm=T)
  affinity=exp(-3*(pmax(thr-sim_mat,0)))
  #affinity=res_cells/thr
  #thr=apply(affinity,1,function(x) {
  #  x=x[order(x,decreasing = T)]
  #  x[min(20,length(x))]
  #})
  #thr=matrix(thr,nrow=nrow(affinity),ncol=ncol(affinity),byrow = F)
  #affinity[affinity<thr]=0
  #affinity=affinity+t(affinity)
  
  #affinity=t(apply(affinity,1,function(x) x/sum(x)))
  
  return(affinity)
  
}

.extra_sconline.MASC=function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL,
                                 verbose = FALSE, save_models = FALSE, save_model_dir = NULL,jackknife=F) {
  # Check inputs
  require(lme4)
  if (is.factor(dataset[[contrast]]) == FALSE) {
    stop("Specified contrast term is not coded as a factor in dataset")
  }
  
  # Convert cluster assignments to string
  cluster <- as.character(cluster)
  # Prepend design matrix generated from cluster assignments
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  # Create output list to hold results
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]
  
  # Create model formulas
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                          paste0("(1|", random_effects, ")", collapse = " + ")),
                        collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else if (!is.null(fixed_effects) && is.null(random_effects)) {
    model_rhs <- paste0(fixed_effects, collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      # For now, do not allow models without mixed effects terms
      stop("No random effects specified")
    }
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else {
    model_rhs <- "1" # only includes intercept
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      stop("No random or fixed effects specified")
    }
  }
  
  # Initialize list to store model objects for each cluster
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]
  
  # Run nested mixed-effects models for each cluster
  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    if (verbose == TRUE) {
      message(paste("Creating logistic mixed models for", test_cluster))
    }
    null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                   model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                   model_rhs), collapse = ""))
    # Run null and full mixed-effects models
    null_model <- lme4::glmer(formula = null_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    full_model <- lme4::glmer(formula = full_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    model_lrt <- anova(null_model, full_model)
    # calculate confidence intervals for contrast term beta
    contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2:length(levels(dataset[[contrast]]))])
    contrast_ci <- confint.merMod(full_model, method = "Wald",
                                  parm = contrast_lvl2)
    # Save model objects to list
    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt <- model_lrt
    cluster_models[[i]]$confint <- contrast_ci
    #do jackknifing
    boot_pvalvec=c()
    boot_coefvec=c()
    boot_stable=1
    if(jackknife){
      if(length(contrast_lvl2)==1){
        for(ibatch in unique(dataset$anno_batch)){
          tmp_dataset=dataset[which(dataset$anno_batch!=ibatch),]
          boot_null_model <- tryCatch({lme4::glmer(formula = null_fm, data = tmp_dataset,
                                                   family = binomial, nAGQ = 1, verbose = 0,
                                                   control = glmerControl(optimizer = "bobyqa"))},error=function(e) {return(F)})
          
          boot_full_model <- tryCatch({lme4::glmer(formula = full_fm, data = tmp_dataset,
                                                   family = binomial, nAGQ = 1, verbose = 0,
                                                   control = glmerControl(optimizer = "bobyqa"))},error=function(e){return(F)})
          if(class(boot_null_model)!=class(T)&class(boot_full_model)!=class(T)){
            boot_model_lrt <- anova(boot_null_model, boot_full_model)
            # calculate confidence intervals for contrast term beta
            boot_pvalvec=c(boot_pvalvec,boot_model_lrt[["Pr(>Chisq)"]][2])
            boot_coefvec=c(boot_coefvec,fixef(boot_full_model)[[contrast_lvl2]])
          } else {
            boot_stable=0
          }
          
          
        }
      } else {
        warning("jackknife for contrasts with more than 2 levels are not yet implemented!")
        boot_pvalvec=(-1)
        boot_coefvec=(-1)
      }
      
    } else {
      boot_pvalvec=(-1)
      boot_coefvec=(-1)
    }
    
    cluster_models[[i]]$boot_pval_median <- median(boot_pvalvec)
    cluster_models[[i]]$boot_pval_mean <- mean(boot_pvalvec)
    cluster_models[[i]]$boot_pval_max <- max(boot_pvalvec)
    
    cluster_models[[i]]$boot_coef_median <- median(boot_coefvec)
    cluster_models[[i]]$boot_coef_mean <- mean(boot_coefvec)
    cluster_models[[i]]$boot_stable <- boot_stable
    cluster_models[[i]]$boot_coef_min <- boot_coefvec[which(abs(boot_coefvec)==min(abs(boot_coefvec)))[1]]
  }
  
  # Organize results into output dataframe
  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
                       size = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
  for(ilvl in 1:length(contrast_lvl2)){
    output[[paste(contrast_lvl2[ilvl], "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2[ilvl]]]))
    output[[paste(contrast_lvl2[ilvl], "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2[ilvl], "2.5 %"]))
    output[[paste(contrast_lvl2[ilvl], "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2[ilvl], "97.5 %"]))
  }
  
  if(jackknife){
    output[[paste(contrast_lvl2,"JK","Min", "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$boot_coef_min))
    output[[paste(contrast_lvl2,"JK","Mean", "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$boot_coef_mean))
    output[[paste(contrast_lvl2,"JK","Median", "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$boot_coef_median))
    output[[paste(contrast_lvl2,"JK","Max", "pvalue", sep = ".")]] <- sapply(cluster_models, function(x) x$boot_pval_max)
    output[[paste(contrast_lvl2,"JK","Mean", "pvalue", sep = ".")]] <- sapply(cluster_models, function(x) x$boot_pval_mean)
    output[[paste(contrast_lvl2,"JK","Median", "pvalue", sep = ".")]] <- sapply(cluster_models, function(x) x$boot_pval_median)
    output[[paste(contrast_lvl2,"JK","Stable", sep = ".")]] <- sapply(cluster_models, function(x) x$boot_stable)
  }
  
  # Return MASC results and save models if specified
  if (save_models == TRUE) {
    saveModelObj(cluster_models, save_dir = save_model_dir)
    return(output)
  } else {
    return(output)
  }
}

.extra_sconline.NetVisFn=function(net,argList,input_pd=NULL,input_umap_centroid=NULL,directional=T,attribute_col=NULL,lable_nodes=T){
  
  if(sum(colnames(net)=="score")==0){
    net$score=1
  }
  if(sum(colnames(net) %in% c("source","target"))!=2){
    stop("source and target cols are required!")
  }
  
  if(is.null(input_pd)){
    input_pd=.sconline.fetch_data("annotation",argList = argList)
  }
  
  if(is.null(input_umap_centroid)){
    input_umap_centroid=.sconline.fetch_data("umap_pseudocells",argList = argList)
  }
  
  net$weight=net$score
  #net=net[net$score>0.8,]
  net$Fnode=net$source
  net$Snode=net$target
  res_net=net
  net=res_net[order(res_net$weight,decreasing = T),]
  #net=net[!duplicated(net$Snode),]
  net$Fnode=gsub("C","",net$Fnode)
  net$Snode=gsub("C","",net$Snode)
  net=net[!duplicated(paste0(net$Fnode,"_",net$Snode)),]
  net = network::network(net[,c("Fnode","Snode")], directed = directional,matrix.type="edgelist")
  
  netVerNames=network::network.vertex.names(net)
  network::set.edge.attribute(net, "weight", ((res_net$weight-min(res_net$weight)+0.05)/(1-min(res_net$weight)+0.05))^3)
  
  centroid_layout=input_umap_centroid[match(gsub("C","",as.character(netVerNames)),as.character(input_umap_centroid$centroid)),-1]
  centroid_layout=as.matrix(centroid_layout)
  colnames(centroid_layout)=c("x","y")
  
  net=ggnetwork:::fortify.network(net,layout = centroid_layout)
  
  pd_summary=input_pd
  
  
  scale_factor=net[!duplicated(net$vertex.names),]
  scale_factor=merge(scale_factor,input_umap_centroid,by.x="vertex.names",by.y="centroid")
  scale_factor1=lm(x~UMAP_1,data=scale_factor)
  pd_summary$UMAP_1=predict(scale_factor1,newdata=pd_summary)
  scale_factor2=lm(y~UMAP_2,data=scale_factor)
  pd_summary$UMAP_2=predict(scale_factor2,newdata=pd_summary)
  
  pd_summary$xend=pd_summary$UMAP_1
  pd_summary$yend=pd_summary$UMAP_2
  pd_summary$vertex.names=""
  pd_summary$color="gray"
  
  #predicting the background color
  library(ggnetwork)
  #centroids_ideal=c("173","192","191","187","127","194")
  if(is.null(attribute_col)){
    p=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend))+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray")
    if(directional){
      p=p+geom_edges( color = "black",aes(size=weight),arrow = arrow(length = unit(6, "pt"), type = "closed"))
    } else {
      p=p+geom_edges( color = "black",aes(size=weight))
    }
    
    if(lable_nodes){
      p=p + geom_label(aes(label=vertex.names))
    } else {
      p=p+geom_point()
    }
      p=p+theme_blank()+scale_size_continuous(range = c(0.06,1))+scale_color_identity()+scale_fill_identity()+theme(legend.position = "none")
  } else {
    p=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend))+geom_point(data=pd_summary,aes_string('UMAP_1','UMAP_2',color=attribute_col),size=0.01)+
      geom_edges(aes(size=weight))+geom_point()
    if(lable_nodes){
      p=p + geom_label(aes(label=vertex.names))
    }
    p=p+theme_blank()+scale_size_continuous(range = c(0.06,1))+theme(legend.position = "none")
    p=p+scale_color_manual(values=hues::iwanthue(length(unique(pd_summary[,attribute_col]))))
  }
  
  
  return(p)
}

varibow <- function(n_colors) {
  sats <- rep_len(c(0.55,0.7,0.85,1),length.out = n_colors)
  vals <- rep_len(c(1,0.8,0.6),length.out = n_colors)
  sub("FF$","",grDevices::rainbow(n_colors, s = sats, v = vals))
}


V<-View

library(stringr)
.mcsaveRDS <- function(object,file,mc.cores=min(parallel::detectCores(),10, na.rm=T)) {
  file = str_replace_all(file, " ", "\\\\ ")
  con <- pipe(paste0("pigz -p",mc.cores," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}
.mcreadRDS <- function(file,mc.cores=min(parallel::detectCores(),10, na.rm=T)) {
  file = str_replace_all(file, " ", "\\\\ ")
  con <- pipe(paste0("pigz -d -c -p",mc.cores," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}
.printDot <- function(){
  cat(".")
  # TODO check if exists first, for jupyter/console
  flush.console()
}
varSizes <- function(){
  sizes <- sapply(ls(envir=.GlobalEnv), function(n) object.size(get(n)), simplify = FALSE);
  print(sapply(sizes[order(as.integer(sizes))], function(s) format(s, unit = 'auto')))
}

.extra_sconline.NetIdFn=function(inputNet,col1="source1",col2="source2"){
  id=paste0(inputNet[,col1],"_",inputNet[,col2])
  id[inputNet[,col2]>inputNet[,col1]]=paste0(inputNet[,col2],"_",inputNet[,col1])[inputNet[,col2]>inputNet[,col1]]
  return(id)
}


.extra_sconline.purityAnalysisFn=function(argList,inputEmbeddings=NULL,inputPhenoData=NULL,expData=NULL,run_harmony=F,organism,batch_variable="anno_batch",addAnno=F,collapse_datasets=T,minCellCountThr=4,analysis_seed=1,extendedMode=F,umap.method='umap-learn',L2Norm=T,mergePseudocells=T,hierarchical_refinement=T,colNormalize=F,merging_strength=0.15){
  #collapse_datasets=F;minCellCountThr=4;analysis_seed=1;extendedMode=F;umap.method='umap-learn'
  require(purrr)
  require(furrr)
  
  options(future.globals.maxSize= 750*1024^4)
  
  res_embeddings=.sconline.embeddingFn(argList,inputEmbeddings=inputEmbeddings,run_harmony=run_harmony,pd=inputPhenoData,inputBatchCol=batch_variable)
  
  res_umap=.sconline.umapFn(argList,umap.method=umap.method)
  
  
  set.seed(analysis_seed)
  supportingFractionThr=argList$DE_supportingFractionThr
  n.adaptiveKernel=argList$DE_n.adaptiveKernel
  nPropIter=argList$DE_nPropIter
  
  #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("pca_centroids",argList=argList))
  
  
  harmony_embeddings=harmony_embeddings[match(row.names(pd),row.names(harmony_embeddings)),,drop=F]
  if(sum(is.na(harmony_embeddings))>0){
    stop("Error in matching Names!")
  }
  
  
  
  if(is.null(expData)){
    tmp=qread(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
    expData=SplitObject(tmp, split.by = "anno_batch")
  } else {
    data=.extra_sconline.exp_creatorFn(argList=argList,inputExpData = expData,batch_variable = batch_variable,organism = organism,addAnno=addAnno)
    expData = data$data
    
    tmp_pd=NULL
    for(i in 1:length(expData)){
      tmp_pd=rbind(tmp_pd,data.frame(sample=colnames(expData[[i]]),anno_batch=expData[[i]]$anno_batch,stringsAsFactors = F))
    }
    #pd=pd[row.names(pd) %in% tmp_pd$sample,]
    tmp_pd=tmp_pd[match(row.names(pd),tmp_pd$sample),]
    if(sum(is.na(tmp_pd$anno_batch))>0){
      stop("Error! phenoData doesn't match with the expression data")
    }
    pd$anno_batch=tmp_pd$anno_batch
  }
  
  pcaList=split(as.data.frame(harmony_embeddings),pd$anno_batch)
  
  
  dataArranged = parallel::mclapply(expData, function(thisExp){
    thisExp=thisExp[,colnames(thisExp) %in% row.names(pd)]
    
    return(list(
      dsName=as.character(thisExp$anno_batch[1]),
      pcaData=pcaList[[as.character(thisExp$anno_batch[1])]]))
  })
  
  if(sum(unlist(lapply(dataArranged,function(x) nrow(x$pcaData))))<sum(unlist(lapply(expData,ncol)))){
    warning("Annotation was found for only a subset of exp data, that subset was only used in the analysis!")
  }
  
  
  #if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="prop_mat_list_step2",uniformImportant=T),argList)))){
  cat("           Constructing the propagation matrices ...\n")
  dataArranged_old = dataArranged
  dataArranged = dataArranged %>% keep(~dim(.$pcaData)[[1]] >= argList$min_ds_size)
  print(paste0("Going from data size ", length(dataArranged_old), " to ",
               length(dataArranged)))
  rm(dataArranged_old)
  
  if(argList$do.split.prop){
    res=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final_v2_split_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=argList$prop.n.neighbors,mc.cores = argList$ncores)
  } else {
    #centroidPCAdata=pca_centroid;exCentroids=NULL;runIndx=1;n.neighbors=argList$prop.n.neighbors;batchPCAdata=harmony_embeddings;n.trees=50;NNmethod="annoy";L2Norm=T;mergePseudocells=T;hierarchical_refinement=T;batch_variable="anno_batch"
    res=.myConcensusDEFn_step2_detail_newprop3_final_v14(dataArranged=dataArranged,centroidPCAdata=pca_centroid,argList=argList,exCentroids=NULL,runIndx=1,batchPCAdata=harmony_embeddings,n.neighbors=argList$prop.n.neighbors,L2Norm=L2Norm,mergePseudocells=T,batch_variable=batch_variable,colNormalize=F,hierarchical_refinement=T,merging_strength=merging_strength)
  }
  
  cat("           Performing purity analysis ...\n")
  #res_prop=res;annoCol="anno_orig_cellState";return_plot=T;min_effective_size=5
  p=.sconline.anno2pseudocell_tmp(res_prop=res$dataArranged,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=collapse_datasets,return_plot=T,min_effective_size=5)
  p$prop_mat=res$prop_mat
  return(p)
}

.myConcensusDEFn_step2=function(argList,expData=NULL,minCellCountThr=4,meta_method="Stouffer",addClusteringModule=F,analysis_seed=1,extendedMode=F,L2Norm=T,mergePseudocells=T,merging_strength=0.3,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2){
  
  #argList=.ArgList;inputEmbeddings=NULL;inputPhenoData=NULL;inputExpData=NULL;organism="Mouse";batch_variable="anno_batch";run_harmony=F;addAnno=F;addClusteringModule=F;L2Norm=T;mergePseudocells=T;generateUMAP=F;extendedMode=F;merging_strength=0.3
  
  #minCellCountThr=4;meta_method="Stouffer";addClusteringModule=F;analysis_seed=1;extendedMode=F;L2Norm=T;mergePseudocells=T;merging_strength=0.3;analysis_seed=1;analysis_seed=1;include.singletons=T;colNormalize=T
  #expData=NULL
  
  
  myMetaAnalysisFn=function(res,argList,prop_mat,secondRun=F,extendedMode=F){
    
    my_pb <- function(total){
      progress::progress_bar$new(
        format = " [:bar] :current/:total (:percent) eta: :eta (already :elapsed)",
        total = total, clear = FALSE, width= 80)
    }
    
    myStoufferFn_overheadissue=function(zscore_array,weight_array,argList){
      weight_array=sqrt(weight_array)
      z_w=apply(zscore_array*weight_array, c(2,3), sum)
      
      w2=matrix(0,nrow=nrow(z_w),ncol=ncol(z_w))
      
      myStoufferWeightFn=function(i,zscore_array,weight_array){
        
        res_stouffer=rep(1,dim(zscore_array)[2])
        tmp_zscore=t(zscore_array[,i,])
        
        empty_ind=which(colSums(tmp_zscore)==0)
        if(length(empty_ind)<(ncol(tmp_zscore))){
          tmp_zscore[tmp_zscore==0]=NA
          tmp_zscore_cor=cor(tmp_zscore,method = "spearman",use="pairwise.complete.obs")
          tmp_zscore_cor[is.na(tmp_zscore_cor)]=0
          tmp_zscore_cor[which(tmp_zscore_cor<=0)]=0
          tmp_w2=unlist(lapply(1:dim(weight_array)[3],function(j) {
            sum(rowSums(Matrix::crossprod(t(weight_array[,i,j]),weight_array[,i,j])*tmp_zscore_cor))
          }))
          res_stouffer=tmp_w2
          
        }
        return(res_stouffer)
      }
      
      #w2=parallel::mclapply(1:dim(zscore_array)[2],myStoufferWeightFn,zscore_array=zscore_array,weight_array=weight_array,mc.cores = argList$ncores)
      w2=lapply(1:dim(zscore_array)[2],myStoufferWeightFn,zscore_array=zscore_array,weight_array=weight_array)
      w2=do.call("rbind",w2)
      w2=sqrt(w2)
      if(sum(w2>0&w2<0.99)>0){
        #warning("Possible issue in the weights")
        w2[which(w2>0&w2<1)]=1
      }
      z_w=z_w/w2
      z_w[which(w2==0)]=0
      return(z_w)
    }
    
    myStoufferFn_org=function(zscore_array,weight_array,argList){
      weight_array=sqrt(weight_array)
      z_w=apply(zscore_array*weight_array, c(2,3), sum)
      
      w2=matrix(0,nrow=nrow(z_w),ncol=ncol(z_w))
      
      myStoufferWeightFn=function(i,zscore_array,weight_array){
        
        res_stouffer=rep(1,dim(zscore_array)[2])
        tmp_zscore=t(zscore_array[,i,])
        
        empty_ind=which(colSums(tmp_zscore)==0)
        if(length(empty_ind)<(ncol(tmp_zscore))){
          tmp_zscore[tmp_zscore==0]=NA
          tmp_zscore_cor=cor(tmp_zscore,method = "spearman",use="pairwise.complete.obs")
          tmp_zscore_cor[is.na(tmp_zscore_cor)]=0
          tmp_zscore_cor[which(tmp_zscore_cor<=0)]=0
          tmp_w2=unlist(lapply(1:dim(weight_array)[3],function(j) {
            sum(rowSums(Matrix::crossprod(t(weight_array[,i,j]),weight_array[,i,j])*tmp_zscore_cor))
          }))
          res_stouffer=tmp_w2
          
        }
        return(res_stouffer)
      }
      
      #w2t=parallel::mclapply(1:dim(zscore_array)[2],myStoufferWeightFn,zscore_array=zscore_array,weight_array=weight_array,mc.cores = argList$ncores)
      w2=list()
      windowsize=100
      if(dim(zscore_array)[2]-max(seq(1,dim(zscore_array)[2],100))<2){
        windowsize=101
      }
      for(i in seq(1,dim(zscore_array)[2],windowsize)){
        tmp=parallel::mclapply(1:length(i:min(i+windowsize-1,dim(zscore_array)[2])),myStoufferWeightFn,zscore_array=zscore_array[,i:min(i+windowsize-1,dim(zscore_array)[2]),],weight_array=weight_array[,i:min(i+windowsize-1,dim(zscore_array)[2]),],mc.cores = argList$ncores)
        w2=c(w2,tmp)
      }
      w2=do.call("rbind",w2)
      w2=sqrt(w2)
      if(sum(w2>0&w2<0.99)>0){
        #warning("Possible issue in the weights")
        w2[which(w2>0&w2<1)]=1
      }
      z_w=z_w/w2
      z_w[which(w2==0)]=0
      return(z_w)
    }
    
    myStoufferFn=function(zscore_array,weight_array,argList){
      
      
      for(i in 1:length(weight_array)){
        weight_array[[i]]=sqrt(weight_array[[i]])
      }
      
      z_w=zscore_array[[1]]*weight_array[[1]]
      for(i in 2:length(zscore_array)){
        z_w=z_w+zscore_array[[i]]*weight_array[[i]]
      }
      
      
      
      myStoufferWeightFn=function(i,zscore_array,weight_array){
        
        res_stouffer=rep(1,dim(zscore_array)[2])
        tmp_zscore=t(zscore_array[,i,])
        
        empty_ind=which(colSums(tmp_zscore)==0)
        if(length(empty_ind)<(ncol(tmp_zscore))){
          tmp_zscore[tmp_zscore==0]=NA
          tmp_zscore_cor=cor(tmp_zscore,method = "spearman",use="pairwise.complete.obs")
          tmp_zscore_cor[is.na(tmp_zscore_cor)]=0
          tmp_zscore_cor[which(tmp_zscore_cor<=0)]=0
          tmp_w2=unlist(lapply(1:dim(weight_array)[3],function(j) {
            sum(rowSums(Matrix::crossprod(t(weight_array[,i,j]),weight_array[,i,j])*tmp_zscore_cor))
          }))
          res_stouffer=tmp_w2
          
        }
        return(res_stouffer)
      }
      
      zscore_array2=array(0,dim = c(length(zscore_array),nrow(zscore_array[[1]]),ncol(zscore_array[[2]])))
      weight_array2=array(0,dim = c(length(weight_array),nrow(weight_array[[1]]),ncol(weight_array[[2]])))
      
      for(itr in 1:length(zscore_array)){
        zscore_array2[itr,,]=as.matrix(zscore_array[[itr]])
      }
      for(itr in 1:length(weight_array)){
        weight_array2[itr,,]=as.matrix(weight_array[[itr]])
      }
      
      rm(zscore_array,weight_array)
      gc()
      
      #w2t=parallel::mclapply(1:dim(zscore_array)[2],myStoufferWeightFn,zscore_array=zscore_array,weight_array=weight_array,mc.cores = argList$ncores)
      w2=list()
      windowsize=100
      if(dim(zscore_array2)[2]-max(seq(1,dim(zscore_array2)[2],100))<2){
        windowsize=101
      }
      for(i in seq(1,dim(zscore_array2)[2],windowsize)){
        
        tmp=parallel::mclapply(1:length(i:min(i+windowsize-1,dim(zscore_array2)[2])),myStoufferWeightFn,zscore_array=zscore_array2[,i:min(i+windowsize-1,dim(zscore_array2)[2]),],weight_array=weight_array2[,i:min(i+windowsize-1,dim(zscore_array2)[2]),],mc.cores = argList$ncores)
        w2=c(w2,tmp)
      }
      w2=do.call("rbind",w2)
      w2=sqrt(w2)
      if(sum(w2>0&w2<0.99)>0){
        #warning("Possible issue in the weights")
        w2[which(w2>0&w2<1)]=1
      }
      z_w=z_w/w2
      z_w[which(w2==0)]=0
      return(z_w)
    }
    
    
    myWeightedMeanFn=function(data_array,weight_array){
      z_w=apply(data_array*weight_array, c(2,3), sum)
      w=apply(weight_array, c(2,3), sum)
      if(sum(w>0&w<0.99)>0){
        #warning("Possible issue in the weights")
        w[which(w>0&w<1)]=1
      }
      z_w=z_w/w
      z_w[which(w==0)]=0
      return(z_w)
    }
    
    myPctAvgfn_org=function(pct_array,centroid_weight_array,effective_size_array){
      myWeightedMeanFn=function(data_array,weight_array){
        z_w=apply(data_array*weight_array, c(2,3), function(x) sum(x,na.rm = T))
        w=apply(weight_array, c(2,3), function(x) sum(x,na.rm = T))
        if(sum(w>0&w<0.99)>0){
          #warning("Possible issue in the weights")
          w[which(w>0&w<1)]=1
        }
        z_w=z_w/w
        z_w[which(w==0)]=0
        return(z_w)
      }
      
      effective_size_array=centroid_weight_array*effective_size_array
      #pct_array[centroid_weight_array<0.9]=NA
      med_pct.1=myWeightedMeanFn(data_array=pct_array,weight_array=effective_size_array)
      row.names(med_pct.1)=dimnames(pct_array)[[2]]
      return(med_pct.1)
      
    }
    
    myPctAvgfn=function(pct_array,effective_size_array){
      
      all(names(pct_array)==names(effective_size_array))
      
      z_w=pct_array[[1]]*effective_size_array[[1]]
      for(itr in 2:length(pct_array)){
        z_w=z_w+pct_array[[itr]]*effective_size_array[[itr]]
      }
      
      w=effective_size_array[[1]]
      for(witr in 2:length(effective_size_array)){
        w=w+effective_size_array[[witr]]
      }
      if(sum(w@x>0&w@x<0.99)>0){
        #warning("Possible issue in the weights")
        w@x[which(w@x>0&w@x<1)]=1
      }
      w_ident=w
      w_ident@x=rep(1,length(w_ident@x))
      z_w=Matrix::drop0(z_w*w_ident)
      z_w@x=z_w@x/w@x
      
      return(z_w)
      
    }
    
    
    myHarmonicPvalFn=function(zscore_array,weight_array,offset=10^(-25)){
      pval=pnorm(zscore_array,lower.tail = F)
      weight_normFactor=apply(weight_array,c(2,3),sum)
      weight_array_n=sweep(weight_array,c(2,3),weight_normFactor,"/")
      
      pval[which(pval<0.1)]=pval[which(pval<0.1)]+offset
      pval[which(pval>0.9)]=pval[which(pval>0.9)]-offset
      resHarmony=apply(weight_array_n/pval,c(2,3),sum)
      resHarmony=1/resHarmony
      resHarmony[is.na(resHarmony)]=0.5
      resHarmony=pmin(resHarmony,1-offset)
      resHarmony=qnorm(resHarmony,lower.tail = F)
      
      return(resHarmony)
    }
    
    library(magrittr)
    library(purrr)
    library(stringr)
    library(future)
    library(furrr)
    
    if(F){
      #Jonah's
      exp_final_fp = file.path(argList$saveDir, "detail_exp_final")
      all_exp_files = list.files(exp_final_fp)
      
      all_exp_files %<>% set_names(., str_remove(., "[.]RDS"))
      
      my_pb <- function(total){
        progress::progress_bar$new(
          format = " [:bar] :current/:total (:percent) eta: :eta (already :elapsed)",
          total = total, clear = FALSE, width= 80)
      }
      this_pb = my_pb(length(all_exp_files))
      res_exp_final = map(all_exp_files, ~{
        this_pb$tick()
        # don't do mc when testing bc doesn't like swap space
        return(.mcreadRDS(file.path(exp_final_fp, .)))
      })
      this_pb$terminate()
      res = res_exp_final
      rm(res_exp_final)
    }
    
    gc()
    plan(sequential)
    
    
    cat("           Arranging z-score results ...\n")
    geneList=list()
    for(gitr in 1:length(res)){
      geneList=c(geneList,list(data.frame(gene=colnames(res[[gitr]]$matWeights),weight=colMaxs(res[[gitr]]$matWeights),stringsAsFactors = F)))
    }
    
    geneList=as.data.frame(data.table::rbindlist(geneList))
    geneList=aggregate(weight~gene,data=geneList,function(x) sum(x>0.9))
    geneList=geneList[which(geneList$weight>(length(res)*argList$DE_supportingFractionThr)),]
    geneList=as.character(geneList$gene)
    
    divby = ceiling(length(res)/50)
    this_pb = my_pb(ceiling(length(res)/divby))
    for(gitr in 1:length(res)){
      if(gitr %% divby == 0){
        this_pb$tick()
      }
      
      set_name_list=setdiff(names(res[[gitr]]),c("dsName","prop_mat","pseudocell_sim_mat","gene_name_list","pseudocell_name_list"))
      if(!secondRun){
        set_name_list=c("pct.1","zscore","matWeights","matEffectiveSize")
      }
      
      for(dsitr in set_name_list){
        tmp=res[[gitr]][[dsitr]]
        tmp_c=setdiff(geneList,colnames(tmp))
        if(length(tmp_c)>0){
          tmp_c=as(matrix(0,nrow=nrow(tmp),ncol=length(tmp_c),dimnames = list(row.names(tmp),tmp_c)),"dgCMatrix")
          tmp=cbind(tmp,tmp_c)
        }
        tmp=tmp[,match(geneList,colnames(tmp))]
        colnames(tmp)=geneList
        
        res[[gitr]][[dsitr]]=as(tmp,"dgCMatrix")
      }
      
      if(gitr %% 10 == 0){
        gc()
      }
    }
    this_pb$terminate()
    
    if(secondRun){
      matWeights=data.frame(centroid=row.names(prop_mat),weight=apply(res[[1]]$matWeights,1,max),effectiveSize=apply(res[[1]]$matEffectiveSize,1,max),dataset=res[[1]]$dsName,stringsAsFactors = F)
      
      if(length(res)>1){
        for(idsScore in 2:length(res)){
          tmp=data.frame(centroid=row.names(prop_mat),weight=apply(res[[idsScore]]$matWeights,1,max),effectiveSize=apply(res[[idsScore]]$matEffectiveSize,1,max),dataset=res[[idsScore]]$dsName,stringsAsFactors = F)
          matWeights=rbind(matWeights,tmp)
        }
      }
      
      matEffectiveSize=reshape2::dcast(centroid~dataset,data = matWeights,value.var = "effectiveSize")
      matWeights=reshape2::dcast(centroid~dataset,data = matWeights,value.var = "weight")
      resDistances=list(matWeights=matWeights,matEffectiveSize=matEffectiveSize)
      save(resDistances,file=do.call('.myFilePathMakerFn',args=c(list(filename="centroid_scores",uniformImportant=T,propImportant=T,qsFormat=T),argList)))
      rm(matEffectiveSize,matWeights,resDistances)
    }
    
    if(extendedMode){
      z_q.8=apply(res_array$zscore,c(2,3),function(x) {
        if(sum(x>3)>0){
          sl_thr=quantile(x[x!=0],0.8)
          sl_ind=which(x>sl_thr)
          x[sl_ind]=sl_thr
        } 
        if(sum(x<(-3))>0){
          sl_thr=quantile(x[x!=0],0.2)
          sl_ind=which(x<sl_thr)
          x[sl_ind]=sl_thr
        }
        x
      })
    }
    
    
    #rearranging the res object
    res_rearranged=NULL
    for(i in setdiff(names(res[[1]]),"dsName")){
      tmp=list()
      for(j in 1:length(res)){
        tmp=c(tmp,list(res[[j]][[i]]))
        names(tmp)[length(tmp)]=res[[j]]$dsName
      }
      res_rearranged=c(res_rearranged,list(tmp))
      names(res_rearranged)[length(res_rearranged)]=i
    }
    
    if(secondRun){
      for(i in names(res_rearranged)){
        tmp=res_rearranged[[i]]
        qsave(tmp,file=.myFilePathMakerFn(paste0("res_array_",i),argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
      }
    }
    
    res_rearranged=res_rearranged[c("zscore","pct.1","matWeights","matEffectiveSize")]
    rm(res)
    
    res_eff_size=list()
    for(i in 1:length(res_rearranged[["matWeights"]])){
      res_eff_size=c(res_eff_size,list(as(res_rearranged$matWeights[[i]]*res_rearranged$matEffectiveSize[[i]],"dgCMatrix")))
      names(res_eff_size)[length(res_eff_size)]=names(res_rearranged$matEffectiveSize)[i]
    }
    
    res_rearranged$matEffectiveSize=res_eff_size
    rm(res_eff_size)
    res_rearranged=res_rearranged[c("zscore","pct.1","matEffectiveSize")]
    
    gc()
    cat("           Calculating meta pct.1 & pct.2 values ...\n")
    med_pct.1=myPctAvgfn(pct_array=res_rearranged$pct.1,effective_size_array=res_rearranged$matEffectiveSize)
    med_pct.2=NULL
    med_logFC=NULL
    res_rearranged$pct.1=NULL
    
    cat("           Calculating meta z-scores ...\n")
    if(meta_method=="Stouffer"){
      meta_z=myStoufferFn(zscore_array = res_rearranged$zscore,weight_array = res_rearranged$matEffectiveSize,argList = argList)
      
      if(extendedMode){
        meta_z.8=myStoufferFn(zscore_array = z_q.8,weight_array = res_array$matEffectiveSize,argList = argList)
      }
      
    } else if(meta_method=="Harmony"){
      meta_z=myHarmonicPvalFn(zscore_array = res_array$zscore,weight_array = res_array$matWeights)
      meta_z[which(meta_z<(-10))]=(-10)
      if(extendedMode){
        meta_z.8=myHarmonicPvalFn(zscore_array = z_q.8,weight_array = res_array$matWeights)
        meta_z.8[which(meta_z.8<(-10))]=(-10)
      }
    } else if(meta_method=="average"){
      meta_z=myPctAvgfn(pct_array=res_array$zscore,centroid_weight_array=res_array$matWeights,effective_size_array=res_array$matEffectiveSize)
    }
    
    return(list(med_pct.1=med_pct.1,med_pct.2=med_pct.2,med_logFC=med_logFC,meta_z=meta_z))
    
  }
  
  require(qs)
  require(purrr)
  require(furrr)
  
  reRunCheck=T
  if(file.exists(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T,propImportant = T))&file.exists(.myFilePathMakerFn("res_DE_wZscore_pathwayAnalysis",argList=argList,uniformImportant=T,propImportant = T))&!argList$newRun){
    reRunCheck=tryCatch({load(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T,propImportant = T));F}, error=function(e) {return(T)})
    if(!reRunCheck){
      reRunCheck=tryCatch({load(.myFilePathMakerFn("res_DE_wZscore_pathwayAnalysis",argList=argList,uniformImportant=T,propImportant = T));F}, error=function(e) {return(T)})
      
    }
    if(!reRunCheck){
      reRunCheck=tryCatch({tst=qs::qread(.myFilePathMakerFn("res_meta",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T));F}, error=function(e) {return(T)})
    }
    
  }
  
  res=""
  if(reRunCheck){
    
    set.seed(analysis_seed)
    supportingFractionThr=argList$DE_supportingFractionThr
    n.adaptiveKernel=argList$DE_n.adaptiveKernel
    nPropIter=argList$DE_nPropIter
    
    #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("pca_centroids",argList=argList))
    
    
    harmony_embeddings=harmony_embeddings[match(row.names(pd),row.names(harmony_embeddings)),,drop=F]
    if(sum(is.na(harmony_embeddings))>0){
      stop("Error in matching Names!")
    }
    
    if(is.null(expData)){
      if(!file.exists(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))){
        stop("Expression data is missing!")
      }
      tmp=qread(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
      expData=SplitObject(tmp, split.by = "anno_batch")
    } else {
      tmp_pd=NULL
      for(i in 1:length(expData)){
        tmp_pd=rbind(tmp_pd,data.frame(sample=colnames(expData[[i]]),anno_batch=expData[[i]]$anno_batch,stringsAsFactors = F))
      }
      #pd=pd[row.names(pd) %in% tmp_pd$sample,]
      tmp_pd=tmp_pd[match(row.names(pd),tmp_pd$sample),]
      if(sum(is.na(tmp_pd$anno_batch))>0){
        stop("Error! phenoData doesn't match with the expression data")
      }
      pd$anno_batch=tmp_pd$anno_batch
    }
    
    pcaList=split(as.data.frame(harmony_embeddings),pd$anno_batch)
    
    if(argList$do.split.prop&F){
      dataArranged = parallel::mclapply(expData, function(thisExp){
        thisExp=thisExp[,colnames(thisExp) %in% row.names(pd)]
        if(class(thisExp)=="SingleCellExperiment"){
          tmp=counts(thisExp)
          #tmp=.extraExport2SeuratFn(thisExp)
          #plan("sequential")
          # plan("multicore", workers = min(parallel::detectCores(), 8, na.rm=T))
          #library(Seurat)
          #normalExp=Seurat::NormalizeData(tmp,verbose=F) # TODO JONAH Maybe not sequential
          # expData[[i]]=tmp
        }else{
          #normalExp = thisExp
          tmp=thisExp@assays$RNA@counts
        }
        
        #tmp=Seurat::NormalizeData(normalExp,normalization.method ="RC",verbose=F)
        
        return(list(
          countData=tmp,
          #logNormData=normalExp@assays$RNA@data,
          #countData=normalExp@assays$RNA@counts,
          #scaleData=tmp@assays$RNA@data,
          dsName=as.character(thisExp$anno_batch[1]),
          pcaData=pcaList[[as.character(thisExp$anno_batch[1])]]))
      })
    } else {
      dataArranged = parallel::mclapply(expData, function(thisExp){
        thisExp=thisExp[,colnames(thisExp) %in% row.names(pd)]
        if(class(thisExp)=="SingleCellExperiment"){
          tmp=counts(thisExp)
          #tmp=.extraExport2SeuratFn(thisExp)
          #plan("sequential")
          # plan("multicore", workers = min(parallel::detectCores(), 8, na.rm=T))
          #library(Seurat)
          #normalExp=Seurat::NormalizeData(tmp,verbose=F) # TODO JONAH Maybe not sequential
          # expData[[i]]=tmp
        }else{
          #normalExp = thisExp
          tmp=thisExp@assays$RNA@counts
        }
        
        #tmp=Seurat::NormalizeData(normalExp,normalization.method ="RC",verbose=F)
        
        return(list(
          countData=tmp,
          dsName=as.character(thisExp$anno_batch[1])))
      })
    }
    
    
    if(sum(unlist(lapply(dataArranged,function(x) ncol(x$countData))))<sum(unlist(lapply(expData,ncol)))){
      warning("Annotation was found for only a subset of exp data, that subset was only used in the analysis!")
    }
    
    res_fd=NULL
    for(i in 1:length(expData)){
      .printDot()
      if(class(expData[[i]])[1]=="SingleCellExperiment"){
        fd=as.data.frame(rowData(expData[[i]]))
      } else {
        fd=as.data.frame(expData[[i]]@assays$RNA@meta.features)
      }
      
      if(sum(colnames(fd)=="ensembl_gene_id")==0){
        fd$ensembl_gene_id=row.names(fd)
      }
      fd$ensembl_gene_id=gsub("_","-",fd$ensembl_gene_id)
      slCols=intersect(colnames(fd),c("gene_name","gene_biotype","symbol","gene_short_name","ensembl_gene_id"))
      fd=fd[,colnames(fd) %in% slCols,drop=F]
      if(sum(is.na(fd$ensembl_gene_id))>0){
        fd$ensembl_gene_id[is.na(fd$ensembl_gene_id)]=row.names(fd)[is.na(fd$ensembl_gene_id)]
      }
      
      if(!is.null(res_fd)){
        fd=fd[!fd$ensembl_gene_id %in% res_fd$ensembl_gene_id,,drop=F]
        if(nrow(fd)>0){
          res_fd=rbind(res_fd,fd)
        }
      } else {
        res_fd=fd
      }
    }
    fd=res_fd
    rm(res_fd,expData,pcaList,harmony_embeddings)
    gc()
    
    #if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="prop_mat_list_step2",uniformImportant=T),argList)))){
    cat("           Constructing the propagation matrices ...\n")
    dataArranged_old = dataArranged
    dataArranged = dataArranged %>% keep(~dim(.$countData)[[2]] >= argList$min_ds_size)
    if(length(dataArranged_old)!=length(dataArranged)){
      cat(paste0("Excluding datasets with less than ", argList$min_ds_size," cells => Retaining ", length(dataArranged), " out of ",
                   length(dataArranged_old)," datasets\n"))
    }
    rm(dataArranged_old)
    
    #centroidPCAdata=pca_centroid;nPropIter=1;n.neighbors=argList$prop.n.neighbors
    
    if(sum(is.na(pca_centroid))>0){
      stop("Error: NA values were found in the pca of pseudocells!")
    }
    
    if(argList$do.split.prop&F){
      res=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final_v2_split_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=argList$prop.n.neighbors,mc.cores = argList$ncores)
    } else {
      #batchPCAdata=NULL;dataArranged=dataArranged;centroidPCAdata=pca_centroid;argList=argList;exCentroids=NULL;runIndx=1;n.neighbors=argList$prop.n.neighbors;n.trees=50;NNmethod="annoy"
      res=.myConcensusDEFn_step2_detail_newprop3_final_v14(dataArranged=dataArranged,centroidPCAdata=pca_centroid,argList=argList,exCentroids=NULL,runIndx=1,batchPCAdata=NULL,n.neighbors=argList$prop.n.neighbors,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength)
    }
    
    #res_limited=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final_v2_limited,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=argList$prop.n.neighbors,mc.cores = argList$ncores)
    
    #res1=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=argList$prop.n.neighbors,mc.cores = argList$ncores)
    
    #Evaluating the performance of the propagation step
    #res_prop=res;argList=argList;annoCol="anno_orig_cellState";collapse_datasets=F;return_plot=T
    #p=.sconline.anno2pseudocell_tmp(res_prop=res,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T,min_effective_size=5)
    #ggsave(plot=p$plot,file="~/myBucket/torm.pdf",width=20,height=20)
    #summary(.myEffSizePropMat(res[[1]]$prop_mat)$effective_sample_size/.myEffSizePropMat(res1[[1]]$prop_mat)$effective_sample_size)
    #tst_res=.myEffSizePropMat(res[[1]]$prop_mat)$effective_sample_size
    #tst_res_limited=.myEffSizePropMat(res_limited[[1]]$prop_mat)$effective_sample_size
    
    if(F){
      p=.sconline.anno2pseudocell_tmp(res_prop=res,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T,min_effective_size=5)
      
      tst=p$results
      tst=apply(tst[,-c(1,ncol(tst))],1,function(x){
        y=colnames(tst)[-c(1,ncol(tst))]
        y=y[which(x==max(x))[1]]
      })
      tst=which(grepl("purkinje",tolower(tst)))
      
      load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
      
      UMAP_centroid=.sconline.fetch_data("umap_centroids",argList=argList)
      
      pd_summary=.sconline.fetch_data("annotation",argList=argList)
      
      pd_summary=pd_summary[grepl("purkinje",tolower(pd_summary$anno_orig_cellState)),]
      
      inputData=p$results[tst,]
      piechart_data=merge(inputData,UMAP_centroid,by.x="pseudocell",by.y="centroid",all=T)
      
      anno_cols=setdiff(colnames(inputData),c("pseudocell","effective_size"))
      pd_summary$anno_orig_cellState=factor(as.character(pd_summary$anno_orig_cellState),levels=anno_cols)
      p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=anno_orig_cellState),size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                cols=anno_cols,r=2) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(length(anno_cols))))+scale_color_manual(values=c("gray",hues::iwanthue(length(anno_cols))))
      
      ggsave(plot=p,file="~/myBucket/torm.pdf",width=12,height = 12)
       
    }
    
    prop_mat=res$prop_mat
    res=res[["dataArranged"]]
    
    tmpValCheck=(unlist(lapply(res,length)))
    if(sum(tmpValCheck==1)>0){
      stop(res[[which(tmpValCheck==1)[1]]])
    }
    rm(tmpValCheck)
    
    ###Checking for outlier centroids
    
    
    
    
    cat("           Calculating dataset specific z-scores ...\n")
    #tst=.myConcensusDEFn_step2_detail_exp_final(res[[4]],argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode)
    #dim(tst$prop_mat)
    
    res2=parallel::mclapply(res,.myConcensusDEFn_step2_detail_exp_final,argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode,mc.preschedule=F,mc.cores = argList$ncores)
    while(sum(unlist(lapply(res2,class))=="NULL")>0){
      gc()
      null_ind=which(unlist(lapply(res2,class))=="NULL")
      cat(paste0("           Redoing the z-score calculations for ",length(null_ind)," datasets ...\n"))
      res3=parallel::mclapply(res[null_ind],.myConcensusDEFn_step2_detail_exp_final,argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode,mc.preschedule=F,mc.cores = 3)
      res2[null_ind]=res3
      rm(res3)
    }
    res=res2
    rm(res2)
    gc()
    
    #merging pseudocells in each batch
    
    
    
    
    tmpValCheck=(unlist(lapply(res,length)))
    if(sum(tmpValCheck==1)>0){
      stop(res[[which(tmpValCheck==1)[1]]])
    }
    rm(tmpValCheck)
    
    #qsave(dataArranged,file="~/torm_dataArranged.qs")
    res_meta=myMetaAnalysisFn(res=res,prop_mat=prop_mat,argList = argList,secondRun=F,extendedMode=F)
    
    med_pct.1=res_meta$med_pct.1
    #med_pct.2=res_meta$med_pct.2
    #med_logFC=res_meta$med_logFC
    meta_z=res_meta$meta_z
    #res_array=res_meta$res_array
    rm(res_meta,res)
    gc()
    
    
    #three_d_image = tf.concat(0, [[r], [g], [b]])
    #three_d_image = tf.stack([r,g,b], axis=2)
    
    
    ##################
    
    #MetaQC: initial implementation
    
    if(F){
      require(GeneMeta)
      meta_z=matrix(0,nrow=dim(res_array$t_adj)[2],ncol=dim(res_array$t_adj)[3])
      for(i in 1:dim(res_array$t_adj)[2]){
        tmp_t=t(res_array$t_adj[,i,])
        tmp_sigmad=t(res_array$sigmad_adj[,i,])
        
        empty_ind=which(colSums(tmp_sigmad)==0|colSums(tmp_t)==0)
        if(length(empty_ind)<(ncol(tmp_t)-1)){
          if(length(empty_ind)>0){
            tmp_t=tmp_t[,-empty_ind]
            tmp_sigmad=tmp_sigmad[,-empty_ind]
          }
          
          my.Q   <- f.Q(tmp_t, tmp_sigmad)
          
          my.tau2.DL<-tau2.DL(my.Q, ncol(tmp_t), my.weights=1/tmp_sigmad)
          #obtain new variances s^2+tau^2
          myvarsDL <- tmp_sigmad + my.tau2.DL
          #compute
          muREM <- mu.tau2(tmp_t, myvarsDL)
          #cumpute mu(tau)
          varREM <- var.tau2(myvarsDL)
          ZREM <- muREM/sqrt(varREM)
          meta_z[i,]=ZREM
        } else if(length(empty_ind)==(ncol(tmp_t)-1)) {
          meta_z[i,]=res_array$zscore[-empty_ind,i,]
        }
        
        
      }
    }
    
    if(F){
      .mycheck.R=function (R, checksym = TRUE, checkna = TRUE, checkpd = FALSE, 
                           nearpd = FALSE, checkcor = FALSE, checkdiag = TRUE, isbase = TRUE, 
                           k, adjust, fun) {
        if (nearpd) 
          checkpd <- TRUE
        if (checkpd) {
          checkna <- TRUE
          checksym <- TRUE
        }
        if (inherits(R, "dpoMatrix")) 
          R <- as.matrix(R)
        if (inherits(R, "data.frame")) 
          R <- as.matrix(R)
        if (checksym && !(is.matrix(R) && isSymmetric(unname(R)))) 
          stop("Argument 'R' must be a symmetric matrix.", call. = FALSE)
        if (checkna && any(is.na(R))) 
          stop("Values in 'R' must not contain NAs.", call. = FALSE)
        if (checkpd && any(eigen(R)$values <= 0)) {
          if (nearpd) {
            warning("Matrix 'R' is not positive definite. Used Matrix::nearPD() to make 'R' positive definite.", 
                    call. = FALSE)
            R <- as.matrix(.find.nonegmat(R))
          }
          else {
            stop("Matrix 'R' is not positive definite.", call. = FALSE)
          }
        }
        if (checkcor && any(abs(R) > 1, na.rm = TRUE)) 
          stop("Argument 'R' must be a correlation matrix, but contains values outside [-1,1].", 
               call. = FALSE)
        if (checkdiag && (any(is.na(diag(R))) || any(diag(R) != 1))) 
          stop("Diagonal values in 'R' must all be equal to 1.", 
               call. = FALSE)
        if (isbase) {
          if (k != nrow(R)) 
            stop("Length of 'p' vector (", k, ") does not match the dimensions of the 'R' matrix (", 
                 nrow(R), ",", ncol(R), ").", call. = FALSE)
          if (adjust == "none") 
            warning("Argument 'R' was specified, but no adjustment method was chosen via the 'adjust' argument.\nTo account for dependence, specify an adjustment method. See help(", 
                    fun, ") for details.", call. = FALSE)
          if (adjust == "user") 
            warning("When 'm' is specified, argument 'R' is irrelevant and ignored.")
        }
        return(R)
      }
      
      .myMefffn=function (R, eigen, method="nyholt", ...) {
        .is.numeric.vector=function (x){
          is.atomic(x) && is.numeric(x) && !is.matrix(x) && !is.null(x)
        }
        
        
        method <- match.arg(method, c("nyholt", "liji", "gao", "galwey"))
        if (missing(eigen)) {
          if (missing(R)) 
            stop("Argument 'R' must be specified.", call. = FALSE)
          R <- .mycheck.R(R, checksym = TRUE, checkna = TRUE, checkpd = FALSE, 
                          nearpd = FALSE, checkcor = TRUE, checkdiag = TRUE, 
                          isbase = FALSE)
          evs <- base::eigen(R)$values
        } else {
          if (!.is.numeric.vector(eigen)) 
            stop("Argument 'eigen' must be a numeric vector.", 
                 call. = FALSE)
          evs <- eigen
        }
        if (any(evs < 0)) 
          warning(paste0("One or more eigenvalues ", ifelse(missing(eigen), 
                                                            "derived from the 'R' matrix ", ""), "are negative."), 
                  call. = FALSE)
        if (method == "nyholt") {
          k <- length(evs)
          m <- 1 + (k - 1) * (1 - var(evs)/k)
        }
        if (method == "liji") {
          abs.evs <- abs(evs) + sqrt(.Machine$double.eps)
          m <- sum(ifelse(abs.evs >= 1, 1, 0) + (abs.evs - floor(abs.evs)))
        }
        if (method == "gao") {
          ddd <- list(...)
          if (!is.null(ddd$C)) {
            C <- ddd$C
          }
          else {
            C <- 0.995
          }
          if (C < 0 || C >= 1) 
            warning("Value of 'C' should be >= 0 and < 1.", call. = FALSE)
          m <- which(cumsum(sort(evs, decreasing = TRUE))/sum(evs) > 
                       C)[1]
        }
        if (method == "galwey") {
          if (any(evs < 0)) {
            warning(paste0("Negative eigenvalues ", ifelse(missing(eigen), 
                                                           "derived from the 'R' matrix ", ""), "were set to 0."), 
                    call. = FALSE)
            evs[evs < 0] <- 0
          }
          m <- sum(sqrt(evs))^2/sum(evs)
        }
        m <- floor(m)
        return(m)
      }
      
      myStoufferFn_org=function(zscore_array,weight_array){
        z_w=apply(zscore_array*weight_array, c(2,3), sum)
        w2=sqrt(apply(weight_array*weight_array, c(2,3), sum))
        if(sum(w2>0&w2<0.99)>0){
          #warning("Possible issue in the weights")
          w2[which(w2>0&w2<1)]=1
        }
        z_w=z_w/w2
        z_w[which(w2==0)]=0
        return(z_w)
      }
    }
    
    
    if(F){
      #alternate slower solution
      myPctAvgfn=function(pct_array,centroid_weight_array,effective_size_array){
        myWeightedMeanFn=function(data_array,weight_array){
          z_w=multiApply::Apply(data_array*weight_array, 1, function(x) sum(x,na.rm = T))[[1]]
          w=multiApply::Apply(weight_array, 1, function(x) sum(x,na.rm = T))[[1]]
          if(sum(w>0&w<0.99)>0){
            #warning("Possible issue in the weights")
            w[which(w>0&w<1)]=1
          }
          z_w=z_w/w
          z_w[which(w==0)]=0
          return(z_w)
        }
        
        pct_array[centroid_weight_array<0.9]=NA
        med_pct.1=myWeightedMeanFn(data_array=pct_array,weight_array=effective_size_array)
        row.names(med_pct.1)=dimnames(pct_array)[[2]]
        return(med_pct.1)
        
      }
      
      tst=myPctAvgfn(pct_array=res_array$pct.1,centroid_weight_array=res_array$matWeights,effective_size_array=res_array$matEffectiveSize)
      
    }
    
    
    
    #tst=list(pct_mat=med_pct.1,argList = argList,meta_z_mat=meta_z,sig1_thr=3,centers=NULL,pct2_thr=0.3,pct_diff_thr=0.2,symmetric=F)
    #med_pct.1=tst$pct_mat;argList = tst$argList;meta_z=tst$meta_z_mat;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
    
    #pct_mat=med_pct.1;meta_z_mat=meta_z;pct_mat_ref=NULL;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
    
    #################################
    #Needs attention
    
    pct_diff_count=.extra_sconline.PctScoreFn_v2(pct_mat=med_pct.1,argList = argList,meta_z_mat=meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F)
    
    #argList$ncores=3
    #pct_diff_count2=.extra_sconline.PctScoreFn(pct_mat=med_pct.1,argList = argList,meta_z_mat=meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F)
    
    
    diff=prop_mat
    diff@x=rep(1,length(diff@x))
    diff=diff %*% t(diff)
    diff=sweep(diff,1,diag(diff),"/")
    diff=as.matrix(diff)
    diff[diff<0.100001]=0.100001
    diff=abs(log10(diff))
    diff=diff+pct_diff_count
    
    diff=hclust(as.dist(diff),method = "complete")
    diff=cutree(diff,h=0.999)
    
    if(F){
      diff=data.frame(pseudocell=names(diff),cluster=diff,stringsAsFactors = F)
      adj2=prop_mat
      load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
      pd=pd[colnames(adj2),]
      
      table(is.na(pd$anno_orig_cellState))
      tmp= adj2 %*% as.matrix(.myOneHotFn(pd$anno_orig_cellState))
      
      tmp=apply(tmp,1,function(x){
        y=which(x==max(x))[1]
        tmp2=data.frame(cluster=colnames(tmp)[y],purity=x[y])
        tmp2
      })
      
      tmp=do.call("rbind",tmp)
      tmp2=aggregate(purity~cluster,data=tmp,median)
      
      table(tmp2[,2]>0.8)
      
      
      
      diff=merge(diff,data.frame(pseudocell=row.names(tmp),tmp,stringsAsFactors = F),by="pseudocell")
      diff2=aggregate(cluster.y~cluster.x,data=diff,function(x){
        y=as.numeric(table(x))
        max(y)/sum(y)
      })
    }
    
    
    if(sum(duplicated(diff))>0){
      cat(paste0("           Merging pseudocells based on their DE count; retaining ",length(unique(diff))," out of ",length(diff),"\n"))
      diff=data.frame(pseudocell=names(diff),cluster=diff,stringsAsFactors = F)
      diff_id=diff[!duplicated(diff$cluster),]
      diff=merge(diff,diff_id,by="cluster")
      diff2=reshape2::dcast(pseudocell.x~pseudocell.y,data=diff,fun.aggregate=length)
      diff2=diff2[match(row.names(prop_mat),diff2[,1]),]
      diff2=t(as.matrix(diff2[,-1]))
      
      diff2=as(diff2,"dgCMatrix")
      diff2=diff2 %*% prop_mat
      diff2 <- .extra_matrix_rowNorm(diff2)#Matrix::Diagonal(x = 1 / (rowSums(diff2)+0.000000000001)) %*% diff2
      
      .diff2=diff2
      
      diff2=.diff2
      
      rowMeans_drop0 <- function (dgCMat) {
        RowInd <- dgCMat@i + 1
        sapply(split(dgCMat@x, RowInd), function(x)quantile(x,0.95))
      }
      
      diff2 = .extra_matrix_rowNorm(input_mat = diff2,rowValues = 1 / (rowMeans_drop0(diff2)))#Matrix::Diagonal(x = 1 / (rowMeans_drop0(diff2))) %*% diff2
      #diff2@x=pmin(diff2@x,1)
      colMax_vals=c()
      for(i in seq(1,ncol(diff2),5000)){
        tmp_max=as.numeric(qlcMatrix::colMax(diff2[,i:min(i+4999,ncol(diff2))]))
        colMax_vals=c(colMax_vals,as.numeric(tmp_max))
      }
      diff2 =  diff2 %*% Matrix::Diagonal(x = 1 / colMax_vals)
      diff2 <- Matrix::Diagonal(x = 1 / rowSums(diff2)) %*% diff2
      
      prop_mat=diff2
      
      matWeights=.myEffSizePropMat(diff2)
      
      matEffectiveSize=matWeights$effective_sample_size
      matWeights=matWeights$centroid_weights
      matWeights=matWeights[match(row.names(diff2),names(matWeights))]
      matEffectiveSize=matEffectiveSize[match(row.names(diff2),names(matEffectiveSize))]
      
      res=dataArranged
      for(i in 1:length(res)){
        tmp_prop_mat=diff2[,match(colnames(res[[i]]$countData),colnames(diff2))]
        #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
        tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
        tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
        tmp_weights[tmp_effsize<4]=0
        tmp=list(prop_mat=tmp_prop_mat,data=res[[i]],matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
        res[[i]]=tmp
      }
      
      gc()
      plan("sequential")
      cat(paste0("           Redoing the propagation after the merging\n"))
      res=parallel::mclapply(res,.myConcensusDEFn_step2_detail_exp_final,argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode,mc.cores = 4)
      
      cat(paste0("           Redoing the meta analysis of zscores\n"))
      res_meta=myMetaAnalysisFn(res=res,prop_mat=prop_mat,argList = argList,secondRun=T,extendedMode=extendedMode)
      
      med_pct.1=res_meta$med_pct.1
      #med_pct.2=res_meta$med_pct.2
      #med_logFC=res_meta$med_logFC
      meta_z=res_meta$meta_z
      #res_array=res_meta$res_array
      rm(res_meta)
      gc()
    }
    
    
    de_pct_res=NULL
    if(addClusteringModule&F){
      
      cosine_dist=.extra_sconline.CosineDistFn(inputMatrix = meta_z,sig_thr = 2)
      
      de_pct_res=.extra_sconline.PctScoreFn(argList = argList,pct_mat=med_pct.1,meta_z_mat=meta_z,pct_mat_ref=NULL,sig1_thr=3,centers=NULL,pct2_thr=0.3,pct_diff_thr=0.2,symmetric=F)
      
    }
    #meta_z_mat=meta_z$meta_z;cosine_dist=cosine_dist;de_dist=de_pct_res;pct_mat=meta_z$med_pct.1;min_marker_thr=20;sig1_thr=3;sig2_thr=1;pct_de_count_thr=1;pct_diff_thr=0.2;pct2_thr=0.3
    #####################
    
    
      
    #####################
    
    
    if(extendedMode){
      res_meta=list(meta_z=meta_z,meta_z.8=meta_z.8,med_pct.1=med_pct.1,med_pct.2=med_pct.2,med_logFC=med_logFC,med_n=med_n,fd=fd)
    } else {
      res_meta=list(meta_z=meta_z,med_pct.1=med_pct.1,fd=fd)
    }
    
    qsave(res_meta,file=.myFilePathMakerFn("res_meta",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
    rm(res_meta)
    
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
    x_z=apply(meta_z,2,max)
    x_z=meta_z[,which(x_z>2)]
    x_z=Seurat::RunPCA(as.matrix(t(x_z)),npcs = min(max(50,argList$nPCs+10),nrow(x_z)-2),verbose = FALSE)
    network=.myConcensusDEFn_step2_FindNeighbors(inputCentroids = x_z@cell.embeddings,argList=argList,verbose = F)
    
    snnNet=as.matrix(network$snn)
    ind=which(snnNet>0,arr.ind = T)
    
    snnNet=data.frame(Fnode=row.names(snnNet)[ind[,1]],Snode=row.names(snnNet)[ind[,2]],weight=snnNet[ind],stringsAsFactors = F)
    snnNet=snnNet[snnNet$Fnode!=snnNet$Snode,]
    snnNet$weight=1-snnNet$weight
    snnNet2=igraph::graph_from_data_frame(snnNet, directed = F, vertices = NULL)
    
    snnNet2=igraph::distances(snnNet2, mode = c("all"), weights = snnNet$weight, algorithm = "dijkstra")
    if(length(which(snnNet2==Inf))>0){
      snnNet2[which(snnNet2==Inf)]=(max(snnNet2[snnNet2!=Inf],na.rm = T)+1)
    }
    
    geneList=colnames(meta_z)
    
    
    cat("           Calculating cor and MI ...\n")
    resCorMI=NULL
    slGenes=apply(meta_z,2,function(x) max(abs(x)))
    slGenes=names(slGenes)[which(slGenes>2)]
    
    myCorMIfn=function(inputGenes,z_score_mat,snnNet2){
      resCorMI=NULL
      for(i in inputGenes){
        tmp=z_score_mat[,which(colnames(z_score_mat)==i)]
        tmpNet=snnNet2[row.names(snnNet2) %in% names(tmp),colnames(snnNet2) %in% names(tmp)]
        tmp=tmp[match(colnames(tmpNet),names(tmp))]
        tmp=data.frame(zscore=tmp,centroid=names(tmp),gene=i,stringsAsFactors = F)
        if(max(abs(tmp$zscore))>2){
          zOrdered=tmp$zscore[order(tmp$zscore,decreasing = T)]
          
          tmpInd=which(row.names(tmpNet) %in% tmp$centroid[which(tmp$zscore>=max(zOrdered[3],1))])
          if(length(tmpInd)>0){
            if(length(tmpInd)>1){
              distMat=tmpNet[tmpInd,]
              distMat=apply(distMat,2,mean)
            } else {
              distMat=tmpNet[tmpInd,]
            }
            
            tmp=data.frame(distance=distMat,heat=tmp$zscore)
            tmpCor=cor(tmp$distance,tmp$heat)
            
            resCorMI=rbind(resCorMI,data.frame(gene=i,cor=tmpCor,stringsAsFactors = F))
            
          } else{
            resCorMI=rbind(resCorMI,data.frame(gene=i,cor=0,stringsAsFactors = F))
          }
        }
      }
      return(resCorMI)
    }
    
    if(argList$ncores>1){
      geneList=split(as.character(slGenes), cut(1:length(as.character(slGenes)), argList$ncores, labels = FALSE)) 
    } else {
      geneList=list(as.character(slGenes))
    }
    
    resCor=parallel::mclapply(geneList,myCorMIfn,z_score_mat=meta_z,snnNet2=snnNet2,mc.cores = argList$ncores)
    resCor=do.call("rbind",resCor)
    
    slInd=which(abs(meta_z)>2,arr.ind = T)
    if(extendedMode){
      if((!all(row.names(meta_z)==row.names(meta_z.8)))|(!all(colnames(meta_z)==colnames(meta_z.8)))){
        stop("Error in matching!")
      }
      if((!all(row.names(meta_z)==row.names(med_n)))|(!all(colnames(meta_z)==colnames(med_n)))){
        stop("Error in matching!")
      }
    }
    
    
    if((!all(row.names(meta_z)==row.names(med_pct.1)))|(!all(colnames(meta_z)==colnames(med_pct.1)))){
      stop("Error in matching!")
    }
    
    
    if(extendedMode){
      res_arranged=data.frame(gene=colnames(meta_z)[slInd[,2]],centroid=row.names(meta_z)[slInd[,1]],zscore=meta_z[slInd],zscore.8=meta_z.8[slInd],pct.1_median=med_pct.1[slInd],pct.2_median=med_pct.2[slInd],logFC_median=med_logFC[slInd],n_median=med_n[slInd],stringsAsFactors = F)
    } else {
      res_arranged=data.frame(gene=colnames(meta_z)[slInd[,2]],centroid=row.names(meta_z)[slInd[,1]],zscore=meta_z[slInd],pct.1_median=med_pct.1[slInd],stringsAsFactors = F)
    }
    
    res_arranged=merge(res_arranged,resCor,by="gene",all.x=T)
    
    #zscore_count measures in how many centroids each gene as a nominal significance
    zscore_counts=aggregate(zscore~gene,data=res_arranged,function(x) sum(x>2))
    colnames(zscore_counts)=c("gene","zscore_count")
    res_arranged=merge(res_arranged,zscore_counts,by="gene",all.x=T)
    zscore_counts=aggregate(zscore~gene,data=res_arranged,function(x) sum(abs(x)>2))
    colnames(zscore_counts)=c("gene","zscore_count_abs")
    res_arranged=merge(res_arranged,zscore_counts,by="gene",all.x=T)
    if(sum(res_arranged$zscore_count>0)>0){
      res_arranged=res_arranged[res_arranged$zscore_count>0,]
    }
    
    res_arranged=res_arranged[order(res_arranged$zscore,decreasing = T),]
    
    save(res_arranged,file=.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T,propImportant = T))
    
    qsave(prop_mat,file=.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    
    
    if(F){
      prop_mat_list=list()
      for(i in 1:length(res)){
        prop_mat_list=c(prop_mat_list,res[[i]]$prop_mat)
        names(prop_mat_list)[length(prop_mat_list)]=res[[i]]$dsName
      }
      
      pseudocell_sim_list=list()
      for(i in 1:length(res)){
        pseudocell_sim_list=c(pseudocell_sim_list,res[[i]]$pseudocell_sim_mat)
        names(pseudocell_sim_list)[length(pseudocell_sim_list)]=res[[i]]$dsName
      }
      
      #res_array$prop_mat=prop_mat_list
      #res_array$pseudocell_sim_list=pseudocell_sim_list
      
    }
    
    
    if(F){
      qsave(res_array,file=.myFilePathMakerFn("res_dataset_array",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
    }
    
    
    
    gc()
    
  }
  
  return("Done")
}

.myConcensusDEFn_step2_fast=function(argList,expData=NULL,input_embedding=NULL,input_pd=NULL,input_prop_mat=NULL,input_pca_centroid=NULL,input_pca_centroids_assignments=NULL,minCellCountThr=4,meta_method="Stouffer",addClusteringModule=F,analysis_seed=1,extendedMode=F,L2Norm=T,mergePseudocells=T,merging_strength=0.3,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,hierarchical_mode=F){
  
  #argList=.ArgList;inputEmbeddings=NULL;inputPhenoData=NULL;inputExpData=NULL;organism="Mouse";batch_variable="anno_batch";run_harmony=F;addAnno=F;addClusteringModule=F;L2Norm=T;mergePseudocells=T;generateUMAP=F;extendedMode=F;merging_strength=0.3
  
  #minCellCountThr=4;meta_method="Stouffer";addClusteringModule=F;analysis_seed=1;extendedMode=F;L2Norm=T;mergePseudocells=T;merging_strength=0.3;analysis_seed=1;analysis_seed=1;include.singletons=T;colNormalize=T;input_embedding=NULL;input_pd=NULL;input_prop_mat=NULL;input_pca_centroid=NULL;input_pca_centroids_assignments=NULL;hierarchical_mode=F
  #expData=NULL
  
  
  require(qs)
  require(purrr)
  require(furrr)
  
  reRunCheck=T
  if(file.exists(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))&!argList$newRun){
    reRunCheck=tryCatch({tst=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T));F}, error=function(e) {return(T)})
    
    if(!reRunCheck){
      reRunCheck=tryCatch({tst=qs::qread(.myFilePathMakerFn("res_meta",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T));F}, error=function(e) {return(T)})
    }
    
  }
  
  res="Done"
  if(reRunCheck|hierarchical_mode){
    
    set.seed(analysis_seed)
    supportingFractionThr=argList$DE_supportingFractionThr
    n.adaptiveKernel=argList$DE_n.adaptiveKernel
    nPropIter=argList$DE_nPropIter
    
    #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    if(is.null(input_pd)){
      load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    } else {
      pd=input_pd
      
    }
    
    if(is.null(input_embedding)){
      load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    } else {
      harmony_embeddings=input_embedding
      
    }
    
    
    harmony_embeddings=harmony_embeddings[match(row.names(pd),row.names(harmony_embeddings)),,drop=F]
    if(sum(is.na(harmony_embeddings))>0){
      stop("Error in matching Names!")
    }
    
    if(is.null(expData)){
      stop("expression data should be provided!")
    }
    
    {
      tmp_pd=data.frame(sample=colnames(expData[[1]]),anno_batch=expData[[1]]$anno_batch,stringsAsFactors = F)
      #pd=pd[row.names(pd) %in% tmp_pd$sample,]
      tmp_pd=tmp_pd[match(row.names(pd),tmp_pd$sample),]
      if(sum(is.na(tmp_pd$anno_batch))>0){
        stop("Error! phenoData doesn't match with the expression data")
      }
     
    }
    
    pcaList=list(as.data.frame(harmony_embeddings))
    
    
    dataArranged = list(list(logNorm=Seurat:::NormalizeData.default(counts(expData[[1]]),normalization.method = "LogNormalize",verbose = F),
                          dsName=as.character("oneDS")))
    
    
    if(ncol(dataArranged[[1]]$logNorm)<ncol(expData[[1]])){
      warning("Annotation was found for only a subset of exp data, that subset was only used in the analysis!")
    }
    
    res_fd=NULL
    for(i in 1:length(expData)){
      .printDot()
      if(class(expData[[i]])[1]=="SingleCellExperiment"){
        fd=as.data.frame(rowData(expData[[i]]))
      } else {
        fd=as.data.frame(expData[[i]]@assays$RNA@meta.features)
      }
      
      if(sum(colnames(fd)=="ensembl_gene_id")==0){
        fd$ensembl_gene_id=row.names(fd)
      }
      fd$ensembl_gene_id=gsub("_","-",fd$ensembl_gene_id)
      slCols=intersect(colnames(fd),c("gene_name","gene_biotype","symbol","gene_short_name","ensembl_gene_id"))
      fd=fd[,colnames(fd) %in% slCols,drop=F]
      if(sum(is.na(fd$ensembl_gene_id))>0){
        fd$ensembl_gene_id[is.na(fd$ensembl_gene_id)]=row.names(fd)[is.na(fd$ensembl_gene_id)]
      }
      
      if(!is.null(res_fd)){
        fd=fd[!fd$ensembl_gene_id %in% res_fd$ensembl_gene_id,,drop=F]
        if(nrow(fd)>0){
          res_fd=rbind(res_fd,fd)
        }
      } else {
        res_fd=fd
      }
    }
    fd=res_fd
    rm(res_fd,expData,pcaList,harmony_embeddings)
    gc()
    
    #centroidPCAdata=pca_centroid;nPropIter=1;n.neighbors=argList$prop.n.neighbors
    
    #centroidPCAdata=pca_centroid;exCentroids=NULL;runIndx=1;batchPCAdata=NULL;n.neighbors=argList$prop.n.neighbors;NNmethod="annoy"
    
    if(is.null(input_prop_mat)){
      if(is.null(input_pca_centroid)){
        load(.myFilePathMakerFn("pca_centroids",argList=argList))
      } else {
        pca_centroid=input_pca_centroid
        rm(input_pca_centroid)
      }
      
      if(sum(is.na(pca_centroid))>0){
        stop("Error: NA values were found in the pca of pseudocells!")
      }
      
      
      #centroidPCAdata=pca_centroid;exCentroids=NULL;runIndx=1;batchPCAdata=input_embedding;n.neighbors=argList$prop.n.neighbors;NNmethod="annoy"
      res=.myConcensusDEFn_step2_detail_newprop3_final_v14(dataArranged=dataArranged,input_pca_centroids_assignments=input_pca_centroids_assignments,centroidPCAdata=pca_centroid,argList=argList,exCentroids=NULL,runIndx=1,batchPCAdata=input_embedding,n.neighbors=argList$prop.n.neighbors,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,hierarchical_mode=hierarchical_mode)
      
      prop_mat=res$prop_mat
      pseudocell_sim_mat=res$pseudocell_sim_mat
      #qsave(pseudocell_sim_mat,file=.myFilePathMakerFn("res_pseudocell_sim",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
      res=res[["dataArranged"]]
    } else {
      prop_mat=input_prop_mat
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      matWeights=.myEffSizePropMat(prop_mat)
      
      matEffectiveSize=matWeights$effective_sample_size
      
      matWeights=matWeights$centroid_weights
      matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
      matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
      
      
      res=dataArranged
      for(i in 1:length(res)){
        tmp_prop_mat=prop_mat[,match(colnames(res[[i]]$logNorm),colnames(prop_mat))]
        #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
        tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
        tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
        tmp_weights[tmp_effsize<4]=0
        tmp=list(prop_mat=tmp_prop_mat,data=res[[i]],matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
        res[[i]]=tmp
      }
    }
    
    
    if(nrow(prop_mat)>100){
      ps_groups=split(row.names(prop_mat),cut(1:length(row.names(prop_mat)),min(max(argList$ncores,2),nrow(prop_mat)/5)))
      
      res2=list()
      for(ips in ps_groups){
        tmp_res=res[1]
        tmp_res[[1]]$prop_mat=tmp_res[[1]]$prop_mat[ips,]
        tmp_res[[1]]$matEffectiveSize=tmp_res[[1]]$matEffectiveSize[match(row.names(tmp_res[[1]]$prop_mat),names(tmp_res[[1]]$matEffectiveSize))]
        tmp_res[[1]]$matWeights=tmp_res[[1]]$matWeights[match(row.names(tmp_res[[1]]$prop_mat),names(tmp_res[[1]]$matWeights))]
        res2=c(res2,tmp_res)
      }
      res=res2
      rm(tmp_res)
    }
    
    
    tmpValCheck=(unlist(lapply(res,length)))
    if(sum(tmpValCheck==1)>0){
      stop(res[[which(tmpValCheck==1)[1]]])
    }
    rm(tmpValCheck)
    
    ###Checking for outlier centroids
    
    
    
    
    cat("           Calculating dataset specific z-scores ...\n")
    #tst=.myConcensusDEFn_step2_detail_exp_final(res[[4]],argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode)
    #dim(tst$prop_mat)
    #qsave(list(res=res,argList=argList,extendedMode=extendedMode),file="~/torm2.qs")
    
    res2=parallel::mclapply(res,.myConcensusDEFn_step2_detail_exp_final,argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode,mc.preschedule=F,mc.cores = argList$ncores)
    while(sum(unlist(lapply(res2,class))!="list")>0){
      gc()
      null_ind=which(unlist(lapply(res2,class))!="list")
      cat(paste0("           Consider increasing RAM! Redoing the z-score calculations for ",length(null_ind)," batch(es) ...\n"))
      res3=parallel::mclapply(res[null_ind],.myConcensusDEFn_step2_detail_exp_final,argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode,mc.preschedule=F,mc.cores = 3)
      res2[null_ind]=res3
      rm(res3)
    }
    res=res2
    
    res_combined=list(zscore=res[[1]]$zscore,pct.1=res[[1]]$pct.1,pct.2=res[[1]]$pct.2,logFC=res[[1]]$logFC,matWeights=res[[1]]$matWeights,exp1=res[[1]]$exp1,gene_name_list=res[[1]]$gene_name_list,pseudocell_name_list=res[[1]]$pseudocell_name_list,dsName=res[[1]]$dsName,matEffectiveSize=res[[1]]$matEffectiveSize)
    
    if(length(res)>1){
      for(ires in 2:length(res)){
        res_combined$zscore=rbind(res_combined$zscore,res[[ires]]$zscore)
        res_combined$pct.1=rbind(res_combined$pct.1,res[[ires]]$pct.1)
        res_combined$pct.2=rbind(res_combined$pct.2,res[[ires]]$pct.2)
        res_combined$logFC=rbind(res_combined$logFC,res[[ires]]$logFC)
        res_combined$exp1=rbind(res_combined$exp1,res[[ires]]$exp1)
        res_combined$matWeights=rbind(res_combined$matWeights,res[[ires]]$matWeights)
        res_combined$matEffectiveSize=rbind(res_combined$matEffectiveSize,res[[ires]]$matEffectiveSize)
        res_combined$pseudocell_name_list=c(res_combined$pseudocell_name_list,res[[ires]]$pseudocell_name_list)
        if(any(res_combined$gene_name_list!=res[[ires]]$gene_name_list)){
          stop("Error")
        }
      }
    }
    
     #qsave(dataArranged,file="~/torm_dataArranged.qs")
    
    med_pct.1=res_combined$pct.1
    #med_pct.2=res_meta$med_pct.2
    #med_logFC=res_meta$med_logFC
    meta_z=res_combined$zscore
    logFC=res_combined$logFC
    #res_array=res_meta$res_array
    rm(res,res2)
    gc()
    
    
    #tst=list(pct_mat=med_pct.1,argList = argList,meta_z_mat=meta_z,sig1_thr=3,centers=NULL,pct2_thr=0.3,pct_diff_thr=0.2,symmetric=F)
    #med_pct.1=tst$pct_mat;argList = tst$argList;meta_z=tst$meta_z_mat;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
    
    #pct_mat=med_pct.1;meta_z_mat=meta_z;pct_mat_ref=NULL;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
    
    #################################
    #Needs attention
    if(argList$ncores>1&(!hierarchical_mode)){
      if(!dir.exists("~/tmp_torm")){
        dir.create("~/tmp_torm",recursive = T)
      }
      if(file.exists("~/tmp_torm/torm.qs")){
        file.remove("~/tmp_torm/torm.qs")
      }
      if(file.exists("~/tmp_torm/torm_res.qs")){
        file.remove("~/tmp_torm/torm_res.qs")
      }
      qsave(list(pct_mat=med_pct.1,argList = argList,meta_z_mat=meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F),"~/tmp_torm/torm.qs")
      system("Rscript ~/serverFiles/Rpctdiffcount.sh")  
      pct_diff_count=qread("~/tmp_torm/torm_res.qs")
      
    } else {
      #pct_mat=med_pct.1;argList = argList;meta_z_mat=meta_z;sig1_thr=sig1_thr;centers=NULL;pct2_thr=pct2_thr;pct_diff_thr=pct_diff_thr;symmetric=F;hierarchical_mode = hierarchical_mode
      pct_diff_count=.extra_sconline.PctScoreFn_v2(pct_mat=med_pct.1,argList = argList,meta_z_mat=meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F,hierarchical_mode = hierarchical_mode)
      
    }
    
    
    #argList$ncores=3
    #pct_diff_count2=.extra_sconline.PctScoreFn(pct_mat=med_pct.1,argList = argList,meta_z_mat=meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F)
    
    
    diff=prop_mat
    diff@x=rep(1,length(diff@x))
    diff=diff %*% t(diff)
    diff=sweep(diff,1,diag(diff),"/")
    diff=as.matrix(diff)
    diff[diff<0.100001]=0.100001
    diff=abs(log10(diff))
    diff=diff+pct_diff_count
    
    diff=hclust(as.dist(diff),method = "complete")
    diff=cutree(diff,h=0.999)
    
    if(F){
      diff=data.frame(pseudocell=names(diff),cluster=diff,stringsAsFactors = F)
      adj2=prop_mat
      load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
      pd=pd[colnames(adj2),]
      
      table(is.na(pd$anno_orig_cellState))
      tmp= adj2 %*% as.matrix(.myOneHotFn(pd$anno_orig_cellState))
      
      tmp=apply(tmp,1,function(x){
        y=which(x==max(x))[1]
        tmp2=data.frame(cluster=colnames(tmp)[y],purity=x[y])
        tmp2
      })
      
      tmp=do.call("rbind",tmp)
      tmp2=aggregate(purity~cluster,data=tmp,median)
      
      table(tmp2[,2]>0.8)
      
      
      
      diff=merge(diff,data.frame(pseudocell=row.names(tmp),tmp,stringsAsFactors = F),by="pseudocell")
      diff2=aggregate(cluster.y~cluster.x,data=diff,function(x){
        y=as.numeric(table(x))
        max(y)/sum(y)
      })
    }
    
    
    if(sum(duplicated(diff))>0&max(res_combined$matEffectiveSize)/ncol(prop_mat)<0.98){
      cat(paste0("           Merging pseudocells based on their DE count; retaining ",length(unique(diff))," out of ",length(diff),"\n"))
      diff=data.frame(pseudocell=names(diff),cluster=diff,stringsAsFactors = F)
      diff_id=diff[!duplicated(diff$cluster),]
      diff=merge(diff,diff_id,by="cluster")
      diff2=reshape2::dcast(pseudocell.x~pseudocell.y,data=diff,fun.aggregate=length)
      diff2=diff2[match(row.names(prop_mat),diff2[,1]),]
      diff2=t(as.matrix(diff2[,-1]))
      
      diff2=as(diff2,"dgCMatrix")
      diff2=diff2 %*% prop_mat
      diff2 <- .extra_matrix_rowNorm(diff2)#Matrix::Diagonal(x = 1 / (rowSums(diff2))) %*% diff2
      
      #.diff2=diff2
      
      rowMeans_drop0 <- function (dgCMat) {
        RowInd <- dgCMat@i + 1
        sapply(split(dgCMat@x, RowInd), function(x)quantile(x,0.95))
      }
      
      diff2 =  .extra_matrix_rowNorm(input_mat=diff2,rowValues=1 / (rowMeans_drop0(diff2)))#Matrix::Diagonal(x = 1 / (rowMeans_drop0(diff2))) %*% diff2
      #diff2@x=pmin(diff2@x,1)
      colMax_vals=c()
      for(i in seq(1,ncol(diff2),5000)){
        tmp_max=as.numeric(qlcMatrix::colMax(diff2[,i:min(i+4999,ncol(diff2))]))
        colMax_vals=c(colMax_vals,as.numeric(tmp_max))
      }
      if(nrow(diff2)>1){
        diff2 =  .extra_matrix_colNorm(input_mat=diff2,colValues=1 / colMax_vals)#diff2 %*% Matrix::Diagonal(x = 1 / colMax_vals)
        diff2 <- .extra_matrix_rowNorm(input_mat=diff2)#Matrix::Diagonal(x = 1 / rowSums(diff2)) %*% diff2
      }
      
      
      
      
      matWeights=.myEffSizePropMat(diff2)
      
      matEffectiveSize=matWeights$effective_sample_size
      matWeights=matWeights$centroid_weights
      
      
      if(max(matEffectiveSize)/ncol(prop_mat)<0.98&nrow(diff2)>1){
        matWeights=matWeights[match(row.names(diff2),names(matWeights))]
        matEffectiveSize=matEffectiveSize[match(row.names(diff2),names(matEffectiveSize))]
        
        prop_mat_org=prop_mat
        prop_mat=diff2
        
        prop_mat_org=prop_mat_org[match(row.names(prop_mat),row.names(prop_mat_org)),]
        changed_ps_list=apply(prop_mat_org - prop_mat,1,function(x) sum(abs(x)))
        changed_ps_list=row.names(prop_mat)[changed_ps_list>10^-10]
        if(length(changed_ps_list)<2){
          changed_ps_list=row.names(prop_mat)
        }
        res_unchanged=NULL
        if(length(setdiff(row.names(prop_mat),changed_ps_list))>0){
          res_unchanged=list(res_combined)
          res_unchanged[[1]]$zscore=res_unchanged[[1]]$zscore[setdiff(row.names(prop_mat),changed_ps_list),,drop=F]
          res_unchanged[[1]]$logFC=res_unchanged[[1]]$logFC[match(row.names(res_unchanged[[1]]$zscore),row.names(res_unchanged[[1]]$logFC)),,drop=F]
          res_unchanged[[1]]$exp1=res_unchanged[[1]]$exp1[match(row.names(res_unchanged[[1]]$zscore),row.names(res_unchanged[[1]]$exp1)),,drop=F]
          res_unchanged[[1]]$pct.1=res_unchanged[[1]]$pct.1[match(row.names(res_unchanged[[1]]$zscore),row.names(res_unchanged[[1]]$pct.1)),,drop=F]
          res_unchanged[[1]]$pct.2=res_unchanged[[1]]$pct.2[match(row.names(res_unchanged[[1]]$zscore),row.names(res_unchanged[[1]]$pct.2)),,drop=F]
          res_unchanged[[1]]$matWeights=res_unchanged[[1]]$matWeights[match(row.names(res_unchanged[[1]]$zscore),row.names(res_unchanged[[1]]$matWeights)),,drop=F]
          res_unchanged[[1]]$matEffectiveSize=res_unchanged[[1]]$matEffectiveSize[match(row.names(res_unchanged[[1]]$zscore),row.names(res_unchanged[[1]]$matEffectiveSize)),,drop=F]
          res_unchanged[[1]]$pseudocell_name_list=row.names(res_unchanged[[1]]$zscore)
        }
        
        
        res=dataArranged
        for(i in 1:length(res)){
          tmp_prop_mat=diff2[,match(colnames(res[[i]]$logNorm),colnames(diff2)),drop=F]
          #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
          tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
          tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
          tmp_weights[tmp_effsize<4]=0
          tmp=list(prop_mat=tmp_prop_mat,data=res[[i]],matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
          res[[i]]=tmp
        }
        
        #changed pseudocells
        
        gc()
        plan("sequential")
        cat(paste0("           Redoing the propagation after the merging\n"))
        
        if(length(changed_ps_list)>0){
          if(length(changed_ps_list)>100){
            ps_groups=split(changed_ps_list,cut(1:length(changed_ps_list),min(max(argList$ncores,2),length(changed_ps_list)/5)))
            
            res2=list()
            for(ips in ps_groups){
              tmp_res=res[1]
              tmp_res[[1]]$prop_mat=tmp_res[[1]]$prop_mat[ips,]
              tmp_res[[1]]$matEffectiveSize=tmp_res[[1]]$matEffectiveSize[match(row.names(tmp_res[[1]]$prop_mat),names(tmp_res[[1]]$matEffectiveSize))]
              tmp_res[[1]]$matWeights=tmp_res[[1]]$matWeights[match(row.names(tmp_res[[1]]$prop_mat),names(tmp_res[[1]]$matWeights))]
              res2=c(res2,tmp_res)
            }
            
          } else {
            res2=list()
            tmp_res=res[1]
            tmp_res[[1]]$prop_mat=tmp_res[[1]]$prop_mat[changed_ps_list,,drop=F]
            tmp_res[[1]]$matEffectiveSize=tmp_res[[1]]$matEffectiveSize[match(row.names(tmp_res[[1]]$prop_mat),names(tmp_res[[1]]$matEffectiveSize))]
            tmp_res[[1]]$matWeights=tmp_res[[1]]$matWeights[match(row.names(tmp_res[[1]]$prop_mat),names(tmp_res[[1]]$matWeights))]
            res2=c(res2,tmp_res)
          }
          
          res3=parallel::mclapply(res2,.myConcensusDEFn_step2_detail_exp_final,argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode,mc.cores = argList$ncores)
          while(sum(unlist(lapply(res3,class))=="NULL")>0){
            gc()
            null_ind=which(unlist(lapply(res3,class))=="NULL")
            cat(paste0("           Consider increasing RAM! Redoing the z-score calculations for ",length(null_ind)," batches ...\n"))
            res4=parallel::mclapply(res2[null_ind],.myConcensusDEFn_step2_detail_exp_final,argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode,mc.preschedule=F,mc.cores = 3)
            res3[null_ind]=res4
            rm(res4)
          }
          res2=res3
          res=c(res2,res_unchanged)
          rm(res2,res_unchanged,res3)
        } else {
          stop("Error!")
        }
        
        
        
        {
          res_combined=list(zscore=res[[1]]$zscore,pct.1=res[[1]]$pct.1,pct.2=res[[1]]$pct.2,logFC=res[[1]]$logFC,matWeights=res[[1]]$matWeights,exp1=res[[1]]$exp1,gene_name_list=res[[1]]$gene_name_list,pseudocell_name_list=res[[1]]$pseudocell_name_list,dsName=res[[1]]$dsName,matEffectiveSize=res[[1]]$matEffectiveSize)
          
          if(length(res)>1){
            for(ires in 2:length(res)){
              if(sum(duplicated(row.names(res_combined$zscore),row.names(res[[ires]]$zscore)))>0){
                stop("Error!!")
              }
              res_combined$zscore=rbind(res_combined$zscore,res[[ires]]$zscore)
              res_combined$pct.1=rbind(res_combined$pct.1,res[[ires]]$pct.1)
              res_combined$pct.2=rbind(res_combined$pct.2,res[[ires]]$pct.2)
              res_combined$logFC=rbind(res_combined$logFC,res[[ires]]$logFC)
              res_combined$exp1=rbind(res_combined$exp1,res[[ires]]$exp1)
              res_combined$matWeights=rbind(res_combined$matWeights,res[[ires]]$matWeights)
              res_combined$matEffectiveSize=rbind(res_combined$matEffectiveSize,res[[ires]]$matEffectiveSize)
              res_combined$pseudocell_name_list=c(res_combined$pseudocell_name_list,res[[ires]]$pseudocell_name_list)
              if(any(res_combined$gene_name_list!=res[[ires]]$gene_name_list)){
                stop("Error")
              }
            }
          }
          
          res_meta=res_combined
        }
        
        med_pct.1=res_meta$pct.1
        #med_pct.2=res_meta$med_pct.2
        #med_logFC=res_meta$med_logFC
        meta_z=res_meta$zscore
        #res_array=res_meta$res_array
        logFC=res_meta$logFC
      }
      
      
      gc()
    }
    
    
    de_pct_res=NULL
    if(addClusteringModule&F){
      
      cosine_dist=.extra_sconline.CosineDistFn(inputMatrix = meta_z,sig_thr = 2)
      
      de_pct_res=.extra_sconline.PctScoreFn(argList = argList,pct_mat=med_pct.1,meta_z_mat=meta_z,pct_mat_ref=NULL,sig1_thr=3,centers=NULL,pct2_thr=0.3,pct_diff_thr=0.2,symmetric=F)
      
    }
    #meta_z_mat=meta_z$meta_z;cosine_dist=cosine_dist;de_dist=de_pct_res;pct_mat=meta_z$med_pct.1;min_marker_thr=20;sig1_thr=3;sig2_thr=1;pct_de_count_thr=1;pct_diff_thr=0.2;pct2_thr=0.3
    #####################
    
    
    
    #####################
    
    
    if(extendedMode){
      res_meta=list(meta_z=meta_z,meta_z.8=meta_z.8,med_pct.1=med_pct.1,med_pct.2=med_pct.2,med_logFC=med_logFC,med_n=med_n,fd=fd)
    } else {
      res_meta=list(meta_z=meta_z,med_pct.1=med_pct.1,logFC=logFC,fd=fd)
    }
    
    if(!hierarchical_mode){
      qsave(res_meta,file=.myFilePathMakerFn("res_meta",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
      qsave(prop_mat,file=.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
      
      #rm(res_meta)
    } else {
      res=list(res_meta=res_meta,prop_mat=prop_mat)
    }
    
    
    
    if(F){
      #archived
      load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
      
      x_z=apply(meta_z,2,max)
      x_z=meta_z[,which(x_z>2),drop=F]
      x_z=Seurat::RunPCA(as.matrix(t(x_z)),npcs = min(max(50,argList$nPCs+10),nrow(x_z)-2),verbose = FALSE)
      network=.myConcensusDEFn_step2_FindNeighbors(inputCentroids = x_z@cell.embeddings,argList=argList,verbose = F)
      
      snnNet=as.matrix(network$snn)
      ind=which(snnNet>0,arr.ind = T)
      
      snnNet=data.frame(Fnode=row.names(snnNet)[ind[,1]],Snode=row.names(snnNet)[ind[,2]],weight=snnNet[ind],stringsAsFactors = F)
      snnNet=snnNet[snnNet$Fnode!=snnNet$Snode,]
      snnNet$weight=1-snnNet$weight
      snnNet2=igraph::graph_from_data_frame(snnNet, directed = F, vertices = NULL)
      
      snnNet2=igraph::distances(snnNet2, mode = c("all"), weights = snnNet$weight, algorithm = "dijkstra")
      if(length(which(snnNet2==Inf))>0){
        snnNet2[which(snnNet2==Inf)]=(max(snnNet2[snnNet2!=Inf],na.rm = T)+1)
      }
      
      geneList=colnames(meta_z)
      
      
      cat("           Calculating cor and MI ...\n")
      resCorMI=NULL
      slGenes=apply(meta_z,2,function(x) max(abs(x)))
      slGenes=names(slGenes)[which(slGenes>2)]
      
      myCorMIfn=function(inputGenes,z_score_mat,snnNet2){
        resCorMI=NULL
        for(i in inputGenes){
          tmp=z_score_mat[,which(colnames(z_score_mat)==i)]
          tmpNet=snnNet2[row.names(snnNet2) %in% names(tmp),colnames(snnNet2) %in% names(tmp)]
          tmp=tmp[match(colnames(tmpNet),names(tmp))]
          tmp=data.frame(zscore=tmp,centroid=names(tmp),gene=i,stringsAsFactors = F)
          if(max(abs(tmp$zscore))>2){
            zOrdered=tmp$zscore[order(tmp$zscore,decreasing = T)]
            
            tmpInd=which(row.names(tmpNet) %in% tmp$centroid[which(tmp$zscore>=max(zOrdered[3],1))])
            if(length(tmpInd)>0){
              if(length(tmpInd)>1){
                distMat=tmpNet[tmpInd,]
                distMat=apply(distMat,2,mean)
              } else {
                distMat=tmpNet[tmpInd,]
              }
              
              tmp=data.frame(distance=distMat,heat=tmp$zscore)
              tmpCor=cor(tmp$distance,tmp$heat)
              
              resCorMI=rbind(resCorMI,data.frame(gene=i,cor=tmpCor,stringsAsFactors = F))
              
            } else{
              resCorMI=rbind(resCorMI,data.frame(gene=i,cor=0,stringsAsFactors = F))
            }
          }
        }
        return(resCorMI)
      }
      
      if(argList$ncores>1){
        geneList=split(as.character(slGenes), cut(1:length(as.character(slGenes)), argList$ncores, labels = FALSE)) 
      } else {
        geneList=list(as.character(slGenes))
      }
      
      resCor=parallel::mclapply(geneList,myCorMIfn,z_score_mat=meta_z,snnNet2=snnNet2,mc.cores = argList$ncores)
      resCor=do.call("rbind",resCor)
      
      slInd=which(abs(meta_z)>2,arr.ind = T)
      if(extendedMode){
        if((!all(row.names(meta_z)==row.names(meta_z.8)))|(!all(colnames(meta_z)==colnames(meta_z.8)))){
          stop("Error in matching!")
        }
        if((!all(row.names(meta_z)==row.names(med_n)))|(!all(colnames(meta_z)==colnames(med_n)))){
          stop("Error in matching!")
        }
      }
      
      
      if((!all(row.names(meta_z)==row.names(med_pct.1)))|(!all(colnames(meta_z)==colnames(med_pct.1)))){
        stop("Error in matching!")
      }
      
      
      if(extendedMode){
        res_arranged=data.frame(gene=colnames(meta_z)[slInd[,2]],centroid=row.names(meta_z)[slInd[,1]],zscore=meta_z[slInd],zscore.8=meta_z.8[slInd],pct.1_median=med_pct.1[slInd],pct.2_median=med_pct.2[slInd],logFC_median=med_logFC[slInd],n_median=med_n[slInd],stringsAsFactors = F)
      } else {
        res_arranged=data.frame(gene=colnames(meta_z)[slInd[,2]],centroid=row.names(meta_z)[slInd[,1]],zscore=meta_z[slInd],pct.1_median=med_pct.1[slInd],stringsAsFactors = F)
      }
      
      res_arranged=merge(res_arranged,resCor,by="gene",all.x=T)
      
      #zscore_count measures in how many centroids each gene as a nominal significance
      zscore_counts=aggregate(zscore~gene,data=res_arranged,function(x) sum(x>2))
      colnames(zscore_counts)=c("gene","zscore_count")
      res_arranged=merge(res_arranged,zscore_counts,by="gene",all.x=T)
      zscore_counts=aggregate(zscore~gene,data=res_arranged,function(x) sum(abs(x)>2))
      colnames(zscore_counts)=c("gene","zscore_count_abs")
      res_arranged=merge(res_arranged,zscore_counts,by="gene",all.x=T)
      if(sum(res_arranged$zscore_count>0)>0){
        res_arranged=res_arranged[res_arranged$zscore_count>0,]
      }
      
      res_arranged=res_arranged[order(res_arranged$zscore,decreasing = T),]
      
      save(res_arranged,file=.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T,propImportant = T))
      
    }
    
    
    
    if(F){
      prop_mat_list=list()
      for(i in 1:length(res)){
        prop_mat_list=c(prop_mat_list,res[[i]]$prop_mat)
        names(prop_mat_list)[length(prop_mat_list)]=res[[i]]$dsName
      }
      
      pseudocell_sim_list=list()
      for(i in 1:length(res)){
        pseudocell_sim_list=c(pseudocell_sim_list,res[[i]]$pseudocell_sim_mat)
        names(pseudocell_sim_list)[length(pseudocell_sim_list)]=res[[i]]$dsName
      }
      
      #res_array$prop_mat=prop_mat_list
      #res_array$pseudocell_sim_list=pseudocell_sim_list
      
    }
    
    
    if(F){
      qsave(res_array,file=.myFilePathMakerFn("res_dataset_array",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
    }
    
    
    
    gc()
    
  }
  
  return(res)
}

.extra_sconline.pseudosim_archive=function(argList,cos_dist=F,binarize=T,n.trees=50,k.param = 20,NNmethod="annoy",L2Norm=T,colNorm=T,rowNorm=T,n.neighbors=20,...){
  #argList=.ArgList;n.trees=50;k.param = 20;NNmethod="annoy";L2Norm=T;binarize=T;colNorm=T;rowNorm=T;n.neighbors=20
  require(Matrix)
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=.extra_matrix_rowNorm(input_mat=inputMat,rowValues=1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  
  load(.myFilePathMakerFn("pca_centroids",argList=argList))
  centroidPCAdata=pca_centroid
  
  
  load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
  batchPCAdata=harmony_embeddings
  
  
  if(sum(argList$singleton.method %in% c("snn","fast"))==0){
    stop("unrecognized singleton.method")
  }
  
  
  
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  
  
  if(NNmethod=="annoy"){
    if(!file.exists(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))){
      idx=Seurat:::AnnoyBuildIndex(data = batchPCAdata, metric = "euclidean", 
                                   n.trees = n.trees)
      nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = batchPCAdata,k=k.param,include.distance = T,search.k = -1)
      qsave(nn.ranked.1,file=.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
    } else {
      nn.ranked.1=qread(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
    }
  } else {
    nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  }
  
  
  
  j <- as.numeric(t(nn.ranked.1$nn.idx))
  i <- ((1:length(j)) - 1)%/%n.neighbors + 1
  
  if(binarize){
    adj = sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(batchPCAdata),nrow(batchPCAdata)))
  } else {
    affinities=.extra_matrix_rowNorm(input_mat=nn.ranked.1$nn.dists,rowValues=1/(nn.ranked.1$nn.dists[,5]+0.000001))#Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,5]+0.000001)) %*% nn.ranked.1$nn.dists
    affinities@x=-1*affinities@x^2
    affinities@x=exp(affinities@x)
    affinities[,1]=affinities[,2]
    x=as.numeric(t(affinities))
    adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(batchPCAdata),nrow(batchPCAdata)))
  }
  
  rownames(adj) <- row.names(batchPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  adj=.extra_matrix_rowNorm(adj)# Matrix::Diagonal(x=1/rowSums(adj)) %*% adj
  
  
  
  load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
  
  adj2=.extra_matrix_rowNorm(input_mat=adj,rowValues=1/qlcMatrix::rowMax(adj)) #Matrix::Diagonal(x=1/qlcMatrix::rowMax(adj)) %*% adj
  adj2=Matrix::drop0(adj2,tol=0.01)
  adj2@x=rep(1,length(adj2@x))
  adj=Matrix::drop0(adj*adj2)
  
  
  adj_t=t(adj)
  adj=adj[sl_pseudo$pseudocell,]
  row.names(adj)=sl_pseudo$cluster
  
  adj=adj %*% adj_t
  
  res1=.extra_sconline.pseudosim(adj=adj,adj_t=adj_t,argList=argList,cos_dist=cos_dist,colNorm=colNorm,rowNorm=rowNorm)
  
  if(F){
    check_adj_sim=T
    counter=0
    adj_sim=adj
    while(check_adj_sim&counter<5){
      counter=counter+1
      adj_sim <- .extra_matrix_rowNorm(adj_sim)#Matrix::Diagonal(x = 1 / rowSums(adj_sim)) %*% adj_sim
      adj_sim <- adj_sim %*% adj_t
      tmp_adj_sim=adj_sim %*% t(adj_sim)
      tmp_adj_sim=.extra_matrix_rowNorm(input_mat=tmp_adj_sim,rowValues=1/qlcMatrix::rowMax(tmp_adj_sim))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(tmp_adj_sim)) %*% tmp_adj_sim
      tmp_adj_sim["ps_100","ps_53"]
      diag(tmp_adj_sim)=0
      tmp_adj_sim=as.numeric(qlcMatrix::rowMax(tmp_adj_sim))
      names(tmp_adj_sim)=row.names(tmp_adj_sim)
      if(mean(tmp_adj_sim)>0.5){
        check_adj_sim=F
      }
      
    }
  }
  
  
  #adj_sim=adj_sim %*% t(adj_sim)
  if(F){
    
    #SNN
    nn.ranked=nn.ranked.1$nn.idx
    snn.matrix <- Seurat:::ComputeSNN(nn_ranked = as.matrix(nn.ranked), prune = 1/16)
    rownames(x = snn.matrix) <- row.names(batchPCAdata)
    colnames(x = snn.matrix) <- row.names(batchPCAdata)
    diff=snn.matrix[sl_pseudo$pseudocell,]
    row.names(diff)=sl_pseudo$cluster
    diff=qlcMatrix::cosSparse(t(diff))
    row.names(diff)=colnames(diff)=sl_pseudo$cluster
    #cosine
    prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=.ArgList,uniformImportant=T,propImportant = T,qsFormat=T))
    prop_mat=.extra_matrix_rowNorm(input_mat=prop_mat,rowValues=1 / (rowSums(prop_mat)+0.000000000001))#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
    
    
    inputData=adj
    inputData=qlcMatrix::cosSparse(t(inputData))
    row.anmes(inputData)=colnames(inputData)=row.names(adj)
    diff=inputData
    
    
    inputData=adj
    inputData=.extra_matrix_rowNorm(input_mat=inputData)#Matrix::Diagonal(x = 1 / (rowSums(inputData)+0.000000000001)) %*% inputData
    inputData=inputData %*% t(adj)
    inputData=.extra_matrix_rowNorm(input_mat=inputData)#Matrix::Diagonal(x = 1 / (rowSums(inputData)+0.000000000001)) %*% inputData
    diff=inputData
    
    
    prop_mat2=myL2normFn(inputMat=adj)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    #c_c_aff=c_c_aff %*%c_c_aff %*% c_c_aff
    
    #c_c_aff=Matrix::drop0(c_c_aff,tol=quantile(c_c_aff@x,0.1))
    c_c_aff <- .extra_matrix_rowNorm(input_mat=c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
    c_c_aff_t=c_c_aff2=c_c_aff
    
    inc_check=T
    inc_count=nrow(prop_mat)
    counter=0
    while((inc_check|counter<10)&counter<15){
      
      counter=counter+1
      
      c_c_aff=c_c_aff %*% c_c_aff_t*0.7+c_c_aff_t *0.3
      
      c_c_aff=Matrix::drop0(c_c_aff,tol=quantile(c_c_aff@x,0.00001))
      c_c_aff <- .extra_matrix_rowNorm(input_mat=c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
      
      tst=as.numeric(qlcMatrix::rowMax(c_c_aff))
      tst2=as.numeric(diag(c_c_aff))
      tst=as.numeric(qlcMatrix::rowMax(c_c_aff))
      tst=.extra_matrix_rowNorm(input_mat=c_c_aff,rowValues=1/tst)#Matrix::Diagonal(x=1/tst) %*% c_c_aff
      tst=which(tst>0.9,arr.ind = T)
      sl_ps_list=unique(tst[,2])
      #print(paste(counter,":",length(sl_ps_list)))
      if(length(sl_ps_list)/inc_count>0.95&length(sl_ps_list)<2000){
        inc_check=F
      }
      inc_count=length(sl_ps_list)
    }
    
    diff=c_c_aff
    
    snn.matrix <- ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(x = snn.matrix) <- rownames(x = object)
    colnames(x = snn.matrix) <- rownames(x = object)
    
  }
  
  
  
  return(res1)
}

.extra_sconline.pseudosim_archive11=function(argList,input_pca_centroids=NULL,input_prop_mat=NULL,input_embeddings=NULL,cos_dist=F,binarize=T,n.trees=50,k.param = 20,NNmethod="annoy",L2Norm=T,colNorm=T,rowNorm=T,n.neighbors=20,hierarchical_mode=F,...){
  #argList=.ArgList;n.trees=50;k.param = 20;NNmethod="annoy";L2Norm=T;binarize=T;colNorm=T;rowNorm=T;n.neighbors=20
  require(Matrix)
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=.extra_matrix_rowNorm(input_mat=inputMat,rowValues=1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  
  if(!file.exists(.myFilePathMakerFn("res_pseudocell_pc_sim",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))|hierarchical_mode){
    if(is.null(input_pca_centroids)){
      load(.myFilePathMakerFn("pca_centroids",argList=argList))
      centroidPCAdata=pca_centroid
    } else {
      centroidPCAdata=input_pca_centroids
    }
    
    
    
    if(is.null(input_embeddings)){
      load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
      batchPCAdata=harmony_embeddings
    } else {
      batchPCAdata=input_embeddings
    }
    
    
    
    if(sum(argList$singleton.method %in% c("snn","fast"))==0){
      stop("unrecognized singleton.method")
    }
    
    
    
    batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
    centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
    
    
    if(NNmethod=="annoy"){
      if(!file.exists(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))){
        idx=Seurat:::AnnoyBuildIndex(data = batchPCAdata, metric = "euclidean", 
                                     n.trees = n.trees)
        nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = batchPCAdata,k=k.param,include.distance = T,search.k = -1)
        if(!hierarchical_mode){
          qsave(nn.ranked.1,file=.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
        }
        
      } else {
        nn.ranked.1=qread(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
      }
    } else {
      nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
    }
    
    
    
    snn.matrix <- Seurat:::ComputeSNN(nn_ranked = nn.ranked.1$nn.idx, prune = 1/15)
    rownames(x = snn.matrix) <- rownames(x = batchPCAdata)
    colnames(x = snn.matrix) <- rownames(x = batchPCAdata)
    
    if(is.null(input_prop_mat)){
      prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=.ArgList,uniformImportant=T,propImportant = T,qsFormat=T))
      
    } else {
      prop_mat=input_prop_mat
    }
    prop_mat=.extra_matrix_rowNorm(input_mat=prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
    
    colMax_vals_m=qlcMatrix::colMax(prop_mat)
    colMax_vals_m=.extra_matrix_colNorm(input_mat=prop_mat,colValues=1/as.numeric(colMax_vals_m))#prop_mat %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=0.9)
    prop_m_hardCluster@x=rep(1,length(prop_m_hardCluster@x))
    
    
    res1=prop_m_hardCluster %*% snn.matrix %*% t(prop_m_hardCluster)
    
    res1=.extra_matrix_rowNorm(input_mat = res1,rowValues = 1/diag(res1))#Matrix::Diagonal(x = 1 / (diag(res1))) %*% res1
    if(!hierarchical_mode){
      qsave(res1,file=.myFilePathMakerFn("res_pseudocell_pc_sim",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
    }
    
  } else {
    res1=qread(.myFilePathMakerFn("res_pseudocell_pc_sim",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
  }
  
  
  
  return(res1)
}


.extra_sconline.pseudosim_archive2=function(argList,cos_dist=F,binarize=T,n.trees=50,k.param = 20,NNmethod="annoy",L2Norm=T,colNorm=T,rowNorm=T,n.neighbors=20,...){
  #argList=.ArgList;n.trees=50;k.param = 20;NNmethod="annoy";L2Norm=T;binarize=T;colNorm=T;rowNorm=T;n.neighbors=20
  require(Matrix)
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  
  load(.myFilePathMakerFn("pca_centroids",argList=argList))
  centroidPCAdata=pca_centroid
  
  
  load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
  batchPCAdata=harmony_embeddings
  
  
  if(sum(argList$singleton.method %in% c("snn","fast"))==0){
    stop("unrecognized singleton.method")
  }
  
  
  
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  
  
  if(NNmethod=="annoy"){
    if(!file.exists(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))){
      idx=Seurat:::AnnoyBuildIndex(data = batchPCAdata, metric = "euclidean", 
                                   n.trees = n.trees)
      nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = batchPCAdata,k=k.param,include.distance = T,search.k = -1)
      qsave(nn.ranked.1,file=.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
    } else {
      nn.ranked.1=qread(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
    }
  } else {
    nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  }
  
  
  
  j <- as.numeric(t(nn.ranked.1$nn.idx))
  i <- ((1:length(j)) - 1)%/%n.neighbors + 1
  
  if(binarize){
    adj = sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(batchPCAdata),nrow(batchPCAdata)))
  } else {
    affinities=.extra_matrix_rowNorm(input_mat = nn.ranked.1$nn.dists,rowValues = 1/(nn.ranked.1$nn.dists[,5]+0.000001))#Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,5]+0.000001)) %*% nn.ranked.1$nn.dists
    affinities@x=-1*affinities@x^2
    affinities@x=exp(affinities@x)
    affinities[,1]=affinities[,2]
    x=as.numeric(t(affinities))
    adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(batchPCAdata),nrow(batchPCAdata)))
  }
  
  rownames(adj) <- row.names(batchPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  adj= .extra_matrix_rowNorm(adj)#Matrix::Diagonal(x=1/rowSums(adj)) %*% adj
  
  
  
  load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
  
  adj2=.extra_matrix_rowNorm(input_mat = adj,rowValues = 1/qlcMatrix::rowMax(adj))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(adj)) %*% adj
  adj2=Matrix::drop0(adj2,tol=0.01)
  adj2@x=rep(1,length(adj2@x))
  adj=Matrix::drop0(adj*adj2)
  
  
  adj_t=t(adj)
  adj=adj[sl_pseudo$pseudocell,]
  row.names(adj)=sl_pseudo$cluster
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=.ArgList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat@x=rep(1,length(prop_mat@x))
  adj=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  adj=adj %*% adj_t
  res1=qlcMatrix::cosSparse(t(adj))
  row.names(res1)=colnames(res1)=row.names(adj)
  #res1=.extra_sconline.pseudosim(adj=adj,adj_t=adj_t,argList=argList,cos_dist=cos_dist,colNorm=colNorm,rowNorm=rowNorm)
  
  return(res1)
}


.extra_sconline.pseudosim=function(adj,adj_t=adj_t,argList,cos_dist=F,colNorm=T,rowNorm=T,...){
  #argList=.ArgList;n.trees=50;k.param = 20;NNmethod="annoy";L2Norm=T
  require(Matrix)
  
  check_adj_sim=T
  counter=0
  adj_sim=adj
  while(check_adj_sim&counter<5){
    counter=counter+1
    adj_sim <- .extra_matrix_rowNorm(adj_sim)#Matrix::Diagonal(x = 1 / rowSums(adj_sim)) %*% adj_sim
    adj_sim <- adj_sim %*% adj_t
    tmp_adj_sim=adj_sim %*% t(adj_sim)
    tmp_adj_sim=.extra_matrix_rowNorm(input_mat = tmp_adj_sim,rowValues = 1/qlcMatrix::rowMax(tmp_adj_sim))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(tmp_adj_sim)) %*% tmp_adj_sim
    
    diag(tmp_adj_sim)=0
    tmp_adj_sim=as.numeric(qlcMatrix::rowMax(tmp_adj_sim))
    if(mean(tmp_adj_sim)>0.1){
      check_adj_sim=F
    }
    
  }
  
  
  #adj_sim=adj_sim %*% t(adj_sim)
  
  
  prop_mat_list=list()
  rowMax_vals=c()
  colMax_vals=c()
  for(i in seq(1,ncol(adj_sim),50000)){
    tmp=adj_sim[,i:min(i+49999,ncol(adj_sim))]
    prop_mat_list=c(prop_mat_list,list(tmp))
    tmp_max=Matrix::drop0(tmp,tol=0.00001)
    tmp_max=as.numeric(qlcMatrix::rowMax(tmp_max))
    if(length(rowMax_vals)>0){
      rowMax_vals=pmax(rowMax_vals,tmp_max)
    } else {
      rowMax_vals=tmp_max
    }
  }
  
  if(sum(rowMax_vals<=0)>0){
    sl_rows=row.names(adj_sim)[which(rowMax_vals==0)]
    tmp_rowMax_vals=c()
    for(i in seq(1,ncol(adj_sim),50000)){
      tmp=c_c_aff[sl_rows,,drop=F] %*% adj_sim[,i:min(i+49999,ncol(adj_sim))]
      tmp_max=as.numeric(qlcMatrix::rowMax(tmp))
      if(length(tmp_rowMax_vals)>0){
        tmp_rowMax_vals=pmax(tmp_rowMax_vals,tmp_max)
      } else {
        tmp_rowMax_vals=tmp_max
      }
    }
    rowMax_vals[match(sl_rows,colnames(c_c_aff))]=tmp_rowMax_vals
  }
  
  for(i in 1:length(prop_mat_list)){
    tmp=prop_mat_list[[i]]
    tmp=Matrix::drop0(tmp,tol=min(rowMax_vals)*0.001)
    #tmp=tmp %*% Matrix::Diagonal(x=1/(qlcMatrix::colMax(tmp)+0.00001))
    if(rowNorm){
      tmp=.extra_matrix_rowNorm(input_mat = tmp,rowValues = 1/rowMax_vals)#Matrix::Diagonal(x = 1 / rowMax_vals) %*% tmp
    }
    
    if(colNorm){
      tmp=.extra_matrix_colNorm(input_mat = tmp,colValues = 1/qlcMatrix::colMax(tmp))#tmp %*% Matrix::Diagonal(x=1/qlcMatrix::colMax(tmp))
    }
    
    #tmp=Matrix::drop0(tmp,tol=0.1)
    prop_mat_list[[i]]=tmp
  }
  
  prop_mat=prop_mat_list[[1]]
  prop_mat=prop_mat[,colSums(prop_mat)>0]
  if(length(prop_mat_list)>1){
    for(i in 2:length(prop_mat_list)){
      tmp=prop_mat_list[[i]]
      tmp=tmp[,colSums(tmp)>0]
      prop_mat=cbind(prop_mat,tmp)
    }
    rm(tmp)
  }
  
  
  if(cos_dist){
    rwnames=row.names(prop_mat)
    prop_mat=qlcMatrix::cosSparse(t(prop_mat))
    row.names(prop_mat)=colnames(prop_mat)=rwnames
  } else {
    prop_mat=prop_mat %*% t(prop_mat)
    diag(prop_mat)=0
    prop_mat=.extra_matrix_rowNorm(input_mat = prop_mat,rowValues = 1/qlcMatrix::rowMax(prop_mat))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(prop_mat)) %*% prop_mat
    diag(prop_mat)=1
  }
  
  return(prop_mat)
}

.myConcensusDEFn_step2_detail_exp_final_org=function(inputData,argList,sd_offset=0.001,reg_sd=T,zscore_cap=10,minCellCountThr=4,extendedMode=F){
  
  #inputData=res[[1]];sd_offset=0.005;reg_sd=T;zscore_cap=10;minCellCountThr=4;extendedMode=F
  if(F){
    #Jonah's
    write.table(data.frame(), file=paste0("/tmp/touch_", as.character(runif(1)*10^20)),
                col.names=FALSE)
    
    library("RhpcBLASctl")
    omp_set_num_threads(4)
    blas_set_num_threads(4)
    
  }
  
  
  #inputData=res[[1]];sd_offset=0.001;reg_sd=T;zscore_cap=10;minCellCountThr=4
  
  #tmp_Propweights=rowSums(as.matrix(inputData$prop_mat))
  
  if(sum(names(argList)=="do.split.prop")==0){
    argList$do.split.prop=T
  }
  
  if(!argList$do.split.prop){
    inputData$prop_mat=.extra_matrix_rowNorm(inputData$prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(inputData$prop_mat)+0.000000000001)) %*% inputData$prop_mat
  }
  
  
  require(Matrix,quietly = T)
  
  zscore_cap=abs(zscore_cap)
  if(zscore_cap<2){
    warning("zscore cap was set to below 2. It was changed to 5.")
    zscore_cap=5
  }
  sd_offset=max(sd_offset,0)
  
  resGeneMeanSd=.extra_gmean(inputData$data$countData)
  resGeneMeanSd=data.frame(gene=names(resGeneMeanSd),geo_mean_count=resGeneMeanSd,stringsAsFactors = F)
  
  
  myWeightedVar=function(inputWeight,inputExp,min_cell_count){
    require(Matrix,quietly = T)
    inputExp=t(inputExp)
    res_binary=(sign(inputWeight)%*% sign(inputExp))
    exp_sq=inputExp^2
    res=(inputWeight%*% exp_sq)
    res=sweep(res,1,rowSums(inputWeight),"*")
    res=res- (inputWeight%*% inputExp)^2
    res=sweep(res,1,rowSums(inputWeight)^2,"/")
    
    singletons=apply(inputWeight,1,function(x) sum(x>0))
    singletons[singletons<=min_cell_count]=0
    singletons[singletons>0]=1
    if(sum(singletons<1)>0){
      res=sweep(res,1,singletons,"*")
    }
    res[res_binary<=min_cell_count]=0
    return(res)
  }
  
  prop_mat=inputData$prop_mat
  
  
  prop_mat_c=prop_mat
  if(quantile(rowSums(prop_mat_c > 0,na.rm = T), 0.25) < (0.85*ncol(prop_mat))){
    prop_mat_c[prop_mat_c>0]=1
  }
  
  prop_mat_c=1-prop_mat_c
  prop_mat_c <- .extra_matrix_rowNorm(prop_mat_c)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat_c))) %*% prop_mat_c
  
  logNormData=Seurat:::NormalizeData.default(inputData$data$countData,normalization.method = "LogNormalize",verbose = F)
  
  
  x_exp=prop_mat %*% t(logNormData)
  if(extendedMode){
    x2_exp=Seurat:::NormalizeData.default(inputData$data$countData,normalization.method = "RC",verbose = F)
    x2_exp=prop_mat %*% t(x2_exp)
  }
  
  y_exp=prop_mat_c %*% t(logNormData)
  
  
  #inputWeight = inputData$prop_mat;inputExp = inputData$data$logNormData
  var_x=myWeightedVar(inputWeight = inputData$prop_mat,inputExp = logNormData,min_cell_count=minCellCountThr)
  var_y=myWeightedVar(inputWeight = prop_mat_c,inputExp = logNormData,min_cell_count=minCellCountThr)
  
  n_x=.myEffSizePropMat(prop_mat = prop_mat)
  n_x=n_x$effective_sample_size
  
  n_y=.myEffSizePropMat(prop_mat = prop_mat_c)
  n_y=n_y$effective_sample_size
  
  vxn=sweep(var_x,1,n_x,"/")
  vyn=sweep(var_y,1,n_y,"/")
  
  vxn2=sweep(var_x,1,n_x-1,"*")
  vyn2=sweep(var_y,1,n_y-1,"*")
  
  
  df=(vxn + vyn)^2
  df=df/(sweep(vxn^2,1,(n_x - 1),"/") + 
           sweep(vyn^2,1,n_y - 1,"/"))
  if(sum(vxn==0,na.rm = T)>0){
    df[Matrix::which(vxn==0)]=matrix(n_y,nrow=nrow(vxn),ncol=ncol(vxn),byrow = F)[Matrix::which(vxn==0)]
  }
  
  vxn[is.na(as.matrix(vxn))]=1
  sxy=sqrt(vxn+vyn)
  
  d_sxy=sqrt(sweep(vxn2+vyn2,1,n_x+n_y-2,"/"))
  rm(vxn2,vyn2)
  
  sd_reg=sxy
  d_sd_reg=d_sxy
  if(reg_sd){
    
    if(F){
      #Vahid's version
      sd_reg=lapply(1:nrow(sxy),function(x){
        y=sxy[x,]
        tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=sxy[x,])
        slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
        if(length(slInd)>0){
          tmp_data=tmp_data[slInd,]
          meanExp=tmp_data$exp
          design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
          fit2=lm.fit(design,log10(tmp_data$sd))
          tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
          
          y[slInd]=tmp_data$fit2
        } else {
          y=pmax(y,1)
        }
        
        return(as.numeric(y))
      })
      sd_reg=do.call("rbind",sd_reg)
    }
    
    
    #Jonah's version
    sd_reg=apply(sxy,1, function(y){
      # browser()
      # y=sxy[x,]
      tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=y)
      slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
      if(length(slInd)>0){
        tmp_data=tmp_data[slInd,]
        meanExp=tmp_data$exp
        design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
        fit2=lm.fit(design,log10(tmp_data$sd))
        tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
        
        y[slInd]=tmp_data$fit2
      } else {
        y=pmax(y,1)
      }
      
      return(as.numeric(y))
    })
    ## browser()
    sd_reg=t(sd_reg)
    
    
    if(F){
      #Vahid's version
      d_sd_reg=lapply(1:nrow(d_sxy),function(x){
        y=d_sxy[x,]
        tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=d_sxy[x,])
        slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
        if(length(slInd)>0){
          tmp_data=tmp_data[slInd,]
          meanExp=tmp_data$exp
          design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
          fit2=lm.fit(design,log10(tmp_data$sd))
          tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
          
          y[slInd]=tmp_data$fit2
        } else {
          y=pmax(y,1)
        }
        
        return(as.numeric(y))
      })
      d_sd_reg=do.call("rbind",d_sd_reg)
    }
    
    
    d_sd_reg=apply(d_sxy, 1, function(y){
      # y=d_sxy[x,]
      tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=y)
      slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
      if(length(slInd)>0){
        tmp_data=tmp_data[slInd,]
        meanExp=tmp_data$exp
        design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
        fit2=lm.fit(design,log10(tmp_data$sd))
        tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
        
        y[slInd]=tmp_data$fit2
      } else {
        y=pmax(y,1)
      }
      
      return(as.numeric(y))
    })
    d_sd_reg = t(d_sd_reg)
    # d_sd_reg=do.call("rbind",d_sd_reg)
    
  }
  rm(d_sxy)
  
  sd_reg=sxy
  sd_offset=0
  
  
  t <- (x_exp - y_exp)/(sd_reg+sd_offset)
  t[which(vxn==0)]=0
  if(F){
    t_adj=sweep(t,1,sqrt((n_x + n_y)/(n_x * n_y)),"*")
    t_adj=t_adj- sweep((3 * t_adj),1,(4 * ((n_x+n_y) - 2) - 1),"/")
    sigmad_adj=1/n_x+1/n_y+sweep((t_adj^2),1,(2*(n_x+n_y)),"/")
  }
  
  if(extendedMode){
    cohen_d=(x_exp - y_exp)/d_sd_reg
    hodge_g=cohen_d*(1-3/(4*(n_x+n_y)-9))
    if(sum(is.na(hodge_g))>0){
      hodge_g[is.na(hodge_g)]=0
    }
    
    se.g=0.5*sweep(hodge_g^2,1,n_x+n_y-3.94,"/")
    se.g <- sqrt(sweep(se.g,1, (n_x+n_y)/(n_x*n_y),"+"))
    if(sum(is.na(se.g))>0){
      se.g[is.na(se.g)]=0
    }
    
  }
  
  
  
  #t2 <- (x_exp - y_exp)/(sd_reg)
  
  t_vec=as.numeric(t)
  zscore=qnorm(pt(as.numeric(abs(t_vec)), as.numeric(df),lower.tail = F,log.p = T),lower.tail = F,log.p = T)
  zscore[is.na(zscore)]=0
  zscore=zscore*sign(t_vec)
  zscore=matrix(zscore,nrow=nrow(t),ncol=ncol(t))
  
  colnames(zscore)=row.names(logNormData)
  row.names(zscore)=row.names(inputData$prop_mat)
  
  exp_binary=t(logNormData)
  exp_binary@x=rep(1,length(exp_binary@x))
  
  pct.1=as.matrix(prop_mat %*% exp_binary)
  pct.2=as.matrix(prop_mat_c %*% exp_binary)
  
  exp_norm=t(expm1(logNormData))
  
  exp.1=log2(prop_mat %*% exp_norm+1)
  exp.2=log2(prop_mat_c %*% exp_norm+1)
  
  logFC=as.matrix(exp.1 - exp.2)
  
  if(extendedMode){
    n=sqrt(sweep(pct.1,1,n_x,"*")*sweep(pct.2,1,n_y,"*"))
  }
  
  
  if(sum(zscore>zscore_cap,na.rm = T)>0){
    zscore[which(zscore>zscore_cap)]=zscore_cap
  }
  
  if(sum(zscore<(-1*zscore_cap),na.rm = T)>0){
    zscore[which(zscore<(-1*zscore_cap))]=(-1*zscore_cap)
  }
  
  
  if(F){
    #moved to the propagation step
    matWeights=.myEffSizePropMat(prop_mat)
    
    matEffectiveSize=matWeights$effective_sample_size
    matWeights=matWeights$centroid_weights
    matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
    matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
    
  }
  
  
  matWeights=matrix(inputData$matWeights,byrow = F,nrow=nrow(zscore),ncol=ncol(zscore))
  matEffectiveSize=matrix(inputData$matEffectiveSize,byrow = F,nrow=nrow(zscore),ncol=ncol(zscore))
  row.names(matWeights)=row.names(matEffectiveSize)=row.names(prop_mat)
  colnames(matWeights)=colnames(matEffectiveSize)=colnames(zscore)
  matWeights[which(pct.1*matEffectiveSize<minCellCountThr&pct.2<0.001)]=0
  matWeights[which(zscore==0|var_x==0)]=0
  
  matEffectiveSize[matWeights<0.0001]=0
  if(extendedMode){
    res=list(zscore=zscore,pct.1=pct.1,pct.2=pct.2,logFC=logFC,se.g=se.g,hodge_g=hodge_g,n=n,matWeights=matWeights,matEffectiveSize=matEffectiveSize,exp1=x_exp,exp2=x2_exp)
  } else {
    res=list(zscore=zscore,pct.1=pct.1,pct.2=pct.2,logFC=logFC,matWeights=matWeights,matEffectiveSize=matEffectiveSize,exp1=x_exp)
  }
  if(F){
    #Jonah's
    output = c(res,list(dsName=inputData$data$dsName,prop_mat=inputData$prop_mat,pseudocell_sim_mat=inputData$pseudocell_sim_mat))
    
    
    res_non_sparse = res
    for(dsitr in c("zscore", "pct.1", "pct.2", "logFC", "n", "matWeights", "matEffectiveSize")){
      res[[dsitr]] = as.sparse(res[[dsitr]])
    }
    gc()
    
    dir.create(file.path(argList$saveDir, "detail_exp_final"), showWarnings = FALSE)
    FN = file.path(argList$saveDir, "detail_exp_final", paste0(inputData$data$dsName, ".RDS"))
    
    .mcsaveRDS <- function(object,file,mc.cores=min(parallel::detectCores(),3, na.rm=T)) {
      con <- pipe(paste0("pigz -p",mc.cores," > ",file),"wb")
      saveRDS(object, file = con)
      close(con)
    }
    print(FN)
    
    .mcsaveRDS(output, FN)
    
    
  }
  
  setdiff(names(inputData) ,c("matWeights","matEffectiveSize"))
  return(c(res,list(dsName=inputData$data$dsName)))
  #return(T)
}

.sconline.RobustFC=function(inputData,batch_col,contrast_col,contrast_1="Abeta",contrast_2="Ctrl",sex_col=NULL,pct_sig_thr=0,logFC_sig_thr=0,ncores=5,groupLevel=F){
  #groupLevel: AD vs all Ctrl
  
  pct_sig_thr=abs(pct_sig_thr)
  logFC_sig_thr=abs(logFC_sig_thr)
  
  mySecondFn=function(j,i,ad_pd,sl_ctrl,data,batch_col){
    tmp_data=data[,which(colData(data)[,batch_col] %in% c(ad_pd[i,batch_col],sl_ctrl[j,batch_col]))]
    tmp_data=suppressWarnings(.extraExport2SeuratFn(inputData = tmp_data))
    tmp_data=suppressWarnings(Seurat::NormalizeData(tmp_data,verbose=F))
    tmp_res=.myEvalMarkers(object=tmp_data, cells.1=colnames(tmp_data)[which(tmp_data@meta.data[,batch_col]==ad_pd[i,batch_col])], cells.2=colnames(tmp_data)[which(tmp_data@meta.data[,batch_col]==sl_ctrl[j,batch_col])], slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cells.1.weight.col=NULL,cluster_name=NULL)
    tmp_res$library_AD=ad_pd$library[i]
    tmp_res$library_Ctrl=sl_ctrl$library[j]
    tmp_res$sbj_AD=ad_pd$anno_batch[i]
    tmp_res$sbj_Ctrl=sl_ctrl$anno_batch[j]
    tmp_res$gene=row.names(tmp_res)
    return(tmp_res)
  }
  
  myThirdFn=function(i,ad_pd,sl_ctrl,data,batch_col){
    tmp_data=data[,which(colData(data)[,batch_col] %in% c(ad_pd[i,batch_col],sl_ctrl[,batch_col]))]
    tmp_data=suppressWarnings(.extraExport2SeuratFn(inputData = tmp_data))
    tmp_data=suppressWarnings(Seurat::NormalizeData(tmp_data,verbose=F))
    tmp_res=.myEvalMarkers(object=tmp_data, cells.1=colnames(tmp_data)[which(tmp_data@meta.data[,batch_col]==ad_pd[i,batch_col])], cells.2=colnames(tmp_data)[which(tmp_data@meta.data[,batch_col] %in% as.character(sl_ctrl[,batch_col]))], slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cells.1.weight.col=NULL,cluster_name=NULL)
    tmp_res$library_AD=ad_pd$library[i]
    tmp_res$sbj_AD=ad_pd$anno_batch[i]
    tmp_res$gene=row.names(tmp_res)
    return(tmp_res)
  }
  
  if(groupLevel){
    stop("Needs to be implemented!")
  }
  
  
  if(is.null(inputData)){
    stop("Expression data is missing!")
  } else if(tolower(class(inputData))=="seurat"){
    inputData=.sconline.convertSeuratToExpressionSet(object=inputData)
  } else if(class(inputData)!="SingleCellExperiment") {
    stop("Unrecognized inputExpression data!")
  }
  
  colData(inputData)[,batch_col]=as.character(colData(inputData)[,batch_col])
  
  tmp_pd=as.data.frame(colData(inputData))
  tmp_pd$status=tmp_pd[,contrast_col]
  tmp_pd=tmp_pd[!duplicated(tmp_pd[,batch_col]),]
  ad_pd=tmp_pd[which(tmp_pd$status==contrast_1),]
  if(nrow(ad_pd)==0){
    stop("contrast_1 was not identified!")
  }
  ctrl_pd=tmp_pd[which(tmp_pd$status==contrast_2),]
  if(nrow(ctrl_pd)==0){
    stop("contrast_2 was not identified!")
  }
  res_robust_FC=list()
  
  for(i in 1:nrow(ad_pd)){
    if(!is.null(sex_col)){
      sl_ctrl=ctrl_pd[which(ctrl_pd$anno_sex==ad_pd$anno_sex[i]),]
    } else {
      sl_ctrl=ctrl_pd
    }
    
    #sl_ctrl=sl_ctrl[which(sl_ctrl$anno_age>=0.9*ad_pd$anno_age[i]&sl_ctrl$anno_age<=1.1*ad_pd$anno_age[i]),]
    if(nrow(sl_ctrl)>0){
      if(groupLevel){
        tmp_res=myThirdFn(i=i,ad_pd=ad_pd,sl_ctrl=sl_ctrl,data=inputData,batch_col=batch_col)
        tmp_res=list(tmp_res)
      } else {
        tmp_res=parallel::mclapply(1:nrow(sl_ctrl),mySecondFn,i=i,ad_pd=ad_pd,sl_ctrl=sl_ctrl,data=inputData,batch_col=batch_col,mc.cores = ncores)
        for(isbj in 1:length(tmp_res)){
          tmp_res[[isbj]]$pct_diff=tmp_res[[isbj]]$pct.1 - tmp_res[[isbj]]$pct.2
        }
        res_count=data.frame(gene=tmp_res[[1]]$gene,cell.count.1=tmp_res[[1]]$cell.count.1,pct_diff_count_up=0,pct_diff_count_down=0,logFC_count_up=0,logFC_count_down=0,ref_count=0)
        if(length(tmp_res)>0){
          for(isbj in 1:length(tmp_res)){
            tmp_count=tmp_res[[isbj]]
            tmp_count=tmp_count[which(tmp_count$cell.count.2>10),]
            tmp_count=tmp_count[match(res_count$gene,tmp_count$gene),]
            tmp_count$pct_diff=tmp_count$pct.1 - tmp_count$pct.2
            tmp_count$pct_diff[is.na(tmp_count$pct_diff)]=0
            tmp_count$avg_logFC[is.na(tmp_count$avg_logFC)]=0
            tmp_count$count=1
            tmp_count$count[is.na(tmp_count$gene)]=0
            res_count$ref_count=res_count$ref_count+tmp_count$count
            res_count$pct_diff_count_up=res_count$pct_diff_count_up+as.numeric(tmp_count$pct_diff>pct_sig_thr)*tmp_count$count
            res_count$pct_diff_count_down=res_count$pct_diff_count_down+as.numeric(tmp_count$pct_diff<(-1*pct_sig_thr))*tmp_count$count
            res_count$logFC_count_up=res_count$logFC_count_up+as.numeric(tmp_count$avg_logFC>logFC_sig_thr)*tmp_count$count
            res_count$logFC_count_down=res_count$logFC_count_down+as.numeric(tmp_count$avg_logFC<(-1*logFC_sig_thr))*tmp_count$count
            
          }
        }
        
        res_count$score_pct=(res_count$pct_diff_count_up - res_count$pct_diff_count_down)/pmax(res_count$ref_count,1)
        res_count$score_logFC=(res_count$logFC_count_up - res_count$logFC_count_down)/pmax(res_count$ref_count,1)
      }
      
      
      res_robust_FC=c(res_robust_FC,list(res_count))
    }
    
  }
  return(res_robust_FC)
}


.myConcensusDEFn_step2_detail_exp_final=function(inputData,argList,sd_offset=0.001,reg_sd=T,zscore_cap=10,minCellCountThr=4,extendedMode=F){
  
  #inputData=res[[1]];sd_offset=0.005;reg_sd=T;zscore_cap=10;minCellCountThr=4;extendedMode=F
  if(F){
    #Jonah's
    write.table(data.frame(), file=paste0("/tmp/touch_", as.character(runif(1)*10^20)),
                col.names=FALSE)
    
    library("RhpcBLASctl")
    omp_set_num_threads(4)
    blas_set_num_threads(4)
    
  }
  #require(bigmemory)
  
  #inputData=res[[1]];sd_offset=0.001;reg_sd=T;zscore_cap=10;minCellCountThr=4
  
  #tmp_Propweights=rowSums(as.matrix(inputData$prop_mat))
  
  if(sum(names(argList)=="do.split.prop")==0){
    argList$do.split.prop=T
  }
  
  inputData$prop_mat=.extra_matrix_rowNorm(inputData$prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(inputData$prop_mat)+0.000000000001)) %*% inputData$prop_mat
  
  
  require(Matrix,quietly = T)
  
  zscore_cap=abs(zscore_cap)
  if(zscore_cap<2){
    warning("zscore cap was set to below 2. It was changed to 5.")
    zscore_cap=5
  }
  sd_offset=max(sd_offset,0)
  
  if(F){
    resGeneMeanSd=.extra_gmean(inputData$data$countData)
    resGeneMeanSd=data.frame(gene=names(resGeneMeanSd),geo_mean_count=resGeneMeanSd,stringsAsFactors = F)
  }
  
  
  
  myWeightedVar=function(inputWeight,inputExp,min_cell_count,mergeSimilars=F){
    #inputWeight = prop_mat_c;inputExp = logNormData;min_cell_count=minCellCountThr
    require(Matrix,quietly = T)
    mapping=data.frame(cluster=row.names(inputWeight),pseudocell=row.names(inputWeight),stringsAsFactors = F)
    inputWeight_reduced=inputWeight
    
    if(mergeSimilars){
      row_sum_counts=inputWeight
      row_sum_counts@x=rep(1,length(row_sum_counts))
      row_sum_counts=rowSums(row_sum_counts)
      
      tmp_cosine=qlcMatrix::cosSparse(t(inputWeight))
      row.names(tmp_cosine)=colnames(tmp_cosine)=row.names(inputWeight)
      tmp_cosine_clust=hclust(as.dist(1-as.matrix(tmp_cosine)),method = "complete")
      tmp_clust=cutree(tmp_cosine_clust,h = 0.99)
      mapping=data.frame(clust_id=tmp_clust,pseudocell=names(tmp_clust),stringsAsFactors = F)
      map_names=mapping[!duplicated(mapping[,1]),]
      colnames(map_names)[2]="cluster"
      mapping=merge(mapping,map_names,by="clust_id")
      if(sum(row_sum_counts<10000)>0){
        mapping$cluster[mapping$pseudocell %in% row.names(inputWeight)[row_sum_counts<10000]]=mapping$pseudocell[mapping$pseudocell %in% row.names(inputWeight)[row_sum_counts<10000]]
      }
      
      inputWeight_reduced=inputWeight[row.names(inputWeight) %in% mapping$cluster,,drop=F]
    }
    #inputExp=t(inputExp)
    res_binary=(sign(inputWeight_reduced)%*% sign(inputExp))
    exp_sq=inputExp^2
    res=(inputWeight_reduced%*% exp_sq)
    res=sweep(res,1,rowSums(inputWeight_reduced),"*")
    res=res- (inputWeight_reduced%*% inputExp)^2
    res=sweep(res,1,rowSums(inputWeight_reduced)^2,"/")
    
    singletons=apply(inputWeight_reduced,1,function(x) sum(x>0))
    singletons[singletons<=min_cell_count]=0
    singletons[singletons>0]=1
    if(sum(singletons<1)>0){
      res=sweep(res,1,singletons,"*")
    }
    res[res_binary<=min_cell_count]=0
    res=res[match(mapping$cluster,row.names(res)),]
    row.names(res)=mapping$pseudocell
    res=res[match(row.names(inputWeight),row.names(res)),]
    return(res)
  }
  
  prop_mat=inputData$prop_mat
  
  
  prop_mat_c=prop_mat
  if(quantile(rowSums(prop_mat_c > 0,na.rm = T), 0.25) < (0.85*ncol(prop_mat))){
    prop_mat_c=Matrix::drop0(prop_mat_c)
    prop_mat_c@x=rep(1,length(prop_mat_c@x))
  }
  
  prop_mat_c=Matrix::drop0(1-prop_mat_c)
  prop_mat_c <- .extra_matrix_rowNorm(prop_mat_c)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat_c))) %*% prop_mat_c
  
  logNormData=t(inputData$data$logNorm)
  logNormData=logNormData[colnames(prop_mat),,drop=F]
  
  x_exp=prop_mat %*% logNormData
  if(extendedMode){
    x2_exp=Seurat:::NormalizeData.default(inputData$data$countData,normalization.method = "RC",verbose = F)
    x2_exp=prop_mat %*% t(x2_exp)
  }
  
  gene_name_list=colnames(logNormData)
  pseudocell_name_list=row.names(prop_mat)
  
  #inputWeight = inputData$prop_mat;inputExp = logNormData;min_cell_count=minCellCountThr
  var_x=myWeightedVar(inputWeight = inputData$prop_mat,inputExp = logNormData,min_cell_count=minCellCountThr)
  var_y=myWeightedVar(inputWeight = prop_mat_c,inputExp = logNormData,min_cell_count=minCellCountThr,mergeSimilars = T)
  
  n_x=.myEffSizePropMat(prop_mat = prop_mat)
  n_x=n_x$effective_sample_size
  
  n_y=.myEffSizePropMat(prop_mat = prop_mat_c)
  n_y=n_y$effective_sample_size
  
  vxn=sweep(var_x,1,n_x,"/")
  vyn=sweep(var_y,1,n_y,"/")
  
  vxn2=sweep(var_x,1,n_x-1,"*")
  vyn2=sweep(var_y,1,n_y-1,"*")
  
  
  df=(vxn + vyn)^2
  df=df/(sweep(vxn^2,1,(n_x - 1),"/") + 
           sweep(vyn^2,1,n_y - 1,"/"))
  if(sum(vxn==0,na.rm = T)>0){
    df[Matrix::which(vxn==0)]=matrix(n_y,nrow=nrow(vxn),ncol=ncol(vxn),byrow = F)[Matrix::which(vxn==0)]
  }
  
  vxn[is.na(as.matrix(vxn))]=1
  sxy=sqrt(vxn+vyn)
  rm(vyn)
  
  d_sxy=sqrt(sweep(vxn2+vyn2,1,n_x+n_y-2,"/"))
  rm(vxn2,vyn2)
  
  sd_reg=sxy
  d_sd_reg=d_sxy
  if(F){
    
    if(reg_sd& extendedMode){
      
      if(F){
        #Vahid's version
        sd_reg=lapply(1:nrow(sxy),function(x){
          y=sxy[x,]
          tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=sxy[x,])
          slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
          if(length(slInd)>0){
            tmp_data=tmp_data[slInd,]
            meanExp=tmp_data$exp
            design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
            fit2=lm.fit(design,log10(tmp_data$sd))
            tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
            
            y[slInd]=tmp_data$fit2
          } else {
            y=pmax(y,1)
          }
          
          return(as.numeric(y))
        })
        sd_reg=do.call("rbind",sd_reg)
      }
      
      
      #Jonah's version
      sd_reg=apply(sxy,1, function(y){
        # browser()
        # y=sxy[x,]
        tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=y)
        slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
        if(length(slInd)>0){
          tmp_data=tmp_data[slInd,]
          meanExp=tmp_data$exp
          design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
          fit2=lm.fit(design,log10(tmp_data$sd))
          tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
          
          y[slInd]=tmp_data$fit2
        } else {
          y=pmax(y,1)
        }
        
        return(as.numeric(y))
      })
      ## browser()
      sd_reg=t(sd_reg)
      
      
      if(F){
        #Vahid's version
        d_sd_reg=lapply(1:nrow(d_sxy),function(x){
          y=d_sxy[x,]
          tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=d_sxy[x,])
          slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
          if(length(slInd)>0){
            tmp_data=tmp_data[slInd,]
            meanExp=tmp_data$exp
            design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
            fit2=lm.fit(design,log10(tmp_data$sd))
            tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
            
            y[slInd]=tmp_data$fit2
          } else {
            y=pmax(y,1)
          }
          
          return(as.numeric(y))
        })
        d_sd_reg=do.call("rbind",d_sd_reg)
      }
      
      
      if(extendedMode){
        d_sd_reg=apply(d_sxy, 1, function(y){
          # y=d_sxy[x,]
          tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=y)
          slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
          if(length(slInd)>0){
            tmp_data=tmp_data[slInd,]
            meanExp=tmp_data$exp
            design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
            fit2=lm.fit(design,log10(tmp_data$sd))
            tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
            
            y[slInd]=tmp_data$fit2
          } else {
            y=pmax(y,1)
          }
          
          return(as.numeric(y))
        })
        d_sd_reg = t(d_sd_reg)
      }
      
      # d_sd_reg=do.call("rbind",d_sd_reg)
      
    }
    
  }
  
  
  rm(sxy,d_sxy)
  sd_offset=0
  
  ##########
  #reduced prop_mat_c
  row_sum_counts=prop_mat_c
  row_sum_counts@x=rep(1,length(row_sum_counts))
  row_sum_counts=rowSums(row_sum_counts)
  
  tmp_cosine=qlcMatrix::cosSparse(t(prop_mat_c))
  row.names(tmp_cosine)=colnames(tmp_cosine)=row.names(prop_mat_c)
  tmp_cosine_clust=hclust(as.dist(1-as.matrix(tmp_cosine)),method = "complete")
  tmp_clust=cutree(tmp_cosine_clust,h = 0.999)
  mapping=data.frame(clust_id=tmp_clust,pseudocell=names(tmp_clust),stringsAsFactors = F)
  map_names=mapping[!duplicated(mapping[,1]),]
  colnames(map_names)[2]="cluster"
  mapping=merge(mapping,map_names,by="clust_id")
  if(sum(row_sum_counts<10000)>0){
    mapping$cluster[mapping$pseudocell %in% row.names(prop_mat_c)[row_sum_counts<10000]]=mapping$pseudocell[mapping$pseudocell %in% row.names(prop_mat_c)[row_sum_counts<10000]]
  }
  prop_mat_c_reduced=prop_mat_c[row.names(prop_mat_c) %in% mapping$cluster,,drop=F]
  
  y_exp=(prop_mat_c_reduced %*% logNormData)
  y_exp=y_exp[match(mapping$cluster,row.names(y_exp)),]
  row.names(y_exp)=row.names(prop_mat_c)
  y_exp=Matrix::drop0(y_exp)
  
  
  
  
  t <- (x_exp - y_exp)/(sd_reg+sd_offset)
  t[which(vxn==0)]=0
  if(F){
    t_adj=sweep(t,1,sqrt((n_x + n_y)/(n_x * n_y)),"*")
    t_adj=t_adj- sweep((3 * t_adj),1,(4 * ((n_x+n_y) - 2) - 1),"/")
    sigmad_adj=1/n_x+1/n_y+sweep((t_adj^2),1,(2*(n_x+n_y)),"/")
  }
  
  
  if(extendedMode){
    cohen_d=(x_exp - y_exp)/d_sd_reg
    hodge_g=cohen_d*(1-3/(4*(n_x+n_y)-9))
    if(sum(is.na(hodge_g))>0){
      hodge_g[is.na(hodge_g)]=0
    }
    
    se.g=0.5*sweep(hodge_g^2,1,n_x+n_y-3.94,"/")
    se.g <- sqrt(sweep(se.g,1, (n_x+n_y)/(n_x*n_y),"+"))
    if(sum(is.na(se.g))>0){
      se.g[is.na(se.g)]=0
    }
    
  }
  rm(vxn,sd_reg,y_exp)
  
  x_exp=as(x_exp, "dgCMatrix")
  gc()
  #t2 <- (x_exp - y_exp)/(sd_reg)
  
  t_vec=as.numeric(t)
  zscore=qnorm(pt(as.numeric(abs(t_vec)), as.numeric(df),lower.tail = F,log.p = T),lower.tail = F,log.p = T)
  zscore[is.na(zscore)]=0
  zscore=zscore*sign(t_vec)
  zscore=as(matrix(zscore,nrow=nrow(t),ncol=ncol(t)), "dgCMatrix")
  #zscore_identity=zscore
  #zscore_identity@x=rep(1,length(zscore_identity@x))
  rm(t_vec,df)
  
  colnames(zscore)=colnames(logNormData)
  row.names(zscore)=row.names(inputData$prop_mat)
  
  exp_binary=logNormData#t(logNormData)
  exp_binary@x=rep(1,length(exp_binary@x))
  
  pct.1=prop_mat %*% exp_binary
  {
    pct.2=(prop_mat_c_reduced %*% exp_binary)
    pct.2=pct.2[match(mapping$cluster,row.names(pct.2)),]
    row.names(pct.2)=row.names(prop_mat_c)
    #pct.2=pct.2*zscore_identity
    pct.2=Matrix::drop0(pct.2)
  }
  
  
  rm(exp_binary)
  
  exp_norm=expm1(logNormData)#t(expm1(logNormData))
  
  logFC=log2(prop_mat_c_reduced %*% exp_norm+1)
  logFC=logFC[match(mapping$cluster,row.names(logFC)),]
  logFC=log2(prop_mat %*% exp_norm+1) - logFC
  rm(exp_norm,logNormData)
  
  if(extendedMode){
    n=sqrt(sweep(pct.1,1,n_x,"*")*sweep(pct.2,1,n_y,"*"))
  }
  rm(n_y)
  
  if(sum(zscore>zscore_cap,na.rm = T)>0){
    zscore[which(zscore>zscore_cap)]=zscore_cap
  }
  
  if(sum(zscore<(-1*zscore_cap),na.rm = T)>0){
    zscore[which(zscore<(-1*zscore_cap))]=(-1*zscore_cap)
  }
  
  
  if(F){
    #moved to the propagation step
    matWeights=.myEffSizePropMat(prop_mat)
    
    matEffectiveSize=matWeights$effective_sample_size
    matWeights=matWeights$centroid_weights
    matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
    matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
    
  }
  
  
  matWeights=matrix(inputData$matWeights,byrow = F,nrow=nrow(zscore),ncol=ncol(zscore))
  matEffectiveSize=matrix(inputData$matEffectiveSize,byrow = F,nrow=nrow(zscore),ncol=ncol(zscore))
  row.names(matWeights)=row.names(matEffectiveSize)=row.names(prop_mat)
  colnames(matWeights)=colnames(matEffectiveSize)=colnames(zscore)
  matWeights[which(pct.1*matEffectiveSize<minCellCountThr&pct.2<0.001)]=0
  matWeights[which(zscore==0|var_x==0)]=0
  
  
  matEffectiveSize[matWeights<0.0001]=0
  if(extendedMode){
    res=list(zscore=as(zscore, "dgCMatrix"),pct.1=as.big.matrix(pct.1),pct.2=as.big.matrix(pct.2),logFC=logFC,se.g=as.big.matrix(se.g),hodge_g=hodge_g,n=n,matWeights=matWeights,matEffectiveSize=matEffectiveSize,exp1=x_exp,exp2=x2_exp)
  } else {
    res=list(zscore=zscore,pct.1=pct.1,pct.2=pct.2,logFC=logFC,matWeights=matWeights,matEffectiveSize=matEffectiveSize,exp1=x_exp,gene_name_list=gene_name_list,pseudocell_name_list=pseudocell_name_list)
  }
  
  gc()
  setdiff(names(inputData) ,c("matWeights","matEffectiveSize"))
  return(c(res,list(dsName=inputData$data$dsName)))
  #return(T)
}


.extra_sconline.argFn=function(runIndx,do.split.prop=T,min_cluster_size,pseudocell_selection_method,saveDir,exNonMicCells=F,prop.n.neighbors = 4,ncores=7,conservativeMapping=F,oldMapping=F,MGIsymbol=F,includeHuman,includeMouse,FinnishSbjBased,DE_supportingFractionThr,DE_n.adaptiveKernel,DE_nPropIter,Leng200=F,uniformZscore=F,dist_zscore_gamma,dist_zscore_norm,dist_zscore_nbinom,regularize,geoMean=F,prefix=NULL,newRun=F,inputDf=NULL,internal_pseudocell_count,pseudocell_size,sensitiveSearch,external_DE_path=NULL,external_DE_name=NULL,do.liger=NULL,singleton.method=singleton.method,include.singletons=T,colNormalize=T){
  
  if(is.null(inputDf)){
    if(includeMouse&(!includeHuman)){
      .dfSetting=.myRunSettingsComplete_mouse()
    } else if((!includeMouse)&(includeHuman)){
      .dfSetting=.myRunSettingsComplete_human()
    } else {
      .dfSetting=.myRunSettingsComplete_mouse()
    }
  } else {
    .dfSetting=inputDf
  }
  
  
  .commonExpressed=.dfSetting$commonExpressed[runIndx]
  .includeHuman=includeHuman
  .includeMouse=includeMouse
  .nPCs=.dfSetting$nPCs[runIndx]
  .HVG_count=.dfSetting$HVG_count[runIndx] #4
  .HVG_method="vst" #"vst" and "sctransform"
  .exNonMicCells=exNonMicCells
  .covariates=.dfSetting$covariates[runIndx]
  .removeHighExp=.dfSetting$removeHighExp[runIndx]
  .UMI_cor_thr=.dfSetting$UMI_cor_thr[runIndx]
  .slMicrogliaClusters=.dfSetting$slMicrogliaClusters[runIndx]
  .breakHammond=T
  if(sum(colnames(.dfSetting)=="breakHammond")>0){
    .dfSetting$breakHammond[runIndx]
  }
  if(.slMicrogliaClusters==""){
    .slMicrogliaClusters=NULL
  }
  
  if(!is.null(.slMicrogliaClusters)){
    .exNonMicCells=T
    .slMicrogliaClusters=unlist(strsplit(.slMicrogliaClusters,","))
  }
  
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
  
  
  .exCellStates=c("CAM","Mnc","BAM","T cells","Neutrophils","cDC","B cells","NK cells","Proliferation","T/NKT cells","cDC2","Monocyte/Mdc","ILC","ydT cells","pDC","cDC1","migDC","MCs","Non-cl. monocytes")
  
  #"depth_per_gene"#'QC_Gene_unique_count'#c('QC_Gene_total_count','QC_top50_pct',)#, "organism","depth_per_gene" 'QC_Gene_total_count','QC_top50_pct','QC_Gene_unique_count'
  .indScaling=.dfSetting$indScaling[runIndx]
  
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
  
  if(is.null(do.liger)){
    do.liger=.dfSetting$do.liger[runIndx]
  }
  
  if(do.liger){
    prefix=paste0(prefix,"_liger")
  }
  
  if(!is.null(external_DE_name)){
    prefix=paste0(prefix,"_",external_DE_name)
  }
  
  if(Leng200){
    prefix=paste0(prefix,"_Leng200")
  }
  
  .saveDir=.mySaveDirMaker(paste0(saveDir,"/",prefix),nPCs = .nPCs,cov=.covariates,exNonMicCells=.exNonMicCells,FinnishSbjBased)
  if(.exNonMicCells){
    .saveDirGlobal=paste0(saveDir,"/",prefix,"-Global-exNonMicCells")
  } else {
    .saveDirGlobal=paste0(saveDir,"/",prefix,"-Global")
  }
  
  if(FinnishSbjBased&includeHuman){
    .saveDirGlobal=paste0(.saveDirGlobal,"-FinnishSbjBased")
  }
  
  argList=list(commonExpressed=.dfSetting$commonExpressed[runIndx],min_cluster_size=min_cluster_size,pseudocell_selection_method=pseudocell_selection_method,do.split.prop=do.split.prop,prop.n.neighbors = prop.n.neighbors,includeHuman=.includeHuman,includeMouse=.includeMouse,conservativeMapping=conservativeMapping,oldMapping=oldMapping,MGIsymbol=MGIsymbol,nPCs=.dfSetting$nPCs[runIndx],HVG_count=.dfSetting$HVG_count[runIndx],HVG_method="vst",include.singletons=include.singletons,colNormalize=colNormalize,exNonMicCells=exNonMicCells,covariates=.covariates,indScaling=.dfSetting$indScaling[runIndx],saveDir=.saveDir,saveDirGlobal=.saveDirGlobal,ncores=ncores,excludeHighExp=.removeHighExp,slMicrogliaClusters=.slMicrogliaClusters,breakHammond=.breakHammond,UMI_cor_thr=.UMI_cor_thr,onlyRankBased=.dfSetting$onlyRankBased[runIndx],varScore.thr=.dfSetting$varScore.thr[runIndx],do.liger=.dfSetting$do.liger[runIndx],allGenesFraction=.dfSetting$allGenesFraction[runIndx],FinnishSbjBased=FinnishSbjBased,uniformZscore=uniformZscore,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=DE_nPropIter,dist_zscore_gamma=dist_zscore_gamma,dist_zscore_norm=dist_zscore_norm,dist_zscore_nbinom=dist_zscore_nbinom,regularize=regularize,geoMean=geoMean,prefix=prefix,internal_pseudocell_count=internal_pseudocell_count,pseudocell_size=pseudocell_size,sensitiveSearch=sensitiveSearch,external_DE_path=external_DE_path,external_DE_name=external_DE_name,singleton.method=singleton.method)
  
  .extraCreateDirFn(argList)
  
  .extraNewRun(argList=argList,newRun=newRun)
  
  return(argList)
  
}

#inputData=res;argList=.ArgList;min_effective_size=5;pd=pd;cell_annoCol=annoCol
.extra_sconline.visPseudocellAnno=function(inputData,argList,min_effective_size=5,point_size=3,pd=NULL,cell_annoCol=NULL,pie_scale=1,piechart.plot=NULL,return_obj=F){
  require(ggplot2)
  require(scatterpie)
  require(hues)
  
  if(F){
    if(sum(colnames(inputData)=="dataset")>0){
      inputData=inputData[,-which(colnames(inputData)=="dataset")]
      absent_pseudocells=apply(inputData[,-which(colnames(inputData) %in% c("pseudocell","effective_size"))],1,sum)
      absent_pseudocells=which(absent_pseudocells<0.2)
      if(length(absent_pseudocells)>0){
        inputData[absent_pseudocells,-which(colnames(inputData)=="pseudocell")]=NA
      }
      
      inputData=aggregate(.~pseudocell,data=inputData,function(x) mean(x,na.rm=T))
    }
  }
  
  
  UMAP_centroid=.sconline.fetch_data("umap_pseudocells",argList=argList)
  
  if(is.null(pd)){
    pd=.sconline.fetch_data("annotation",argList=argList)
  }
  
  if(is.null(cell_annoCol)){
    pd_summary=.extra_sconline.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"),n=300)
  } else {
    pd_summary=pd
  }
  
  
  if(sum(colnames(inputData)=="dataset")>0){
    inputData=inputData[which(inputData$effective_size>min_effective_size),]
    
  }
  piechart_data=merge(inputData,UMAP_centroid,by.x="pseudocell",by.y="centroid",all.x=T)
  
  anno_cols=setdiff(colnames(inputData),c("pseudocell","effective_size","dataset"))
  
  if(is.null(piechart.plot)){
    do.scatterPlot=F
    if(length(anno_cols)>1){
      do.scatterPlot=T
    }
  } else {
    do.scatterPlot=(!piechart.plot)
  }
  
  
  
  
  if(is.null(cell_annoCol)){
    if(do.scatterPlot){
      color_num=length(anno_cols)
      
      
      if(sum(colnames(inputData)=="dataset")>0){
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                  cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(color_num)))+facet_wrap(~dataset)
        
        
      } else {
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                  cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(length(anno_cols))))
        
      }
    } else {
      if(sum(colnames(inputData)=="dataset")>0){
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_point(aes_string(x='UMAP_1', y='UMAP_2', color=anno_cols), data=piechart_data) + coord_equal()+theme_classic()+scale_fill_gradient(low="white",high="red")+facet_wrap(~dataset)
        
        
      } else {
        
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_point(aes_string(x='UMAP_1', y='UMAP_2', fill=anno_cols),color="black", data=piechart_data,shape=21,size=point_size) + coord_equal()+theme_classic()+scale_fill_gradient(low="white",high="red")
        
      }
    }
    
  } else {
    pd$cluster_anno_res=pd[,cell_annoCol]
    if(do.scatterPlot){
      color_num=max(length(unique(pd$cluster_anno_res)),length(anno_cols))
      if(length(setdiff(pd$cluster_anno_res,anno_cols))==0){
        pd$cluster_anno_res=factor(as.character(pd$cluster_anno_res),levels = anno_cols)
      } else {
        pd$cluster_anno_res=factor(as.character(pd$cluster_anno_res))
      }
      
      if(sum(colnames(inputData)=="dataset")>0){
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=cluster_anno_res),size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                  cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c(hues::iwanthue(color_num)))+scale_color_manual(values=c(hues::iwanthue(length(anno_cols))))+facet_wrap(~dataset)
        
        
      } else {
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=cluster_anno_res),size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                  cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c(hues::iwanthue(length(unique(pd$cluster_anno_res)))))+scale_color_manual(values=c(hues::iwanthue(length(anno_cols))))
        
      }
    } else {
      if(sum(colnames(inputData)=="dataset")>0){
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=cluster_anno_res),size=0.3) + geom_point(aes(x=UMAP_1, y=UMAP_2, fill=cluster),color="black",shape=21)+
                                                                                                             coord_equal()+theme_classic()+scale_fill_manual(values=c(hues::iwanthue(length(unique(pd$cluster_anno_res)))))+facet_wrap(~dataset)+scale_color_manual(values=c(hues::iwanthue(length(unique(pd$cluster_anno_res)))))
        
        
      } else {
        color_num=max(length(piechart_data$cluster),length(unique(pd_summary$cluster_anno_res)))
        p=.myDimPlotFn(object=pd_summary, dimCols = c("UMAP_1", "UMAP_2"), attCol = "cluster_anno_res",set_col=F) + geom_point(data=piechart_data,aes(x=UMAP_1, y=UMAP_2, fill=cluster),color="black",shape=21,size=point_size)+ 
                                                                                                             coord_equal()+theme_classic()+scale_fill_gradient(low="white",high="red")+scale_color_manual(values=c(hues::iwanthue(length(unique(pd$cluster_anno_res)))))
        
      }
    }
  }
  
  if(return_obj){
    return(list(pseudocell=piechart_data,cell=pd_summary))
  } else {
    return(p)
  }
  
}

.extra_sconline.visPseudocellAnno_cluster=function(inputData,argList,min_effective_size=5,pd=NULL,cell_annoCol=NULL,pie_scale=1,piechart.plot=NULL,return_obj=F){
  
  #pd=NULL;cell_annoCol=NULL;pie_scale=1;piechart.plot=NULL;return_obj=F
  require(ggplot2)
  require(scatterpie)
  require(hues)
  
  if(F){
    if(sum(colnames(inputData)=="dataset")>0){
      inputData=inputData[,-which(colnames(inputData)=="dataset")]
      absent_pseudocells=apply(inputData[,-which(colnames(inputData) %in% c("pseudocell","effective_size"))],1,sum)
      absent_pseudocells=which(absent_pseudocells<0.2)
      if(length(absent_pseudocells)>0){
        inputData[absent_pseudocells,-which(colnames(inputData)=="pseudocell")]=NA
      }
      
      inputData=aggregate(.~pseudocell,data=inputData,function(x) mean(x,na.rm=T))
    }
  }
  
  
  UMAP_centroid=.sconline.fetch_data("umap_pseudocells",argList=argList)
  
  if(is.null(pd)){
    pd=.sconline.fetch_data("annotation",argList=argList)
  }
  
  if(is.null(cell_annoCol)){
    pd_summary=.extra_sconline.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"),n=300)
  } else {
    pd_summary=pd
  }
  
  
  if(sum(colnames(inputData)=="dataset")>0){
    inputData=inputData[which(inputData$effective_size>min_effective_size),]
    
  }
  piechart_data=merge(inputData,UMAP_centroid,by.x="pseudocell",by.y="centroid",all.x=T)
  
  anno_cols=setdiff(colnames(inputData),c("pseudocell","effective_size","dataset"))
  
  if(is.null(piechart.plot)){
    do.scatterPlot=F
    if(length(anno_cols)>1){
      do.scatterPlot=T
    }
  } else {
    do.scatterPlot=(!piechart.plot)
  }
  
  
  
  
  if(is.null(cell_annoCol)){
    if(do.scatterPlot){
      color_num=length(anno_cols)
      
      
      if(sum(colnames(inputData)=="dataset")>0){
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                  cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(color_num)))+facet_wrap(~dataset)
        
        
      } else {
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                  cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(length(anno_cols))))
        
      }
    } else {
      if(sum(colnames(inputData)=="dataset")>0){
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=cluster), data=piechart_data) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(length(unique(pd$cluster)))))+facet_wrap(~dataset)+scale_color_manual(values=c("gray",hues::iwanthue(length(unique(pd$cluster)))))
        
        
      } else {
        
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_point(aes(x=UMAP_1, y=UMAP_2, fill=cluster),color="black", data=piechart_data) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(length(unique(pd$cluster)))))+scale_color_manual(values=c("gray",hues::iwanthue(length(unique(pd$cluster)))))
        
      }
    }
    
  } else {
    pd$cluster_anno_res=pd[,cell_annoCol]
    if(do.scatterPlot){
      color_num=max(length(unique(pd$cluster_anno_res)),length(anno_cols))
      if(length(setdiff(pd$cluster_anno_res,anno_cols))==0){
        pd$cluster_anno_res=factor(as.character(pd$cluster_anno_res),levels = anno_cols)
      } else {
        pd$cluster_anno_res=factor(as.character(pd$cluster_anno_res))
      }
      
      if(sum(colnames(inputData)=="dataset")>0){
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=cluster_anno_res),size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                       cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c(hues::iwanthue(color_num)))+scale_color_manual(values=c(hues::iwanthue(length(anno_cols))))+facet_wrap(~dataset)
        
        
      } else {
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=cluster_anno_res),size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                       cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c(hues::iwanthue(length(unique(pd$cluster_anno_res)))))+scale_color_manual(values=c(hues::iwanthue(length(anno_cols))))
        
      }
    } else {
      if(sum(colnames(inputData)=="dataset")>0){
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=cluster_anno_res),size=0.3) + geom_point(aes(x=UMAP_1, y=UMAP_2, fill=cluster),color="black",shape=21)+
          coord_equal()+theme_classic()+scale_fill_manual(values=c(hues::iwanthue(length(unique(pd$cluster_anno_res)))))+facet_wrap(~dataset)+scale_color_manual(values=c(hues::iwanthue(length(unique(pd$cluster_anno_res)))))
        
        
      } else {
        color_num=max(length(piechart_data$cluster),length(unique(pd_summary$cluster_anno_res)))
        p=.myDimPlotFn(object=pd_summary, dimCols = c("UMAP_1", "UMAP_2"), attCol = "cluster_anno_res",set_col=F) + geom_point(data=piechart_data,aes(x=UMAP_1, y=UMAP_2, fill=cluster),color="black",shape=21)+ 
          coord_equal()+theme_classic()+scale_fill_manual(values=c(hues::iwanthue(color_num)))+scale_color_manual(values=c(hues::iwanthue(length(unique(pd$cluster_anno_res)))))
        
      }
    }
  }
  
  if(return_obj){
    return(list(pseudocell=piechart_data,cell=pd_summary))
  } else {
    return(p)
  }
  
}


.extra_sconline.CosineDistFn=function(inputMatrix,sig_thr=3){
  
  mymakeCosine=function(inputMatrix,wJaccardDistance=F,sig_only=T,sig_thr=3){
    
    res_mat=NULL
    if(wJaccardDistance){
      ##########################
      res_mat1=res_mat2=inputMatrix
      res_mat1[which(inputMatrix<3)]=0
      res_mat1[which(res_mat1>0)]=1
      
      res_mat2[which(inputMatrix<2)]=0
      res_mat2[which(res_mat2>0)]=1
      
      res_mat=res_mat1 %*% t(res_mat2)
      res_mat=sweep(res_mat,1,pmax(rowSums(res_mat1),1),"/")
      row.names(res_mat)=colnames(res_mat)=row.names(inputMatrix)
      
      res_mat2=t(res_mat)
      res_mat[res_mat<res_mat2]=res_mat2[res_mat<res_mat2]
      rm(res_mat1,res_mat2)
      
      ##########################
    }
    
    if(sig_only){
      inputMatrix=mymakeCosine2(inputMatrix = inputMatrix,sig_thr = sig_thr)
    } else {
      inputMatrix_cosine=inputMatrix*inputMatrix
      inputMatrix_cosine=sqrt(rowSums(inputMatrix_cosine))
      inputMatrix=sweep(inputMatrix,1,inputMatrix_cosine,"/")
      inputMatrix=inputMatrix
      inputMatrix=inputMatrix %*% t(inputMatrix)
    }
    
    
    if(wJaccardDistance){
      inputMatrix=res_mat
    }
    
    return(inputMatrix)
  }
  
  
  res_mat=NULL
  
  inputMatrix_org=inputMatrix
  inputMatrix_cosine=inputMatrix*inputMatrix
  inputMatrix_cosine=sqrt(rowSums(inputMatrix_cosine))
  inputMatrix=sweep(inputMatrix,1,inputMatrix_cosine,"/")
  
  
  inputMatrix2=matrix(0,nrow=nrow(inputMatrix),ncol=nrow(inputMatrix))
  
  for(i in 1:nrow(inputMatrix)){
    #print(i)
    for(j in i:nrow(inputMatrix)){
      sl_ind=apply(inputMatrix_org[c(i,j),],2,max)
      sl_ind=which(sl_ind>=sig_thr)
      if(length(sl_ind)>10){
        sl_ind=mymakeCosine(inputMatrix_org[c(i,j),sl_ind],wJaccardDistance = F,sig_only = F,sig_thr = sig_thr)
        inputMatrix2[i,j]=sl_ind[1,2]
      } else {
        inputMatrix2[i,j]=NA
      }
      
    }
  }
  inputMatrix2=inputMatrix2+t(inputMatrix2)
  diag(inputMatrix2)=1
  
  row.names(inputMatrix2)=colnames(inputMatrix2)=row.names(inputMatrix)
  
  return(inputMatrix2)
}

.extra_sconline.PctScoreFn_slow=function(argList,pct_mat,meta_z_mat,sig1_thr=3,centers=NULL,pct2_thr=0.3,pct_diff_thr=0.2,symmetric=F){
  
  
  if(is.null(centers)){
    centers=row.names(pct_mat)
  }
  
  meta_z_mat=apply(meta_z_mat,2,function(x) as.numeric(x>sig1_thr))
  
  pctDiffCountFn_org=function(i, pct_mat,pct_mat_ref,centers,pct_diff_thr,pct2_thr,meta_z_mat,sig1_thr){
    res_counts=unlist(lapply(1:length(centers),function(j){
      sum((.extra_sconline.PctDiffscoreFn(pct_mat[i,]-pct_mat_ref[jList[j],],pct_diff_thr = pct_diff_thr)*.extra_sconline.Pct2scoreFn(pct_mat_ref[jList[j],],pct2_thr = pct2_thr)>0.99)&meta_z_mat[i,]>sig1_thr)
    }))
    return(res_counts)
  }
  
  
  pctDiffCountFn=function(ilist, pct_mat2,pct_diff_thr,pct2_thr,meta_z_mat2){
    res=matrix(0,nrow=length(ilist),ncol=nrow(pct_mat2))
    counter=0
    for(i in ilist){
      counter=counter+1
      p1=sweep(pct_mat2,2,pct_mat2[i,],"-")*(-1)
      #p1=.extra_sconline.PctDiffscoreFn(p1,pct_diff_thr = pct_diff_thr)
      #p2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2,pct2_thr = pct2_thr)
      p1=unlist(lapply(1:nrow(p1),function(x) {
        x1=.extra_sconline.PctDiffscoreFn(p1[x,],pct_diff_thr = pct_diff_thr)
        x2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2[x,],pct2_thr = pct2_thr)
        sum(as.numeric(x1*x2>0.99)*meta_z_mat2[i,])
      }))
      #p1=do.call("rbind",p1)
      #p1=Matrix::rowSums(p1)
      res[counter,]=p1
    }
    row.names(res)=ilist
    
    return(res)
  }
  
  tmp_groups=split(1:nrow(pct_mat),ceiling(1:nrow(pct_mat)/(nrow(pct_mat)/(argList$ncores*3))))
  
  tstFn=function(tmp_groups){
    res_mat=parallel::mclapply(tmp_groups,pctDiffCountFn,pct_mat2=pct_mat2,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat2=meta_z_mat2,mc.cores = argList$ncores)
    res_mat=do.call("rbind",res_mat)
    return(res_mat)
  }
  
  #pct_mat2=pct_mat;meta_z_mat2=meta_z_mat
  my.env <- new.env()
  my.env$meta_z_mat2 <- meta_z_mat
  my.env$pct_mat2 <- pct_mat
  my.env$pct_diff_thr <- pct_diff_thr
  my.env$pct2_thr <- pct2_thr
  my.env$argList <- argList
  
  with_env <- function(f, e=parent.frame()) {
    stopifnot(is.function(f))
    environment(f) <- e
    f
  }
  
  res_mat=with_env(tstFn,my.env)(tmp_groups)
  
  #res_mat=lapply(1:nrow(prop_mat),pctDiffCountFn_org,centers=centers,pct_mat_ref=pct_mat,pct_mat=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat=meta_z,sig1_thr=sig1_thr)
  #tst=pctDiffCountFn(ilist=tmp_groups[[1]],pct_mat2=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat2=meta_z_mat)
  #res_mat=do.call("rbind",res_mat)
  gc()
  if(symmetric){
    res_mat2=res_mat+t(res_mat)
  } else {
    res_mat2=res_mat
  }
  
  row.names(res_mat2)=row.names(pct_mat)
  colnames(res_mat2)=row.names(pct_mat)
  
  return(res_mat2)
}


.extra_sconline.PctScoreFn=function(argList,pct_mat,meta_z_mat,sig1_thr=3,centers=NULL,pct2_thr=0.3,pct_diff_thr=0.2,symmetric=F,cell_count_diff_thr=0,hierarchical_mode=F){
  
  
  if(is.null(centers)){
    centers=row.names(pct_mat)
  }
  
  meta_z_mat=apply(meta_z_mat,2,function(x) as.numeric(x>sig1_thr))
  
  p2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat,pct2_thr = pct2_thr)
  
  res_mat=matrix(0,nrow=nrow(pct_mat),ncol=length(centers))
  pctDiffCountFn_archive=function(i, pct_mat,pct_mat_ref,centers,pct_diff_thr,pct2_thr,meta_z_mat,sig1_thr){
    res_counts=unlist(lapply(1:length(centers),function(j){
      sum((.extra_sconline.PctDiffscoreFn(pct_mat[i,]-pct_mat_ref[jList[j],],pct_diff_thr = pct_diff_thr)*.extra_sconline.Pct2scoreFn(pct_mat_ref[jList[j],],pct2_thr = pct2_thr)>0.99)&meta_z_mat[i,]>sig1_thr)
    }))
    return(res_counts)
  }
  
  
  pctDiffCountFn_org=function(ilist, pct_mat, pct_mat2,pct_diff_thr,meta_z_mat){
    res=matrix(0,nrow=length(ilist),ncol=ncol(pct_mat))
    counter=0
    
    for(i in ilist){
      counter=counter+1
      #print(counter)
      p1=Matrix::drop0(pct_mat[,i] - pct_mat)
      #p1=sweep(pct_mat,2,pct_mat[i,],"-")*(-1)
      p1=.extra_sconline.PctDiffscoreFn(p1,pct_diff_thr = pct_diff_thr)
      #p2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2,pct2_thr = pct2_thr)
      p1=colSums((p1*pct_mat2>0.99)*meta_z_mat[,i])
      #p1=do.call("rbind",p1)
      #p1=Matrix::rowSums(p1)
      res[counter,]=p1
    }
    row.names(res)=ilist
    
    return(res)
  }
  
  pctDiffCountFn=function(ilist, pct_mat, pct_mat2,pct_diff_thr,meta_z_mat,cell_count_diff_thr,cellCountMat){
    res=matrix(0,nrow=length(ilist),ncol=ncol(pct_mat))
    counter=0
    
    for(i in ilist){
      counter=counter+1
      #print(counter)
      p1=Matrix::drop0(pct_mat[,i] - pct_mat)
      if(cell_count_diff_thr>0){
        p_count=Matrix::drop0(cellCountMat[,i] - cellCountMat,tol = cell_count_diff_thr)
        p_count@x=rep(1,length(p_count@x))
        p1=Matrix::drop0(p1*p_count)
        p1@x[p1@x<0]=0
        p1=Matrix::drop0(p1)
      }
      #p1=sweep(pct_mat,2,pct_mat[i,],"-")*(-1)
      p1=.extra_sconline.PctDiffscoreFn(p1,pct_diff_thr = pct_diff_thr)
      #p2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2,pct2_thr = pct2_thr)
      p1=colSums((p1*pct_mat2>0.99)*meta_z_mat[,i])
      #p1=do.call("rbind",p1)
      #p1=Matrix::rowSums(p1)
      res[counter,]=p1
    }
    row.names(res)=ilist
    
    return(res)
  }
  
  cell_counts=NULL
  if(cell_count_diff_thr>0){
    prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
    matWeights=.myEffSizePropMat(prop_mat)
    matEffectiveSize=matWeights$effective_sample_size
    matEffectiveSize=matEffectiveSize[row.names(pct_mat)]
    cell_counts=.extra_matrix_rowNorm(input_mat = pct_mat,rowValues = matEffectiveSize)#Matrix::Diagonal(x = matEffectiveSize) %*% pct_mat
    
  }
  
  
  tmp_groups=split(1:nrow(pct_mat),ceiling(1:nrow(pct_mat)/(nrow(pct_mat)/(argList$ncores*3))))
  
  
  tstFn=function(tmp_groups){
    res_mat=parallel::mclapply(tmp_groups,pctDiffCountFn,pct_mat=p_mat,pct_mat2=p_mat2,pct_diff_thr=pct_diff_thr,meta_z_mat=meta_z_mat,cellCountMat=cellCountMat,cell_count_diff_thr=cell_count_diff_thr,mc.cores = argList$ncores)
    res_mat=do.call("rbind",res_mat)
    return(res_mat)
  }
  
  if(hierarchical_mode){
    argList$ncores=2
  }
  
  #pct_mat2=pct_mat;meta_z_mat2=meta_z_mat
  my.env <- new.env()
  my.env$meta_z_mat <- t(meta_z_mat)
  my.env$p_mat2 <- t(p2)
  my.env$p_mat <- t(pct_mat)
  my.env$pct_diff_thr <- pct_diff_thr
  my.env$argList <- argList
  if(!is.null(cell_counts)){
    my.env$cellCountMat <- t(cell_counts)
  } else {
    my.env$cellCountMat <- NULL
  }
  
  my.env$cell_count_diff_thr=cell_count_diff_thr
  
  with_env <- function(f, e=parent.frame()) {
    stopifnot(is.function(f))
    environment(f) <- e
    f
  }
  
  res_mat=with_env(tstFn,my.env)(tmp_groups)
  
  #res_mat=lapply(1:nrow(prop_mat),pctDiffCountFn_org,centers=centers,pct_mat_ref=pct_mat,pct_mat=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat=meta_z,sig1_thr=sig1_thr)
  #tst=pctDiffCountFn(ilist=tmp_groups[[1]],pct_mat2=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat2=meta_z_mat)
  #res_mat=do.call("rbind",res_mat)
  gc()
  if(symmetric){
    res_mat2=res_mat+t(res_mat)
  } else {
    res_mat2=res_mat
  }
  
  row.names(res_mat2)=row.names(pct_mat)
  colnames(res_mat2)=row.names(pct_mat)
  
  return(res_mat2)
}

.extra_sconline.PctScoreFn_v2=function(argList,pct_mat,meta_z_mat,sig1_thr=3,centers=NULL,pct2_thr=0.3,pct_diff_thr=0.2,symmetric=F,hierarchical_mode=F){
  
  if(is.null(centers)){
    centers=row.names(pct_mat)
  }
  
  #meta_z_mat=apply(meta_z_mat,2,function(x) as.numeric(x>sig1_thr))
  
  meta_z_mat=meta_z_mat>sig1_thr
  meta_z_mat=t(apply(meta_z_mat,1,as.numeric))
  p2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat,pct2_thr = pct2_thr)
  
  res_mat=matrix(0,nrow=nrow(pct_mat),ncol=length(centers))
  pctDiffCountFn_org=function(i, pct_mat,pct_mat_ref,centers,pct_diff_thr,pct2_thr,meta_z_mat,sig1_thr){
    res_counts=unlist(lapply(1:length(centers),function(j){
      sum((.extra_sconline.PctDiffscoreFn(pct_mat[i,]-pct_mat_ref[jList[j],],pct_diff_thr = pct_diff_thr)*.extra_sconline.Pct2scoreFn(pct_mat_ref[jList[j],],pct2_thr = pct2_thr)>0.99)&meta_z_mat[i,]>sig1_thr)
    }))
    return(res_counts)
  }
  
  
  pctDiffCountFn=function(ilist, pct_mat, pct_mat2,pct_diff_thr,meta_z_mat){
    res=matrix(0,nrow=length(ilist),ncol=ncol(pct_mat))
    counter=0
    
    for(i in ilist){
      counter=counter+1
      #print(counter)
      p1=Matrix::drop0(pct_mat[,i] - pct_mat)
      #p1=sweep(pct_mat,2,pct_mat[i,],"-")*(-1)
      p1=.extra_sconline.PctDiffscoreFn(p1,pct_diff_thr = pct_diff_thr)
      #p2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2,pct2_thr = pct2_thr)
      p1=colSums((p1*pct_mat2>0.99)*meta_z_mat[,i])
      #p1=do.call("rbind",p1)
      #p1=Matrix::rowSums(p1)
      res[counter,]=p1
    }
    row.names(res)=ilist
    
    return(res)
  }
  
  pctDiffCountFn2=function(ilist, argList,hierarchical_mode=F,tmp){
    
    if(is.null(tmp)){
      tmp=qread(.myFilePathMakerFn("tmp_file",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
    }
    
    
    pct_mat=tmp$pct_mat
    pct_mat2=tmp$pct_mat2
    pct_diff_thr=tmp$pct_diff_thr
    meta_z_mat=tmp$meta_z_mat
    rm(tmp)
    
    res=matrix(0,nrow=length(ilist),ncol=ncol(pct_mat))
    counter=0
    
    for(i in ilist){
      counter=counter+1
      #print(counter)
      p1=Matrix::drop0(pct_mat[,i] - pct_mat)
      #p1=sweep(pct_mat,2,pct_mat[i,],"-")*(-1)
      p1=.extra_sconline.PctDiffscoreFn(p1,pct_diff_thr = pct_diff_thr)
      #p2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2,pct2_thr = pct2_thr)
      p1=colSums((p1*pct_mat2>0.99)*meta_z_mat[,i])
      #p1=do.call("rbind",p1)
      #p1=Matrix::rowSums(p1)
      res[counter,]=p1
    }
    row.names(res)=ilist
    
    return(res)
  }
  
  tmp_groups=split(1:nrow(pct_mat),ceiling(1:nrow(pct_mat)/(nrow(pct_mat)/(argList$ncores*3))))
  
  tmp_res=NULL
  if(!hierarchical_mode){
    qsave(list(meta_z_mat=t(meta_z_mat),pct_mat2=t(p2),pct_mat=t(pct_mat),pct_diff_thr=pct_diff_thr),file=.myFilePathMakerFn("tmp_file",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
  } else {
    tmp_res=list(meta_z_mat=t(meta_z_mat),pct_mat2=t(p2),pct_mat=t(pct_mat),pct_diff_thr=pct_diff_thr)
  }
  
  
  if(F){
    #pct_mat2=pct_mat;meta_z_mat2=meta_z_mat
    my.env <- new.env()
    my.env$meta_z_mat <- t(meta_z_mat)
    my.env$pct_mat2 <- t(p2)
    my.env$pct_mat <- t(pct_mat)
    my.env$pct_diff_thr <- pct_diff_thr
    my.env$argList <- argList
    
    with_env <- function(f, e=parent.frame()) {
      stopifnot(is.function(f))
      environment(f) <- e
      f
    }
    
    res_mat=with_env(tstFn,my.env)(tmp_groups)
  } else {
    res_mat=parallel::mclapply(tmp_groups,pctDiffCountFn2,argList=argList,tmp=tmp_res,mc.cores = 2)
    res_mat=do.call("rbind",res_mat)
    
  }
  
  
  #res_mat=lapply(1:nrow(prop_mat),pctDiffCountFn_org,centers=centers,pct_mat_ref=pct_mat,pct_mat=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat=meta_z,sig1_thr=sig1_thr)
  #tst=pctDiffCountFn(ilist=tmp_groups[[1]],pct_mat2=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat2=meta_z_mat)
  #res_mat=do.call("rbind",res_mat)
  gc()
  if(symmetric){
    res_mat2=res_mat+t(res_mat)
  } else {
    res_mat2=res_mat
  }
  
  row.names(res_mat2)=row.names(pct_mat)
  colnames(res_mat2)=row.names(pct_mat)
  
  return(res_mat2)
}

.extra_sconline.PctScoreFn2=function(){
  
  
  if(is.null(centers)){
    centers=row.names(pct_mat)
  }
  
  meta_z_mat=apply(meta_z_mat,2,function(x) as.numeric(x>sig1_thr))
  
  res_mat=matrix(0,nrow=nrow(pct_mat),ncol=length(centers))
  pctDiffCountFn_org=function(i, pct_mat,pct_mat_ref,centers,pct_diff_thr,pct2_thr,meta_z_mat,sig1_thr){
    res_counts=unlist(lapply(1:length(centers),function(j){
      sum((.extra_sconline.PctDiffscoreFn(pct_mat[i,]-pct_mat_ref[jList[j],],pct_diff_thr = pct_diff_thr)*.extra_sconline.Pct2scoreFn(pct_mat_ref[jList[j],],pct2_thr = pct2_thr)>0.99)&meta_z_mat[i,]>sig1_thr)
    }))
    return(res_counts)
  }
  
  
  pctDiffCountFn=function(ilist, pct_mat2,pct_diff_thr,pct2_thr,meta_z_mat2){
    res=matrix(0,nrow=length(ilist),ncol=nrow(pct_mat2))
    counter=0
    for(i in ilist){
      counter=counter+1
      p1=sweep(pct_mat2,2,pct_mat2[i,],"-")*(-1)
      #p1=.extra_sconline.PctDiffscoreFn(p1,pct_diff_thr = pct_diff_thr)
      #p2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2,pct2_thr = pct2_thr)
      p1=unlist(lapply(1:nrow(p1),function(x) {
        x1=.extra_sconline.PctDiffscoreFn(p1[x,],pct_diff_thr = pct_diff_thr)
        x2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2[x,],pct2_thr = pct2_thr)
        sum(as.numeric(x1*x2>0.99)*meta_z_mat2[i,])
      }))
      #p1=do.call("rbind",p1)
      #p1=Matrix::rowSums(p1)
      res[counter,]=p1
    }
    row.names(res)=ilist
    
    return(res)
  }
  
  tmp_groups=split(1:nrow(pct_mat),ceiling(1:nrow(pct_mat)/(nrow(pct_mat)/(argList$ncores*3))))
  
  tstFn=function(tmp_groups){
    res_mat=parallel::mclapply(tmp_groups,pctDiffCountFn,pct_mat2=pct_mat2,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat2=meta_z_mat2,mc.cores = argList$ncores)
    res_mat=do.call("rbind",res_mat)
    return(res_mat)
  }
  
  #pct_mat2=pct_mat;meta_z_mat2=meta_z_mat
  my.env <- new.env()
  my.env$meta_z_mat2 <- meta_z_mat
  my.env$pct_mat2 <- pct_mat
  my.env$pct_diff_thr <- pct_diff_thr
  my.env$pct2_thr <- pct2_thr
  my.env$argList <- argList
  
  with_env <- function(f, e=parent.frame()) {
    stopifnot(is.function(f))
    environment(f) <- e
    f
  }
  
  res_mat=with_env(tstFn,my.env)(tmp_groups)
  
  #res_mat=lapply(1:nrow(prop_mat),pctDiffCountFn_org,centers=centers,pct_mat_ref=pct_mat,pct_mat=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat=meta_z,sig1_thr=sig1_thr)
  #tst=pctDiffCountFn(ilist=tmp_groups[[1]],pct_mat2=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat2=meta_z_mat)
  #res_mat=do.call("rbind",res_mat)
  gc()
  if(symmetric){
    res_mat2=res_mat+t(res_mat)
  } else {
    res_mat2=res_mat
  }
  
  row.names(res_mat2)=row.names(pct_mat)
  colnames(res_mat2)=row.names(pct_mat)
  
  return(res_mat2)
}

.extra_sconline.PctDiffscoreFn=function(pct_diff_value,pct_diff_thr){
  #pct_diff_value=1-pct_diff_value
  #pct_diff_thr=1-pct_diff_thr
  pct_diff_value=pct_diff_thr - pct_diff_value# - pct_diff_thr
  pct_diff_value@x=(-1)*pmax(pct_diff_value@x,0)*20 - pmin(pct_diff_value@x,0)*2
  pct_diff_value@x=1-exp(pct_diff_value@x)
  pct_diff_value=1-as.matrix(pct_diff_value)
  #pct_diff=exp(-pmax(pct_diff_value,0)*20-pmin(pct_diff_value,0)*2)
  return(pct_diff_value)
}

.extra_sconline.Pct2scoreFn=function(pct2_value,pct2_thr){
  pct2=round(exp(-(pmax((pct2_value-pct2_thr),0)/(0.1))^2),3)
  return(pct2)
}

.extra_sconline.FixedSizekmeansFn=function(harmony_embeddings,nPCs,pseudocell_count){
  
  myPseudoAffinityMakerFn=function(harmony_embeddings,k.param=20,prune.SNN=1/15,n.trees = 50){
    #nn.ranked.1 <- RANN::nn2(harmony_embeddings, k = 10, eps = 0)
    
    idx=Seurat:::AnnoyBuildIndex(data = harmony_embeddings, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = harmony_embeddings,k=k.param,include.distance = T,search.k = -1)
    
    
    nn.ranked=nn.ranked.1$nn.idx
    graph= Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(graph) <- rownames(x = harmony_embeddings)
    colnames(graph) <- rownames(x = harmony_embeddings)
    
    
    #graph=Seurat::FindNeighbors(harmony_embeddings,compute.SNN=T)
    #graph=as(graph[["snn"]], "dgCMatrix")
    
    #j <- as.numeric(t(nn.ranked.1$nn.idx))
    #i <- ((1:length(j)) - 1)%/%ncol(nn.ranked.1$nn.idx) + 1
    #k=1#as.numeric(t(nn.ranked.1$nn.dists))
    #graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(harmony_embeddings), nrow(harmony_embeddings)))
    
    #graph=graph+t(graph)
    #graph@x=rep(1,length(graph@x))
    
    if(F){
      graph=graph %*% t(graph)
      diag(graph)=0
      
      
      listCols_sparse<-function(X){
        #converts a sparse Matrix into a list of its columns
        #each list item contains only the nonzero elements of the column
        X<-as(X,"CsparseMatrix")
        res<-split(X@x, findInterval(seq_len(nnzero(X)), X@p, left.open=TRUE))
        names(res)<-colnames(X)
        res
      }
      
      colapply_sparse_nonzero<-function(X,FUN,...,mc.cores=1){
        #apply a function FUN to NONZERO elements of each column of sparse Matrix X
        #for an alternative that operates on all values, see colapply_full
        #mc: should parallel processing be used? Only recommended if FUN is slow
        #... additional args passed to mclapply or to FUN
        #this function always returns a list of length ncol(X)
        if(mc.cores>1){
          res=mclapply(listCols_sparse(X),FUN,...,mc.cores=mc.cores)
        } else {
          res=lapply(listCols_sparse(X),FUN,...)
        }
        res=unlist(res)
        X@x=res
        return(X)
      }
      
      affinities=colapply_sparse_nonzero(X=t(graph),FUN=function(x) exp((-3)*((max(x)-x)/(max(x)+1))^2),mc.cores=argList$ncores)
      affinities=sqrt(t(affinities)*affinities)
      diag(affinities)=0
    } else {
      affinities=graph
    }
    
    return(affinities)
  }
  
  set.seed(1)
  doClustering=T
  itrClustering=0
  while(doClustering&itrClustering<10){
    itrClustering=itrClustering+1
    res_clust=kmeans(harmony_embeddings[,1:nPCs,drop=F],pseudocell_count,iter.max = 1000) #,algorithm = "Lloyd")
    
    
    
    if(sum(is.na(res_clust$centers))==0){
      doClustering=F
    }
  }
  
  if(sum(is.na(res_clust$centers))>0){
    stop("Error in identification of the pseudocells")
  }
  
  res_clusters=data.frame(sample=names(res_clust$cluster),cluster_id=res_clust$cluster,stringsAsFactors = F)
  
  pseudo_affinity=myPseudoAffinityMakerFn(harmony_embeddings = harmony_embeddings)
  row.names(pseudo_affinity)=colnames(pseudo_affinity)=row.names(harmony_embeddings)
  
  sl_pseudo=NULL
  for(x in unique(res_clusters$cluster_id)){
    if(sum(res_clusters$cluster_id==x)==1){
      x2=data.frame(cluster=x,pseudocell=res_clusters$sample[res_clusters$cluster_id==x],stringsAsFactors = F)
    } else {
      tmp=res_clusters$sample[res_clusters$cluster_id==x]
      tmp=rowSums(as.matrix(pseudo_affinity[tmp,tmp]))
      tmp=tmp[order(tmp,decreasing = T)]
      tmp=names(tmp)[1]
      x2=data.frame(cluster=x,pseudocell=tmp,stringsAsFactors = F)
    }
    sl_pseudo=rbind(sl_pseudo,x2)
  }
  
  #pca_centroid=res_clust$centers
  pca_centroid=harmony_embeddings[sl_pseudo$pseudocell,,drop=F]
  row.names(pca_centroid)=sl_pseudo$cluster
  pca_centroid=pca_centroid[row.names(res_clust$centers)[row.names(res_clust$centers) %in% row.names(pca_centroid)],]
  
  output=sl_pseudo$pseudocell
  return(output)
}

.extra_sconline.FixedSizeFn=function(inputExp,inputEmbedding=NULL,pseudocell_size=40,n.neighbors=20,n.trees=50,nPCs=30,k.param = 20,include.outlier.score=F,rand_pseudobulk_mod=T,cols_to_sum=NULL){
  #inputExp=datalist[[i]];inputEmbedding=NULL;pseudocell_size=50;nPCs=30
  #n.neighbors=20;n.trees=50;k.param = 20
  
  if(ncol(inputExp)==1){
    inputExp$pseudocell_size=1
    return(inputExp)
  }
  
  pseudocell_count=round(ncol(inputExp)/pseudocell_size)
  if(pseudocell_count<=1){
    expData=matrix(as.numeric(rowSums(counts(inputExp))),ncol=1)
    row.names(expData)=row.names(inputExp)
    colnames(expData)=colnames(inputExp)[1]
    
    pd=as.data.frame(colData(inputExp))
    pd=pd[1,]
    pd$pseudocell_outlier_score=NA
    pd$pseudocell_size=ncol(inputExp)
    
    if(!is.null(cols_to_sum)){
      for(icoltosum in cols_to_sum){
        pd[1,icoltosum]=sum(as.data.frame(colData(inputExp))[,icoltosum])
      }
    }
    
    fd=as.data.frame(rowData(inputExp))
    
    res=SingleCellExperiment(assays = list(counts = expData),colData = pd,rowData=fd)
    return(res)
  }
  
  if(rand_pseudobulk_mod){
    prop_m_hardCluster=sample(colnames(inputExp))
    prop_m_hardCluster=split(prop_m_hardCluster,floor(ecdf(seq_along(prop_m_hardCluster))(seq_along(prop_m_hardCluster))*pseudocell_count-0.001))
    prop_m_hardCluster=lapply(prop_m_hardCluster,function(x){
      x=data.frame(pseudocell=x[1],cells=x,stringsAsFactors = F)
      return(x)
    })
    prop_m_hardCluster=do.call("rbind",prop_m_hardCluster)
    
    prop_m_hardCluster$name=prop_m_hardCluster[,1]
    prop_m_hardCluster$name=gsub(" ",".",as.character(prop_m_hardCluster$name))
    
    prop_m_hardCluster$name=gsub("[[:punct:]]+", ".", prop_m_hardCluster$name)
    
    m2=.myOneHotFn(inputVector=prop_m_hardCluster$name)
    row.names(m2)=prop_m_hardCluster[,2]
    m2_colname=prop_m_hardCluster[!duplicated(prop_m_hardCluster[,1]),]
    m2_colname=m2_colname[match(colnames(m2),m2_colname$name),]
    colnames(m2)=m2_colname[,1]
    prop_m_hardCluster=t(as.matrix(m2))
    rm(m2)
  } else {
    if(is.null(inputEmbedding)){
      tmpData=suppressWarnings(.extraExport2SeuratFn(inputExp))
      tmpData = NormalizeData(tmpData,verbose =F)
      tmpData = FindVariableFeatures(tmpData, selection.method = "vst", nfeatures = 2000,verbose =F)
      tmpData <- ScaleData(tmpData,verbose =F)
      tmpData <- RunPCA(tmpData,npcs =nPCs,verbose =F)
      inputEmbedding=tmpData@reductions$pca@cell.embeddings[,1:nPCs]
    } else {
      inputEmbedding=inputEmbedding[colnames(inputExp),1:nPCs]
    }
    
    harmony_embeddings=inputEmbedding
    
    
    pseudocell_names=.extra_sconline.FixedSizekmeansFn(harmony_embeddings=harmony_embeddings,nPCs = nPCs,pseudocell_count = pseudocell_count)#,kmeansMethod=kmeans_method)
    
    
    #pca_centroid=res_clust$centers
    
    pca_centroid=harmony_embeddings[pseudocell_names,,drop=F]
    row.names(pca_centroid)=paste0("ps_",1:nrow(pca_centroid))
    sl_pseudo=data.frame(cluster=paste0("ps_",1:nrow(pca_centroid)),pseudocell=pseudocell_names,stringsAsFactors = F)
    
    
    idx=Seurat:::AnnoyBuildIndex(data = harmony_embeddings, metric = "euclidean", n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = harmony_embeddings,k=k.param,include.distance = T,search.k = -1)
    
    affinities=.extra_matrix_rowNorm(input_mat = nn.ranked.1$nn.dists,rowValues = 1/(nn.ranked.1$nn.dists[,2]+0.000001))#Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,2]+0.000001)) %*% nn.ranked.1$nn.dists
    affinities@x=-1*affinities@x^2
    affinities@x=exp(affinities@x)
    affinities[,1]=affinities[,2]
    j <- as.numeric(t(nn.ranked.1$nn.idx))
    i <- ((1:length(j)) - 1)%/%n.neighbors + 1
    x=as.numeric(t(affinities))
    adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(inputEmbedding),nrow(inputEmbedding)))
    rownames(adj) <- row.names(inputEmbedding)
    colnames(adj)=c(row.names(inputEmbedding))
    adj= .extra_matrix_rowNorm(adj)#Matrix::Diagonal(x=1/rowSums(adj)) %*% adj
    
    
    adj_t=t(adj)
    adj_t=.extra_matrix_rowNorm(adj_t)#Matrix::Diagonal(x=1/rowSums(adj_t)) %*% adj_t
    adj=adj[sl_pseudo$pseudocell,]
    #row.names(adj)=sl_pseudo$cluster
    adj=adj %*% adj_t
    
    colMax_vals_m=qlcMatrix::colMax(adj)
    colMax_vals_m=.extra_matrix_colNorm(input_mat = adj,colValues = 1/as.numeric(colMax_vals_m))#adj %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    itol=0.5
    tst_mat=Matrix::drop0(colMax_vals_m,tol=itol)
    tst_mat@x=rep(1,length(tst_mat@x))
    while(all(rowSums(tst_mat)>=pseudocell_size)&itol<0.95){
      itol=itol+0.05
      tst_mat=Matrix::drop0(colMax_vals_m,tol=itol)
      tst_mat@x=rep(1,length(tst_mat@x))
    }
    prop_m_hardCluster=Matrix::drop0(colMax_vals_m,tol=itol)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    
    o=table(prop_m_hardCluster$i)
    o=o[order(as.numeric(o),decreasing = F)]
    
    prop_m_hardCluster$groups=factor(as.character(prop_m_hardCluster$i),levels=names(o))
    prop_m_hardCluster=prop_m_hardCluster[order(prop_m_hardCluster$groups,decreasing = F),]
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
    
    prop_m_hardCluster=sparseMatrix(i = prop_m_hardCluster$i, j = prop_m_hardCluster$j, x = rep(1,nrow(prop_m_hardCluster)), dims = c(nrow(adj), ncol(adj)))
    row.names(prop_m_hardCluster)=row.names(adj)
    colnames(prop_m_hardCluster)=row.names(inputEmbedding)
    prop_m_hardCluster=prop_m_hardCluster[which(rowSums(prop_m_hardCluster)>0),]
    
  }
  expData=t(counts(inputExp))
  expData=t(prop_m_hardCluster %*% expData[colnames(prop_m_hardCluster),])
  
  #inputData
  #sl_pseudo=sl_pseudo[match(row.names(prop_m_hardCluster),sl_pseudo$pseudocell),]
  
  if(is.null(cols_to_sum)&sum(colnames(colData(inputExp)) %in% cols_to_sum)==0){
    pd=as.data.frame(colData(inputExp))
    pd=pd[match(row.names(prop_m_hardCluster),colnames(inputExp)),]
    pd$pseudocell_size=rowSums(prop_m_hardCluster)
    
    
  } else {
    if(length(setdiff(cols_to_sum,colnames(colData(inputExp))))>0){
      warning("some of provided cols_to_sum cols were not identified in the dataset")
    }
    
    pd=as.data.frame(colData(inputExp))[,!colnames(colData(inputExp)) %in% cols_to_sum]
    pd=pd[match(row.names(prop_m_hardCluster),colnames(inputExp)),]
    pd$pseudocell_size=rowSums(prop_m_hardCluster)
    
    
    sums_res=prop_m_hardCluster %*% as.matrix(as.data.frame(colData(inputExp)[,cols_to_sum]))
    if(any(row.names(sums_res)!=row.names(pd),na.rm = F)){
      print("Error in summing the cols!")
    }
    pd=cbind(pd,sums_res)
  }
  
  
  
  fd=as.data.frame(rowData(inputExp))
  
  pd$pseudocell_outlier_score=NA
  
  if(include.outlier.score&ncol(expData)>3){
      bkg_genes=counts(inputExp)
      bkg_genes=rowSums(bkg_genes>0)/max(ncol(bkg_genes),10)
      if(sum(bkg_genes>0.1)>100){
        bkg_genes=row.names(expData)[bkg_genes>0.1]
        
        logCPM=edgeR::cpm(expData[bkg_genes,],normalized.lib.sizes = F, log = TRUE, prior.count = 1)
        
        tocheck=require(WGCNA,quietly = T)
        if(!tocheck){
          stop("WGCNA package is missing!")
        }
        
        normadj <- (0.5+0.5*bicor(logCPM, use='pairwise.complete.obs'))^2
        netsummary <- fundamentalNetworkConcepts(normadj)
        outlier_score=scale(netsummary$Connectivity)
        pd$pseudocell_outlier_score=as.numeric(outlier_score)
      }
  }
  
  res=SingleCellExperiment(assays = list(counts = expData),colData = pd,rowData=fd)
  
  return(res)
}

.extra_sconline.duplicateCorrelation=function (object, design = NULL, ndups = 2, spacing = 1, block = NULL, 
                                               trim = 0.15, weights = NULL) {
  require(limma)
  y <- limma:::getEAWP(object)
  M <- y$exprs
  ngenes <- nrow(M)
  narrays <- ncol(M)
  if (is.null(design)) 
    design <- y$design
  if (is.null(design)) 
    design <- matrix(1, ncol(y$exprs), 1)
  else {
    design <- as.matrix(design)
    if (mode(design) != "numeric") 
      stop("design must be a numeric matrix")
  }
  if (nrow(design) != narrays) 
    stop("Number of rows of design matrix does not match number of arrays")
  ne <- limma:::nonEstimable(design)
  if (!is.null(ne)) 
    cat("Coefficients not estimable:", paste(ne, collapse = " "), 
        "\n")
  nbeta <- ncol(design)
  if (missing(ndups) && !is.null(y$printer$ndups)) 
    ndups <- y$printer$ndups
  if (missing(spacing) && !is.null(y$printer$spacing)) 
    spacing <- y$printer$spacing
  if (missing(weights) && !is.null(y$weights)) 
    weights <- y$weights
  if (!is.null(weights)) {
    weights <- asMatrixWeights(weights, dim(M))
    weights[weights <= 0] <- NA
    M[!is.finite(weights)] <- NA
  }
  if (is.null(block)) {
    if (ndups < 2) {
      warning("No duplicates: correlation between duplicates not estimable")
      return(list(cor = NA, cor.genes = rep(NA, nrow(M))))
    }
    if (is.character(spacing)) {
      if (spacing == "columns") 
        spacing <- 1
      if (spacing == "rows") 
        spacing <- object$printer$nspot.c
      if (spacing == "topbottom") 
        spacing <- nrow(M)/2
    }
    Array <- rep(1:narrays, rep(ndups, narrays))
  }
  else {
    ndups <- 1
    nspacing <- 1
    Array <- block
  }
  if (is.null(block)) {
    M <- limma:::unwrapdups(M, ndups = ndups, spacing = spacing)
    ngenes <- nrow(M)
    if (!is.null(weights)) 
      weights <- limma:::unwrapdups(weights, ndups = ndups, spacing = spacing)
    design <- design %x% rep(1, ndups)
  }
  if (!requireNamespace("statmod", quietly = TRUE)) 
    stop("statmod package required but is not installed")
  rho <- rep(NA, ngenes)
  nafun <- function(e) NA
  for (i in 1:ngenes) {
    y <- drop(M[i, ])
    o <- is.finite(y)
    A <- factor(Array[o])
    nobs <- sum(o)
    nblocks <- length(levels(A))
    if (nobs > (nbeta + 2) && nblocks > 1 && nblocks < nobs - 1) {
      y <- y[o]
      X <- design[o, , drop = FALSE]
      Z <- model.matrix(~0 + A)
      if (!is.null(weights)) {
        w <- drop(weights[i, ])[o]
        s <- tryCatch(statmod::mixedModel2Fit(y, X, Z, 
                                              w, only.varcomp = TRUE, maxit = 20)$varcomp, error = nafun)
      } else{
        s <- tryCatch(statmod::mixedModel2Fit(y, X, 
                                              Z, only.varcomp = TRUE, maxit = 20)$varcomp, 
                      error = nafun)
      } 
      if (!is.na(s[1])) 
        rho[i] <- s[2]/sum(s)
    }
  }
  arho <- atanh(pmax(-1, rho))
  mrho <- tanh(mean(arho, trim = trim, na.rm = TRUE))
  list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho,value.list=arho)
}


######################
#main sconline functions

#Performs UMAP on the PC space and maps pseudocells to it
.sconline.umapFn_org=function(argList,umap.method='umap-learn'){
  reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F));F}, error=function(e) {return(T)})
  
  if(!reRunCheck){
    reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_centroid",argList=argList));F}, error=function(e) {return(T)})
  }
  
  if(reRunCheck|argList$newRun){
    cat("Running UMAP\n")
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    
    load(.myFilePathMakerFn("pca_centroids",argList=argList))
    
    pca_centroid=pca_centroid[,1:argList$nPCs]
    
    
    
    tst=.reductionUMAPFn(harmony_embeddings,umap.method=umap.method,testPCAembeddings=pca_centroid)
    resUMAP=tst$embedding
    .mySaveFn(resUMAP,file=.myFilePathMakerFn("UMAP_res",argList=argList))
    
    x=as.data.frame(resUMAP)
    x=cbind(x,pd)
    row.names(x)=row.names(pd)
    pd=x
    .mySaveFn(pd,file=.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
    tst=tst[["test"]]
    UMAP_centroid=data.frame(centroid=row.names(pca_centroid),UMAP_1=tst[,1],UMAP_2=tst[,2],stringsAsFactors = F)
    .mySaveFn(UMAP_centroid,file=.myFilePathMakerFn("UMAP_centroid",argList=argList))
    
    
  }
  return("Done!")
}

#argList=.ArgList;inputExpData=data;inputGenes="DDR3"
.sconline.markerPlot_archive=function(argList,inputExpData,inputGenes){
  library(dplyr)
  library(purrr)
  library(cowplot)
  library(patchwork)
  
  load(.myFilePathMakerFn("consensusDE_markers",argList=argList,uniformImportant=T))
  load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
  
  cat("           Generating marker figs (cell based) ...\n")
  load(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  tmp=NULL
  if(is.null(inputExpData)&file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),argList)))){
    load(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),argList)))
  } else if(!is.null(inputExpData)){
    tmp=inputExpData
  }
  
  
  if(!is.null(tmp)){
    tmp=tmp[,colnames(tmp) %in% row.names(pd)]
    pd2=pd[match(colnames(tmp),row.names(pd)),]
  }
  
  
  
  plotFn_integrated=function(index,geneNameList,expData,annotationFile,input_cell_bkg_umap,prefix){
    
    genes=unlist(geneNameList[[index]]$gene_short_name)
    
    figWidth=49
    if(length(genes)<21){
      figWidth=49/3*ceiling(length(genes)/7)
    }
    
    figHieght=35
    if(length(genes)<7){
      figHieght=35/7*length(genes)
    }
    
    if(!is.null(expData)){
      expData=expData[match(c(geneNameList[[index]]$gene,setdiff(row.names(expData),geneNameList[[index]]$gene)),row.names(expData)),]
    }
    
    p=list()
    for(i in seq(1,length(genes),7)){
      p1=NULL
      if(!is.null(expData)){
        p1=.myFeaturePlot(inputSeurat = expData,inputDimData = annotationFile,inputGenes = genes[i:min((i+6),length(genes))],combine_figs=F)
      }
      p2=.my2dPlot_counts(inputPCA=annotationFile,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=genes[i:min((i+6),length(genes))],geneNameCol="gene_short_name",expData=expData,combine_figs=F)
      
      tmp_UMAP=list()
      for(ij in i:min((i+6),length(genes))){
        tmp=UMAP_centroid
        tmp$gene_id=geneNameList[[index]]$gene[which(geneNameList[[index]]$gene_short_name==genes[ij])]
        tmp$gene=genes[ij]
        tmp_res=res_arranged[which(res_arranged$gene==unique(tmp$gene_id)),]
        tmp_res=tmp_res[match(as.character(tmp$centroid),as.character(tmp_res$centroid)),]
        tmp$zscore=tmp_res$zscore
        tmp$effectiveSize=tmp_res$effective_size
        tmp$size=tmp_res$count
        if(sum(is.na(tmp$zscore))>0){
          tmp$zscore[is.na(tmp$zscore)]=0
          tmp$effectiveSize[is.na(tmp$effectiveSize)]=0
          tmp$size[is.na(tmp$size)]=0
        }
        tmp$effectiveSize=tmp$effectiveSize+1
        tmp_UMAP=c(tmp_UMAP,list(tmp))
      }
      tmp_UMAP=do.call("rbind",tmp_UMAP)
      
      tmp_UMAP$gene=factor(as.character(tmp_UMAP$gene),levels=genes[i:min((i+6),length(genes))])
      
      tmp_UMAP$effectiveSize=as.numeric(as.character(tmp_UMAP$effectiveSize))
      tmp_UMAP$effectiveSize[is.na(tmp_UMAP$effectiveSize)]=min(min(tmp_UMAP$effectiveSize,na.rm = T),1)
      if(sum(tmp_UMAP$zscore>5)>0){
        tmp_UMAP$zscore[which(tmp_UMAP$zscore>5)]=5
      }
      if(sum(tmp_UMAP$zscore<(-5))>0){
        tmp_UMAP$zscore[which(tmp_UMAP$zscore<(-5))]=(-5)
      }
      p3=tmp_UMAP %>% 
        group_split(gene) %>% 
        map(
          ~ggplot(.,aes(UMAP_1,UMAP_2,color=zscore,size=effectiveSize))+geom_point(data=input_cell_bkg_umap,aes(UMAP_1,UMAP_2),color="lightgrey",size=0.1)+geom_point()+theme_cowplot()+ theme(plot.title = element_text(hjust = 0.5))+ggtitle(toupper(unique(.$gene)))+scale_color_gradientn(colors = c("lightblue","lightblue","black","yellow","orange","red"),breaks=c(-5,-2,0,2,3,5),limits=c(-5,5))+scale_size_continuous(limits=c(0,max(tmp_UMAP$effectiveSize)),breaks=seq(0,max(tmp_UMAP$effectiveSize),5))
        ) #%>% 
      #plot_grid(plotlist = ., align = 'hv',ncol=1)
      
      if(length(genes)<(i+6)&length(genes)>7){
        dfNull=data.frame()
        pNull=ggplot(dfNull)+geom_point()
        for(inull in 1:(i+6-length(genes))){
          p1=c(p1,list(pNull))
          p2=c(p2,list(pNull))
          p3=c(p3,list(pNull))
        }
      }
      
      p=c(p,p1,p3,p2)
      
    }
    p=wrap_plots(p,nrow = min(length(genes),7),byrow = F)
    if(is.null(prefix)){
      ggsave(plot=p,file=.myFilePathMakerFn(filename = paste0("marker_",names(geneNameList)[index]),argList = argList,pdf = T,uniformImportant=T,makePlotDir=T),width = figWidth,height = figHieght)
    } else {
      ggsave(plot=p,file=.myFilePathMakerFn(filename = paste0(prefix,"_marker_",names(geneNameList)[index]),argList = argList,pdf = T,uniformImportant=T,makePlotDir=T),width = figWidth,height = figHieght)
    }
    
    return("Done")
  }
  
  if(!is.null(inputGenes)){
    geneEffectSizes=geneEffectSizes[which(toupper(geneEffectSizes$gene) %in% toupper(inputGenes)),]
  }
  slGenes=geneEffectSizes[order(geneEffectSizes$effectSize,decreasing = T),]
  id=rep(0,nrow(slGenes))
  for(ic in unique(slGenes$cluster)){
    tmpId=which(slGenes$cluster==ic)
    tmpId=split(tmpId,ceiling(seq_along(tmpId)/21))
    counter=0
    for(ij in 1:length(tmpId)){
      counter=counter+1
      slGenes$cluster[tmpId[[ij]]]=paste0("cluster",slGenes$cluster[tmpId[[ij]]],"_",counter)
    }
  }
  slGenes=split(slGenes, slGenes$cluster)
  
  input_cell_bkg_umap=.scatterPlot_summary2d(object=pd2,reductionCols=c("UMAP_1","UMAP_2"))
  
  plotMI=parallel::mclapply(1:length(slGenes),plotFn_integrated,geneNameList=slGenes,expData=tmp[unique(c(1:100,which(row.names(tmp) %in% geneEffectSizes$gene))),],annotationFile=pd2,input_cell_bkg_umap=input_cell_bkg_umap,prefix=prefix,mc.cores = argList$ncores)
  save(plotMI,file="final_errors.rda")
  
  
  return("Done!")
}


.extra_sconline.Fit_LimmaTrendFn_archived=function(sl_data,covariates,blocking_var="anno_batch"){
  require(edgeR)
  require(limma)
  
  #Note: don't foget to scale numeric variables when needed. one example of such covariate is QC_Gene_unique_count
  #sl_data$anno_age=as.numeric(sl_data$anno_age) it can be log scale if it makes more sense!
  #sl_data$anno_age_scale=scale(sl_data$anno_age)
  
  cov_str="~0"
  for(icov in covariates){
    cov_str=paste0(cov_str,"+",icov)
  }
  model = model.matrix(as.formula(cov_str),data=colData(sl_data))
  dge <- DGEList(counts=counts(sl_data))
  keep <- filterByExpr(counts(sl_data),min.count = 0.9, 
                       min.total.count = max(0.02*ncol(sl_data),10))
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  #dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=T, prior.count=3)
  
  
  blocked_analysis=F
  if(!is.null(blocking_var)){
    sl_data$anno_batch=colData(sl_data)[,blocking_var]
    dc <- .extra_sconline.duplicateCorrelation(logCPM,design=model, block=sl_data[,blocking_var])
    if(!is.nan(dc$consensus.correlation)){
      if(abs(dc$consensus.correlation)<0.9){
        fit <- lmFit(logCPM, model,block = sl_data[,blocking_var], correlation=dc$consensus.correlation)
        blocked_analysis=T
      } else {
        fit <- lmFit(logCPM, model)
      }
    } else {
      fit <- lmFit(logCPM, model)
    }
    
  } else {
    fit <- lmFit(logCPM, model)
  }
  
  
  return(list(fit=fit,dc=dc,model=model,blocked_analysis=blocked_analysis))
}

.extra_sconline.BinarizePropMatrix=function(inputPorpMat=NULL,argList=NULL,tol=0.95){
  #either inputPorpMat or argList should be provided
  #if inputPorpMat is provided, argList argument will be ignored
  
  require(Matrix)
  
  if(is.null(inputPorpMat)){
    inputPorpMat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    inputPorpMat=.extra_matrix_rowNorm(inputPorpMat)#Matrix::Diagonal(x = 1 / (rowSums(inputPorpMat)+0.000000000001)) %*% inputPorpMat
  }
  
  if(is.null(row.names(inputPorpMat))){
    row.names(inputPorpMat)=1:nrow(inputPorpMat)
  }
  if(is.null(colnames(inputPorpMat))){
    colnames(inputPorpMat)=1:ncol(inputPorpMat)
  }
  
  colMax_vals_m=as.numeric(qlcMatrix::colMax(inputPorpMat))
  colMax_vals_m[which(colMax_vals_m==0)]=1
  colMax_vals_m=.extra_matrix_colNorm(input_mat = inputPorpMat,colValues = 1/as.numeric(colMax_vals_m))#inputPorpMat %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=Matrix::drop0(colMax_vals_m,tol=tol)
  prop_m_hardCluster@x=rep(1,length(prop_m_hardCluster@x))
  row.names(prop_m_hardCluster)=row.names(inputPorpMat)
  colnames(prop_m_hardCluster)=colnames(inputPorpMat)
  
  binary_df=as.data.frame(summary(prop_m_hardCluster))
  binary_df$i=row.names(prop_m_hardCluster)[binary_df$i]
  binary_df$j=colnames(prop_m_hardCluster)[binary_df$j]
  binary_df=binary_df[,c("i","j")]
  
  return(list(binary_mat=prop_m_hardCluster,binary_df=binary_df))
}


.extra_sconline.NormVST=function(inputCountData,fitType="parametric"){
  #fitType: used for the estimatation of dispersions. acceptable values: parametric, local, and mean
  require(DESeq2)
  
  keep <- rowSums(counts(inputCountData)>0) >= ncol(inputCountData)*0.1
  
  tmp=DESeq2::varianceStabilizingTransformation(as.matrix(counts(inputCountData)), blind = TRUE, fitType = fitType)
  
  counts(inputCountData)=.matrixExtraction(tmp)
  
  if(length(keep)==nrow(inputCountData)){
    inputCountData2=inputCountData[keep,]
  } else {
    print("error")
  }
  
  return(list(originalData=inputCountData,filteredData=inputCountData2))
}


#random pseudocell creation

.myRNAseqNormVSTfn=function(inputCountData,fitType="parametric"){
  #fitType: used for the estimatation of dispersions. acceptable values: parametric, local, and mean
  require(DESeq2)
  
  keep <- rowSums(counts(inputCountData)>0) >= ncol(inputCountData)*0.1
  
  tmp=DESeq2::varianceStabilizingTransformation(as.matrix(counts(inputCountData)), blind = TRUE, fitType = fitType)
  
  counts(inputCountData)=.matrixExtraction(tmp)
  
  if(length(keep)==nrow(inputCountData)){
    inputCountData2=inputCountData[keep,]
  } else {
    print("error")
  }
  
  return(list(originalData=inputCountData,filteredData=inputCountData2))
}


.extra_sconline.Fit_LimmaTrendFn=function(sl_data,model,random_effect=NULL,TMMnorm=F,VSTnorm=F,prior.count=1,quantile.norm=F,bkg_genes=NULL,no_normalization=F,dc.object=NULL,include.malat1.as.covariate=F){
  require(edgeR)
  require(limma)
  library(Matrix)
  
  if(VSTnorm){
    logCPM=.myRNAseqNormVSTfn(inputCountData=sl_data,fitType="parametric")
    logCPM=logCPM$originalData
    logCPM=counts(logCPM)
    if(!is.null(bkg_genes)){
      keep=row.names(sl_data) %in% bkg_genes
    } else {
      keep <- rowSums(logCPM>3)>(0.05*ncol(logCPM))
    }
    
    logCPM <- logCPM[keep,]
  } else if(no_normalization){
    logCPM=counts(sl_data)
    if(quantile.norm){
      logCPM=limma::normalizeQuantiles(logCPM)
    }
  } else {
    if(!is.null(bkg_genes)){
      keep=row.names(sl_data) %in% bkg_genes
    } else {
      tmpCount2=apply(counts(sl_data),1,function(x) sum(x>0))
      tmpCount=rowSums(edgeR::cpm(as.matrix(counts(sl_data))))
      keep=tmpCount>max(0.01*ncol(sl_data),min(15,ncol(sl_data)/3))
      keep=keep & tmpCount2>max(0.01*ncol(sl_data),min(10,ncol(sl_data)/3))
    }
    
    #tmpCount=apply(counts(sl_data),1,function(x) sum(x>0))
    #keep=rowData(sl_data)$exp.pct>0.03
    #keep=filterByExpr(counts(sl_data),min.count = 0.9, 
    #                  min.total.count = max(0.02*ncol(sl_data),10))
    print(paste("Number of expressed genes:",sum(keep)))
    tmpexp=counts(sl_data)[keep,]
    dge <- DGEList(tmpexp)
    
    if(TMMnorm){
      dge <- calcNormFactors(dge)
    }
    
    logCPM <- new("EList")
    logCPM$E <- edgeR::cpm(dge,normalized.lib.sizes = TMMnorm, log = TRUE, prior.count = prior.count)
    if(quantile.norm){
      logCPM$E=limma::normalizeQuantiles(logCPM$E)
    }
  }
  
  
  if(include.malat1.as.covariate){
    malat_gene=row.names(logCPM)[grepl("malat1",tolower(row.names(logCPM)))]
    if(class(logCPM)[1]=="EList"){
      model=cbind(model,malat=as.numeric(logCPM$E[malat_gene,]))
    } else {
      model=cbind(model,malat=as.numeric(logCPM[malat_gene,]))
    }
    
  }
  
  dc=NULL
  if(!is.null(random_effect)){
    if(is.null(dc.object)){
      dc <- .extra_sconline.duplicateCorrelation(logCPM,design=model, block=colData(sl_data)[,random_effect])
    } else {
      dc=dc.object
    }
    
  }
  
  #
  blocked_analysis=F
  if(!is.null(dc)){
    if(!is.nan(dc$consensus.correlation)){
      if(abs(dc$consensus.correlation)<0.9){
        fit <- lmFit(logCPM, model,block = as.character(colData(sl_data)[,random_effect]), correlation=dc$consensus.correlation)
        blocked_analysis=T
      } else {
        fit <- lmFit(logCPM, model)
      }
    } else {
      fit <- lmFit(logCPM, model)
    }
  } else {
    fit <- lmFit(logCPM, model)
  }
  
  
  return(list(fit=fit,dc=dc,model=model,normData=logCPM,blocked_analysis=blocked_analysis))
}


.extra_sconline.Fit_LimmaDreamFn=function(sl_data,pd,model,random_effect=NULL,TMMnorm=F,VSTnorm=F,prior.count=1,quantile.norm=F,bkg_genes=NULL,no_normalization=F,dc.object=NULL,include.malat1.as.covariate=F,ncores=4){
  require(edgeR)
  require(limma)
  require(variancePartition)
  
  param = SnowParam(ncores, "SOCK", progressbar=TRUE)
  
  # estimate weights using linear mixed model of dream
  
  sl_data=sl_data[rowSds(as.matrix(counts(sl_data))) > 0, ]
  
  dge <- DGEList(counts=counts(sl_data))
  if(is.null(bkg_genes)){
    tmpCount=apply(counts(sl_data),1,function(x) sum(x>0))
    keep=tmpCount>10
    
  } else {
    keep=row.names(sl_data) %in% bkg_genes
  }
  
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  if(TMMnorm){
    dge <- calcNormFactors(dge)
  }
  
 
  
  # estimate weights using linear mixed model of dream
  #model=as.formula("~status + anno_sex +nUMI_scaled  + pseudocell_size_scale +(1 | anno_batch)")
  
  if(quantile.norm){
    vobjDream = voomWithDreamWeights( dge, model,as.data.frame(pd), BPPARAM=param, normalize.method="quantile" )
  } else {
    vobjDream = voomWithDreamWeights( dge, model,as.data.frame(pd), BPPARAM=param )
  }
  
  
  
  
  fit = dream( vobjDream, model, as.data.frame(pd) )
  
  
  return(list(fit=fit,dc=NULL,model=model,normData=NULL,blocked_analysis=NULL))
}


.extra_sconline.Fit_LimmaVoomFn=function(sl_data,model,random_effect=NULL,quantile.norm=F,sample.weights=F,TMMnorm=F,bkg_genes=NULL,dc.object=NULL){
  require(edgeR)
  require(limma)
  
  
  dge <- DGEList(counts=counts(sl_data))
  if(is.null(bkg_genes)){
    tmpCount=apply(counts(sl_data),1,function(x) sum(x>0))
    keep=tmpCount>10
    
  } else {
    keep=row.names(sl_data) %in% bkg_genes
  }
  
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  if(TMMnorm){
    dge <- calcNormFactors(dge)
  }
  
  if(sample.weights){
    if(quantile.norm){
      logCPM <- voomWithQualityWeights(dge, model, normalize.method="quantile",plot = T,save.plot=T)
    } else {
      logCPM <- voomWithQualityWeights(dge, model,plot = T,save.plot=T)
    }
  } else {
    if(quantile.norm){
      logCPM <- voom(dge, model, normalize.method="quantile",plot = T,save.plot=T)
    } else {
      logCPM <- voom(dge, model,plot = T,save.plot=T)
    }
  }
  
  dc=NULL
  if(!is.null(random_effect)){
    if(is.null(dc.object)){
      dc <- .extra_sconline.duplicateCorrelation(logCPM,design=model, block=colData(sl_data)[,random_effect])
    } else {
      dc=dc.object
    }
  }
  
  
  blocked_analysis=F
  if(!is.null(dc)){
    if(!is.nan(dc$consensus.correlation)){
      if(abs(dc$consensus.correlation)<0.9){
        fit <- lmFit(logCPM, model,block = colData(sl_data)[,random_effect], correlation=dc$consensus.correlation)
        blocked_analysis=T
      } else {
        fit <- lmFit(logCPM, model)
      }
    } else {
      fit <- lmFit(logCPM, model)
    }
  } else {
    fit <- lmFit(logCPM, model)
  }
  
  
  #logCPM=preprocessCore::normalize.quantiles(logCPM)
  
  
  return(list(fit=fit,dc=dc,model=model,blocked_analysis=blocked_analysis,normData=logCPM))
}

#inputExpData=x;covariates=c("genotype","scaleQC_GeneTotalCount","scalePseudocellSize");randomEffect="orig.ident"
#DEmethod="Trend";normalization="CPM";quantile.norm=F;VST_fitType="parametric";prior.count=3


.sconline.fitLimmaFn=function(inputExpData,covariates,randomEffect,DEmethod="Trend",normalization="CPM",quantile.norm=F,bkg_genes=NULL,VST_fitType="parametric",prior.count=1,include.malat1.as.covariate=F,dc.object=NULL,dream_ncores=4){
  
  #for covariates use the scaled versions of nUMI (), nGene(QC_Gene_unique_count_scale), and pseudocell_size as appropriate
  #include cellType as covariate as appropriate
  
  #DEmethod="Trend";normalization="CPM";quantile.norm=F;VST_fitType="parametric";prior.count=1
  
  normalization=match.arg(normalization,c("CPM", "TMM", "VST", "rmTop50","none"))
  DEmethod=match.arg(DEmethod,c("Trend", "Voom", "VoomSampleWeights","Dream"))
  
  if(DEmethod=="Voom"&normalization=="VST"){
    stop("normalization method can be either Voom or VST but not both")
  }
  
  norm.tmm=normalization=="TMM"
  norm.rmTop50=normalization=="rmTop50"
  norm.vst=normalization=="VST"
  no_normalization=normalization=="none"
  
  if(is.null(inputExpData)){
    stop("Expression data is missing!")
  } else if(tolower(class(inputExpData))=="seurat"){
    inputExpData=.sconline.convertSeuratToExpressionSet(object=inputExpData)
  } else if(class(inputExpData)!="SingleCellExperiment") {
    stop("Unrecognized inputExpression data!")
  }
  
  if(sum(colnames(colData(inputExpData)) %in% covariates)<length(covariates)){
    stop(paste0("Covariates ",paste(setdiff(covariates,colnames(colData(inputExpData))),collapse = ", ")," were not identified in the inputExpData!"))
  }
  
  
  
  if(DEmethod!="Dream"){
    model_matrix="~0"
    for(icov in covariates){
      if(length(unique(colData(inputExpData)[,icov]))>1){
        if(class(colData(inputExpData)[,icov])==class(factor)){
          colData(inputExpData)[,icov]=as.character(colData(inputExpData)[,icov])
        }
        model_matrix=paste0(model_matrix,"+",icov)
      } else {
        warning(paste0("Excluding ",icov," covariate as it has only one level!"))
      }
    }
    model_matrix = model.matrix(as.formula(model_matrix),data=colData(inputExpData))
  } else {
    model_matrix="~"
    for(icov in covariates){
      if(length(unique(colData(inputExpData)[,icov]))>1){
        if(class(colData(inputExpData)[,icov])==class(factor)){
          colData(inputExpData)[,icov]=as.character(colData(inputExpData)[,icov])
        }
        if(model_matrix=="~"){
          model_matrix=paste0(model_matrix,icov)
        } else {
          model_matrix=paste0(model_matrix,"+",icov)
        }
        
      } else {
        warning(paste0("Excluding ",icov," covariate as it has only one level!"))
      }
    }
    model_matrix = as.formula(paste0(model_matrix," + (1|",randomEffect,")"))
  }
  
  if(norm.rmTop50){
    if(sum(colnames(rowData(inputExpData))=="QC_top50_expressed")>0){
      inputExpData=inputExpData[which(rowData(inputExpData)$QC_top50_expressed=="No"),]
    } else {
      warning("Column QC_top50_expressed was not identified in the gene attribute dataframe. skipping the removal of top 50 expressed genes!")
    }
  }
  
  #sl_data=inputExpData;model=model_matrix;random_effect=randomEffect;quantile.norm = quantile.norm;TMMnorm=norm.tmm;VSTnorm=norm.vst;prior.count=prior.count;pd=colData(inputExpData)
  switch(DEmethod,
         Trend=.extra_sconline.Fit_LimmaTrendFn(sl_data=inputExpData,model=model_matrix,random_effect=randomEffect,quantile.norm = quantile.norm,TMMnorm=norm.tmm,VSTnorm=norm.vst,prior.count=prior.count,bkg_genes=bkg_genes,no_normalization=no_normalization,dc.object=dc.object,include.malat1.as.covariate),
         Dream=.extra_sconline.Fit_LimmaDreamFn(sl_data=inputExpData,model=model_matrix,pd=as.data.frame(colData(inputExpData)),random_effect=randomEffect,quantile.norm = quantile.norm,TMMnorm=norm.tmm,VSTnorm=norm.vst,prior.count=prior.count,bkg_genes=bkg_genes,no_normalization=no_normalization,dc.object=dc.object,include.malat1.as.covariate,ncores = dream_ncores),
         Voom=.extra_sconline.Fit_LimmaVoomFn(sl_data=inputExpData,model=model_matrix,random_effect=randomEffect,quantile.norm=quantile.norm,sample.weights=F,TMMnorm=norm.tmm,bkg_genes=bkg_genes,dc.object=dc.object),
         VoomSampleWeights=.extra_sconline.Fit_LimmaVoomFn(sl_data=inputExpData,model=model_matrix,random_effect=randomEffect,quantile.norm=quantile.norm,sample.weights=T,TMMnorm=norm.tmm,bkg_genes=bkg_genes,dc.object=dc.object))
  
}

#argList=.ArgList;umap.method='uwot';generateUMAP=T;saveFiles=T;input_UMAP_embedding=NULL
.sconline.umapFn=function(argList,umap.method='umap-learn',generateUMAP=T,saveFiles=T,input_UMAP_embedding=NULL){
  reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F));F}, error=function(e) {return(T)})
  
  if((reRunCheck|argList$newRun)&generateUMAP&is.null(input_UMAP_embedding)){
    cat("Running UMAP\n")
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    
    #load(.myFilePathMakerFn("pca_centroids",argList=argList))
    #pca_centroid=pca_centroid[,1:argList$nPCs]
    
    
    
    tst=.reductionUMAPFn(harmony_embeddings,umap.method=umap.method,testPCAembeddings=NULL,n.neighbors =40)
    resUMAP=tst$embedding
    .mySaveFn(resUMAP,file=.myFilePathMakerFn("UMAP_res",argList=argList,pseudoImportant = F))
    
    x=as.data.frame(resUMAP)
    x=cbind(x,pd)
    row.names(x)=row.names(pd)
    pd=x
    .mySaveFn(pd,file=.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
  } else if(!is.null(input_UMAP_embedding)){
    load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    if(sum(row.names(pd) %in% row.names(input_UMAP_embedding))<nrow(pd)){
      stop("some of the cells don't have UMAP embedding. row.names of the input_UMAP_embedding should contain all cell names")
    }
    input_UMAP_embedding=input_UMAP_embedding[match(row.names(pd),row.names(input_UMAP_embedding)),]
    pd$UMAP_1=input_UMAP_embedding[,"UMAP_1"]
    pd$UMAP_2=input_UMAP_embedding[,"UMAP_2"]
    .mySaveFn(pd,file=.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else if(!generateUMAP){
    load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    pd$UMAP_1=1
    pd$UMAP_2=2
    .mySaveFn(pd,file=.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  }
  
  reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_centroid",argList=argList));F}, error=function(e) {return(T)})
  if((reRunCheck|argList$newRun)|(!is.null(input_UMAP_embedding))){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
    UMAP_centroid=data.frame(centroid=sl_pseudo$cluster,pd[sl_pseudo$pseudocell,c("UMAP_1","UMAP_2")],stringsAsFactors = F)
    row.names(UMAP_centroid)=UMAP_centroid$centroid
    .mySaveFn(UMAP_centroid,file=.myFilePathMakerFn("UMAP_centroid",argList=argList))
  }
  
  
  return("Done!")
}

.extra_sconline.JaccardDistanceFn=function(meta_z_mat,meta_z_mat_bkg=NULL,sig1_thr=3,sig2_thr=1,mode="ratio"){
  res_mat1=meta_z_mat
  if(is.null(meta_z_mat_bkg)){
    meta_z_mat_bkg=meta_z_mat
  }
  res_mat2=meta_z_mat_bkg
  res_mat1[which(meta_z_mat<sig1_thr)]=0
  res_mat1[which(res_mat1>0)]=1
  
  res_mat2[which(meta_z_mat_bkg<sig2_thr)]=0
  res_mat2[which(res_mat2>0)]=1
  
  res_mat=res_mat1 %*% t(res_mat2)
  if(mode=="ratio"){
    res_mat=sweep(res_mat,1,pmax(rowSums(res_mat1),1),"/")
  } else if(mode=="difference"){
    res_mat=sweep(res_mat,1,pmax(rowSums(res_mat1),1),"-")*(-1)
  }
  
  row.names(res_mat)=colnames(res_mat)=row.names(meta_z_mat)
  return(res_mat)
}

.extra_sconline.densityPeakClustering=function(meta_z_mat,cosine_dist,de_dist,pct_de_count_thr=0,sig1_thr,min_marker_thr,pseudocells=NULL){
  
  if(!is.null(pseudocells)){
    cosine_dist=cosine_dist[pseudocells,pseudocells]
    de_dist=de_dist[pseudocells,pseudocells]
    if(length(pseudocells)<3){
      return(list(D=pseudocells[1],cluster_assignments=data.frame(local_density=0,centroid=pseudocells,sigma=0,source_name=pseudocells[1])))
    }
  }
  
  local_density=cosine_dist#de_dist
  local_density[is.na(local_density)]=0
  local_density[t(de_dist)>pct_de_count_thr]=0
  local_density=.extra_sconline.AffinityFn(sim_mat = local_density)
  local_density_thr=apply(local_density,1,function(x) quantile(x,0.95))
  local_density2=apply(local_density,1,function(x) sum(x[x>=quantile(local_density_thr,0.5)],na.rm = T))
  local_density=apply(local_density,1,function(x) sum(x>=quantile(local_density_thr,0.75)))
  
  local_density=local_density[order(local_density,local_density2,decreasing = T)]
  #ggplot(data=data$input_umap_centroid,aes(UMAP_1,UMAP_2,label=centroid))+geom_point()+geom_label(data=data$input_umap_centroid[data$input_umap_centroid$centroid %in% gsub("C","",names(local_density[1:40])),])
  
  sig_thr=apply(meta_z_mat,1,function(x) sum(x>=sig1_thr))
  local_density=local_density[names(local_density) %in% row.names(meta_z_mat)[sig_thr>=min_marker_thr]]
  
  de_dist_sym=de_dist+t(de_dist)
  
  res_min_dist=0
  res_source_name=names(local_density)[1]
  for(i in 2:length(local_density)){
    source_name_list=names(local_density)[1:(i-1)]
    tmp_score=de_dist_sym[names(local_density[i]),source_name_list]
    sl_id=which(tmp_score==min(tmp_score))[1]
    if(length(sl_id)>1){
      stop("Here")
    }
    res_source_name=c(res_source_name,source_name_list[sl_id])
    res_min_dist=c(res_min_dist,tmp_score[sl_id])
  }
  res_min_dist[1]=max(res_min_dist)
  
  
  res=data.frame(local_density=local_density,centroid=names(local_density),sigma=res_min_dist,source_name=res_source_name,stringsAsFactors = F)
  res=res[order(res$sigma,decreasing = T),]
  
  #ggplot(data=data$input_umap_centroid,aes(UMAP_1,UMAP_2,label=centroid))+geom_point()+geom_label(data=data$input_umap_centroid[data$input_umap_centroid$centroid %in% gsub("C","",res$centroid[1:3]),])
  
  #.de_dist=de_dist
  #de_dist=de_dist+t(de_dist)
  checkFlag=T
  counter=1
  D=res$centroid[1]
  while(checkFlag&counter<nrow(de_dist_sym)){
    counter=counter+1
    #print(counter)
    
    tmp=de_dist_sym[res$centroid[counter],D]
    
    if(min(tmp)>=pct_de_count_thr&res$sigma[counter]>=pct_de_count_thr){
      D=c(D,res$centroid[counter])
      
    } else {
      checkFlag=F
    }
  }
  
  for(center in D){
    res$source_name[res$centroid==center]=center
  }
  
  for(center in D){
    res$source_name[res$centroid==center]=center
    centerList=res$centroid[res$source_name==center]
    preCenters=""
    while(!identical(centerList,preCenters)){
      preCenters=centerList
      res$source_name[res$source_name %in% centerList]=center
      centerList=res$centroid[res$source_name==center]
    }
  }
  
  return(list(D=D,cluster_assignments=res))
}

.extra_sconline.AdjancyMatConstructor=function(meta_z_mat,pct_mat,de_dist,cosine_dist,pct_diff_thr,centers=NULL,symmetric=F){
  
  if(is.null(centers)){
    centers=row.names(pct_mat)
  }
  
  if(length(centers)==1){
    res=matrix(1,nrow=nrow(pct_mat),ncol=1)
    row.names(res)=row.names(pct_mat)
    colnames(res)=centers
    return(res)
  }
  
  diag(de_dist)=max(as.numeric(de_dist),na.rm = T)
  res_mat2=matrix(de_dist[,centers],ncol=length(centers))
  row.names(res_mat2)=row.names(de_dist)
  colnames(res_mat2)=centers
  
  res_mat2=t(apply(res_mat2,1,function(x) {
    x[x>min(x)]=NA
    x
  }))
  res_mat2[!is.na(res_mat2)]=1
  res_mat2[is.na(res_mat2)]=0
  if(sum(rowSums(res_mat2,na.rm = T)>1)>0){
    for(i in which(rowSums(res_mat2,na.rm = T)>1)){
      sl_ind=which(res_mat2[i,]>0)
      sl_scores=cosine_dist[i,sl_ind]
      sl_ind2=sl_ind[which(sl_scores==max(sl_scores))]
      res_mat2[i,sl_ind]=0
      res_mat2[i,sl_ind2]=1
    }
  }
  
  res_mat2_binary=matrix(pct_mat[centers,],nrow=length(centers))
  row.names(res_mat2_binary)=centers
  colnames(res_mat2_binary)=colnames(meta_z_mat)
  res_mat_residual=res_mat2 %*% res_mat2_binary
  
  res_mat1=pct_mat
  res_mat1[which(.extra_sconline.PctDiffscoreFn(res_mat1-res_mat_residual,pct_diff_thr = pct_diff_thr)*.sconline_extra_Pct2scoreFn(res_mat_residual,pct2_thr = pct2_thr)>0.99)]=0
  
  de_dist2=.sconline_extra_PctScoreFn(meta_z_mat = meta_z_mat,pct_mat=res_mat1,pct_mat_ref = pct_mat,pct2_thr = pct2_thr,pct_diff_thr = pct_diff_thr,centers=centers,symmetric = symmetric)
  de_dist2=de_dist2+max(as.numeric(de_dist2))*res_mat2
  de_dist2[cbind(match(centers,row.names(de_dist2)),match(centers,colnames(de_dist2)))]=max(as.numeric(de_dist2))
  thr=apply(de_dist2,1,min)
  res_mat3=t(apply(de_dist2,1,function(x) {
    x[x>min(x)]=NA
    x
  }))
  res_mat3[!is.na(res_mat3)]=1
  res_mat3[is.na(res_mat3)]=0
  if(sum(rowSums(res_mat3,na.rm = T)>1)>0){
    for(i in which(rowSums(res_mat3,na.rm = T)>1)){
      sl_ind=which(res_mat3[i,]>0)
      sl_scores=cosine_dist[i,sl_ind]
      sl_ind2=sl_ind[which(sl_scores==max(sl_scores))]
      res_mat3[i,sl_ind]=0
      res_mat3[i,sl_ind2]=1
    }
  }
  
  res_mat_net=res_mat2+res_mat3
  res_mat_net[cbind(match(centers,row.names(res_mat_net)),match(centers,colnames(res_mat_net)))]=0
  final_residual=apply(de_dist[,centers],1,min)
  thr=thr[match(names(final_residual),names(thr))]
  final_residual=data.frame(fdeg=final_residual,sdeg=thr,stringsAsFactors = F)
  return(list(res_mat_net=res_mat_net,residual=final_residual))
}

.sconline_ClusteringFn_archive=function(argList,min_marker_thr=20,sig1_thr=3,sig2_thr=1,pct_de_count_thr=1,pct_diff_thr=0.2,pct2_thr=0.3){
  
  #min_marker_thr=50;sig1_thr=3;sig2_thr=1;pct_de_count_thr=2;pct_diff_thr=0.2;pct2_thr=0.3
  
  meta_z_list=.sconline.fetch_data("meta_z",argList = argList)
  
  meta_z_mat=meta_z_list$meta_z
  cosine_dist=meta_z_list$cosine_dist
  de_dist=meta_z_list$de_pct_dist
  pct_mat=meta_z_list$med_pct.1
  
  j_dist=.extra_sconline.JaccardDistanceFn(meta_z_mat=meta_z_mat,meta_z_mat_bkg=NULL,sig1_thr=sig1_thr,sig2_thr=sig2_thr,mode="difference")
  j_dist=j_dist[row.names(de_dist),colnames(de_dist)]
  
  de_dist[j_dist<min_marker_thr]=0
  
  
  net_components=which(de_dist<=pct_de_count_thr,arr.ind = T)
  net_components=data.frame(from=row.names(de_dist)[net_components[,1]],to=colnames(de_dist)[net_components[,2]],stringsAsFactors = F)
  net_components=igraph::graph_from_data_frame(net_components,directed = F)
  net_components=igraph::components(net_components)
  net_components=data.frame(pseudocell=names(net_components$membership),membership=net_components$membership,stringsAsFactors = F)
  D=c()
  cluster_assignments=NULL
  for(icomponent in unique(net_components$membership)){
    sl_pseudocell=net_components$pseudocell[net_components$membership==icomponent]
    res_d=.extra_sconline.densityPeakClustering(meta_z_mat=meta_z_mat,cosine_dist=cosine_dist,de_dist=de_dist,pseudocells=sl_pseudocell,pct_de_count_thr = pct_de_count_thr,sig1_thr=sig1_thr,min_marker_thr=min_marker_thr)
    D=c(D,res_d$D)
    cluster_assignments=rbind(cluster_assignments,res_d$cluster_assignments[,c("centroid","source_name")])
  }
  
  
  #tst=.extra_sconline.AdjancyMatConstructor(meta_z_mat = meta_z_mat,pct_mat=pct_mat,de_dist=de_dist+t(de_dist),cosine_dist=cosine_dist,pct_diff_thr=pct_diff_thr,centers=D,symmetric=F)
  
  res_intermediates=NULL
  for(i in setdiff(row.names(pct_mat),D)){
    m2 =glmnet::glmnet(t(pct_mat[setdiff(D,i),]), y=pct_mat[i,], lower.limits = 0,lambda = 0,alpha=1,family="gaussian", intercept = FALSE)
    m2=coef(m2)
    m2=m2[row.names(m2)!="(Intercept)",]
    m2=m2[order(m2,decreasing = T)]
    
    m2 =glmnet::glmnet(t(meta_z_mat[names(m2)[1:2],]), y=meta_z_mat[i,], lower.limits = 0,lambda = 0,alpha=1,family="gaussian", intercept = FALSE)
    m2=coef(m2)
    m2=m2[row.names(m2)!="(Intercept)",]
    m2=m2[order(m2,decreasing = T)]
    
    
    res_intermediates=rbind(res_intermediates,data.frame(source1=names(m2)[1],source2=names(m2)[2],target=i,value1=m2[1],value2=m2[2],ratio=m2[1]/m2[2],cosine_count=sum(cosine_dist[i,names(m2)[1:2]]>0),stringsAsFactors = F))
  }
  
  res_intermediates=res_intermediates[which(res_intermediates$ratio<7&res_intermediates$cosine_count>1),]
  res_intermediates=res_intermediates[order(res_intermediates$ratio,decreasing = F),]
  
  res_intermediates=res_intermediates[!duplicated(.extra_sconline.NetIdFn(res_intermediates)),]
  tormList=c()
  for(i in 1:nrow(res_intermediates)){
    if(sum(de_dist[res_intermediates$target[i],cluster_assignments$centroid[cluster_assignments$source_name==res_intermediates$source1[i]]]>0)==0){
      tormList=c(tormList,i)
    }
    if(sum(de_dist[res_intermediates$target[i],cluster_assignments$centroid[cluster_assignments$source_name==res_intermediates$source2[i]]]>0)==0){
      tormList=c(tormList,i)
    }
  }
  if(length(tormList)>0){
    tormList=unique(tormList)
    res_intermediates=res_intermediates[-tormList,]
  }
  
  net_directional=data.frame(source=c(res_intermediates$source1,res_intermediates$source2),target=c(res_intermediates$target,res_intermediates$target),score=1,stringsAsFactors = F)
  
  
  #exNodes=row.names(tst$residual)[which(tst$residual[,1]==0)]
  
  #ggplot(data=data$input_umap_centroid,aes(UMAP_1,UMAP_2,label=centroid))+geom_point()+geom_label(data=data$input_umap_centroid[data$input_umap_centroid$centroid %in% gsub("C","",D),])
  
  gc()
  
  
  
  #cosine_dist=mymakeCosine(meta_z$meta_z,wJaccardDistance = F)
  #cosine_dist[D,D]
  
  #min_marker_thr=10
  
  
  
  #res_mat_net=.sconline_extra_AdjancyMatConstructor(pct_mat = pct_mat[D,],de_dist = de_dist[D,D],cosine_dist = cosine_dist[D,D])
  #res_mat_net=res_mat_net$res_mat_net
  #res_mat_net_indx=which(res_mat_net>0,arr.ind = T)
  #net=data.frame(target=row.names(res_mat_net)[res_mat_net_indx[,1]],source=row.names(res_mat_net)[res_mat_net_indx[,2]],score=1)
  #.extra_sconline.NetVisFn(net[!duplicated(myNetIdFn(net)),])
  
  #D=c("C16","C13","C126")
  
  colnames(cluster_assignments)=c("pseudocell","cluster")
  
  
  
  #cluster_assignments=.sconline.extra_clustering_assignments(res_mat=cosine_dist,centroids=D,not_sig_centroids_org=NULL,binary_thr=0.8,return_affinities=F)
  #cluster_assignments$cluster[cluster_assignments$pseudocell=="C160"]="C_14"
  #cluster_assignments$cluster[cluster_assignments$cluster=="C_3"]="C_1"
  #cluster_assignments$cluster[cluster_assignments$pseudocell %in% c("C26","C157")]="C_4"
  
  if(length(unique(cluster_assignments$cluster))<12&length(unique(cluster_assignments$cluster))>3){
    cluster_assignments$color=rcartocolor::carto_pal(length(unique(cluster_assignments$cluster)), "Vivid")[as.numeric(factor(cluster_assignments$cluster))]
  } else {
    cluster_assignments$color=hues::iwanthue(length(unique(cluster_assignments$cluster)))[as.numeric(factor(cluster_assignments$cluster))]
  }
  
  cluster_assignments$color[which(cluster_assignments$cluster=="C_0")]="gray"
  
  if(sum(is.na(cluster_assignments$color))>0){
    colPallette=cluster_assignments[cluster_assignments$pseudocell %in% D,]
    clust_affinities=.sconline.extra_clustering_assignments(res_mat=cosine_dist,centroids=D,not_sig_centroids_org=NULL,binary_thr=0.6,return_affinities=T)
    clust_affinities=clust_affinities[which(row.names(clust_affinities) %in% cluster_assignments$pseudocell[is.na(cluster_assignments$cluster)]),]
    
    col_set=.sconline_extra_clustering_color_gradient(affinity_mat = clust_affinities,col_rgb_pallette = colPallette)
    col_set=data.frame(pseudocell=row.names(clust_affinities),color2=col_set,stringsAsFactors = F)
    col_set=col_set[match(cluster_assignments$pseudocell,col_set$pseudocell),]
    cluster_assignments$color[is.na(cluster_assignments$color)]=col_set$color2[is.na(cluster_assignments$color)]
  }
  
  
  p2=""
  p1=""
  {
    
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
    input_umap_centroid=UMAP_centroid
    
    require(ggnetwork)
    res_net=which(upper.tri(de_dist)&de_dist<1&cosine_dist>0.8,arr.ind = T)
    res_net=data.frame(Fnode=row.names(de_dist)[res_net[,1]],Snode=colnames(de_dist)[res_net[,2]],score=de_dist[res_net],stringsAsFactors = F)
    res_net$score=(2.5*(res_net$score-0.5))^3
    
    
    #res_net=which(res_mat>=0.5,arr.ind = T)
    #res_net=data.frame(Fnode=row.names(res_mat)[res_net[,1]],Snode=colnames(res_mat)[res_net[,2]],score=res_mat[res_net],stringsAsFactors = F)
    #res_net$score=(2.5*(res_net$score-0.5))^3
    
    #id=paste0(res_net$Fnode,"_",res_net$Snode)
    #id[res_net$Fnode>res_net$Snode]=paste0(res_net$Snode,"_",res_net$Fnode)[res_net$Fnode>res_net$Snode]
    #res_net=res_net[!duplicated(id),]
    
    #res_net=sg2
    #res_net$Fnode=res_net$from
    #res_net$Snode=res_net$to
    #res_net$score=res_net$weight
    #res_net=res_net[,c("Fnode","Snode","score")]
    
    if(length(setdiff(colnames(de_dist),c(res_net$Fnode,res_net$Snode)))>0){
      res_net=rbind(res_net,data.frame(Fnode=setdiff(colnames(de_dist),c(res_net$Fnode,res_net$Snode)),Snode=setdiff(colnames(de_dist),c(res_net$Fnode,res_net$Snode)),score=0,stringsAsFactors = F))
    }
    net = network::network(res_net, directed = FALSE,matrix.type="edgelist")
    
    
    library(gplots)
    library(devtools)
    library(ggnetwork)
    require(network)
    require(sna)
    require(ggplot2)
    require(Matrix)
    
    netVerNames=network::network.vertex.names(net)
    resClustering=cluster_assignments
    #resClustering$cluster[which(resClustering$cluster=="C_0")]=NA
    net %v% "cluster" = resClustering$cluster[match(netVerNames,as.character(resClustering$pseudocell))]
    network::set.edge.attribute(net, "weight", res_net$score)
    
    
    centroid_layout=input_umap_centroid[match(gsub("C","",as.character(netVerNames)),as.character(input_umap_centroid$centroid)),-1]
    centroid_layout=as.matrix(centroid_layout)
    colnames(centroid_layout)=c("x","y")
    
    net=ggnetwork:::fortify.network(net,layout = centroid_layout)
    
    clusterCol=resClustering[!duplicated(resClustering$cluster),]
    
    pd_summary=.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"))
    
    p1=.extra_sconline.NetVisFn(net=net_directional,input_pd=pd,input_umap_centroid=input_umap_centroid)
    
    resClustering$pseudocell=gsub("C","",resClustering$pseudocell)
    lda.training.data=merge(resClustering,input_umap_centroid,by.x="pseudocell",by.y="centroid")
    #lda.fit=MASS::lda(cluster~UMAP_1+UMAP_2,data=lda.training.data)
    #lda.pred=predict(lda.fit,pd_summary)
    
    knet=RANN::nn2(query=pd_summary[,c("UMAP_1","UMAP_2")],data = lda.training.data[,c("UMAP_1","UMAP_2")],k=12,eps=0)
    affinity_mat=matrix(0,nrow=nrow(knet$nn.idx),ncol=length(unique(lda.training.data$pseudocell)))
    colnames(affinity_mat)=unique(lda.training.data$pseudocell)
    
    tmp_affinities=t(apply(knet$nn.dists,1,function(x) exp((-1)*(x/x[2])^2)))
    tmp_affinities=t(apply(tmp_affinities,1,function(x) x/sum(x)))
    for(itr in 1:ncol(knet$nn.idx)){
      tst=.myOneHotFn(inputVector=factor(lda.training.data$pseudocell[knet$nn.idx[,itr]],levels=unique(lda.training.data$pseudocell)))
      tst=tst[,match(colnames(affinity_mat),colnames(tst))]
      tst=as.matrix(tst)
      tst=tmp_affinities[,itr]*tst
      affinity_mat=affinity_mat+tst
      #print(paste(itr,":",affinity_mat[6151,2]))
      rm(tst)
    }
    
    
    col_rgb_pallette=resClustering
    col_rgb_pallette=col_rgb_pallette[match(colnames(affinity_mat),col_rgb_pallette$pseudocell),]
    col_rgb1=matrix(grDevices::col2rgb(col_rgb_pallette$color)[1,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    col_rgb2=matrix(grDevices::col2rgb(col_rgb_pallette$color)[2,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    col_rgb3=matrix(grDevices::col2rgb(col_rgb_pallette$color)[3,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    
    colnames(col_rgb1)=col_rgb_pallette$cluster
    colnames(col_rgb2)=col_rgb_pallette$cluster
    colnames(col_rgb3)=col_rgb_pallette$cluster
    
    col_rgb1=affinity_mat*col_rgb1
    col_rgb2=affinity_mat*col_rgb2
    col_rgb3=affinity_mat*col_rgb3
    
    col_rgb1=rowSums(col_rgb1)
    col_rgb2=rowSums(col_rgb2)
    col_rgb3=rowSums(col_rgb3)
    
    pd_summary$color=rgb(red=col_rgb1,green=col_rgb2,blue=col_rgb3,maxColorValue = 255)
    pd_summary_main=pd_summary
    
    scale_factor=net[!duplicated(net$vertex.names),]
    scale_factor$vertex.names=gsub("C","",scale_factor$vertex.names)
    scale_factor=merge(scale_factor,input_umap_centroid,by.x="vertex.names",by.y="centroid")
    scale_factor1=lm(x~UMAP_1,data=scale_factor)
    pd_summary$UMAP_1=predict(scale_factor1,newdata=pd_summary)
    scale_factor2=lm(y~UMAP_2,data=scale_factor)
    pd_summary$UMAP_2=predict(scale_factor2,newdata=pd_summary)
    
    pd_summary$xend=pd_summary$UMAP_1
    pd_summary$yend=pd_summary$UMAP_2
    
    #predicting the background color
    net=merge(net,cluster_assignments[,c("pseudocell","color")],by.x="vertex.names",by.y="pseudocell")
    #net$pct.sig="#FFFFFF"
    #net$pct.sig[net$vertex.names %in% cluster_assignments$pseudocell[cluster_assignments$cluster %in% cluster_assignments$cluster[cluster_assignments$pseudocell %in% res_clustering$centroids2]]]="#FCFAD2"
    
    
    if(F){
      p=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) + geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=color))+
        geom_edges( color = "grey50",aes(size=weight)) +
        geom_point(data=net[!duplicated(net$vertex.names),],size=2,aes(fill=color),shape=21)+
        geom_nodelabel(data=net,aes(color = color, label = as.character(cluster),fill=pct.sig),fontface = "bold")+
        theme_blank()+scale_size_continuous(range = c(0.04,1.3))+scale_color_identity()+scale_fill_identity()+theme(legend.position = "none")
      
    }
    
    p2=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend))+ geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=color))+
      geom_edges( color = "black") +
      #geom_point(data=net[!duplicated(net$vertex.names)&is.na(net$cluster),],aes(fill=color,size=2),shape=21)+
      geom_point(data=net[!duplicated(net$vertex.names)&!is.na(net$cluster),],aes(fill=color,size=5),shape=21,stroke = 1)+
      #geom_label(aes(label=cluster))+
      theme_blank()+scale_size_continuous(range = c(0.3,2))+scale_fill_identity()+scale_color_identity()+theme(legend.position = "none")
    
    
  }
  
  #p
  #ggsave("~/Desktop/FullCB_clusters.png",width = 8,height=8)
  
  #net=data.frame(target=row.names(res_mat_net)[res_mat_net_indx[,1]],source=row.names(res_mat_net)[res_mat_net_indx[,2]],score=1)
  #.extra_sconline.NetVisFn(net[!duplicated(myNetIdFn(net)),])
  
  #net=data.frame(source=c("C16","C9"),target=c("C9","C126"),score=1)
  
  #myNetVisFn2(net,data_pd=pd_summary_main,convertUMAP = T)
  #ggsave("~/Desktop/UBC_pattern.png",width = 8,height=8)
  
  return(list(net=net,plot1=p1,plot2=p2,cluster_assignments=cluster_assignments,centers=D))
}

.dev.pseudocell_assignment=function(argList,anno_col="anno_cellState"){
  #argList=.ArgList;anno_col="anno_cellState"
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  pd=.sconline.fetch_data("annotation",argList = argList)
  
  ps_anno=as.matrix(.myOneHotFn(inputVector=pd[,anno_col]))
  ps_anno=prop_mat %*% ps_anno
  ps_anno=lapply(1:nrow(ps_anno),function(x) {
    y=which(ps_anno[x,]==max(ps_anno[x,]))[1]
    y=colnames(ps_anno)[y]
    return(data.frame(ps=row.names(ps_anno)[x],anno=y,purity=max(ps_anno[x,]),stringsAsFactors = F))
  })
  ps_anno=do.call("rbind",ps_anno)
  ps_anno$anno[is.na(ps_anno$anno)]="NA"
  res=matrix(0,nrow=nrow(ps_anno),ncol=nrow(ps_anno))
  row.names(res)=colnames(res)=ps_anno$ps
  
  for(i in unique(ps_anno$anno)){
    res[ps_anno$ps[which(ps_anno$anno==i)],ps_anno$ps[which(ps_anno$anno==i)]]=1
  }
  
  return(list(adj_mat=res,anno=ps_anno))
  
}

.dev.metric_value=function(argList,metric_matrix,anno_col="anno_cellState"){
  #metric_matrix=sim_mat
  
  ps_anno=.dev.pseudocell_assignment(argList = argList,anno_col = anno_col)
  adj_mat=ps_anno$adj_mat
  adj_mat=adj_mat[row.names(metric_matrix),colnames(metric_matrix)]
  ps_anno=ps_anno$anno
  ps_anno=ps_anno[match(row.names(metric_matrix),ps_anno$ps),]
  sim_mat=metric_matrix
  sim_mat[adj_mat==0]=NA
  sim_mat=apply(sim_mat,1,function(x){
    x=x[!is.na(x)]
    x=x[order(x,decreasing = T)]
    x[min(3,length(x))]
  })
  
  diff_mat=metric_matrix
  diff_mat[adj_mat==1]=NA
  diff_mat=apply(diff_mat,1,function(x){
    x=x[!is.na(x)]
    x=x[order(x,decreasing = T)]
    x[min(3,length(x))]
  })
  res=data.frame(pseudocell=row.names(metric_matrix),sim=sim_mat,diff=diff_mat,cluster=ps_anno$anno,stringsAsFactors = F)
  return(res)
}

.dev.pseudocell_metric=function(argList,cell_anno_col="anno_cellState",tol_level=0.9,inputPd){
  ps=.dev.pseudocell_assignment(argList,anno_col=cell_anno_col)
  ps=ps$anno
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  ps=ps[match(row.names(prop_mat),ps$ps),]
  
  if(is.null(inputPd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPd
  }
  
  
  
  pd=pd[!is.na(pd[,cell_anno_col]),]
  prop_mat=prop_mat[,colnames(prop_mat) %in% row.names(pd)]
  pd=pd[match(row.names(pd),colnames(prop_mat)),]
  
  prop_merged=ps$anno
  names(prop_merged)=ps$ps
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=(prop_merged))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
  prop_merged=.extra_matrix_rowNorm(input_mat = prop_merged,rowValues = 1/ qlcMatrix::rowMax(prop_merged))#Matrix::Diagonal(x=1/ qlcMatrix::rowMax(prop_merged)) %*% prop_merged
  #prop anno
  
  colMax_vals_m=qlcMatrix::colMax(prop_merged)
  colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
  prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
  #prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
  prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
  prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
  pd$anno_cluster_res=paste0("C",as.character(prop_m_hardCluster$i))
  
  
  
  ari = mclust::adjustedRandIndex(pd$anno_cluster_res, pd[,cell_anno_col])
  nmi=aricode::NMI(pd$anno_cluster_res, pd[,cell_anno_col],variant="sum")
  ami=aricode::AMI(pd$anno_cluster_res, pd[,cell_anno_col])
  ariScores=data.frame(cluster_count=n_clusters,ari=ari,nmi=nmi,ami=ami,stringsAsFactors = F)
  
  
  tmp_pd=pd
  
  tmp_clust_size=as.data.frame(table(tmp_pd[,cell_anno_col]))
  tmp_clust_size=tmp_clust_size[scale(tmp_clust_size[,2])<2,]
  sl_ind=tmp_pd[,cell_anno_col] %in% as.character(tmp_clust_size[,1])
  tmp_pd=tmp_pd[sl_ind,]
  
  
  ari = mclust::adjustedRandIndex(tmp_pd$anno_cluster_res, tmp_pd[,cell_anno_col])
  nmi=aricode::NMI(tmp_pd$anno_cluster_res, tmp_pd[,cell_anno_col],variant="sum")
  ami=aricode::AMI(tmp_pd$anno_cluster_res, tmp_pd[,cell_anno_col])
  ariScores_balanced=data.frame(cluster_count=n_clusters,ari=ari,nmi=nmi,ami=ami,stringsAsFactors = F)
  
  return(list(all_cellTypes=ariScores,rm_LargeCellTypes=ariScores_balanced))
  
}


.dev.cluster.matrix.expression.distances=function(counts, groups=NULL, useVariablegenes=F, variableGenes=NULL, dist='cor',
                                                     use.scaled.data=FALSE, min.cluster.size=1, max.n.cells=Inf) {
  require(abind)
  require(sccore)
  if(is.null(groups)) {
    stop('no groups specified')
  } else {
    groups <- as.factor(groups)
  }
  
  valid.dists <- c('JS','cor');
  if(!dist %in% valid.dists) stop(paste('only the following distance types are supported:',paste(valid.dists,collapse=', ')))
  
  cl <- factor(groups[match(rownames(counts), names(groups))],levels=levels(groups))
  tt <- table(cl)
  empty <- tt< min.cluster.size
  
  # aggregate clusters
  if(any(tt>max.n.cells)) { # need subsampling
    scn <- unlist(tapply(names(cl), cl, function(x) sample(x,min(max.n.cells,length(x)))))
    cl[!(names(cl) %in% scn)] <- NA; tt <- table(cl);
  }
  
  
  
  {
    tc <- sccore:::colSumByFactor(counts, cl);
    tc <- tc[-1,,drop=F]  # omit NA cells
    # correlation distance
    if(dist=='JS') {
      tc <- t(tc/pmax(1,rowSums(tc)))
      tcd <- pagoda2:::jsDist(tc); dimnames(tcd) <- list(colnames(tc),colnames(tc));
    } else { # correlation distance
      if(use.scaled.data){
        tc <- t(tc)
      } else{
        tc <- log10(t(tc/pmax(1,rowSums(tc)))*1e3+1)
      }
      
      if(useVariablegenes){
        if(is.null(variableGenes)){
          stop("no variable genesets provided")
        } else{
          tc=tc[row.names(tc) %in% variableGenes,]
        }
      }
      tcd <- 1-cor(tc)
    }
  }
  return(tcd)
}


.sconline.cluster_archive=function(argList,inputExpData=NULL,clustering_version=NULL,scaling_factor=10,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,affinity_param=NULL,n_clusters=NULL,inputPd=NULL,clustering_method="average",new_DEcount_run=F,tol_level=0.9,cell_count_diff_thr=0){

  #clustering_method="average"
  #argList=.ArgList;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;affinity_param=5;scaling_factor=10;inputPd=NULL;n_clusters=115;new_DEcount_run=F;cell_count_diff_thr=25
  #clustering_version="v22"
  #inputExpData=data
  
  .myaff=function (diff, K = 3, sigma = 0.5) 
  {
    #K = 3; sigma = 0.5
    diff=as.matrix(diff)
    N <- nrow(diff)
    diff <- (diff + t(diff))/2
    diag(diff) <- 0
    sortedColumns <- as.matrix(t(apply(diff, 2, sort)))
    finiteMean <- function(x) {
      return(mean(x[is.finite(x)]))
    }
    means <- apply(sortedColumns[, 1:K + 1], 1, finiteMean) + 
      .Machine$double.eps
    avg <- function(x, y) {
      return((x + y)/2)
    }
    
    Sig <- outer(means, means, avg)/3 * 2 + diff/3 + .Machine$double.eps
    Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
    
    Sig=quantile(Sig,0.1)
    densities <- dnorm(as.matrix(diff), 0, sigma * Sig, log = FALSE)
    W <- (densities + t(densities))/2
    return(W)
  }
  
  .mynetAff=function(inputData,rev=T,do.scale=T){
    if(do.scale){
      aff=t(apply(log2(inputData+1),1,scale))
      aff[is.na(aff)]=0
      aff=(aff+t(aff))/sqrt(2)
      row.names(aff)=colnames(aff)=row.names(inputData)
    } else {
      aff=inputData
    }
    
    
    
    #aff=pnorm(aff,lower.tail = F)
    if(T){
      
      if(rev){
        diff=aff
        if(do.scale){
          diff[aff>0|t(aff)>0]=0
        }
        
        cor_map=diff
        cor_map=cor_map*(-1)
      } else {
        diff=aff
        if(do.scale){
          diff[aff<0|t(aff)<0]=0
        }
        
        cor_map=diff
      }
      
      net=NULL
      nn.rank=matrix(0,nrow=nrow(cor_map),ncol=5)
      row.names(nn.rank)=row.names(cor_map)
      for(i in 1:nrow(cor_map)){
        tmp=cor_map[i,]
        names(tmp)=colnames(cor_map)
        tmp=tmp[order(tmp,decreasing = T)]
        
        #sl_neighbors=names(tmp)[which(tmp>0.1)]
        sl_neighbors=names(tmp)[tmp>=(0.98*tmp[4])]
        if(length(sl_neighbors)>0){
          net=rbind(net,data.frame(source=row.names(cor_map)[i],target=sl_neighbors,score=cor_map[row.names(cor_map)[i],sl_neighbors],stringsAsFactors = F))
        }
        
      }
      net$name=paste0(net$source,"_",net$target)
      net$name[net$source>net$target]=paste0(net$target,"_",net$source)[net$source>net$target]
      net=net[!duplicated(net$name),]
      #net$score=1
      net=reshape2::dcast(source~target,data=net,value.var = "score")
      row.names(net)=net[,1]
      net=net[,-1]
      net_c_names=setdiff(row.names(aff),colnames(net))
      net_c=matrix(0,nrow=nrow(net),ncol=length(net_c_names))
      colnames(net_c)=net_c_names
      net=cbind(net,net_c)
      net=net[match(row.names(aff),row.names(net)),colnames(aff)]
      net[is.na(net)]=0
      row.names(net)=row.names(aff)
      net=net+t(net)
      net=as.matrix(net)/2
      #diag(net)=1
      net=as(net,"dgCMatrix")
      net=.extra_matrix_rowNorm(net)#Matrix::Diagonal(x=1/rowSums(net)) %*% net
      net2=net
      for(iii in 1:20){
        net2=(net2) %*% t(net)
        net2=.extra_matrix_rowNorm(net2)#Matrix::Diagonal(x=1/rowSums(net2)) %*% net2
        net2=net2+net
        net2=net2/2
      }
      net=net2
      if(F){
        diag(net)=0
      }
      
      net=.extra_matrix_rowNorm(input_mat = net,rowValues = 1/qlcMatrix::rowMax(net))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(net)) %*% net
      #diag(net)=1
      
    }
    return(net)
  }
  
  
  
  if(is.null(clustering_version)){
    clustering_version="v58"
  }
  
  match.arg(clustering_version,c(paste0("v",1:60),"random"))
  
  clustering_method=tolower(clustering_method)
  if(sum(clustering_method %in% c("complete", "average", "mcquitty"))==0){
    stop("Acceptable clustering methods: complete, average, mcquitty")
  }
  
  if(is.null(inputPd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPd
  }
  
  
  if(is.null(n_clusters)){
    
    stop("n_clusters argument should be provided!")
  }
  
  #res_array=.sconline.fetch_data("sconline_arrays",argList = argList)
  meta_data=.sconline.fetch_data("meta_z",argList = argList)
  de_prefix=paste0("de_sig",sig1_thr,"_pctdiff",pct_diff_thr,"_pct2",pct2_thr)
  if(cell_count_diff_thr>0){
    de_prefix=paste0(de_prefix,"_cellCount",cell_count_diff_thr)
  }
  #pct_mat=meta_data$med_pct.1;argList = argList;meta_z_mat=meta_data$meta_z;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
  if(!file.exists(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))|new_DEcount_run){
    de_pct_res=.extra_sconline.PctScoreFn(pct_mat=meta_data$med_pct.1,argList = argList,meta_z_mat=meta_data$meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F,cell_count_diff_thr=cell_count_diff_thr)
    
    #de_pct_res=.extra_sconline.PctScoreFn(pct_mat=.data$pct_mat,argList = .data$argList,meta_z_mat=.data$meta_z_mat,sig1_thr=.data$sig1_thr,centers=NULL,pct2_thr=.data$pct2_thr,pct_diff_thr=.data$pct_diff_thr,symmetric=F)
    
    qsave(de_pct_res,file=.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    de_pct_res=qread(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  }
  
  if(F){
    z_cosine=meta_data$med_pct.1
    z_rm=qlcMatrix::colMax(z_cosine)
    z_cosine=z_cosine[,as.numeric(z_rm)>0.2]
    z_cosine=qlcMatrix::cosSparse(t(de_pct_res))
  } else {
    z_cosine3=qlcMatrix::corSparse(t(meta_data$meta_z))
    row.names(z_cosine3)=colnames(z_cosine3)=row.names(meta_data$meta_z)
    z_cosine3[is.na(z_cosine3)]=0
    row.names(z_cosine3)=row.names(de_pct_res)
    colnames(z_cosine3)=colnames(de_pct_res)
    
    if(!clustering_version %in% c("v23","v29","v27")){
      if(clustering_version=="v21"){
        de_pct_res[z_cosine3<0]=100
        z_cosine=qlcMatrix::corSparse(t(log2(de_pct_res+1)))
        clustering_version="v20"
      } else {
        z_cosine=qlcMatrix::corSparse(t(de_pct_res+1))
      }
      
      
      
      z_cosine[is.na(z_cosine)]=0
      row.names(z_cosine)=row.names(de_pct_res)
      colnames(z_cosine)=colnames(de_pct_res)
      z_cosine2=z_cosine
      z_cosine[z_cosine<0.05]=0
      z_cosine[de_pct_res>20&t(de_pct_res)>20]=0
    }
    
  }
  
  
  
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(n_clusters>nrow(prop_mat)){
    warning(paste("Number of clusters can't be larger than the number of pseudocells! setting n_clusters to",nrow(prop_mat)))
    n_clusters=nrow(prop_mat)
  }
  
  if(F){
    ps_anno=as.matrix(.myOneHotFn(inputVector=pd$anno_orig_cellState))
    ps_anno=prop_mat %*% ps_anno
    ps_anno=lapply(1:nrow(ps_anno),function(x) {
      y=which(ps_anno[x,]==max(ps_anno[x,]))[1]
      y=colnames(ps_anno)[y]
      return(data.frame(ps=row.names(ps_anno)[x],anno=y,stringsAsFactors = F))
    })
    ps_anno=do.call("rbind",ps_anno)
    
  }
  
  {
    myL2normFn=function(inputMat){
      prop_mat2=rowSums(inputMat^2)
      prop_mat2=sqrt(prop_mat2)
      res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
      return(res)
    }
    
    #tst=prop_mat
    #tst@x=rep(1,length(tst@x))
    if(clustering_version %in% paste0("v",1:21)){
      tst=myL2normFn(inputMat = prop_mat)
      diff=t(tst)
      diff=tst %*% diff
      #diff@x=2*(exp(diff@x/max(quantile(diff@x,0.95),0.1))/(1+exp(diff@x/max(quantile(diff@x,0.95),0.1)))-0.5)
      diff=(1-diff)
    }
    
    if(clustering_version=="v1"){
      diff=diff+(1-exp(-1*(de_pct_res)/affinity_param))+pmax(1-as.matrix(z_cosine),0)^0.25
      diff=scaling_factor*diff
    } else if(clustering_version=="v2"){
      diff=diff+(1-exp(-1*(de_pct_res)/affinity_param))+pmax(1-as.matrix(z_cosine2),0)^0.5
      diff=scaling_factor*diff
    } else if(clustering_version=="v3"){
      diff=scaling_factor*(diff+(1-exp(-1*(de_pct_res)/affinity_param)))+pmax(1-as.matrix(z_cosine2),0)
    } else if(clustering_version %in% c("v22","v28")){
      
      sim_mat=.extra_sconline.pseudosim_archive(argList=argList,binarize = T)
      if(F){
        sim_mat=sim_mat=.extra_matrix_rowNorm(input_mat = sim_mat,rowValues = 1/qlcMatrix::rowMax(sim_mat))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(sim_mat)) %*% sim_mat
        sim_mat=sim_mat %*% Matrix::Diagonal(x=1/qlcMatrix::colMax(sim_mat))
        sim_mat=sim_mat %*% t(sim_mat)
        diag(sim_mat)=0
        sim_mat=.extra_matrix_rowNorm(input_mat = sim_mat,rowValues = 1/qlcMatrix::rowMax(sim_mat))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(sim_mat)) %*% sim_mat
      }
      
      
      if(F){
        sim_mat=prop_mat
        sim_mat@x=rep(1,length(sim_mat@x))
        sim_mat=sim_mat %*% t(sim_mat)
        sim_mat=.extra_matrix_rowNorm(input_mat = sim_mat,rowValues = 1/qlcMatrix::rowMax(sim_mat))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(sim_mat)) %*% sim_mat
      }
      
      
      #diag(sim_mat)=0
      #sim_mat=Matrix::Diagonal(x=1/qlcMatrix::rowMax(sim_mat)) %*% sim_mat
      sim_mat=exp(-3*(1-sim_mat))-0.04979139
      sim_mat=sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      z_cosine_binary=z_cosine3
      z_cosine_binary[z_cosine_binary<0]=0
      z_cosine_binary=Matrix::drop0(z_cosine_binary)
      z_cosine_binary@x=rep(1,length(z_cosine_binary))
      sim_mat=sim_mat*z_cosine_binary
      
      
      sim_mat=1-sim_mat
      diag(sim_mat)=0
      diff=scaling_factor*(sim_mat+(1-exp(-1*(de_pct_res)/affinity_param))+(1-exp(-1*(de_pct_res)/20)))+(pmax(1-as.matrix(z_cosine3),0))
    }else if(clustering_version=="v49"){
      
      sim_mat=.extra_sconline.pseudosim_archive(argList=argList,binarize = ncol(prop_mat)>100000)
      
      
      #diag(sim_mat)=0
      #sim_mat=Matrix::Diagonal(x=1/qlcMatrix::rowMax(sim_mat)) %*% sim_mat
      
      sim_mat=sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      sim_mat=.mynetAff(inputData = sim_mat,rev = F)
      
      
      
      
      z_cosine_binary=z_cosine3
      z_cosine_binary[z_cosine_binary<0]=0
      z_cosine_binary=Matrix::drop0(z_cosine_binary)
      z_cosine_binary@x=rep(1,length(z_cosine_binary))
      sim_mat=sim_mat*z_cosine_binary
      
      de2=.mynetAff(de_pct_res,rev=T)
      #diff=scaling_factor*(sim_mat+(1-exp(-1*(de_pct_res)/affinity_param))+(1-exp(-1*(de_pct_res)/20)))+(pmax(1-as.matrix(z_cosine3),0))+(1-cosine3)
      
      de2=cor(as.matrix(de2),as.matrix(sim_mat))
      de2[is.na(de2)]=(-1)
      diff=scaling_factor*((1-de2)*(1-exp(-1*(de_pct_res)/20))+(1-exp(-1*(de_pct_res)/affinity_param)))
      diff=pmax(diff,0)
      diff=diff+t(diff)
      diag(diff)=0
    } else if(clustering_version=="v51"){
      
      sim_mat=.extra_sconline.pseudosim_archive(argList=argList,binarize = ncol(prop_mat)<100000,cos_dist = T)
      #sim_mat2=.mynetAff(inputData = sim_mat,rev = F)
      
      #diag(sim_mat)=0
      #sim_mat=Matrix::Diagonal(x=1/qlcMatrix::rowMax(sim_mat)) %*% sim_mat
      
      sim_mat=sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      
      z_cosine_binary=.mynetAff(inputData = z_cosine3,rev = F)
      z_cosine_binary[z_cosine_binary<0.00001]=0
      z_cosine_binary=Matrix::drop0(z_cosine_binary)
      z_cosine_binary@x=rep(1,length(z_cosine_binary))
      sim_mat=sim_mat*z_cosine_binary
      
      
      
      
      diff=1-sim_mat
      diff=pmax(diff,0)
      diff=diff+t(diff)
      diag(diff)=0
    } else if(clustering_version=="v52"){
      colMax_vals_m=qlcMatrix::colMax(prop_mat)
      colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/as.numeric(colMax_vals_m))#prop_mat %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
      prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=0.9)
      prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
      prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
      prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
      prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
      groups=row.names(prop_mat)[prop_m_hardCluster$i]
      names(groups)=colnames(prop_mat)
      inputExpData=inputExpData[,names(groups)]
      
      exp_sim_mat=1-.dev.cluster.matrix.expression.distances(counts = t(counts(inputExpData)),groups = groups,dist="cor", useVariablegenes=FALSE,  use.scaled.data=FALSE)
      diff=1-exp_sim_mat
      diff=pmax(diff,0)
      diff=diff+t(diff)
      diag(diff)=0
    } else if(clustering_version=="v58"){
      
      
      affinity_param=min(affinity_param,quantile(as.numeric(de_pct_res),0.1))
      affinity_param=max(affinity_param,2)
      ps_sim_mat=.extra_sconline.pseudosim_archive11(argList=argList,binarize = ncol(prop_mat)<100000,cos_dist = T)
      #sim_mat2=.mynetAff(inputData = sim_mat,rev = F)
      
      #diag(sim_mat)=0
      #sim_mat=Matrix::Diagonal(x=1/qlcMatrix::rowMax(sim_mat)) %*% sim_mat
      
      ps_sim_mat=ps_sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      
      #sim_mat=cor(as.matrix(exp_sim_mat),as.matrix(ps_sim_mat))
      
      if(!is.null(affinity_param)){
        sl_param=affinity_param
      } else {
        sl_param=NULL
        sl_score=0
        for(affinity_param in 2:20){
          de_sim_mat=exp(-1*(de_pct_res)/affinity_param)
          
          tst=cor(as.matrix(de_sim_mat),as.matrix(ps_sim_mat))
          diag(tst)=0
          tmp_score=median(apply(tst,1,max))
          if(tmp_score>sl_score){
            sl_score=tmp_score
            sl_param=affinity_param
          }
        }
      }
      
      
      de_sim_mat=exp(-1*(de_pct_res)/sl_param)
      #diff=diff=scaling_factor*((1-ps_sim_mat)^0.5+(1-exp(-1*(de_pct_res)/affinity_param)))
      
      diag(ps_sim_mat)=0
      diag(de_sim_mat)=0
      ps_max_vals=qlcMatrix::rowMax(ps_sim_mat)
      de_max_vals=qlcMatrix::rowMax(de_sim_mat)
      ps_sim_mat=.extra_matrix_rowNorm(input_mat = ps_sim_mat,rowValues = 1/(ps_max_vals+0.0000001)) #Matrix::Diagonal(x=1/(ps_max_vals+0.0000001)) %*% ps_sim_mat
      de_sim_mat=.extra_matrix_rowNorm(input_mat = de_sim_mat,rowValues = 1/(de_max_vals+0.0000001))#Matrix::Diagonal(x=1/(de_max_vals+0.0000001)) %*% de_sim_mat
      ps_sim_mat=sweep(as.matrix(ps_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
      de_sim_mat=sweep(as.matrix(de_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
      diag(ps_sim_mat)=1
      diag(de_sim_mat)=1
      diff=diff=scaling_factor*((1-ps_sim_mat)+(1-de_sim_mat))
      
      diff=pmax(diff,0)
      diff=diff+t(diff)
      diag(diff)=0
      
    } else if(clustering_version=="v59"){
      
      
      ps_sim_mat=.extra_sconline.pseudosim_archive11(argList=argList,binarize = ncol(prop_mat)<100000,cos_dist = T)
      #sim_mat2=.mynetAff(inputData = sim_mat,rev = F)
      
      #diag(sim_mat)=0
      #sim_mat=Matrix::Diagonal(x=1/qlcMatrix::rowMax(sim_mat)) %*% sim_mat
      
      ps_sim_mat=ps_sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      
      #sim_mat=cor(as.matrix(exp_sim_mat),as.matrix(ps_sim_mat))
      
      de_sim_mat=exp(-1*(de_pct_res)/affinity_param)
      
      #diff=diff=scaling_factor*((1-ps_sim_mat)^0.5+(1-exp(-1*(de_pct_res)/affinity_param)))
      
      diag(ps_sim_mat)=0
      diag(de_sim_mat)=0
      ps_max_vals=qlcMatrix::rowMax(ps_sim_mat)
      de_max_vals=qlcMatrix::rowMax(de_sim_mat)
      ps_sim_mat=.extra_matrix_rowNorm(input_mat = ps_sim_mat,rowValues = 1/(ps_max_vals+0.000001))#Matrix::Diagonal(x=1/(ps_max_vals+0.0000001)) %*% ps_sim_mat
      de_sim_mat=.extra_matrix_rowNorm(input_mat = de_sim_mat,rowValues = 1/(de_max_vals+0.0000001))#Matrix::Diagonal(x=1/(de_max_vals+0.0000001)) %*% de_sim_mat
      ps_sim_mat=sweep(as.matrix(ps_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
      de_sim_mat=sweep(as.matrix(de_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
      diag(ps_sim_mat)=1
      diag(de_sim_mat)=1
      diff=diff=scaling_factor*((1-ps_sim_mat)+(1-de_sim_mat))
      
      diff=pmax(diff,0)
      diff=diff+t(diff)
      diag(diff)=0
      
    } else if(clustering_version=="v50"){
      
      colMax_vals_m=qlcMatrix::colMax(prop_mat)
      colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/as.numeric(colMax_vals_m))#prop_mat %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
      prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=0.99)
      prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
      prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
      prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
      prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
      groups=row.names(prop_mat)[prop_m_hardCluster$i]
      names(groups)=colnames(prop_mat)
      inputExpData=inputExpData[,names(groups)]
      
      exp_sim_mat=1-.dev.cluster.matrix.expression.distances(counts = t(counts(inputExpData)),groups = groups,dist="cor", useVariablegenes=FALSE,  use.scaled.data=FALSE)
      
      pc_sim_mat=.extra_sconline.pseudosim_archive(argList=argList,binarize = ncol(prop_mat)<100000,cos_dist = T)
      
      de2=.mynetAff(de_pct_res,rev=T)
      z_cosine_binary=de2
      z_cosine_binary=quantile(as.numeric(z_cosine_binary[z_cosine3<0]),0.95)
      de3=Matrix::drop0(de2,tol=z_cosine_binary)
      #diff=scaling_factor*(sim_mat+(1-exp(-1*(de_pct_res)/affinity_param))+(1-exp(-1*(de_pct_res)/20)))+(pmax(1-as.matrix(z_cosine3),0))+(1-cosine3)
      de_sim_mat=((de2+de3)/2)
      
      
      pc_sim_mat=pc_sim_mat+t(pc_sim_mat)
      de_sim_mat=de_sim_mat+t(de_sim_mat)
      exp_sim_mat=exp_sim_mat+t(exp_sim_mat)
      
      exp_sim_mat=exp_sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      pc_sim_mat=pc_sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      de_sim_mat=de_sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      
      
      de_distance=(de_pct_res+t(de_pct_res))/2
      de_quantile_thr=max(quantile(as.numeric(de_distance),0.05),20)
      df=data.frame(de=1-as.numeric(de_sim_mat),exp=1-as.numeric(exp_sim_mat),pc=1-as.numeric(pc_sim_mat),target=as.numeric(1-exp(-1*(de_distance)/de_quantile_thr)))
      df$weights=1
      df$weights[df$target<0.5]=sum(df$target>0.5)/sum(df$target<0.5)
      library(glmnet)
      res <- glmnet(df[,c("de","exp","pc")], df$target, lambda = 0, lower.limits = 0, intercept = FALSE,weights=df$weights)
      coef_res=as.numeric(coef(res))
      names(coef_res)=row.names(coef(res))
      
      sim_mat=de_sim_mat*coef_res["de"]+exp_sim_mat*coef_res["exp"]+pc_sim_mat*coef_res["pc"]
      sim_mat=sim_mat=.extra_matrix_rowNorm(input_mat = sim_mat,rowValues = 1/diag(sim_mat))#Matrix::Diagonal(x=1/diag(sim_mat)) %*% sim_mat
      sim_mat[sim_mat>1]=1
      sim_mat=sim_mat^3
      
      
      #sim_mat=.mynetAff(sim_mat,rev=F)
      #sim_mat2=.mynetAff(inputData = sim_mat,rev = F)
      
      #diag(sim_mat)=0
      #sim_mat=Matrix::Diagonal(x=1/qlcMatrix::rowMax(sim_mat)) %*% sim_mat
      
      
      
      #z_cosine_binary=.mynetAff(inputData = z_cosine3,rev = F)
      #z_cosine_binary[z_cosine_binary<0.00001]=0
      #z_cosine_binary=Matrix::drop0(z_cosine_binary)
      #z_cosine_binary@x=rep(1,length(z_cosine_binary))
      #sim_mat=sim_mat*z_cosine_binary
      
      diff=1-(sim_mat)
      diff=pmax(diff,0)
      diff=(diff+t(diff))/2
      diag(diff)=0
      
    } else if(clustering_version=="v55"){
      
      colMax_vals_m=qlcMatrix::colMax(prop_mat)
      colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/as.numeric(colMax_vals_m))#prop_mat %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
      prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=0.99)
      prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
      prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
      prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
      prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
      groups=row.names(prop_mat)[prop_m_hardCluster$i]
      names(groups)=colnames(prop_mat)
      inputExpData=inputExpData[,names(groups)]
      
      exp_sim_mat=1-.dev.cluster.matrix.expression.distances(counts = t(counts(inputExpData)),groups = groups,dist="cor", useVariablegenes=FALSE,  use.scaled.data=FALSE)
      
      pc_sim_mat=.extra_sconline.pseudosim_archive(argList=argList,binarize = ncol(prop_mat)<100000,cos_dist = T)
      
      de2=.mynetAff(de_pct_res,rev=T)
      z_cosine_binary=de2
      z_cosine_binary=quantile(as.numeric(z_cosine_binary[z_cosine3<0]),0.95)
      de3=Matrix::drop0(de2,tol=z_cosine_binary)
      #diff=scaling_factor*(sim_mat+(1-exp(-1*(de_pct_res)/affinity_param))+(1-exp(-1*(de_pct_res)/20)))+(pmax(1-as.matrix(z_cosine3),0))+(1-cosine3)
      de_sim_mat=((de2+de3)/2)
      
      
      pc_sim_mat=pc_sim_mat+t(pc_sim_mat)
      de_sim_mat=de_sim_mat+t(de_sim_mat)
      exp_sim_mat=exp_sim_mat+t(exp_sim_mat)
      
      exp_sim_mat=exp_sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      pc_sim_mat=pc_sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      de_sim_mat=de_sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      
      
      de_distance=(de_pct_res+t(de_pct_res))/2
      de_quantile_thr=max(quantile(as.numeric(de_distance),0.05),5)
      df=data.frame(de=1-as.numeric(de_sim_mat[upper.tri(de_sim_mat)]),exp=1-as.numeric(exp_sim_mat[upper.tri(de_sim_mat)]),pc=1-as.numeric(pc_sim_mat[upper.tri(de_sim_mat)]),target=as.numeric(1-exp(-1*(de_distance)/de_quantile_thr)[upper.tri(de_sim_mat)]))
      df$weights=1
      df$weights[df$target<0.5]=sum(df$target>0.5)/sum(df$target<0.5)
      library(glmnet)
      res <- glmnet(df[,c("de","exp","pc")], df$target, lambda = 0, lower.limits = 0, intercept = FALSE,weights=df$weights)
      coef_res=as.numeric(coef(res))
      names(coef_res)=row.names(coef(res))
      
      sim_mat=de_sim_mat*coef_res["de"]+exp_sim_mat*coef_res["exp"]+pc_sim_mat*coef_res["pc"]
      
      sim_mat=sim_mat=.extra_matrix_rowNorm(input_mat = sim_mat,rowValues = 1/diag(sim_mat))#Matrix::Diagonal(x=1/diag(sim_mat)) %*% sim_mat
      sim_mat[sim_mat>1]=1
      
      
      
      #sim_mat=.mynetAff(sim_mat,rev=F)
      #sim_mat2=.mynetAff(inputData = sim_mat,rev = F)
      
      #diag(sim_mat)=0
      #sim_mat=Matrix::Diagonal(x=1/qlcMatrix::rowMax(sim_mat)) %*% sim_mat
      
      
      
      z_cosine_binary=.mynetAff(inputData = z_cosine3,rev = F)
      z_cosine_binary[z_cosine_binary<0.00001]=0
      z_cosine_binary=Matrix::drop0(z_cosine_binary)
      z_cosine_binary@x=rep(1,length(z_cosine_binary))
      sim_mat=sim_mat*z_cosine_binary
      
      
      
      
      diff=cor(as.matrix(sim_mat),as.matrix(exp(-1*(de_pct_res)/de_quantile_thr)))
      diff=1- .extra_matrix_rowNorm(prop_mat=diff,rowValues = 1/qlcMatrix::rowMax(diff))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(diff)) %*% diff
      
      diff=pmax(diff,0)
      diff=(diff+t(diff))/2
      diag(diff)=0
      
    }  else if(clustering_version=="v56"){
      
      de2=.mynetAff(de_pct_res,rev=T)
      z_cosine_binary=de2
      z_cosine_binary=quantile(as.numeric(z_cosine_binary[z_cosine3<0]),0.95)
      de3=Matrix::drop0(de2,tol=z_cosine_binary)
      #diff=scaling_factor*(sim_mat+(1-exp(-1*(de_pct_res)/affinity_param))+(1-exp(-1*(de_pct_res)/20)))+(pmax(1-as.matrix(z_cosine3),0))+(1-cosine3)
      
      
      diff=(1-(de2+de3)/2)
      
      
      diff=pmax(diff,0)
      diff=diff+t(diff)
      diag(diff)=0
      
    } else if(clustering_version=="v23"){
      
      if(!file.exists(.myFilePathMakerFn("res_pseudocell_sim",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))){
        sim_mat=.extra_sconline.pseudosim_archive(argList=argList,binarize = nrow(prop_mat)<100000)
        qsave(sim_mat,file=.myFilePathMakerFn("res_pseudocell_sim",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
      } else {
        sim_mat=qread(.myFilePathMakerFn("res_pseudocell_sim",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
      }
      
      
      
      #diag(sim_mat)=0
      #sim_mat=Matrix::Diagonal(x=1/qlcMatrix::rowMax(sim_mat)) %*% sim_mat
      
      if(T){
        sim_mat=exp(-3*(1-sim_mat))-0.04979139
        sim_mat=sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      } else {
        sim_mat=.myaff(1-sim_mat)
        sim_mat=sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
      }
      
      
      
      z_cosine_binary=z_cosine3
      z_cosine_binary[z_cosine_binary<0]=0
      z_cosine_binary=Matrix::drop0(z_cosine_binary)
      z_cosine_binary@x=rep(1,length(z_cosine_binary))
      sim_mat=sim_mat*z_cosine_binary
      
      sim_mat=1-sim_mat
      diag(sim_mat)=0
      
      if(F){
        tst1=(1-exp(-1*(de_pct_res)/affinity_param))+(1-exp(-1*(de_pct_res)/20))
        tst1=tst1+t(tst1)
        tst2=(pmax(1-as.matrix(z_cosine3),0))
        tst3=sim_mat
        tst4=cor(t(as.matrix(meta_data$logFC)))
        tst4=1-tst4[row.names(tst1),colnames(tst1)]
        
        tst1+(exp((tst2))/exp(1.6))*2+tst3*4
      }
      
      cosine3=.myaff(1-z_cosine3,K=2)
      
      diff=scaling_factor*(sim_mat+(1-exp(-1*(de_pct_res)/affinity_param))+(1-exp(-1*(de_pct_res)/20)))+(pmax(1-as.matrix(z_cosine3),0))+(1-cosine3)
      diff=pmax(diff,0)
      diff=diff+t(diff)
    } else if(clustering_version=="random"){
      diff=diff+(1-exp(-1*(de_pct_res)/affinity_param))+pmax(1-as.matrix(z_cosine),0)^0.25
      diff=scaling_factor*diff
      diff2=matrix(sample(as.numeric(as.matrix(diff))),nrow=nrow(diff),ncol=ncol(diff))
      row.names(diff2)=row.names(diff)
      colnames(diff2)=colnames(diff)
      diff=diff2
    }
    
    
    
  }
  
  if(!clustering_version %in% c("v23","v27","v28","v29","v40","v")|F){
    diff=diff + t(diff)
    diff_clust=hclust(as.dist(diff),method = clustering_method)
    
  } else {
    diff_clust=hclust(as.dist(diff),method = clustering_method,members = .myEffSizePropMat(prop_mat)$effective_sample_size)
    
  }
  
  
  
  
  
    d_conMat=cutree(diff_clust,k=n_clusters)
    prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
    prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),row.names(pd)]
    prop_merged=.extra_matrix_rowNorm(prop_merged) #Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
    #prop anno
    
    colMax_vals_m=qlcMatrix::colMax(prop_merged)
    colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
    prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
    pd2=data.frame(cluster=paste0("C",as.character(prop_m_hardCluster$i)))
    row.names(pd2)=row.names(pd)
    pseudocell_cluster_assignments=data.frame(pseudocell=names(d_conMat),cluster=paste0("C",d_conMat),stringsAsFactors=F)
  
  return(list(cluster_object=diff_clust,distance_matrix=diff,n_clusters=n_clusters,cell_cluster_assignments=pd2,pseudocell_cluster_assignments=pseudocell_cluster_assignments))
}


.sconline.cluster=function(argList,sig1_thr=3,pct2_thr=0.3,input_meta_z=NULL,input_prop_mat=NULL,hierarchical_mode=F,clustering_version=NULL,pct_diff_thr=0.2,affinity_param=NULL,prune_clusters=F,n_clusters=NULL,inputPd=NULL,clustering_method="average",new_DEcount_run=F,tol_level=0.9,cell_count_diff_thr=0,combinatorial_pct_tol=1,inputExpData=NULL,forgiveness_factor=1,input_pca_centroids=NULL,input_embeddings=NULL,...){
  
  #clustering_method="average"
  #argList=.ArgList;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;affinity_param=5;inputPd=pd;n_clusters=5;new_DEcount_run=F;cell_count_diff_thr=0
  #input_meta_z=inputData$meta_data;input_prop_mat=inputData$prop_mat;hierarchical_mode=T
  
  #clustering_version="new";affinity_param=NULL;cell_count_diff_thr=0
  
  if(is.null(clustering_version)){
    clustering_version="new"
  }
  
  if(F){
    if(is.null(inputExpData)&prune_clusters){
      stop("inputExpData needs to be provided for cluster pruning")
    }
  }
  
  
  .myaff=function (diff, K = 3, sigma = 0.5) 
  {
    #K = 3; sigma = 0.5
    diff=as.matrix(diff)
    N <- nrow(diff)
    diff <- (diff + t(diff))/2
    diag(diff) <- 0
    sortedColumns <- as.matrix(t(apply(diff, 2, sort)))
    finiteMean <- function(x) {
      return(mean(x[is.finite(x)]))
    }
    means <- apply(sortedColumns[, 1:K + 1], 1, finiteMean) + 
      .Machine$double.eps
    avg <- function(x, y) {
      return((x + y)/2)
    }
    
    Sig <- outer(means, means, avg)/3 * 2 + diff/3 + .Machine$double.eps
    Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
    
    Sig=quantile(Sig,0.1)
    densities <- dnorm(as.matrix(diff), 0, sigma * Sig, log = FALSE)
    W <- (densities + t(densities))/2
    return(W)
  }
  
  .mynetAff=function(inputData,rev=T,do.scale=T){
    if(do.scale){
      aff=t(apply(log2(inputData+1),1,scale))
      aff[is.na(aff)]=0
      aff=(aff+t(aff))/sqrt(2)
      row.names(aff)=colnames(aff)=row.names(inputData)
    } else {
      aff=inputData
    }
    
    
    
    #aff=pnorm(aff,lower.tail = F)
    if(T){
      
      if(rev){
        diff=aff
        if(do.scale){
          diff[aff>0|t(aff)>0]=0
        }
        
        cor_map=diff
        cor_map=cor_map*(-1)
      } else {
        diff=aff
        if(do.scale){
          diff[aff<0|t(aff)<0]=0
        }
        
        cor_map=diff
      }
      
      net=NULL
      nn.rank=matrix(0,nrow=nrow(cor_map),ncol=5)
      row.names(nn.rank)=row.names(cor_map)
      for(i in 1:nrow(cor_map)){
        tmp=cor_map[i,]
        names(tmp)=colnames(cor_map)
        tmp=tmp[order(tmp,decreasing = T)]
        
        #sl_neighbors=names(tmp)[which(tmp>0.1)]
        sl_neighbors=names(tmp)[tmp>=(0.98*tmp[4])]
        if(length(sl_neighbors)>0){
          net=rbind(net,data.frame(source=row.names(cor_map)[i],target=sl_neighbors,score=cor_map[row.names(cor_map)[i],sl_neighbors],stringsAsFactors = F))
        }
        
      }
      net$name=paste0(net$source,"_",net$target)
      net$name[net$source>net$target]=paste0(net$target,"_",net$source)[net$source>net$target]
      net=net[!duplicated(net$name),]
      #net$score=1
      net=reshape2::dcast(source~target,data=net,value.var = "score")
      row.names(net)=net[,1]
      net=net[,-1]
      net_c_names=setdiff(row.names(aff),colnames(net))
      net_c=matrix(0,nrow=nrow(net),ncol=length(net_c_names))
      colnames(net_c)=net_c_names
      net=cbind(net,net_c)
      net=net[match(row.names(aff),row.names(net)),colnames(aff)]
      net[is.na(net)]=0
      row.names(net)=row.names(aff)
      net=net+t(net)
      net=as.matrix(net)/2
      #diag(net)=1
      net=as(net,"dgCMatrix")
      net=.extra_matrix_rowNorm(net)#Matrix::Diagonal(x=1/rowSums(net)) %*% net
      net2=net
      for(iii in 1:20){
        net2=(net2) %*% t(net)
        net2=.extra_matrix_rowNorm(net2)#Matrix::Diagonal(x=1/rowSums(net2)) %*% net2
        net2=net2+net
        net2=net2/2
      }
      net=net2
      if(F){
        diag(net)=0
      }
      
      net=.extra_matrix_rowNorm(input_mat = net,rowValues = 1/qlcMatrix::rowMax(net))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(net)) %*% net
      #diag(net)=1
      
    }
    return(net)
  }
  
  
  if(sum(clustering_method %in% c("complete", "average", "mcquitty"))==0){
    stop("Acceptable clustering methods: complete, average, mcquitty")
  }
  
  if(is.null(inputPd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPd
  }
  
  
  if(is.null(n_clusters)){
    stop("n_clusters argument should be provided!")
  }
  
  #res_array=.sconline.fetch_data("sconline_arrays",argList = argList)
  if(!is.null(input_meta_z)){
    meta_data=input_meta_z
  } else {
    meta_data=.sconline.fetch_data("meta_z",argList = argList)
  }
  
  if(!hierarchical_mode){
    de_prefix=paste0("de_sig",sig1_thr,"_pctdiff",pct_diff_thr,"_pct2",pct2_thr)
    if(cell_count_diff_thr>0){
      de_prefix=paste0(de_prefix,"_cellCount",cell_count_diff_thr)
    }
    
    #pct_mat=meta_data$med_pct.1;argList = argList;meta_z_mat=meta_data$meta_z;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
    if(!file.exists(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))|new_DEcount_run){
      
      de_pct_res=.extra_sconline.PctScoreFn(pct_mat=meta_data$med_pct.1,argList = argList,meta_z_mat=meta_data$meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F,cell_count_diff_thr=cell_count_diff_thr)
      
      #de_pct_res=.extra_sconline.PctScoreFn(pct_mat=.data$pct_mat,argList = .data$argList,meta_z_mat=.data$meta_z_mat,sig1_thr=.data$sig1_thr,centers=NULL,pct2_thr=.data$pct2_thr,pct_diff_thr=.data$pct_diff_thr,symmetric=F)
      
      qsave(de_pct_res,file=.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    } else {
      
      de_pct_res=qread(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    }
  } else {
    de_pct_res=.extra_sconline.PctScoreFn(pct_mat=meta_data$med_pct.1,argList = argList,meta_z_mat=meta_data$meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F,cell_count_diff_thr=cell_count_diff_thr,hierarchical_mode=hierarchical_mode)
    
  }
  
  if(is.null(input_prop_mat)){
    prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    prop_mat=input_prop_mat
  }
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(n_clusters>nrow(prop_mat)){
    warning(paste("Number of clusters can't be larger than the number of pseudocells! setting n_clusters to",nrow(prop_mat)))
    n_clusters=nrow(prop_mat)
  }
  
  tst_cluster=hclust(as.dist(de_pct_res))
  tst_cluster2=cutree(tst_cluster,h=0.01)
  if(length(unique(tst_cluster2))<=n_clusters){
    n_clusters=length(unique(tst_cluster2))-1
  }
  
  if(n_clusters==1&nrow(prop_mat)==2){
    if(de_pct_res[1,2]<5){
      n_clusters=1
    } else {
      n_clusters=2
    }
    de_pct_res=1-exp(-1*(de_pct_res)/5)
    diff_clust=hclust(as.dist(de_pct_res))
    if(n_clusters==1){
      cell_cluster_assignments=data.frame(cluster="C1",cell=colnames(prop_mat))
      row.names(cell_cluster_assignments)=colnames(prop_mat)
      pseudocell_cluster_assignments=data.frame(pseudocell=row.names(prop_mat),cluster="C1",stringsAsFactors = F)
    } else {
      prop_merged=prop_mat
      colMax_vals_m=qlcMatrix::colMax(prop_merged)
      colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
      prop_m_hardCluster=Matrix::drop0(colMax_vals_m,tol=tol_level)
      prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
      prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
      prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
      prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
      pd2=data.frame(cluster=paste0("C",as.character(prop_m_hardCluster$i)))
      row.names(pd2)=row.names(pd)
      pseudocell_cluster_assignments=data.frame(pseudocell=row.names(prop_mat),cluster=paste0("C",1:2),stringsAsFactors=F)
      
    }
    
    return(list(cluster_object=diff_clust,distance_matrix=de_pct_res,n_clusters=n_clusters,cell_cluster_assignments=cell_cluster_assignments,pseudocell_cluster_assignments=pseudocell_cluster_assignments))
  }
  
  
  if(F){
    ps_anno=as.matrix(.myOneHotFn(inputVector=pd$anno_orig_cellState))
    ps_anno=prop_mat %*% ps_anno
    ps_anno=lapply(1:nrow(ps_anno),function(x) {
      y=which(ps_anno[x,]==max(ps_anno[x,]))[1]
      y=colnames(ps_anno)[y]
      return(data.frame(ps=row.names(ps_anno)[x],anno=y,stringsAsFactors = F))
    })
    ps_anno=do.call("rbind",ps_anno)
    
  }
  
  {
    #affinity_param=min(affinity_param,quantile(as.numeric(de_pct_res),0.1))
    #affinity_param=max(affinity_param,2)
    if(!is.null(input_embeddings)&!is.null(input_prop_mat)){
      input_embeddings=input_embeddings[match(colnames(input_prop_mat),row.names(input_embeddings)),]
    }
    
    ps_sim_mat=.extra_sconline.pseudosim_archive11(argList=argList,input_prop_mat=input_prop_mat,hierarchical_mode=hierarchical_mode,binarize = ncol(prop_mat)<100000,cos_dist = T,input_pca_centroids=input_pca_centroids,input_embeddings=input_embeddings)
    
    ps_sim_mat=ps_sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
    
    
    
    if(!is.null(affinity_param)){
      sl_param=affinity_param
    } else {
      if(clustering_version=="old"){
        sl_param=NULL
        sl_score=0
        for(affinity_param in 2:20){
          de_sim_mat=exp(-1*(de_pct_res)/affinity_param)
          
          tst=cor(as.matrix(de_sim_mat),as.matrix(ps_sim_mat))
          
          diag(tst)=0
          tmp_score=median(apply(tst,1,function(x) max(x,na.rm = T)),na.rm = T)
          if(tmp_score>sl_score){
            sl_score=tmp_score
            sl_param=affinity_param
          }
        }
      } else if(clustering_version=="new") {
        
        sl_param=NULL
        sl_score=0
        for(affinity_param in 2:20){
          de_sim_mat=exp(-1*(de_pct_res)/affinity_param)
          
          tst=cor(as.matrix(de_sim_mat),as.matrix(ps_sim_mat))
          diag(tst)=0
          tmp_score=median(apply(tst,1,max),na.rm = T)
          if(tmp_score>sl_score){
            sl_score=tmp_score
            sl_param=affinity_param
          }
        }
        sl_param1=sl_param
        
        tst=de_pct_res
        diag(tst)=100
        sl_param=apply(as.matrix(tst),1,function(x) {#quantile(x,0.01);
          x=x[order(x,decreasing = F)]
          (x[2]+quantile(x,0.01))/2})
        sl_param=median(sl_param)
        sl_param=max(sl_param,2)
        if(!is.null(sl_param1)){
          sl_param=(sl_param+sl_param1)/2
        }
        
      }
      
      print(paste("Selected affinity_param:",sl_param))
    }
    
    #
    de_sim_mat=exp(-1*(de_pct_res)/sl_param)
    
    diag(ps_sim_mat)=0
    diag(de_sim_mat)=0
    ps_max_vals=qlcMatrix::rowMax(ps_sim_mat)
    de_max_vals=qlcMatrix::rowMax(de_sim_mat)
    ps_sim_mat=.extra_matrix_rowNorm(input_mat = ps_sim_mat,rowValues = 1/(ps_max_vals+0.0000001))#Matrix::Diagonal(x=1/(ps_max_vals+0.0000001)) %*% ps_sim_mat
    de_sim_mat=.extra_matrix_rowNorm(input_mat = de_sim_mat,rowValues = 1/(de_max_vals+0.0000001))#Matrix::Diagonal(x=1/(de_max_vals+0.0000001)) %*% de_sim_mat
    ps_sim_mat=sweep(as.matrix(ps_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
    de_sim_mat=sweep(as.matrix(de_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
    diag(ps_sim_mat)=1
    diag(de_sim_mat)=1
    diff=((1-ps_sim_mat)+(1-de_sim_mat))
    
    diff=pmax(diff,0)
    diff=diff+t(diff)
    diag(diff)=0
    
  }
  
  res_eff_size=.myEffSizePropMat(prop_mat)$effective_sample_size
  names(res_eff_size)=row.names(prop_mat)
  diff_clust=hclust(as.dist(diff),method = clustering_method,members = res_eff_size[colnames(diff)])
  
  
  d_conMat=cutree(diff_clust,k=max(n_clusters,1))
  
  if(prune_clusters){
    
    #.tst=list(cluster_assignments=d_conMat,clust_obj=diff_clust,
    #          inputExpData=inputExpData,
    #          argList=argList,
    #          combinatorial_pct_tol=combinatorial_pct_tol,marker_sig1_thr=sig1_thr,marker_pct2_thr=pct2_thr,marker_pct_diff_thr=pct_diff_thr,forgiveness_factor=forgiveness_factor,input_meta_data=input_meta_z,inputPhenoData=pd,input_prop_mat=prop_mat)
    
    if(!is.null(prop_mat)){
      pd=pd[match(colnames(prop_mat),row.names(pd)),]
    }
    d_conMat=.sconline.cluster_pruning(cluster_assignments=d_conMat,clust_obj=diff_clust,
                              inputExpData=inputExpData,
                              argList=argList,
                              combinatorial_pct_tol=combinatorial_pct_tol,marker_sig1_thr=sig1_thr,marker_pct2_thr=pct2_thr,marker_pct_diff_thr=pct_diff_thr,forgiveness_factor=forgiveness_factor,input_meta_data=input_meta_z,inputPhenoData=pd,input_prop_mat=prop_mat)
  }
  
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),row.names(pd)]
  prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  #prop anno
  
  colMax_vals_m=qlcMatrix::colMax(prop_merged)
  colMax_vals_m=.extra_matrix_colNorm(prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=Matrix::drop0(colMax_vals_m,tol=tol_level)
  prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
  prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
  prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
  prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
  pd2=data.frame(cluster=paste0("C",as.character(prop_m_hardCluster$i)))
  row.names(pd2)=row.names(pd)
  pseudocell_cluster_assignments=data.frame(pseudocell=names(d_conMat),cluster=paste0("C",d_conMat),stringsAsFactors=F)
  
  return(list(cluster_object=diff_clust,distance_matrix=diff,n_clusters=n_clusters,cell_cluster_assignments=pd2,pseudocell_cluster_assignments=pseudocell_cluster_assignments))
}

.sconline.cluster.spectral_clustering=function(argList,clustering_version=NULL,scaling_factor=10,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,affinity_param=5,n_clusters=NULL,inputPd=NULL,clustering_method="average",new_DEcount_run=F,tol_level=0.9,cell_count_diff_thr=25){
  
  #clustering_method="average"
  #argList=.ArgList;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;affinity_param=5;scaling_factor=10;inputPd=NULL;n_clusters=115;new_DEcount_run=F;cell_count_diff_thr=25
  #clustering_version="v22"
  
  .myaff=function (diff, K = 3, sigma = 0.5) 
  {
    #K = 3; sigma = 0.5
    diff=as.matrix(diff)
    N <- nrow(diff)
    diff <- (diff + t(diff))/2
    diag(diff) <- 0
    sortedColumns <- as.matrix(t(apply(diff, 2, sort)))
    finiteMean <- function(x) {
      return(mean(x[is.finite(x)]))
    }
    means <- apply(sortedColumns[, 1:K + 1], 1, finiteMean) + 
      .Machine$double.eps
    avg <- function(x, y) {
      return((x + y)/2)
    }
    
    Sig <- outer(means, means, avg)/3 * 2 + diff/3 + .Machine$double.eps
    Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
    
    Sig=quantile(Sig,0.1)
    densities <- dnorm(as.matrix(diff), 0, sigma * Sig, log = FALSE)
    W <- (densities + t(densities))/2
    return(W)
  }
  
  if(is.null(clustering_version)){
    clustering_version=30
  }
  
  
  clustering_method=tolower(clustering_method)
  if(sum(clustering_method %in% c("complete", "average", "mcquitty"))==0){
    stop("Acceptable clustering methods: complete, average, mcquitty")
  }
  
  if(is.null(inputPd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPd
  }
  
  
  if(is.null(n_clusters)){
    
    stop("n_clusters argument should be provided!")
  }
  
  #res_array=.sconline.fetch_data("sconline_arrays",argList = argList)
  meta_data=.sconline.fetch_data("meta_z",argList = argList)
  de_prefix=paste0("de_sig",sig1_thr,"_pctdiff",pct_diff_thr,"_pct2",pct2_thr)
  if(cell_count_diff_thr>0){
    de_prefix=paste0(de_prefix,"_cellCount",cell_count_diff_thr)
  }
  #pct_mat=meta_data$med_pct.1;argList = argList;meta_z_mat=meta_data$meta_z;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
  if(!file.exists(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))|new_DEcount_run){
    de_pct_res=.extra_sconline.PctScoreFn(pct_mat=meta_data$med_pct.1,argList = argList,meta_z_mat=meta_data$meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F,cell_count_diff_thr=cell_count_diff_thr)
    
    #de_pct_res=.extra_sconline.PctScoreFn(pct_mat=.data$pct_mat,argList = .data$argList,meta_z_mat=.data$meta_z_mat,sig1_thr=.data$sig1_thr,centers=NULL,pct2_thr=.data$pct2_thr,pct_diff_thr=.data$pct_diff_thr,symmetric=F)
    
    qsave(de_pct_res,file=.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    de_pct_res=qread(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  }
  
  z_cosine3=qlcMatrix::corSparse(t(meta_data$meta_z))
  z_cosine3[is.na(z_cosine3)]=0
  row.names(z_cosine3)=row.names(de_pct_res)
  colnames(z_cosine3)=colnames(de_pct_res)
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(n_clusters>nrow(prop_mat)){
    warning(paste("Number of clusters can't be larger than the number of pseudocells! setting n_clusters to",nrow(prop_mat)))
    n_clusters=nrow(prop_mat)
  }
  
    #tst=prop_mat
    #tst@x=rep(1,length(tst@x))
    sim_mat=.extra_sconline.pseudosim_archive(argList=argList,binarize = T)
    
    
    sim_mat=exp(-3*(1-sim_mat))-0.04979139
    sim_mat=sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
    
    
    
    
    z_cosine_binary=z_cosine3
    z_cosine_binary[z_cosine_binary<0]=0
    z_cosine_binary=Matrix::drop0(z_cosine_binary)
    z_cosine_binary@x=rep(1,length(z_cosine_binary))
    sim_mat=sim_mat*z_cosine_binary
    
    sim_mat=1-sim_mat
    diag(sim_mat)=0
    
    cosine3=.myaff(1-z_cosine3,K=2)
    
    de2=t(apply(log2(de_pct_res+1),1,scale))
    de2[is.na(de2)]=0
    de2=(de2+t(de2))/sqrt(2)
    row.names(de2)=colnames(de2)=row.names(de_pct_res)
    
    
    
    #de2=pnorm(de2,lower.tail = F)
    if(T){
      
      diff=de2
      diff[de2>0|t(de2)>0]=0
      cor_map=diff+t(diff)
      cor_map=cor_map*(-1)
      net=NULL
      nn.rank=matrix(0,nrow=nrow(cor_map),ncol=5)
      row.names(nn.rank)=row.names(cor_map)
      for(i in 1:nrow(cor_map)){
        tmp=cor_map[i,]
        names(tmp)=colnames(cor_map)
        tmp=tmp[order(tmp,decreasing = T)]
        
        #sl_neighbors=names(tmp)[which(tmp>0.1)]
        sl_neighbors=names(tmp)[tmp>=(0.98*tmp[4])]
        if(length(sl_neighbors)>0){
          net=rbind(net,data.frame(source=row.names(cor_map)[i],target=sl_neighbors,score=cor_map[row.names(cor_map)[i],sl_neighbors],stringsAsFactors = F))
        }
        
      }
      net$name=paste0(net$source,"_",net$target)
      net$name[net$source>net$target]=paste0(net$target,"_",net$source)[net$source>net$target]
      net=net[!duplicated(net$name),]
      #net$score=1
      net=reshape2::dcast(source~target,data=net,value.var = "score")
      row.names(net)=net[,1]
      net=net[,-1]
      net_c_names=setdiff(row.names(de2),colnames(net))
      net_c=matrix(0,nrow=nrow(net),ncol=length(net_c_names))
      colnames(net_c)=net_c_names
      net=cbind(net,net_c)
      net=net[match(row.names(de2),row.names(net)),colnames(de2)]
      net[is.na(net)]=0
      row.names(net)=row.names(de2)
      net=net+t(net)
      net=as.matrix(net)/2
      #diag(net)=1
      net=as(net,"dgCMatrix")
      net=.extra_matrix_rowNorm(net)#Matrix::Diagonal(x=1/rowSums(net)) %*% net
      if(T){
        net=(net) %*% t(net)
      }
      
      if(F){
        diag(net)=0
      }
      
      net=.extra_matrix_rowNorm(input_mat = net,rowValues = 1/qlcMatrix::rowMax(net))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(net)) %*% net
      #diag(net)=1
      de2=1-net
      
    }
    #diff=scaling_factor*(sim_mat+(1-exp(-1*(de_pct_res)/affinity_param))+(1-exp(-1*(de_pct_res)/20)))+(pmax(1-as.matrix(z_cosine3),0))+(1-cosine3)
    
    
    
    de2=cor(as.matrix(de2),as.matrix(sim_mat))
    de2[is.na(de2)]=(-1)
    #diff=(de2*(exp(-1*(de_pct_res)/20))+(exp(-1*(de_pct_res)/affinity_param)))+((pmin(as.matrix(z_cosine3),1))+(cosine3))/scaling_factor
    if(clustering_version=="v30"){
      diff=(de2)+(exp(-1*(de_pct_res)/affinity_param))
    } else if(clustering_version=="v31"){
      diff=(de2)
    }
    #+((pmin(as.matrix(z_cosine3),1))+(cosine3))/scaling_factor
    
    #diff=diff/10
    #diff=(de2)+(exp(-1*(de_pct_res)/affinity_param))/2
    
    LW <- length(diff)
    normalize <- function(X) {
      row.sum.mdiag <- rowSums(X) - diag(X)
      row.sum.mdiag[row.sum.mdiag == 0] <- 1
      X <- X/(2 * (row.sum.mdiag))
      diag(X) <- 0.5
      return(X)
    }
    
    W <- diff/LW
    W[W<0]=0
    W <- normalize(W)
    W <- (W + t(W))/2
    
  d_conMat=SNFtool::spectralClustering(W, K=n_clusters, type = 3)
  names(d_conMat)=row.names(diff)
  
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),row.names(pd)]
  prop_merged=.extra_matrix_rowNorm(prop_merged)# Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  #prop anno
  
  colMax_vals_m=qlcMatrix::colMax(prop_merged)
  colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
  prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
  prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
  prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
  prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
  pd2=data.frame(cluster=paste0("C",as.character(prop_m_hardCluster$i)))
  row.names(pd2)=row.names(pd)
  pseudocell_cluster_assignments=data.frame(pseudocell=names(d_conMat),cluster=paste0("C",d_conMat),stringsAsFactors=F)
  
  return(list(affinity=W,n_clusters=n_clusters,cell_cluster_assignments=pd2,pseudocell_cluster_assignments=pseudocell_cluster_assignments))
}

.sconline.pseudocellMapfn.DEmap.count=function(argList,cluster_obj=NULL,attribute_col=NULL,include_directional_changes=T,n_clusters=NULL,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,DE_quantile=0.05,DE_count_thr=NULL,new_DEcount_run=F,label_nodes=F,cell_count_diff_thr=25){
  require(hues)
  #argList=.ArgList;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;inputPd=NULL;DE_quantile=0.05;new_DEcount_run=F;label_nodes=F
  #cluster_obj=clust1;attribute_col=NULL;n_clusters=40;DE_count_thr=NULL
  
  
  input_umap_centroid=NULL
  
  if(DE_quantile>1|DE_quantile<0){
    stop("DE_quantile is expected to be between 0 and 1")
  } else if(DE_quantile>0.5){
    warning("set DE_quantile value is too high!")
  }
  
  if(!is.null(cluster_obj)){
    clust_data=.sconline.cluster.Vis(argList=argList,cluster_obj=cluster_obj,attribute_col=attribute_col,n_clusters=n_clusters,return_obj=T)
    pd=clust_data$cell
    attribute_col="cluster_anno_res"
    input_umap_centroid=clust_data$pseudocell
    input_umap_centroid$centroid=input_umap_centroid$pseudocell
  } else {
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  }
  
  
  
  #res_array=.sconline.fetch_data("sconline_arrays",argList = argList)
  meta_data=.sconline.fetch_data("meta_z",argList = argList)
  de_prefix=paste0("de_sig",sig1_thr,"_pctdiff",pct_diff_thr,"_pct2",pct2_thr)
  if(cell_count_diff_thr>0){
    de_prefix=paste0(de_prefix,"_cellCount",cell_count_diff_thr)
  }
  #pct_mat=meta_data$med_pct.1;argList = argList;meta_z_mat=meta_data$meta_z;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
  if(!file.exists(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))|new_DEcount_run){
    de_pct_res=.extra_sconline.PctScoreFn(pct_mat=meta_data$med_pct.1,argList = argList,meta_z_mat=meta_data$meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F)
    
    #de_pct_res=.extra_sconline.PctScoreFn(pct_mat=.data$pct_mat,argList = .data$argList,meta_z_mat=.data$meta_z_mat,sig1_thr=.data$sig1_thr,centers=NULL,pct2_thr=.data$pct2_thr,pct_diff_thr=.data$pct_diff_thr,symmetric=F)
    
    qsave(de_pct_res,file=.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    de_pct_res=qread(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  }
  
  if(is.null(input_umap_centroid)){
    input_umap_centroid=.sconline.fetch_data("umap_pseudocells",argList)
  }
  
  input_umap_centroid=input_umap_centroid[input_umap_centroid$centroid %in% row.names(meta_data$meta_z),]
  
  sl_net=(de_pct_res+t(de_pct_res))
  DE_tol=de_pct_res[upper.tri(de_pct_res,diag = F)]
  if(!is.null(DE_count_thr)){
    message("DE_count_thr is provided; ignoring DE_quantile value")
    DE_tol=DE_count_thr
  } else {
    DE_tol=quantile(DE_tol,DE_quantile)
    message(paste("set DE count threshold:", DE_tol))
  }
  
  sl_net=which(sl_net<=DE_tol&upper.tri(sl_net),arr.ind = T)
  net=data.frame(Fnode=row.names(de_pct_res)[sl_net[,1]],Snode=row.names(de_pct_res)[sl_net[,2]],stringsAsFactors = F)#,weight=(1-de_pct_res[sl_net]/(DE_tol+1))^2
  net$directed=F
  if(include_directional_changes){
    myL2normFn=function(inputMat){
      prop_mat2=rowSums(inputMat^2)
      prop_mat2=sqrt(prop_mat2)
      res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
      return(res)
    }
    
    prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
    
    
    
    #tst=prop_mat
    #tst@x=rep(1,length(tst@x))
    tst=myL2normFn(inputMat = prop_mat)
    diff=t(tst)
    diff=tst %*% diff
    diff=1-as.matrix(diff)
    #diff@x=2*(exp(diff@x/max(quantile(diff@x,0.95),0.1))/(1+exp(diff@x/max(quantile(diff@x,0.95),0.1)))-0.5)
    
    dir_net=which(de_pct_res<=DE_tol&t(de_pct_res)>DE_tol,arr.ind = T)
    dir_net=data.frame(Fnode=row.names(de_pct_res)[dir_net[,"row"]],Snode=row.names(de_pct_res)[dir_net[,"col"]],stringsAsFactors = F)#,weight=(1-de_pct_res[sl_net]/(DE_tol+1))^2
    dir_net$DEcount=de_pct_res[cbind(dir_net[,"Snode"],dir_net[,"Fnode"])]
    dir_net$diff=as.matrix(diff)[cbind(dir_net[,"Snode"],dir_net[,"Fnode"])]
    dir_net=dir_net[order(dir_net$diff,dir_net$DEcount,decreasing = F),]
    dir_net=dir_net[!duplicated(dir_net$Fnode),]
    dir_net$directed=T
    
    net=rbind.fill(net,dir_net)
  }
  
  res_net=net
  
  
  #net=res_net[order(res_net$weight,decreasing = T),]
  
  net = network::network(net, directed = F,matrix.type="edgelist")
  
  netVerNames=network::network.vertex.names(net)
  #network::set.edge.attribute(net, "weight", res_net$weight)
  
  centroid_layout=input_umap_centroid[match(as.character(netVerNames),as.character(input_umap_centroid$centroid)),]
  if(sum(colnames(input_umap_centroid)=="cluster")>0){
    network::set.vertex.attribute(net,"cluster",centroid_layout$cluster)
  }
  centroid_layout=centroid_layout[,c("UMAP_1","UMAP_2")]
  centroid_layout=as.matrix(centroid_layout)
  colnames(centroid_layout)=c("x","y")
  
  net=ggnetwork:::fortify.network(net,layout = centroid_layout)
  
  
  
  pd_summary=pd
  
  scale_factor=net[!is.na(net$directed),]
  scale_factor=scale_factor[!duplicated(scale_factor$vertex.names),]
  scale_factor=merge(scale_factor,input_umap_centroid,by.x="vertex.names",by.y="centroid")
  scale_factor1=lm(x~UMAP_1,data=scale_factor)
  pd_summary$UMAP_1=predict(scale_factor1,newdata=pd_summary)
  input_umap_centroid$UMAP_1=predict(scale_factor1,newdata=input_umap_centroid)
  scale_factor2=lm(y~UMAP_2,data=scale_factor)
  pd_summary$UMAP_2=predict(scale_factor2,newdata=pd_summary)
  input_umap_centroid$UMAP_2=predict(scale_factor2,newdata=input_umap_centroid)
  
  pd_summary$xend=pd_summary$UMAP_1
  pd_summary$yend=pd_summary$UMAP_2
  pd_summary$vertex.names=""
  if(!is.null(attribute_col)){
    pd_summary$color=pd_summary[,attribute_col]
  } else {
    pd_summary$color="gray"
  }
  
  input_umap_centroid$xend=input_umap_centroid$UMAP_1
  input_umap_centroid$yend=input_umap_centroid$UMAP_2
  input_umap_centroid$vertex.names=""
  if(sum(colnames(input_umap_centroid)=="cluster")>0){
    colnames(input_umap_centroid)[which(colnames(input_umap_centroid)=="cluster")]=attribute_col
  }
  input_umap_centroid$color="gray"
  
  tmp=input_umap_centroid#[!input_umap_centroid$centroid %in% net$vertex.names,]
  
  
  #predicting the background color
  library(ggnetwork)
  #centroids_ideal=c("173","192","191","187","127","194")
  
  if(is.null(attribute_col)){
    p=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend))+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.01)+
      geom_point(data=tmp,aes(UMAP_1,UMAP_2),color="black",fill="red",shape=21)+
      geom_edges(data=net[(net$directed==F),], color = "#0B1D87")
    if(include_directional_changes){
      p=p+geom_edges(data=net[which(net$directed==T),],color = "#4F0C0C",arrow = arrow(length = unit(3, "pt"),type = "closed"))
      }
    #,aes(size=weight)
      
    p=p+theme_blank()+scale_size_continuous(range = c(0.06,1))+theme(legend.position = "none")
    p=p+scale_color_identity()+scale_fill_identity()
  } else {
    p=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend))+geom_point(data=pd_summary,aes_string('UMAP_1','UMAP_2',color=attribute_col),size=0.01)+
      geom_point(data=tmp,aes_string('UMAP_1','UMAP_2',fill=attribute_col),color="black",shape=21)+
      geom_edges(data=net[(net$directed==F),], color = "#0B1D87")
    if(include_directional_changes){
      p=p+geom_edges(data=net[which(net$directed==T),],color = "#4F0C0C",arrow = arrow(length = unit(3, "pt"),type = "closed"))
    }
    p=p+theme_blank()+scale_size_continuous(range = c(0.06,1))+theme(legend.position = "none")
    p=p+scale_fill_manual(values=hues::iwanthue(length(unique(tmp[,attribute_col]))))+scale_color_manual(values=hues::iwanthue(length(unique(pd_summary[,attribute_col]))))
  }
  if(label_nodes){
    p=p+geom_label(aes(label=vertex.names))
  }
  
  
  return(p)
}

.sconline.pseudocellMapfn.integrated=function(argList,cluster_obj=NULL,attribute_col=NULL,include_directional_changes=T,n_clusters=NULL,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,new_DEcount_run=F,label_nodes=F,cell_count_diff_thr=25,scaling_factor=10,affinity_param=5){
  require(hues)
  #argList=.ArgList;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;inputPd=NULL;new_DEcount_run=F;label_nodes=F
  #cluster_obj=NULL;attribute_col="anno_cellState";n_clusters=40;scaling_factor=10;affinity_param=5
  
  
  input_umap_centroid=NULL
  
  if(!is.null(cluster_obj)){
    clust_data=.sconline.cluster.Vis(argList=argList,cluster_obj=cluster_obj,attribute_col=attribute_col,n_clusters=n_clusters,return_obj=T)
    pd=clust_data$cell
    attribute_col="cluster_anno_res"
    input_umap_centroid=clust_data$pseudocell
    input_umap_centroid$centroid=input_umap_centroid$pseudocell
  } else {
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  }
  
  
  
  #res_array=.sconline.fetch_data("sconline_arrays",argList = argList)
  meta_data=.sconline.fetch_data("meta_z",argList = argList)
  de_prefix=paste0("de_sig",sig1_thr,"_pctdiff",pct_diff_thr,"_pct2",pct2_thr)
  if(cell_count_diff_thr>0){
    de_prefix=paste0(de_prefix,"_cellCount",cell_count_diff_thr)
  }
  #pct_mat=meta_data$med_pct.1;argList = argList;meta_z_mat=meta_data$meta_z;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
  if(!file.exists(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))|new_DEcount_run){
    de_pct_res=.extra_sconline.PctScoreFn(pct_mat=meta_data$med_pct.1,argList = argList,meta_z_mat=meta_data$meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F)
    
    #de_pct_res=.extra_sconline.PctScoreFn(pct_mat=.data$pct_mat,argList = .data$argList,meta_z_mat=.data$meta_z_mat,sig1_thr=.data$sig1_thr,centers=NULL,pct2_thr=.data$pct2_thr,pct_diff_thr=.data$pct_diff_thr,symmetric=F)
    
    qsave(de_pct_res,file=.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    de_pct_res=qread(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  }
  
  if(is.null(input_umap_centroid)){
    input_umap_centroid=.sconline.fetch_data("umap_pseudocells",argList)
  }
  
  .myaff=function (diff, K = 3, sigma = 0.5) {
    #K = 3; sigma = 0.5
    diff=as.matrix(diff)
    N <- nrow(diff)
    diff <- (diff + t(diff))/2
    diag(diff) <- 0
    sortedColumns <- as.matrix(t(apply(diff, 2, sort)))
    finiteMean <- function(x) {
      return(mean(x[is.finite(x)]))
    }
    means <- apply(sortedColumns[, 1:K + 1], 1, finiteMean) + 
      .Machine$double.eps
    avg <- function(x, y) {
      return((x + y)/2)
    }
    
    Sig <- outer(means, means, avg)/3 * 2 + diff/3 + .Machine$double.eps
    Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
    
    Sig=quantile(Sig,0.1)
    densities <- dnorm(as.matrix(diff), 0, sigma * Sig, log = FALSE)
    W <- (densities + t(densities))/2
    return(W)
  }
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(!file.exists(.myFilePathMakerFn("res_pseudocell_sim",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))){
    sim_mat=.extra_sconline.pseudosim_archive(argList=argList,binarize = nrow(prop_mat)<100000)
    qsave(sim_mat,file=.myFilePathMakerFn("res_pseudocell_sim",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
  } else {
    sim_mat=qread(.myFilePathMakerFn("res_pseudocell_sim",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
  }
  
  
  sim_mat=exp(-3*(1-sim_mat))-0.04979139
  sim_mat=sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
  
  z_cosine3=qlcMatrix::corSparse(t(meta_data$meta_z))
  z_cosine3[is.na(z_cosine3)]=0
  row.names(z_cosine3)=row.names(de_pct_res)
  colnames(z_cosine3)=colnames(de_pct_res)
  
  z_cosine_binary=z_cosine3
  z_cosine_binary[z_cosine_binary<0]=0
  z_cosine_binary=Matrix::drop0(z_cosine_binary)
  z_cosine_binary@x=rep(1,length(z_cosine_binary))
  sim_mat=sim_mat*z_cosine_binary
  
  sim_mat=1-sim_mat
  diag(sim_mat)=0
  
  cosine3=.myaff(1-z_cosine3,K=2)
  
  diff=scaling_factor*(sim_mat+(1-exp(-1*(de_pct_res)/affinity_param))+(1-exp(-1*(de_pct_res)/20)))+(pmax(1-as.matrix(z_cosine3),0))+(1-cosine3)
  diff=pmax(diff,0)
  diff=diff+t(diff)
  
  
  for(i in 1:nrow(diff)){
    tmp=diff[i,]
    names(tmp)=colnames(diff)
    tmp_thr=tmp[order(tmp,decreasing = F)]
    tmp_thr=tmp_thr[3]
    #sl_neighbors=names(tmp)[1:min(3,length(tmp)/2)]
    diff[i,]=exp(-1*(tmp/min(tmp_thr,quantile(as.numeric(diff),0.01)))^2)
  }
  
  
  net=NULL
  for(i in 1:nrow(diff)){
    tmp=diff[i,]
    names(tmp)=colnames(diff)
    tmp=tmp[order(tmp,decreasing = F)]
    tmp=tmp[-1]
    #sl_neighbors=names(tmp)[1:min(3,length(tmp)/2)]
    tmp_score=exp(-1*(tmp/min(tmp[2],quantile(as.numeric(diff),0.01)))^2)
    sl_neighbors=names(tmp)[tmp_score>0.2]
    
    if(length(sl_neighbors)>0){
      net=rbind(net,data.frame(source=row.names(diff)[i],target=sl_neighbors,distance=tmp[tmp_score>0.2],stringsAsFactors = F))
    }
    
  }
  net$name=paste0(net$source,"_",net$target)
  net$name[net$source>net$target]=paste0(net$target,"_",net$source)[net$source>net$target]
  net=net[!duplicated(net$name),]
  #net$score=exp(-1*net$distance/quantile(net$distance,0.25))
  p=.extra_sconline.NetVisFn(net,argList=argList,directional=F)
  ggsave(plot=p,file="~/myBucket/torm.pdf",width=40,height=40)
  
  res_net=net
  
  
  #net=res_net[order(res_net$weight,decreasing = T),]
  
  net = network::network(net, directed = F,matrix.type="edgelist")
  
  netVerNames=network::network.vertex.names(net)
  #network::set.edge.attribute(net, "weight", res_net$weight)
  
  centroid_layout=input_umap_centroid[match(as.character(netVerNames),as.character(input_umap_centroid$centroid)),]
  if(sum(colnames(input_umap_centroid)=="cluster")>0){
    network::set.vertex.attribute(net,"cluster",centroid_layout$cluster)
  }
  centroid_layout=centroid_layout[,c("UMAP_1","UMAP_2")]
  centroid_layout=as.matrix(centroid_layout)
  colnames(centroid_layout)=c("x","y")
  
  net=ggnetwork:::fortify.network(net,layout = centroid_layout)
  
  
  
  pd_summary=pd
  
  scale_factor=net[!is.na(net$directed),]
  scale_factor=scale_factor[!duplicated(scale_factor$vertex.names),]
  scale_factor=merge(scale_factor,input_umap_centroid,by.x="vertex.names",by.y="centroid")
  scale_factor1=lm(x~UMAP_1,data=scale_factor)
  pd_summary$UMAP_1=predict(scale_factor1,newdata=pd_summary)
  input_umap_centroid$UMAP_1=predict(scale_factor1,newdata=input_umap_centroid)
  scale_factor2=lm(y~UMAP_2,data=scale_factor)
  pd_summary$UMAP_2=predict(scale_factor2,newdata=pd_summary)
  input_umap_centroid$UMAP_2=predict(scale_factor2,newdata=input_umap_centroid)
  
  pd_summary$xend=pd_summary$UMAP_1
  pd_summary$yend=pd_summary$UMAP_2
  pd_summary$vertex.names=""
  if(!is.null(attribute_col)){
    pd_summary$color=pd_summary[,attribute_col]
  } else {
    pd_summary$color="gray"
  }
  
  input_umap_centroid$xend=input_umap_centroid$UMAP_1
  input_umap_centroid$yend=input_umap_centroid$UMAP_2
  input_umap_centroid$vertex.names=""
  if(sum(colnames(input_umap_centroid)=="cluster")>0){
    colnames(input_umap_centroid)[which(colnames(input_umap_centroid)=="cluster")]=attribute_col
  }
  input_umap_centroid$color="gray"
  
  tmp=input_umap_centroid#[!input_umap_centroid$centroid %in% net$vertex.names,]
  
  
  #predicting the background color
  library(ggnetwork)
  #centroids_ideal=c("173","192","191","187","127","194")
  
  if(is.null(attribute_col)){
    p=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend))+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.01)+
      geom_point(data=tmp,aes(UMAP_1,UMAP_2),color="black",fill="red",shape=21)+
      geom_edges(data=net[(net$directed==F),], color = "#0B1D87")
    if(include_directional_changes){
      p=p+geom_edges(data=net[which(net$directed==T),],color = "#4F0C0C",arrow = arrow(length = unit(3, "pt"),type = "closed"))
    }
    #,aes(size=weight)
    
    p=p+theme_blank()+scale_size_continuous(range = c(0.06,1))+theme(legend.position = "none")
    p=p+scale_color_identity()+scale_fill_identity()
  } else {
    p=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend))+geom_point(data=pd_summary,aes_string('UMAP_1','UMAP_2',color=attribute_col),size=0.01)+
      geom_point(data=tmp,aes_string('UMAP_1','UMAP_2',fill=attribute_col),color="black",shape=21)+
      geom_edges(data=net[(net$directed==F),], color = "#0B1D87")
    if(include_directional_changes){
      p=p+geom_edges(data=net[which(net$directed==T),],color = "#4F0C0C",arrow = arrow(length = unit(3, "pt"),type = "closed"))
    }
    p=p+theme_blank()+scale_size_continuous(range = c(0.06,1))+theme(legend.position = "none")
    p=p+scale_fill_manual(values=hues::iwanthue(length(unique(tmp[,attribute_col]))))+scale_color_manual(values=hues::iwanthue(length(unique(pd_summary[,attribute_col]))))
  }
  if(label_nodes){
    p=p+geom_label(aes(label=vertex.names))
  }
  
  
  return(p)
}

.sconline.pseudocellMapfn=function(argList,pseudo_sim_mode="integrated",label_nodes=F,cluster_obj=NULL,attribute_col=NULL,include_directional_changes=T,n_clusters=NULL,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,new_DEcount_run=F,cell_count_diff_thr=0,affinity_param=NULL,prune_clusters=F,combinatorial_pct_tol=1,inputExpData=NULL){
  require(hues)
  #argList=.ArgList;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;inputPd=NULL;new_DEcount_run=F;label_nodes=F
  #cluster_obj=cluster_res_avg;attribute_col=NULL;n_clusters=34;affinity_param=NULL;pseudo_sim_mode="integrated";cell_count_diff_thr=0
  
  match.arg(pseudo_sim_mode,c("integrated","DE","PC"))
  
  input_umap_centroid=NULL
  
  if(!is.null(cluster_obj)){
    if(is.null(inputExpData)&prune_clusters){
      stop("inputExpData needs to be provided for cluster pruning")
    }
    clust_data=.sconline.cluster.Vis(argList=argList,cluster_obj=cluster_obj,prune_clusters = prune_clusters,attribute_col=attribute_col,n_clusters=n_clusters,return_obj=T,combinatorial_pct_tol=combinatorial_pct_tol,inputExpData=inputExpData)
    pd=clust_data$cell
    attribute_col="cluster_anno_res"
    input_umap_centroid=clust_data$pseudocell
    input_umap_centroid$centroid=input_umap_centroid$pseudocell
  } else {
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  }
  
  
  
 
  
  
  if(is.null(input_umap_centroid)){
    input_umap_centroid=.sconline.fetch_data("umap_pseudocells",argList)
  }
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  ps_sim_mat=.extra_sconline.pseudosim_archive11(argList=argList,binarize = ncol(prop_mat)<100000,cos_dist = T)
  
  meta_data=.sconline.fetch_data("meta_z",argList = argList)
  de_prefix=paste0("de_sig",sig1_thr,"_pctdiff",pct_diff_thr,"_pct2",pct2_thr)
  if(cell_count_diff_thr>0){
    de_prefix=paste0(de_prefix,"_cellCount",cell_count_diff_thr)
  }
  #pct_mat=meta_data$med_pct.1;argList = argList;meta_z_mat=meta_data$meta_z;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
  if(!file.exists(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))|new_DEcount_run){
    de_pct_res=.extra_sconline.PctScoreFn(pct_mat=meta_data$med_pct.1,argList = argList,meta_z_mat=meta_data$meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F)
    
    #de_pct_res=.extra_sconline.PctScoreFn(pct_mat=.data$pct_mat,argList = .data$argList,meta_z_mat=.data$meta_z_mat,sig1_thr=.data$sig1_thr,centers=NULL,pct2_thr=.data$pct2_thr,pct_diff_thr=.data$pct_diff_thr,symmetric=F)
    
    qsave(de_pct_res,file=.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    de_pct_res=qread(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  }
  ps_sim_mat=ps_sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
  
  if(!is.null(affinity_param)){
    sl_param=affinity_param
  } else {
    
    sl_param=NULL
    sl_score=0
    for(affinity_param in 2:20){
      de_sim_mat=exp(-1*(de_pct_res)/affinity_param)
      
      tst=cor(as.matrix(de_sim_mat),as.matrix(ps_sim_mat))
      diag(tst)=0
      tmp_score=median(apply(tst,1,max))
      if(tmp_score>sl_score){
        sl_score=tmp_score
        sl_param=affinity_param
      }
    }
    sl_param1=sl_param
    
    tst=de_pct_res
    diag(tst)=100
    sl_param=apply(as.matrix(tst),1,function(x) {#quantile(x,0.01);
      x=x[order(x,decreasing = F)]
      (x[2]+quantile(x,0.01))/2})
    sl_param=median(sl_param)
    sl_param=max(sl_param,2)
    sl_param=(sl_param+sl_param1)/2
  }
  
  
  if(pseudo_sim_mode %in% c("integrated","DE")){
    
      print(paste("Selected affinity_param:",sl_param))
      
      if(pseudo_sim_mode=="integrated"){
        diag(ps_sim_mat)=0
        diag(de_sim_mat)=0
        ps_max_vals=qlcMatrix::rowMax(ps_sim_mat)
        de_max_vals=qlcMatrix::rowMax(de_sim_mat)
        ps_sim_mat=.extra_matrix_rowNorm(input_mat = ps_sim_mat,rowValues = 1/(ps_max_vals+0.0000001))#Matrix::Diagonal(x=1/(ps_max_vals+0.0000001)) %*% ps_sim_mat
        de_sim_mat=.extra_matrix_rowNorm(input_mat = de_sim_mat,rowValues = 1/(de_max_vals+0.0000001))#Matrix::Diagonal(x=1/(de_max_vals+0.0000001)) %*% de_sim_mat
        ps_sim_mat=sweep(as.matrix(ps_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
        de_sim_mat=sweep(as.matrix(de_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
        diag(ps_sim_mat)=1
        diag(de_sim_mat)=1
        diff=((1-ps_sim_mat)+(1-de_sim_mat))
        
        diff=pmax(diff,0)
        diff=diff+t(diff)
        diag(diff)=0
      } else {
        diag(ps_sim_mat)=0
        diag(de_sim_mat)=0
        ps_max_vals=qlcMatrix::rowMax(ps_sim_mat)
        de_max_vals=qlcMatrix::rowMax(de_sim_mat)
        ps_sim_mat=.extra_matrix_rowNorm(input_mat = ps_sim_mat,rowValues = 1/(ps_max_vals+0.0000001))#Matrix::Diagonal(x=1/(ps_max_vals+0.0000001)) %*% ps_sim_mat
        de_sim_mat=.extra_matrix_rowNorm(input_mat = de_sim_mat,rowValues = 1/(de_max_vals+0.0000001))#Matrix::Diagonal(x=1/(de_max_vals+0.0000001)) %*% de_sim_mat
        ps_sim_mat=sweep(as.matrix(ps_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
        de_sim_mat=sweep(as.matrix(de_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
        diag(ps_sim_mat)=1
        diag(de_sim_mat)=1
        diff=(1-de_sim_mat)
        
        diff=pmax(diff,0)
        diff=diff+t(diff)
        diag(diff)=0
      }
      
  } else {
    diag(ps_sim_mat)=0
    diag(de_sim_mat)=0
    ps_max_vals=qlcMatrix::rowMax(ps_sim_mat)
    de_max_vals=qlcMatrix::rowMax(de_sim_mat)
    ps_sim_mat=.extra_matrix_rowNorm(input_mat = ps_sim_mat,rowValues = 1/(ps_max_vals+0.0000001))#Matrix::Diagonal(x=1/(ps_max_vals+0.0000001)) %*% ps_sim_mat
    de_sim_mat=.extra_matrix_rowNorm(input_mat = de_sim_mat,rowValues = 1/(de_max_vals+0.0000001))#Matrix::Diagonal(x=1/(de_max_vals+0.0000001)) %*% de_sim_mat
    
    ps_sim_mat=sweep(as.matrix(ps_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
    de_sim_mat=sweep(as.matrix(de_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
    diag(ps_sim_mat)=1
    diag(de_sim_mat)=1
    diff=(1-ps_sim_mat)
    
    diff=pmax(diff,0)
    diff=diff+t(diff)
    diag(diff)=0
  }
  
  
  .mynetAff=function(inputData,rev=T,do.scale=T){
    if(do.scale){
      aff=t(apply(log2(inputData+1),1,scale))
      aff[is.na(aff)]=0
      aff=(aff+t(aff))/sqrt(2)
      row.names(aff)=colnames(aff)=row.names(inputData)
    } else {
      aff=inputData
    }
    
    
    
    #aff=pnorm(aff,lower.tail = F)
    if(T){
      
      if(rev){
        diff=aff
        if(do.scale){
          diff[aff>0|t(aff)>0]=0
        }
        
        cor_map=diff
        cor_map=cor_map*(-1)
      } else {
        diff=aff
        if(do.scale){
          diff[aff<0|t(aff)<0]=0
        }
        
        cor_map=diff
      }
      
      net=NULL
      nn.rank=matrix(0,nrow=nrow(cor_map),ncol=5)
      row.names(nn.rank)=row.names(cor_map)
      for(i in 1:nrow(cor_map)){
        tmp=cor_map[i,]
        names(tmp)=colnames(cor_map)
        tmp=tmp[order(tmp,decreasing = T)]
        
        #sl_neighbors=names(tmp)[which(tmp>0.1)]
        sl_neighbors=names(tmp)[tmp>=(0.98*tmp[4])]
        if(length(sl_neighbors)>0){
          net=rbind(net,data.frame(source=row.names(cor_map)[i],target=sl_neighbors,score=cor_map[row.names(cor_map)[i],sl_neighbors],stringsAsFactors = F))
        }
        
      }
      net$name=paste0(net$source,"_",net$target)
      net$name[net$source>net$target]=paste0(net$target,"_",net$source)[net$source>net$target]
      net=net[!duplicated(net$name),]
      #net$score=1
      net=reshape2::dcast(source~target,data=net,value.var = "score")
      row.names(net)=net[,1]
      net=net[,-1]
      net_c_names=setdiff(row.names(aff),colnames(net))
      net_c=matrix(0,nrow=nrow(net),ncol=length(net_c_names))
      colnames(net_c)=net_c_names
      net=cbind(net,net_c)
      net=net[match(row.names(aff),row.names(net)),colnames(aff)]
      net[is.na(net)]=0
      row.names(net)=row.names(aff)
      net=net+t(net)
      net=as.matrix(net)/2
      #diag(net)=1
      net=as(net,"dgCMatrix")
      net=.extra_matrix_rowNorm(net)#Matrix::Diagonal(x=1/rowSums(net)) %*% net
      net2=net
      for(iii in 1:20){
        net2=(net2) %*% t(net)
        net2=.extra_matrix_rowNorm(net2)#Matrix::Diagonal(x=1/rowSums(net2)) %*% net2
        net2=net2+net
        net2=net2/2
      }
      net=net2
      if(F){
        diag(net)=0
      }
      
      net=.extra_matrix_rowNorm(input_mat = net,rowValues = 1/qlcMatrix::rowMax(net))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(net)) %*% net
      #diag(net)=1
      
    }
    return(net)
  }
  
  sim_map=.mynetAff(diff)
  
  sim_net=as.data.frame(which(upper.tri(sim_map),arr.ind = T))
  sim_net$score=sim_map[cbind(sim_net[,1],sim_net[,2])]
  sim_net$source=row.names(sim_map)[sim_net[,1]]
  sim_net$target=row.names(sim_map)[sim_net[,2]]
  net=sim_net[sim_net$score>0.1,]
  
  p=.extra_sconline.NetVisFn(net,argList=argList,input_pd = pd,directional=F,attribute_col=attribute_col,lable_nodes = label_nodes)
  
  
  return(p)
}

.sconline.createShinyObj=function(argList,inputPd=NULL,fileName,cell_count_diff_thr=25){
  
  if(is.null(inputPd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPd
  }
  
  
  #res_array=.sconline.fetch_data("sconline_arrays",argList = argList)
  meta_data=.sconline.fetch_data("meta_z",argList = argList)
  de_prefix=paste0("de_sig",sig1_thr,"_pctdiff",pct_diff_thr,"_pct2",pct2_thr)
  if(cell_count_diff_thr>0){
    de_prefix=paste0(de_prefix,"_cellCount",cell_count_diff_thr)
  }
  #pct_mat=meta_data$med_pct.1;argList = argList;meta_z_mat=meta_data$meta_z;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
  if(!file.exists(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))|new_DEcount_run){
    de_pct_res=.extra_sconline.PctScoreFn(pct_mat=meta_data$med_pct.1,argList = argList,meta_z_mat=meta_data$meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F)
    qsave(de_pct_res,file=.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    de_pct_res=qread(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  }
  
  
  input_umap_centroid=.sconline.fetch_data("umap_pseudocells",argList)
  input_umap_centroid=input_umap_centroid[input_umap_centroid$centroid %in% row.names(meta_data$meta_z),]
  
  #net=net[!duplicated(net$Snode),]
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(F){
    ps_anno=as.matrix(.myOneHotFn(inputVector=pd$anno_orig_cellState))
    ps_anno=prop_mat %*% ps_anno
    ps_anno=lapply(1:nrow(ps_anno),function(x) {
      y=which(ps_anno[x,]==max(ps_anno[x,]))[1]
      y=colnames(ps_anno)[y]
      return(data.frame(ps=row.names(ps_anno)[x],anno=y,stringsAsFactors = F))
    })
    ps_anno=do.call("rbind",ps_anno)
    
  }
  
  if(F){
    diff=prop_mat
    diff@x=rep(1,length(diff@x))
    diff=diff %*% t(diff)
    diff=sweep(diff,1,diag(diff),"/")
    diff=as.matrix(diff)
    diff[diff<1/scaling_factor]=1/scaling_factor
    diff=abs(log10(diff))
    diff=diff+scaling_factor*(1-exp(-1*de_pct_res/affinity_param))
  } else {
    myL2normFn=function(inputMat){
      prop_mat2=rowSums(inputMat^2)
      prop_mat2=sqrt(prop_mat2)
      res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
      return(res)
    }
    
    #tst=prop_mat
    #tst@x=rep(1,length(tst@x))
    tst=myL2normFn(inputMat = prop_mat)
    diff=t(tst)
    diff=tst %*% diff
    #diff@x=2*(exp(diff@x/max(quantile(diff@x,0.95),0.1))/(1+exp(diff@x/max(quantile(diff@x,0.95),0.1)))-0.5)
    diff=scaling_factor*(1-diff)
    diff=diff+scaling_factor*(1-exp(-1*(de_pct_res)/affinity_param))
    
  }
  
  if(F){
    matWeights=.myEffSizePropMat(prop_mat)
    matEffectiveSize=matWeights$effective_sample_size
  }
  
  
  diff=diff + t(diff)
  
  local_save_dir_path="~/shinyData"
  if(!dir.exists(local_save_dir_path)){
    dir.create(local_save_dir_path,recursive = T)
  }
  
  meta_z_pct=meta_data$meta_z
  meta_z_pct@x[which(meta_z_pct@x<3)]=0
  meta_z_pct=Matrix::drop0(meta_z_pct)
  meta_z_pct=meta_z_pct*meta_data$med_pct.1
  
  sl_cols=apply(meta_z_pct,2,max)
  meta_z_pct=meta_z_pct[,sl_cols>0]
  
  fd=meta_data$fd[match(colnames(meta_z_pct),row.names(meta_data$fd)),]
  rm_dup_genes=which(!duplicated(fd$gene_short_name))
  meta_z_pct=meta_z_pct[,rm_dup_genes]
  fd=fd[rm_dup_genes,]
  
  rm_dup_genes=which(!is.na(fd$gene_short_name))
  meta_z_pct=meta_z_pct[,rm_dup_genes]
  fd=fd[rm_dup_genes,]
  colnames(meta_z_pct)=fd$gene_short_name
  meta_z_pct=t(meta_z_pct)
  
  qsave(list(clust_data=diff,de_count=de_pct_res,centroid_ump=input_umap_centroid,cell_pd=pd,meta_z_pct=meta_z_pct),file=file.path(local_save_dir_path,paste0("data_shiny",fileName,".qs")))
  
  system(paste0("gsutil -m cp ",file.path(local_save_dir_path,paste0("data_shiny",fileName,".qs"))," gs://macosko_data/vgazesta/shinyData"))
  
  
  return(paste0("transferred to gs://macosko_data/vgazesta/shinyData/",paste0("data_shiny",fileName,".qs")))
}

.sconline.clusteringFn_modified=function(argList,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,affinity_param=5,scaling_factor=10,n_clusters=NULL,inputPd=NULL,clustering_method="average",new_DEcount_run=F,tol_level=0.9,cell_count_diff_thr=25){
  
  #argList=.ArgList;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;affinity_param=5;scaling_factor=10;createPlot=T;anno_col='anno_orig_cellState';inputPd=NULL;cluster_selection_regularExp=NULL;n_clusters=115;new_DEcount_run=F
  
  clustering_method=tolower(clustering_method)
  if(sum(clustering_method %in% c("complete", "average", "mcquitty"))==0){
    stop("Acceptable clustering methods: complete, average, mcquitty")
  }
  
  if(is.null(inputPd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPd
  }
  
  
  if(is.null(n_clusters)){
    
    stop("n_clusters argument should be provided!")
  }
  
  #res_array=.sconline.fetch_data("sconline_arrays",argList = argList)
  meta_data=.sconline.fetch_data("meta_z",argList = argList)
  de_prefix=paste0("de_sig",sig1_thr,"_pctdiff",pct_diff_thr,"_pct2",pct2_thr)
  if(cell_count_diff_thr>0){
    de_prefix=paste0(de_prefix,"_cellCount",cell_count_diff_thr)
  }
  #pct_mat=meta_data$med_pct.1;argList = argList;meta_z_mat=meta_data$meta_z;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
  if(!file.exists(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))|new_DEcount_run){
    de_pct_res=.extra_sconline.PctScoreFn(pct_mat=meta_data$med_pct.1,argList = argList,meta_z_mat=meta_data$meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F)
    qsave(de_pct_res,file=.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    de_pct_res=qread(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  }
  
  
  if(F){
    umap_pseudocells=.sconline.fetch_data("umap_pseudocells",argList)
    p=ggplot(umap_pseudocells,aes(UMAP_1,UMAP_2))+geom_point(data=pd,aes(color=anno_orig_cellState))+geom_label(aes(label=centroid))+theme_bw()
    ggsave(plot=p,file="~/myBucket/torm1.pdf",width=20,height = 13)
    
    sl_net=(de_pct_res+t(de_pct_res))
    sl_net=which(sl_net<5&upper.tri(sl_net),arr.ind = T)
    sl_net=data.frame(Fnode=row.names(de_pct_res)[sl_net[,1]],Snode=row.names(de_pct_res)[sl_net[,2]],de_count=de_pct_res[sl_net],stringsAsFactors = F)
    
    
    ###############
    require(ggnetwork)
    res_net=sl_net
    res_net$score=(1-res_net$de_count/6)^3
    
    net = network::network(res_net, directed = FALSE,matrix.type="edgelist")
    
    
    library(gplots)
    library(devtools)
    library(ggnetwork)
    require(network)
    require(sna)
    require(ggplot2)
    require(Matrix)
    
    netVerNames=network::network.vertex.names(net)
    if(F){
      resClustering=cluster_assignments
      #resClustering$cluster[which(resClustering$cluster=="C_0")]=NA
      net %v% "cluster" = resClustering$cluster[match(netVerNames,as.character(resClustering$pseudocell))]
    }
    
    network::set.edge.attribute(net, "weight", res_net$score)
    
    
    centroid_layout=umap_pseudocells[match(as.character(netVerNames),as.character(umap_pseudocells$centroid)),-1]
    centroid_layout=as.matrix(centroid_layout)
    colnames(centroid_layout)=c("x","y")
    
    net=ggnetwork:::fortify.network(net,layout = centroid_layout)
    
    #clusterCol=resClustering[!duplicated(resClustering$cluster),]
    
    pd_summary=.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"))
    
    sl_net$score=1-sl_net$de_count/6
    sl_net$source=sl_net$Fnode
    sl_net$target=sl_net$Snode
    p1=.extra_sconline.NetVisFn(net=sl_net,input_pd=pd,input_umap_centroid=umap_pseudocells)
    
    
  }
  
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(F){
    ps_anno=as.matrix(.myOneHotFn(inputVector=pd$anno_orig_cellState))
    ps_anno=prop_mat %*% ps_anno
    ps_anno=lapply(1:nrow(ps_anno),function(x) {
      y=which(ps_anno[x,]==max(ps_anno[x,]))[1]
      y=colnames(ps_anno)[y]
      return(data.frame(ps=row.names(ps_anno)[x],anno=y,stringsAsFactors = F))
    })
    ps_anno=do.call("rbind",ps_anno)
    
  }
  
  if(F){
    diff=prop_mat
    diff@x=rep(1,length(diff@x))
    diff=diff %*% t(diff)
    diff=sweep(diff,1,diag(diff),"/")
    diff=as.matrix(diff)
    diff[diff<1/scaling_factor]=1/scaling_factor
    diff=abs(log10(diff))
    diff=diff+scaling_factor*(1-exp(-1*de_pct_res/affinity_param))
  } else {
    myL2normFn=function(inputMat){
      prop_mat2=rowSums(inputMat^2)
      prop_mat2=sqrt(prop_mat2)
      res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
      return(res)
    }
    
    #tst=prop_mat
    #tst@x=rep(1,length(tst@x))
    tst=myL2normFn(inputMat = prop_mat)
    diff=t(tst)
    diff=tst %*% diff
    #diff@x=2*(exp(diff@x/max(quantile(diff@x,0.95),0.1))/(1+exp(diff@x/max(quantile(diff@x,0.95),0.1)))-0.5)
    diff=scaling_factor*(1-diff)
    diff=diff+scaling_factor*(1-exp(-1*(de_pct_res)/affinity_param))
    
  }
  
  if(F){
    matWeights=.myEffSizePropMat(prop_mat)
    matEffectiveSize=matWeights$effective_sample_size
  }
  
  
  diff=diff + t(diff)
  
  diff_clust=hclust(as.dist(diff),method = clustering_method)
  
  
  d_conMat=cutree(diff_clust,k=n_clusters)
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),row.names(pd)]
  prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  #prop anno
  
  colMax_vals_m=qlcMatrix::colMax(prop_merged)
  colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
  prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
  prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
  prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
  prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
  pd2=data.frame(cluster=paste0("C",as.character(prop_m_hardCluster$i)))
  row.names(pd2)=row.names(pd)
  pseudocell_cluster_assignments=data.frame(pseudocell=names(d_conMat),cluster=paste0("C",d_conMat),stringsAsFactors=F)
  
  return(list(cluster_object=diff_clust,distance_matrix=diff,n_clusters=n_clusters,cell_cluster_assignments=pd2,pseudocell_cluster_assignments=pseudocell_cluster_assignments))
}

.extra_sconline.clusteringValidationFn=function(argList,clust_obj,anno_col="anno_orig_cellState",inputPd=NULL,n_clusters=NULL,cluster_selection_regularExp=NULL,flip_xy=F){
  
  require(ggraph)
  library(igraph)
  require(tidyverse) 
  require(ape)
  require(seriation)
  require(qs)
  
  #load("full_CB/sconline/full_cb470-nPCs60/UMAP_anno_vst_HVG3_indScaling_commonExp_varRankThr5000_UMIcor_0.3.rda")
  #prop_mat=qread("full_CB/sconline/full_cb470-nPCs60/res_prop_mat_merged_vst_HVG3_noSplitPorp_regularized_propIter4_sensitiveSearch1_indScaling_commonExp_varRankThr5000_UMIcor_0.3_pseudocell250.qs")
  #meta_data=qread("full_CB/sconline/full_cb470-nPCs60/res_meta_vst_HVG3_noSplitPorp_regularized_propIter4_sensitiveSearch1_indScaling_commonExp_varRankThr5000_UMIcor_0.3_pseudocell250.qs")
  
  
  if(is.null(n_clusters)){
    n_clusters=clust_obj$n_clusters
  }
  
  if(is.null(inputPd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPd
  }
  
  if(sum(colnames(pd)==anno_col)==0){
    stop("Provided anno_col was not identified!")
  }
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  adj2=prop_mat
  
  #adj2=diff2
  
  
  pd=pd[colnames(adj2),]
  pd[,anno_col]=gsub(" ",".",as.character(pd[,anno_col]))
  
  pd[,anno_col]=gsub("[[:punct:]]+", ".", pd[,anno_col])
  tmp= adj2 %*% as.matrix(.myOneHotFn(pd[,anno_col]))
  
  tmp=apply(tmp,1,function(x){
    y=which(x==max(x))[1]
    tmp2=data.frame(cluster=colnames(tmp)[y],purity=x[y])
    tmp2
  })
  
  tmp=do.call("rbind",tmp)
  
  if(!is.null(cluster_selection_regularExp)){
    sl_ind=which(grepl(tolower(cluster_selection_regularExp),tolower(tmp$cluster)))
    diff=diff[sl_ind,sl_ind]
    tmp=tmp[sl_ind,]
    
    clust_obj$cluster_object=hclust(as.dist(diff),method = "average")
    clust_obj$cluster_object=seriation:::reorder.hclust(x=clust_obj$cluster_object, dist=as.dist(diff), method = "OLO")
    pd=pd[which(grepl(cluster_selection_regularExp,tolower(pd[,anno_col]))),]
    
  }
  #tmp2=aggregate(purity~cluster,data=tmp,median)
  
  
  
  graph_tree = as.phylo(clust_obj$cluster_object)
  edges = graph_tree$edge
  node_names=graph_tree$tip.label[edges[,2]]
  node_names[is.na(node_names)]=edges[is.na(node_names),2]
  edges[,2]=node_names
  edges=as.data.frame(edges)
  colnames(edges)=c("from","to")
  # Let's add a column with the group of each name. It will be useful later to color points
  vertices = data.frame(
    name = unique(c(as.character(edges$from), as.character(edges$to)))
  ) 
  vertices$group=tmp$cluster[match(vertices$name,row.names(tmp))]
  #Let's add information concerning the label we are going to add: angle, horizontal adjustement and potential flip
  #calculate the ANGLE of the labels
  vertices$id=NA
  myleaves=which(!is.na(vertices$group))
  nleaves=length(myleaves)
  vertices$id[ myleaves ] = 1:nleaves
  vertices$angle= 90 - 360 * vertices$id / nleaves
  
  # calculate the alignment of labels: right or left
  # If I am on the left part of the plot, my labels have currently an angle < -90
  vertices$hjust<-ifelse( vertices$angle < -90, 1, 0)
  
  # flip angle BY to make them readable
  vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)
  
  # Create a graph object
  vertices$group=as.factor(vertices$group)
  mygraph <- igraph::graph_from_data_frame( edges, vertices=vertices )
  
  # Make the plot
  p_dendogram=ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
    geom_edge_elbow(color="gray") +
    scale_edge_colour_distiller(palette = "RdPu") +
    geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=group, angle = angle, hjust=hjust, colour=group), size=2.7, alpha=1) +
    geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, alpha=0.2)) +
    scale_colour_manual(values= hues::iwanthue(length(unique(vertices$group)))) +
    scale_size_continuous( range = c(0.1,10) ) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0),"cm"),
    ) +
    expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
  #ggsave(plot=p_dendogram,file="~/myBucket/torm.pdf")
  
  #Create confusion matrix
  
  
  d_conMat=cutree(clust_obj$cluster_object,k=n_clusters)
  
  
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
  
  colMax_vals=c()
  for(i in seq(1,ncol(prop_merged),5000)){
    tmp_max=as.numeric(qlcMatrix::colMax(prop_merged[,i:min(i+4999,ncol(prop_merged))]))
    colMax_vals=c(colMax_vals,as.numeric(tmp_max))
  }
  prop_merged = .extra_matrix_colNorm(input_mat = prop_merged,colValues = 1 / (colMax_vals+0.00000000001)) #prop_merged %*% Matrix::Diagonal(x = 1 / (colMax_vals+0.00000000001))
  #prop_mat=Matrix::drop0(prop_mat,0.1)
  prop_merged <- .extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x = 1 / rowSums(prop_merged)) %*% prop_merged
  
  matWeights=.myEffSizePropMat(prop_merged)
  matEffectiveSize=matWeights$effective_sample_size
  
  anno=prop_merged %*% as.matrix(.myOneHotFn(inputVector=pd[,anno_col]))
  matWeights=.myEffSizePropMat(prop_merged)
  matEffectiveSize=matWeights$effective_sample_size
  tmp_coverage=apply(anno,1,function(x) colnames(anno)[which(x==max(x))])
  tmp_coverage=data.frame(sconline=names(tmp_coverage),eff_size=matEffectiveSize,anno=tmp_coverage,stringsAsFactors = F)
  tmp_coverage=merge(tmp_coverage,as.data.frame(table(pd[,anno_col])),by.x="anno",by.y="Var1")
  p_coverage=ggplot(tmp_coverage,aes(anno,eff_size/Freq))+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.8,hjust=1,color="black"))
  
  if(!flip_xy){
    anno=apply(anno,2,function(x) x/sum(x))
  }
  
  
  anno=as.data.frame(anno)
  anno$ps=row.names(prop_merged)
  anno=reshape::melt(anno,id.vars="ps")
  anno=anno[order(anno$value,decreasing = T),]
  anno$ps=factor(as.character(anno$ps),levels=unique(anno$ps))
  anno=anno[order(anno$ps),]
  tst_lbl=as.character(anno$variable)[anno$value!=0]
  tst_lbl=tst_lbl[!duplicated(tst_lbl)]
  tst_lbl=c(tst_lbl,setdiff(as.character(anno$variable),tst_lbl))
  anno$variable=factor(as.character(anno$variable),levels=tst_lbl)
  recovered_clusters=anno[order(anno$value,decreasing = T),]
  recovered_clusters=recovered_clusters[!duplicated(recovered_clusters$ps),]
  
  tmp_anno=anno[order(anno$value,decreasing = T),]
  tmp_anno=tmp_anno[!duplicated(tmp_anno$ps),]
  
  
  
  tmp_recall=unlist(lapply(1:nrow(tmp_anno),function(x){
    length(unique(tmp_anno$variable[1:x]))/length(unique(anno$variable))
  }))
  tmp_recall=data.frame(purity=tmp_anno$value,Recall=tmp_recall,stringsAsFactors = F)
  
  
  colMax_vals_m=qlcMatrix::colMax(prop_merged)
  colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=0.99)
  anno_m=colMax_vals_m %*% as.matrix(.myOneHotFn(inputVector=pd[,anno_col]))
  anno_m=.extra_matrix_rowNorm(as.matrix(anno_m))#as.matrix(Matrix::Diagonal(x=1/rowSums(anno_m)) %*% anno_m)
  anno_m=as.data.frame(anno_m)
  anno_m$cluster=row.names(prop_merged)
  anno_m=reshape::melt(anno_m,id.vars="cluster")
  anno_m=anno_m[order(anno_m$value,decreasing = T),]
  anno_m=anno_m[!duplicated(anno_m$cluster),]
  
  cdf_binary=unlist(lapply(1:nrow(anno_m), function(x) length(unique(anno_m$variable[1:x]))))
  cdf_binary=data.frame(purity=anno_m$value,Recall=cdf_binary/length(unique(anno$variable)),stringsAsFactors = F)
  cdf_binary=cdf_binary[order(cdf_binary$Recall,decreasing = T),]
  cdf_binary=cdf_binary[!duplicated(cdf_binary$purity),]
  cdf_binary=cdf_binary[order(cdf_binary$purity,decreasing = T),]
  cdf_binary$status="Binary"
  tmp_recall$status="Continuous"
  cdf_analysis=rbind(tmp_recall,cdf_binary)
  
  if(flip_xy){
    p_anno_celllevel=ggplot(data=anno,aes(variable,ps,fill=value))+geom_tile(color="black")+scale_fill_gradient(low="white",high="red")+ggtitle(paste("Total # recoved clusters:",nrow(recovered_clusters)))
  } else {
    p_anno_celllevel=ggplot(data=anno,aes(ps,variable,fill=value))+geom_tile(color="black")+scale_fill_gradient(low="white",high="red")+ggtitle(paste("Total # recoved clusters:",nrow(recovered_clusters)))
  }
  
  #ggsave(plot=p_anno_celllevel,file="~/myBucket/torm.pdf")
  
  d_conMat=cutree(clust_obj$cluster_object,k=n_clusters)
  d_conMat=data.frame(tmp,cluster_hclust=d_conMat,stringsAsFactors = F)
  d_conMat=reshape2::dcast(cluster~cluster_hclust,data=d_conMat,value.var = "purity",fun.aggregate =length)
  d_conMat=d_conMat[match(unique(pd[,anno_col]),d_conMat[,1]),]
  row.names(d_conMat)=unique(pd[,anno_col])
  d_conMat=d_conMat[,-1]
  d_conMat[is.na(d_conMat)]=(0)
  d_conMat=sweep(d_conMat,2,colSums(d_conMat),"/")
  d_conMat$cluster=row.names(d_conMat)
  d_conMat=reshape::melt(d_conMat,id.vars="cluster")
  colnames(d_conMat)[colnames(d_conMat)=="variable"]="sconline_cluster"
  d_conMat=d_conMat[order(d_conMat$value,as.character(d_conMat$cluster),decreasing = T),]
  d_conMat$sconline_cluster=factor(as.character(d_conMat$sconline_cluster),levels=unique(d_conMat$sconline_cluster))
  d_conMat=d_conMat[order(d_conMat$sconline_cluster),]
  tst_lbl=d_conMat$cluster[d_conMat$value!=0]
  tst_lbl=tst_lbl[!duplicated(tst_lbl)]
  tst_lbl=c(tst_lbl,setdiff(as.character(d_conMat$cluster),tst_lbl))
  d_conMat$cluster=factor(as.character(d_conMat$cluster),levels=tst_lbl[!duplicated(tst_lbl)])
  if(flip_xy){
    p_confusionMatrix=ggplot(data=d_conMat,aes(cluster,sconline_cluster,fill=value))+geom_tile(color="black")+scale_fill_gradient(low="white",high="red")+ylab("")
  } else {
    p_confusionMatrix=ggplot(data=d_conMat,aes(sconline_cluster,cluster,fill=value))+geom_tile(color="black")+scale_fill_gradient(low="white",high="red")+ylab("")
  }
  
  
  return(list(plot_dendrogram=p_dendogram,plot_confusionMatrix=p_confusionMatrix,p_coverage=p_coverage,p_anno_celllevel=p_anno_celllevel,cdf_analysis=cdf_analysis))
}

.sconline.recoveredClusters=function(inputPhenoData,cluster_col,cellType_col){
  purity_anno=aggregate(as.formula(paste0(cellType_col,"~",cluster_col)),data=inputPhenoData,function(x) {x=table(x); y=x/sum(x); x=names(x)[order(y,decreasing = T)];x[1]})
  colnames(purity_anno)=c(cluster_col,"cellType")
  purity_value=aggregate(as.formula(paste0(cellType_col,"~",cluster_col)),data=inputPhenoData,function(x) {x=table(x); max(x)/sum(x)})
  colnames(purity_value)=c(cluster_col,"purity")
  purity_res=merge(purity_anno,purity_value,by=cluster_col)
  return(purity_res)
}

.sconline.recoveredClusters_resiprocal=function(inputPhenoData,cluster_col,cellType_col){
  purity_anno=aggregate(as.formula(paste0(cluster_col,"~",cellType_col)),data=inputPhenoData,function(x) {x=table(x); max(x)/sum(x)})
  colnames(purity_anno)=c(cellType_col,"purity_ctype")
  
  purity_link=aggregate(as.formula(paste0(cellType_col,"~",cluster_col)),data=inputPhenoData,function(x) {x=table(x); y=x/sum(x); x=names(x)[order(y,decreasing = T)];x[1]})
  colnames(purity_link)=c(cluster_col,cellType_col)
  
  purity_cluster=aggregate(as.formula(paste0(cellType_col,"~",cluster_col)),data=inputPhenoData,function(x) {x=table(x); max(x)/sum(x)})
  colnames(purity_cluster)=c(cluster_col,"purity_cluster")
  
  x=merge(purity_anno,purity_link,by=cellType_col)
  x=merge(x,purity_cluster,cluster_col)
  
  return(x)
}


.sconline.mastDEanalysis=function(argList,inputExpData,clust_obj,maineffect="genotype",covList=c("dissection"),randEffect_var="anno_batch",contrast_df,n_clusters=NULL,inputPd=NULL,verbose = TRUE,ncores=10,tol_level=0.95){
  #argList=.ArgList;inputExpData=xpo_data;clust_obj=clust_obj;maineffect="genotype";covList=c("dissection");randEffect_var="anno_batch";n_clusters=NULL;inputPd=NULL;verbose = TRUE;ncores=10;tol_level=0.95
  #contrast_df=data.frame(g1=c("WT","WT"),g2=c("HET","KO"),stringsAsFactors = F)
  #formula_str="~genotype+dissection+(1|anno_batch)"
  
  require(scater)
  require(MAST)
  
  #contrast_df=data.frame(g1=c("WT","WT"),g2=c("HET","KO"),stringsAsFactors = F)
  #formula_str="~genotype+dissection+(1|anno_batch)"
  
  
  
  if(ncores>1){
    options("mc.cores"=ncores)
  }
  
  if(is.null(n_clusters)){
    n_clusters=clust_obj$n_clusters
  }
  
  if(is.null(inputPd)){
    inputPd=.sconline.fetch_data("annotation",argList)
  }
  
  min_obs_size=length(unique(inputPd[,maineffect]))
  if(!is.null(covList)){
    for(icov in covList){
      min_obs_size=max(min_obs_size,length(unique(inputPd[,icov])))
    }
  }
  
  if(!is.null(randEffect_var)){
    min_obs_size=max(min_obs_size,length(unique(inputPd[,randEffect_var])))
  }
  
  min_obs_size=max(min_obs_size,10)
  
  d_conMat=cutree(clust_obj$cluster_object,k=n_clusters)
  
  inputExpData=inputExpData[,colnames(inputExpData) %in% row.names(inputPd)]
  inputPd=inputPd[match(colnames(inputExpData),row.names(inputPd)),]
  
  {
    for(icov in 1:ncol(colData(inputExpData))){
      colData(inputExpData)[,icov]=gsub(" ",".",as.character(colData(inputExpData)[,icov]))
    }
    
  }
  
  
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=prop_mat[,row.names(inputPd)]
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),row.names(inputPd)]
  prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  #prop anno
  
  colMax_vals_m=qlcMatrix::colMax(prop_merged)
  colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_merged,colValues = 1/as.numeric(colMax_vals_m))# prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
  prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
  prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
  prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
  prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(inputPd),prop_m_hardCluster$j),]
  inputExpData$cluster_anno_res=inputPd$cluster_anno_res=paste0("C",as.character(prop_m_hardCluster$i))
  
  inputExpData=.mySplitObject(inputExpData,"cluster_anno_res")
  
  #formula_str=paste0(varList,collapse = "+")
  formula_str=paste0("~",maineffect)
  randomEffectModel=F
  if(!is.null(randEffect_var)){
    randomEffectModel=T
    randEffect_var=paste0("(1|",randEffect_var,")")
  }
  
  
  
  de_res=list()
  
  for(iclust in names(inputExpData)){
    for(icontrast in 1:nrow(contrast_df)){
      
      tmpexp=inputExpData[[iclust]]
      tmpexp=tmpexp[,colData(tmpexp)[,maineffect] %in% c(contrast_df[icontrast,1],contrast_df[icontrast,2])]
      if(length(unique(colData(tmpexp)[,maineffect]))>1){
        colData(tmpexp)[,maineffect]=factor(as.character(colData(tmpexp)[,maineffect]),levels=c(contrast_df[icontrast,1],contrast_df[icontrast,2]))
        tmpgroup=inputPd[match(colnames(tmpexp),row.names(inputPd)),]
        tmp_formula_str=formula_str
        for(icov in covList){
          if(length(unique(tmpgroup[,icov]))>1){
            tmp_formula_str=paste0(tmp_formula_str,"+",icov)
          }
        }
        
        if(randomEffectModel){
          tmp_formula_str=paste0(tmp_formula_str,"+",randEffect_var)
        }
        
        tmpexp_binary=counts(inputExpData[[iclust]])
        
        tmpexp_binary@x=rep(1,length(tmpexp_binary@x))
        tmpCount=rowSums(tmpexp_binary)
        tmpexp=tmpexp[tmpCount>min_obs_size,]
        
        
        tmpexp = scater::logNormCounts(tmpexp)
        zlm.res=list()
        for(igene in seq(1,nrow(tmpexp),5000)){
          sca = MAST::SceToSingleCellAssay(tmpexp[igene:min(igene+4999,nrow(tmpexp)),])
          
          if(randomEffectModel){
            #,fitArgsD=list(nAGQ=0)
            zlmCond <- MAST::zlm(as.formula(tmp_formula_str), sca = sca,method='glmer', ebayes=FALSE,parallel = T)
          } else {
            zlmCond <- MAST::zlm(as.formula(tmp_formula_str), sca = sca,parallel = T)
          }
          
          zlm.lr <- MAST::lrTest(zlmCond, maineffect)
          zlm.lr=zlm.lr[,,"Pr(>Chisq)"]
          zlm.lr=as.data.frame(zlm.lr)
          zlm.res=c(zlm.res,list(zlm.lr))
        }
        zlm.res=do.call("rbind",zlm.res)
        
        tmp_seurat=.extraExport2SeuratFn(tmpexp)
        tmp_seurat=Seurat::NormalizeData(tmp_seurat)
        tmp_fc_data=.myEvalMarkers(object=tmp_seurat, cells.1=colnames(tmp_seurat)[tmp_seurat@meta.data[,maineffect]==contrast_df[icontrast,2]], cells.2=colnames(tmp_seurat)[tmp_seurat@meta.data[,maineffect]==contrast_df[icontrast,1]], slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cluster_name=iclust)
        tmp_arranged=merge(data.frame(gene=row.names(tmp_fc_data),tmp_fc_data,stringsAsFactors = F),data.frame(gene=row.names(zlm.res),zlm.res,stringsAsFactors = F),by="gene",all=T)
        tmp_arranged$cluster=iclust
        tmp_arranged$contrast=paste0(contrast_df[icontrast,2],"_vs_",contrast_df[icontrast,1])
        de_res=c(de_res,list(tmp_arranged))
      }
      
    }
  }
  de_res=do.call("rbind",de_res)
  
  return(de_res)
}


.sconline.cluster.markers=function(argList,input_meta_data=NULL,inputPhenoData=NULL,cluster_obj,input_prop_mat=NULL,pseudocell_assignments=NULL,inputExpData=NULL,n_clusters=NULL,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,only.pos.logFC=T,logFC_thr=0.01,return_pct_mat=F,hierarchical_mode=F,return_marker_data=F){
  
  #argList=.ArgList;cluster_obj=.tmp;n_clusters=NULL;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;inputExpData=.xpo_data
  
  #cluster_obj: output of .sconline.clusteringFn()
  #n_clusters: number of desired clusters
  #sig1_thr: z-score threshold
  #pct2_thr: marker pct.2 threshold (% expressed in other clusters)
  #pct_diff_thr: pct.1 - pct.2 threshold
  #only.pos.logFC: report only upregulated genes in the cluster
  
  #require(ggraph)
  #library(igraph)
  #require(tidyverse)
  #require(ape)
  #require(seriation) 
  
  
  if(!is.null(pseudocell_assignments)){
    d_conMat=pseudocell_assignments
  } else {
    if(is.null(n_clusters)){
      n_clusters=cluster_obj$n_clusters
    }
    if(is.null(n_clusters)){
      stop("Number of desired clusters should be provided!")
    }
    
    diff_clust=cluster_obj$cluster_object
    if(is.null(cluster_obj$cluster_object)){
      stop("A proper cluster object needs to be provided")
    }
    d_conMat=cutree(diff_clust,k=n_clusters)
  }
  
  
  if(is.null(input_prop_mat)){
    prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    prop_mat=input_prop_mat
  }
  
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(F){
    marker_data=NULL
    if(!file.exists(.myFilePathMakerFn("res_marker_data",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))|hierarchical_mode){
      prop_mat_c=prop_mat
      if(quantile(rowSums(prop_mat_c > 0,na.rm = T), 0.25) < (0.85*ncol(prop_mat))){
        prop_mat_c=Matrix::drop0(prop_mat_c)
        prop_mat_c@x=rep(1,length(prop_mat_c@x))
      }
      
      prop_mat_c=Matrix::drop0(1-prop_mat_c)
      prop_mat_c <- .extra_matrix_rowNorm(prop_mat_c)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat_c))) %*% prop_mat_c
      
      row_sum_counts=prop_mat_c
      row_sum_counts@x=rep(1,length(row_sum_counts))
      row_sum_counts=rowSums(row_sum_counts)
      
      tmp_cosine=qlcMatrix::cosSparse(t(prop_mat_c))
      row.names(tmp_cosine)=colnames(tmp_cosine)=row.names(prop_mat_c)
      tmp_cosine_clust=hclust(as.dist(1-as.matrix(tmp_cosine)),method = "complete")
      tmp_clust=cutree(tmp_cosine_clust,h = 0.999)
      mapping=data.frame(clust_id=tmp_clust,pseudocell=names(tmp_clust),stringsAsFactors = F)
      map_names=mapping[!duplicated(mapping[,1]),]
      colnames(map_names)[2]="cluster"
      mapping=merge(mapping,map_names,by="clust_id")
      if(sum(row_sum_counts<10000)>0){
        mapping$cluster[mapping$pseudocell %in% row.names(prop_mat_c)[row_sum_counts<10000]]=mapping$pseudocell[mapping$pseudocell %in% row.names(prop_mat_c)[row_sum_counts<10000]]
      }
      prop_mat_c_reduced=prop_mat_c[row.names(prop_mat_c) %in% mapping$cluster,,drop=F]
      
      
      
      logNormData=t(Seurat:::NormalizeData.default(counts(inputExpData)[,colnames(prop_mat)],normalization.method = "RC",verbose = F))
      exp_binary=logNormData
      exp_binary@x=rep(1,length(exp_binary@x))
      
      fc_1=prop_mat %*% logNormData
      batch_size=100
      if(nrow(prop_mat_c_reduced) %% 100 ==0){
        batch_size=99
      }
      
      if(nrow(prop_mat_c_reduced)>batch_size){
        fc_2=parallel::mclapply(seq(1,nrow(prop_mat_c_reduced),batch_size),function(x) {prop_mat_c_reduced[x:min(x+batch_size-1,nrow(prop_mat_c_reduced)),] %*% logNormData},mc.cores = argList$ncores)
        fc_2=do.call("rbind",fc_2)
      } else {
        fc_2=prop_mat_c_reduced %*% logNormData
      }
      
      fc_2=fc_2[match(mapping$cluster,row.names(fc_2)),]
      row.names(fc_2)=mapping$pseudocell
      
      fc_1@x=log2(fc_1@x+1)
      fc_2@x=log2(fc_2@x+1)
      
      logFC_res=fc_1 - fc_2#log2(prop_mat %*% logNormData+1) - log2(prop_mat_c %*% logNormData+1)
      
      
      pct.1=prop_mat %*% exp_binary
      
      sl_ind=which(pct.1<logFC_thr)
      if(length(sl_ind)>0){
        logFC_res[sl_ind]=0
      }
      
      if(!hierarchical_mode){
        qsave(list(pct.1=pct.1,logFC=logFC_res),file=.myFilePathMakerFn("res_marker_data",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
      }
      
      marker_data=list(pct.1=pct.1,logFC=logFC_res)
      
    } else {
      marker_data=qread(.myFilePathMakerFn("res_marker_data",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
      pct.1=marker_data$pct.1
      logFC_res=marker_data$logFC
    }
  }
  
  if(is.null(input_meta_data)){
    meta_data=.sconline.fetch_data("meta_z",argList = argList)
  } else {
    meta_data=input_meta_data
  }
  
  pct.1=meta_data$med_pct.1
  logFC_res=meta_data$logFC
  
  if(is.null(inputPhenoData)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPhenoData
  }
  
  
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
  prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  
  
  z_mat=meta_data$meta_z
  matWeights=.myEffSizePropMat(prop_mat)
  matEffectiveSize=matWeights$effective_sample_size
  z_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  z_merged_c=Matrix::drop0(1-z_merged)
  z_merged=sweep(z_merged,2,matEffectiveSize,"*")
  z_merged_c=sweep(z_merged_c,2,matEffectiveSize,"*")
  z_merged=.extra_matrix_rowNorm(z_merged)#Matrix::Diagonal(x=1/rowSums(z_merged)) %*% z_merged
  z_merged_c=.extra_matrix_rowNorm(z_merged_c)#Matrix::Diagonal(x=1/rowSums(z_merged_c)) %*% z_merged_c
  pct.2=z_merged_c %*% pct.1[,colnames(meta_data$meta_z)]
  pct.1=z_merged %*% pct.1[,colnames(meta_data$meta_z)]
  
  logFC_res=z_merged %*% logFC_res[colnames(z_merged),colnames(meta_data$meta_z)]
  z_merged= z_merged %*% meta_data$meta_z[colnames(z_merged),]
  
  pct.1_mat=NULL
  z_mat=NULL
  if(return_pct_mat){
    pct.1_mat=pct.1
    z_mat=z_merged
  }
  
  pct_diff=pct.1 - pct.2
  pct_diff[pct_diff<0]=0
  pct_diff=Matrix::drop0(pct_diff,tol=pct_diff_thr)
  
  
  
  pct_diff=pct_diff[,colnames(z_merged)]
  pct.1=pct.1[,colnames(z_merged)]
  pct.2=pct.2[,colnames(z_merged)]
  logFC_res=logFC_res[,colnames(z_merged)]
  z_merged=Matrix::drop0(z_merged,tol=sig1_thr)
  
  if(only.pos.logFC){
    z_merged@x[z_merged@x<0]=0
    z_merged=Matrix::drop0(z_merged,tol=sig1_thr)
  }
  
  pct_diff@x=rep(1,length(pct_diff@x))
  pct.2_binary=as.matrix(pct.2)
  pct.2_binary[pct.2_binary>pct2_thr]=(-1)
  pct.2_binary[pct.2_binary>=0]=1
  pct.2_binary[pct.2_binary<0]=0
  z_merged=z_merged*pct_diff*pct.2_binary
  
  z_binary=z_merged
  z_binary[z_binary!=0]=1
  pct.1=Matrix::drop0(pct.1*z_binary)
  pct.2=Matrix::drop0(pct.2*z_binary)
  logFC_res=Matrix::drop0(logFC_res*z_binary)
  pct_diff=Matrix::drop0(pct.1 - pct.2)
  
  my_sparseMat_summary=function(inputMat,val_col_name=NULL){
    res=as.data.frame(summary(inputMat))
    res$i=row.names(inputMat)[res$i]
    res$j=colnames(inputMat)[res$j]
    if(!is.null(val_col_name)){
      colnames(res)[colnames(res)=="x"]=val_col_name
    }
    return(res)
  }
  
  z_merged=my_sparseMat_summary(z_merged,val_col_name = "zscore")
  pct.1=my_sparseMat_summary(pct.1,val_col_name = "pct.1")
  pct.2=my_sparseMat_summary(pct.2,val_col_name = "pct.2")
  logFC_res=my_sparseMat_summary(logFC_res,val_col_name = "avg_log2FC")
  
  res=merge(z_merged,pct.1,by=c("i","j"))
  res=merge(res,pct.2,by=c("i","j"))
  res=merge(res,logFC_res,by=c("i","j"))
  colnames(res)[colnames(res)=="i"]="cluster"
  colnames(res)[colnames(res)=="j"]="gene"
  res=merge(data.frame(rowname=row.names(meta_data$fd),meta_data$fd,stringsAsFactors = F),res,by.y="gene",by.x="rowname",all.y=T)
  #res=res[,colnames(res)!="rowname"]
  
  res=res[which(res$zscore*res$avg_log2FC>0),]
  
  if(return_pct_mat){
    res=list(res=res)
    res=c(res,list(pct.1=pct.1_mat,z=z_mat))
  }
  
  return(res)
  
}


.sconline.cluster_pruning=function(cluster_assignments,clust_obj, inputExpData, argList,input_meta_data=NULL,inputPhenoData=NULL,input_prop_mat=NULL,
                                   combinatorial_pct_tol=1,marker_sig1_thr=3,marker_pct2_thr=0.3,marker_pct_diff_thr=0.2,forgiveness_factor=1){
  n_clusters=length(unique(cluster_assignments))
  
  combinatorial_counts=.sconline.marker.combinatorial_count(argList=argList,inputExpData=inputExpData,combinatorial_pct_tol=combinatorial_pct_tol,marker_sig1_thr=marker_sig1_thr,marker_pct2_thr=marker_pct2_thr,marker_pct_diff_thr=marker_pct_diff_thr,cluster_assignments=cluster_assignments,fast_mode=T,input_meta_data=input_meta_data,inputPhenoData=inputPhenoData,input_prop_mat=input_prop_mat)
  clust_list=data.frame(pseudocell=names(cluster_assignments),cluster=cluster_assignments,stringsAsFactors = F)
  for(iclust in (n_clusters):2){
    tmp=cutree(clust_obj,k=iclust)
    tmp=tmp[match(clust_list$pseudocell,names(tmp))]
    clust_list$new=paste0(iclust,"_",tmp)
    colnames(clust_list)[colnames(clust_list)=="new"]=paste0("x",iclust)
  }
  
  merged_net=cluster_assignments
  
  while(sum(combinatorial_counts[,2]==0)>forgiveness_factor){
    if(length(unique(merged_net))==2){
      merged_net=rep(1,length(cluster_assignments))
      names(merged_net)=names(cluster_assignments)
      break;
    } else {
      i=combinatorial_counts[which(combinatorial_counts[,2]==0),1]
      {
        rank_list=unlist(lapply(i,function(x){
          tmp=clust_list[which(clust_list$cluster==x),]
          tmp=unlist(lapply(colnames(tmp)[grepl("^x",colnames(tmp))],function(y) {
            sum(clust_list[,y]==unique(tmp[,y])[1])>nrow(tmp)&length(unique(tmp[,y]))==1
          }))
          which(tmp)[1]
        }))
        i=i[order(rank_list,decreasing = F)]
        i=i[1]
        i_ind=rank_list[order(rank_list,decreasing = F)]
        i_ind=colnames(clust_list)[grepl("^x",colnames(clust_list))][i_ind][1]
        if(is.na(i_ind)){
          merged_net=rep(1,length(cluster_assignments))
          names(merged_net)=names(cluster_assignments)
          break;
        }
      }
      
      #sl=which(clust_list[,i_ind]==unique(clust_list[clust_list$cluster==i,i_ind]))
      #table(clust_list$cluster[sl])
      
      clust_list$cluster[which(clust_list[,i_ind]==unique(clust_list[clust_list$cluster==i,i_ind]))]=unique(clust_list[clust_list$cluster==i,i_ind])
      
      merged_net=clust_list$cluster
      names(merged_net)=clust_list$pseudocell
      combinatorial_counts=.sconline.marker.combinatorial_count(argList=argList,inputExpData=inputExpData,cluster_obj=clust_obj,n_clusters=n_clusters,cluster_assignments=merged_net,combinatorial_pct_tol=combinatorial_pct_tol,marker_sig1_thr=marker_sig1_thr,marker_pct2_thr=marker_pct2_thr,marker_pct_diff_thr=marker_pct_diff_thr,fast_mode=T,input_meta_data=input_meta_data,inputPhenoData=inputPhenoData,input_prop_mat=input_prop_mat)
    }
    
    
    
  }
  
  
  return(merged_net)
  
}

.sconline.marker.combinatorial_count=function(argList,inputExpData,cluster_obj,n_clusters,combinatorial_pct_tol=1,marker_sig1_thr=3,marker_pct2_thr=0.3,marker_pct_diff_thr=0.2,cluster_assignments=NULL,fast_mode=F,input_meta_data=NULL,inputPhenoData=NULL,input_prop_mat=NULL){
  
  #argList=.ArgList;inputExpData=data;cluster_obj=cluster_res_avg;n_clusters=5
  #combinatorial_pct_tol=0.95;marker_sig1_thr=3;marker_pct2_thr=0.3;marker_pct_diff_thr=0.2;cluster_assignments=NULL;fast_mode=F;input_meta_data=NULL;inputPhenoData=NULL;input_prop_mat=NULL;
  markers=.sconline.cluster.markers(argList=argList,cluster_obj=cluster_obj,inputExpData=inputExpData,n_clusters=n_clusters,sig1_thr=marker_sig1_thr,pct2_thr=marker_pct2_thr,pct_diff_thr=marker_pct_diff_thr,only.pos.logFC=T,pseudocell_assignments=cluster_assignments,input_meta_data=input_meta_data,inputPhenoData=inputPhenoData,input_prop_mat=input_prop_mat)
  
  if(is.null(input_meta_data)){
    marker_data=.sconline.fetch_data("meta_z",argList = argList)
  } else {
    marker_data=input_meta_data
  }
  
  
  if(is.null(cluster_assignments)){
    if(is.null(n_clusters)){
      n_clusters=cluster_obj$n_clusters
    }
    if(is.null(n_clusters)){
      stop("Number of desired clusters should be provided!")
    }
    
    diff_clust=cluster_obj$cluster_object
    if(is.null(cluster_obj$cluster_object)){
      stop("A proper cluster object needs to be provided")
    }
    d_conMat=cutree(diff_clust,k=n_clusters)
  } else {
    d_conMat=cluster_assignments
  }
  
  if(combinatorial_pct_tol>1|combinatorial_pct_tol<0){
    stop("combinatorial_pct_tol is expected to be between 0 and 1")
  }
  
  if(is.null(input_prop_mat)){
    prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    prop_mat=input_prop_mat
  }
  
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  pct.1=marker_data$med_pct.1
  
  matWeights=.myEffSizePropMat(prop_mat)
  matEffectiveSize=matWeights$effective_sample_size
  p_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  p_merged=sweep(p_merged,2,matEffectiveSize,"*")
  p_merged=.extra_matrix_rowNorm(p_merged)#Matrix::Diagonal(x=1/rowSums(p_merged)) %*% p_merged
  pct.1=p_merged %*% pct.1[colnames(p_merged),]
  
  res_combination_counts=NULL
  for(i in row.names(pct.1)){
    combinations=c()
    if(sum(markers$cluster==i)>0){
      tmp_pct.1=pct.1[-which(row.names(pct.1)==i),colnames(pct.1) %in% markers$rowname[markers$cluster==i],drop=F]
      tmp_pct.1_colmin=qlcMatrix::colMin(tmp_pct.1)
      pct.diff=(-1)*t(sweep(as.matrix(tmp_pct.1),2,as.numeric(pct.1[which(row.names(pct.1)==i),colnames(pct.1) %in% markers$rowname[markers$cluster==i]]),"-"))
      tmp_pct.1=t(sweep(as.matrix(tmp_pct.1),2,tmp_pct.1_colmin,"-"))
      tmp_pct.1[which(pct.diff<marker_pct_diff_thr)]=2
      tmp_marker_count=unlist(lapply(1:nrow(tmp_pct.1),function(x) sum(tmp_pct.1[x,]<combinatorial_pct_tol)))
      if(sum(tmp_marker_count==ncol(tmp_pct.1))>0){
        combinations=c(combinations,rep(row.names(tmp_pct.1),sum(tmp_marker_count==ncol(tmp_pct.1))))
        tmp_pct.1=tmp_pct.1[-which(tmp_marker_count==ncol(tmp_pct.1)),,drop=F]
      }
      
      if(fast_mode&length(combinations)>0){
        res_combination_counts=rbind(res_combination_counts,data.frame(cluster=i,combination_count=length(combinations),stringsAsFactors = F))
        next;
      }
      if(nrow(tmp_pct.1)>0){
        if(fast_mode){
          for(j in row.names(tmp_pct.1)){
            tmp=tmp_pct.1[-which(row.names(tmp_pct.1)==j),which(tmp_pct.1[j,]>=combinatorial_pct_tol),drop=F]
            tmp2=apply(tmp,1,function(x) sum(x<combinatorial_pct_tol))
            if(sum(tmp2==ncol(tmp))>0){
              combinations=c(combinations,row.names(tmp)[which(tmp2==ncol(tmp))])
              break;
            }
          }
        } else {
          for(j in row.names(tmp_pct.1)){
            tmp=tmp_pct.1[-which(row.names(tmp_pct.1)==j),which(tmp_pct.1[j,]>=combinatorial_pct_tol),drop=F]
            tmp2=apply(tmp,1,function(x) sum(x<combinatorial_pct_tol))
            if(sum(tmp2==ncol(tmp))>0){
              #stop("here")
              combinations=c(combinations,row.names(tmp)[which(tmp2==ncol(tmp))])
            }
          }
        }
      }
      
    }
    res_combination_counts=rbind(res_combination_counts,data.frame(cluster=i,combination_count=length(combinations),stringsAsFactors = F))
  }
  return(res_combination_counts)
}

.sconline.cluter.marker_exploration=function(argList,inputExpData,clust_obj,cluster_count_range,combinatorial_pct_tol=1,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,only.pos.logFC=T){
  
  DE_w_dup=NULL
  DE_wo_dup=NULL
  DE_wo_dup_combinatorial=NULL
  #singleton_count=NULL
  cluster_size_dist=NULL
  
  for(iclust in cluster_count_range){
    markers=.sconline.cluster.markers(argList=argList,cluster_obj=clust_obj,inputExpData=inputExpData,n_clusters=iclust,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,only.pos.logFC=T)
    markers$cluster=factor(as.character(markers$cluster),levels=1:iclust)
    
    w_dup_count=as.data.frame(table(markers$cluster))
    wo_dup_count=markers$rowname[duplicated(markers$rowname)]
    wo_dup_count=markers[!markers$rowname %in% wo_dup_count,]
    wo_dup_count=as.data.frame(table(wo_dup_count$cluster))
    colnames(w_dup_count)[1]=colnames(wo_dup_count)[1]="cluster"
    
    cluster_size=.sconline.cluster(argList,n_clusters=iclust,clustering_method="average")
    #singleton_ps_count=as.data.frame(table(cluster_size$pseudocell_cluster_assignments[,2]))
    #singleton_ps_count=sum(singleton_ps_count[,2]==1)
    cluster_size=as.data.frame(table(cluster_size$cell_cluster_assignments[,1]))
    cluster_size$cluster_count=w_dup_count$cluster_count=wo_dup_count$cluster_count=iclust
    #singleton_ps_count=data.frame(cluster_count=iclust,count=singleton_ps_count)
    DE_w_dup=rbind(DE_w_dup,w_dup_count)
    DE_wo_dup=rbind(DE_wo_dup,wo_dup_count)
    cluster_size_dist=rbind(cluster_size_dist,cluster_size)
    #singleton_count=rbind(singleton_count,singleton_ps_count)
    
    #double marker count
    combinatorial_counts=.sconline.marker.combinatorial_count(argList=argList,inputExpData=inputExpData,cluster_obj=clust_obj,n_clusters=iclust,combinatorial_pct_tol=combinatorial_pct_tol,marker_sig1_thr=sig1_thr,marker_pct2_thr=pct2_thr,marker_pct_diff_thr=pct_diff_thr)
    DE_wo_dup_combinatorial=rbind(DE_wo_dup_combinatorial,data.frame(cluster_count=iclust,pct_wo_markers=sum(combinatorial_counts[,2]!=0)/nrow(combinatorial_counts)*100))
  }
  
  
  
  DE_w_dup$cluster_count=as.factor(DE_w_dup$cluster_count)
  DE_wo_dup$cluster_count=as.factor(DE_wo_dup$cluster_count)
  #singleton_count$cluster_count=as.factor(singleton_count$cluster_count)
  cluster_size_dist$cluster_count=as.factor(cluster_size_dist$cluster_count)
  DE_wo_dup_combinatorial$cluster_count=as.factor(DE_wo_dup_combinatorial$cluster_count)
  
  p1=ggplot(DE_w_dup,aes(cluster_count,Freq,fill=as.numeric(cluster_count)))+geom_violin()+geom_hline(yintercept = 0)+
    cowplot::theme_cowplot()+ylab("# Markers (can be duplicated between clusters)")+
    scale_fill_gradient(low="#E1AD9D",high="steelblue")+theme(legend.position = "none")
  
  p2=ggplot(DE_wo_dup,aes(cluster_count,Freq,fill=as.numeric(cluster_count)))+geom_violin()+geom_hline(yintercept = 0)+
    cowplot::theme_cowplot()+ylab("# Unique Markers (unique to each cluster)")+
    scale_fill_gradient(low="#E1AD9D",high="steelblue")+theme(legend.position = "none")
  
  
  p3=ggplot(DE_wo_dup_combinatorial,aes(cluster_count,pct_wo_markers))+geom_point()+
    cowplot::theme_cowplot()+ylab("% clusters with unique combinatorial markers")
  
  
  p4=ggplot(cluster_size_dist,aes(cluster_count,Freq,fill=as.numeric(cluster_count)))+geom_violin()+
    cowplot::theme_cowplot()+ylab("Cluster size distribution")+scale_y_log10()+
    scale_fill_gradient(low="#E1AD9D",high="steelblue")+theme(legend.position = "none")
  
  p=p1+p2+p3+p4
  return(p)
}

.sconline.eval_markers=function(argList,cluster_obj,inputExpData,n_clusters=NULL,tol_level=0.9){
  
  #argList=.ArgList;cluster_obj=.tmp;n_clusters=NULL;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;inputExpData=.xpo_data
  
  #cluster_obj: output of .sconline.clusteringFn()
  #n_clusters: number of desired clusters
  #sig1_thr: z-score threshold
  #pct2_thr: marker pct.2 threshold (% expressed in other clusters)
  #pct_diff_thr: pct.1 - pct.2 threshold
  #only.pos.logFC: report only upregulated genes in the cluster
  
  {
    require(ggraph)
    library(igraph)
    require(tidyverse)
    require(ape)
    require(seriation)
  } 
  
  if(is.null(n_clusters)){
    n_clusters=cluster_obj$n_clusters
  }
  if(is.null(n_clusters)){
    stop("Number of desired clusters should be provided!")
  }
  
  diff_clust=cluster_obj$cluster_object
  d_conMat=cutree(diff_clust,k=n_clusters)
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
  prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  #prop anno
  
  colMax_vals_m=qlcMatrix::colMax(prop_merged)
  colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
  prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
  prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
  prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
  prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
  
  prop_m_hardCluster$cluster=paste0("C",prop_m_hardCluster$i)
  prop_m_hardCluster$cell=row.names(pd)[prop_m_hardCluster$j]
  prop_m_hardCluster=prop_m_hardCluster[,c("cluster","cell")]
  
  if(class(inputExpData)!="Seurat"){
    inputExpData=.extraExport2SeuratFn(inputExpData)
  }
  inputExpData=Seurat::NormalizeData(inputExpData,verbose=F)
  inputExpData=inputExpData[,colnames(inputExpData) %in% prop_m_hardCluster$cell]
  if(length(setdiff(prop_m_hardCluster$cell,colnames(inputExpData)))>0){
    warning("Expression data for some of the cells were not found in the inputExpData!")
  }
  
  res=list()
  for(iclust in unique(prop_m_hardCluster$cluster)){
    tmp=.myEvalMarkers(object=inputExpData, cells.1=prop_m_hardCluster$cell[which(prop_m_hardCluster$cluster==iclust)], cells.2=prop_m_hardCluster$cell[which(prop_m_hardCluster$cluster!=iclust)], slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cluster_name=iclust,cells.1.weight.col=NULL)
    tmp$gene=row.names(tmp)
    res=c(res,list(tmp))
  }
  
  res=do.call("rbind",res)
  
  return(list(markers=res,clusters=prop_m_hardCluster))
}

.extra_sconline.PseudobulkFn_archive=function(inputData,colName,mode="sum",min_size_limit=20,ncores=1){
  
  #colName: the column in the meta/pheno data that specifies the pseudobulk level.
  #colName: usually a column that is a combination of subjectId/libraryId + cellType
  #min_size_limit: the minimum acceptable size of the pseudobulk data.
  #min_size_limit: pseudobulks with cells less than this threshold are excluded from the analysis 
  
  pdlist=split(as.data.frame(colData(inputData)),colData(inputData)[,colName])
  res=matrix(0,nrow=nrow(inputData),ncol=length(pdlist))
  colnames(res)=names(pdlist)
  if(mode=="sum"){
    
    res=parallel::mclapply(1:length(pdlist),function(x,inputData,pdlist) {
      tmp=inputData[,colnames(inputData) %in% row.names(pdlist[[x]])]
      rowSums(counts(tmp))
    },inputData=inputData,pdlist=pdlist,mc.cores = ncores)
    res=do.call("cbind",res)
    colnames(res)=names(pdlist)
  } else if (mode=="mean"){
    for(i in 1:length(pdlist)){
      tmp=inputData[,colnames(inputData) %in% row.names(pdlist[[i]])]
      res[,i]=rowMeans(counts(tmp))
    }
  }
  
  
  res_pd=as.data.frame(colData(inputData))
  res_pd=res_pd[!duplicated(res_pd[,colName]),]
  row.names(res_pd)=res_pd[,colName]
  res_pd=res_pd[match(colnames(res),row.names(res_pd)),]
  res=SingleCellExperiment(assays = list(counts = res),colData = res_pd,rowData=as.data.frame(rowData(inputData)))
  
  res$QC_Gene_total_count=apply(counts(res),2,sum)
  res$QC_Gene_unique_count=apply(counts(res),2,function(x) sum(x>0))
  lib_sizes=as.data.frame(table(colData(inputData)[,colName]))
  lib_sizes=lib_sizes[match(colData(res)[,colName],lib_sizes[,1]),]
  res$pseudocell_size=lib_sizes[,2]
  if(!is.null(min_size_limit)){
    res=res[,which(res$pseudocell_size>=min_size_limit)]
  }
  
  return(res)
  
}

.extra_sconline.PseudobulkFn=function(inputData,colName,mode="sum",min_size_limit=20,ncores=1,cols_to_sum=NULL){
  
  #colName: the column in the meta/pheno data that specifies the pseudobulk level.
  #colName: usually a column that is a combination of subjectId/libraryId + cellType
  #min_size_limit: the minimum acceptable size of the pseudobulk data.
  #min_size_limit: pseudobulks with cells less than this threshold are excluded from the analysis 
  
  #pdlist=split(as.data.frame(colData(inputData)),colData(inputData)[,colName])
  #res=matrix(0,nrow=nrow(inputData),ncol=length(pdlist))
  #colnames(res)=names(pdlist)
  
  require(Matrix)
  
  design_mat=as.matrix(.myOneHotFn(colData(inputData)[,colName]))
  if(!is.null(min_size_limit)){
    design_mat=design_mat[,colSums(design_mat)>=min_size_limit,drop=F]
  }
  
  design_mat=t(design_mat)
  design_mat=as(design_mat,"dgCMatrix")
  
  if(mode=="sum"){
    agg_mat=design_mat %*% t(counts(inputData))
    agg_mat=t(agg_mat)
    
  } else if (mode=="mean"){
    design_mat=Matrix::diag(x=1/rowSums(design_mat)) %*% design_mat
    agg_mat=design_mat %*% t(counts(inputData))
    agg_mat=t(agg_mat)
  }
  
  if(is.null(cols_to_sum)&sum(colnames(colData(inputData)) %in% cols_to_sum)==0){
    res_pd=as.data.frame(colData(inputData))
    res_pd=res_pd[!duplicated(res_pd[,colName]),]
    row.names(res_pd)=res_pd[,colName]
    res_pd=res_pd[match(colnames(agg_mat),row.names(res_pd)),]
    
  } else {
    if(length(setdiff(cols_to_sum,colnames(colData(inputData))))>0){
      warning("some of provided cols_to_sum cols were not identified in the dataset")
    }
    res_pd=as.data.frame(colData(inputData)[,!colnames(colData(inputData)) %in% cols_to_sum])
    res_pd=res_pd[!duplicated(res_pd[,colName]),]
    row.names(res_pd)=res_pd[,colName]
    res_pd=res_pd[match(colnames(agg_mat),row.names(res_pd)),]
    
    sums_res=design_mat %*% as.matrix(as.data.frame(colData(inputData)[,cols_to_sum]))
    if(any(row.names(sums_res)!=row.names(res_pd),na.rm = F)){
      print("Error in summing the cols!")
    }
    res_pd=cbind(res_pd,sums_res)
  }
  
  
  res=SingleCellExperiment(assays = list(counts = agg_mat),colData = res_pd,rowData=as.data.frame(rowData(inputData)))
  
  res$QC_Gene_total_count=apply(counts(res),2,sum)
  res$QC_Gene_unique_count=apply(counts(res),2,function(x) sum(x>0))
  lib_sizes=as.data.frame(table(colData(inputData)[,colName]))
  lib_sizes=lib_sizes[match(colData(res)[,colName],lib_sizes[,1]),]
  res$pseudocell_size=lib_sizes[,2]
  
  return(res)
  
}

.sconline.Vis_pseuodcell_piechart=function(argList,annoCol,pie_scale=1,tol_level=0.9,return_obj=F){
  
  #argList=.ArgList;cluster_obj=.tmp;n_clusters=NULL
  
  {
    require(ggraph)
    library(igraph)
    require(tidyverse)
    require(ape)
    require(seriation)
  } 
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  if(is.numeric(pd[,annoCol])){
    x=matrix(pd[,annoCol],ncol=1)
  } else {
    x=as.matrix(.myOneHotFn(inputVector=pd[,annoCol]))
  }
  
  if(sum(is.na(pd$UMAP_1))>0|sum(colnames(pd)=="UMAP_1")==0){
    warning("UMAP coordinates are missing for some annotations. there might an issue in the data!")
  }
  res=prop_mat %*% x
  
  res=as.data.frame(res)
  
  res$pseudocell=row.names(prop_mat)
  pd$cluster_anno_res=pd[,annoCol]
  
  
  #inputData=res;argList=argList;min_effective_size=5;pd=pd;cell_annoCol="cluster_anno_res";pie_scale=1
  p=.extra_sconline.visPseudocellAnno_cluster(inputData=res,argList=argList,min_effective_size=5,pd=pd,cell_annoCol="cluster_anno_res",pie_scale=pie_scale,return_obj=return_obj)
  
  
  
  return(p)
}

.sconline.Vis_attribute=function(argList,attribute_col,inputPd=NULL,pie_scale=1,tol_level=0.9){
  
  #argList=.ArgList;cluster_obj=.tmp;n_clusters=NULL
  require(ggplot2)
  require(scatterpie)
  require(hues)
  
  UMAP_centroid=.sconline.fetch_data("umap_pseudocells",argList=argList)
  if(is.null(inputPd)){
    pd=.sconline.fetch_data("annotation",argList=argList)
  } else {
    pd=inputPd
  }
  
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  
  piechart.plot=T
  if(is.numeric(pd[,attribute_col])){
    x=matrix(pd[,attribute_col],ncol=1)
    colnames(x)=attribute_col
    piechart.plot=F
  } else {
    x=as.matrix(.myOneHotFn(inputVector=pd[,attribute_col]))
    pd$cluster_anno_res=pd[,attribute_col]
  }
  
  if(sum(is.na(pd$UMAP_1))>0|sum(colnames(pd)=="UMAP_1")==0){
    warning("UMAP coordinates are missing for some annotations. there might an issue in the data!")
  }
  res=prop_mat %*% x
  
  res=as.data.frame(res)
  
  res$pseudocell=row.names(prop_mat)
  
  pd_summary=pd
  
  
  piechart_data=merge(res,UMAP_centroid,by.x="pseudocell",by.y="centroid",all.x=T)
  
  anno_cols=setdiff(colnames(res),c("pseudocell","effective_size","dataset"))
  
  if(piechart.plot){
    color_num=length(anno_cols)
    p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                              cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(length(anno_cols))))
  } else {
    p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_point(aes_string("UMAP_1", "UMAP_2", fill=anno_cols),color="black", data=piechart_data,shape=21) + coord_equal()+theme_classic()+scale_fill_gradient(low = "white",high="red")
  }
  
  return(p)
}

.sconline.cluster.Vis=function(argList,cluster_obj,attribute_col=NULL,n_clusters=NULL,tol_level=0.9,return_obj=F,prune_clusters=F,inputExpData=NULL,combinatorial_pct_tol=1,marker_sig1_thr=3,marker_pct2_thr=0.3,marker_pct_diff_thr=0.2,forgiveness_factor=1){
  
  #argList=.ArgList;cluster_obj=.tmp;n_clusters=NULL;
  #argList=.ArgList;cluster_obj=cluster_res_avg;attribute_col=NULL;n_clusters=10;tol_level=0.9;prune_clusters=T;inputExpData=data;combinatorial_pct_tol=0.1;marker_sig1_thr=3;marker_pct2_thr=0.3;marker_pct_diff_thr=0.2
  
  require(ggraph)
  library(igraph)
  require(tidyverse)
  require(ape)
  require(seriation)
  require(cowplot)
  
  
  if(is.null(n_clusters)){
    n_clusters=cluster_obj$n_clusters
  }
  if(is.null(n_clusters)){
    stop("Number of desired clusters should be provided!")
  }
  
  if(is.null(inputExpData)&prune_clusters){
    stop("inputExpData needs to be provided for cluster pruning")
  }
  
  diff_clust=cluster_obj$cluster_object
  
  d_conMat=cutree(diff_clust,k=n_clusters)
  
  if(prune_clusters){
    #cluster_assignments=d_conMat;clust_obj=diff_clust;inputExpData=inputExpData;argList=argList;combinatorial_pct_tol=combinatorial_pct_tol;marker_sig1_thr=marker_sig1_thr;marker_pct2_thr=marker_pct2_thr;marker_pct_diff_thr=marker_pct_diff_thr
    d_conMat=.sconline.cluster_pruning(cluster_assignments=d_conMat,clust_obj=diff_clust,
                                       inputExpData=inputExpData,
                                       argList=argList,
                                       combinatorial_pct_tol=combinatorial_pct_tol,marker_sig1_thr=marker_sig1_thr,marker_pct2_thr=marker_pct2_thr,marker_pct_diff_thr=marker_pct_diff_thr,forgiveness_factor=forgiveness_factor)
  }
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  {
    prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
    if(sum(colnames(pd)=="sample")==1&length(setdiff(pd$sample,colnames(prop_mat)))==0){
      prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample,drop=F]
    } else {
      prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),row.names(pd)]
    }
    
    prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
    #prop anno
    
    colMax_vals_m=qlcMatrix::colMax(prop_merged)
    colMax_vals_m=.extra_matrix_colNorm(prop_merged,colValues=1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),,drop=F]
    prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),,drop=F]
    if(!is.null(attribute_col)){
      if(sum(colnames(pd)==attribute_col)==1){
        pd$cluster_anno_res=pd[,attribute_col]
      } else {
        warning(paste("Issue with",attribute_col,"attribute column in meta data -- ignoring it!"))
        pd$cluster_anno_res=paste0("C",as.character(prop_m_hardCluster$i))
      }
      
    } else {
      pd$cluster_anno_res=paste0("C",as.character(prop_m_hardCluster$i))
    }
    
    res=data.frame(pseudocell=names(d_conMat),cluster=paste0("C",d_conMat),stringsAsFactors=F)
  }
  
  #inputData=res;argList=argList;min_effective_size=5;pd=pd;cell_annoCol="cluster_anno_res";pie_scale=1
  p=.extra_sconline.visPseudocellAnno_cluster(inputData=res,argList=argList,min_effective_size=5,pd=pd,cell_annoCol="cluster_anno_res",return_obj=return_obj)
  
  
  
  return(p)
}

.sconline.cluster.Vis2=function(argList,cluster_obj,attribute_col=NULL,n_clusters=NULL,tol_level=0.9,return_obj=F){
  
  #argList=.ArgList;cluster_obj=cluster_res_avg;n_clusters=7;attribute_col=NULL;tol_level=0.9;return_obj=F
  
  require(ggraph)
  library(igraph)
  require(tidyverse)
  require(ape)
  require(seriation)
  require(cowplot)
  
  
  if(is.null(n_clusters)){
    n_clusters=cluster_obj$n_clusters
  }
  if(is.null(n_clusters)){
    stop("Number of desired clusters should be provided!")
  }
  
  
  d_conMat=SNFtool::spectralClustering(cluster_obj$affinity, K=n_clusters, type = 3)
  names(d_conMat)=row.names(cluster_obj$affinity)
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  {
    prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
    if(sum(colnames(pd)=="sample")==1&length(setdiff(pd$sample,colnames(prop_mat)))==0){
      prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
    } else {
      prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),row.names(pd)]
    }
    
    prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
    #prop anno
    
    colMax_vals_m=qlcMatrix::colMax(prop_merged)
    colMax_vals_m=.extra_matrix_colNorm(prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
    prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
    if(!is.null(attribute_col)){
      if(sum(colnames(pd)==attribute_col)==1){
        pd$cluster_anno_res=pd[,attribute_col]
      } else {
        warning(paste("Issue with",attribute_col,"attribute column in meta data -- ignoring it!"))
        pd$cluster_anno_res=paste0("C",as.character(prop_m_hardCluster$i))
      }
      
    } else {
      pd$cluster_anno_res=paste0("C",as.character(prop_m_hardCluster$i))
    }
    
    res=data.frame(pseudocell=names(d_conMat),cluster=paste0("C",d_conMat),stringsAsFactors=F)
  }
  
  #inputData=res;argList=argList;min_effective_size=5;pd=pd;cell_annoCol="cluster_anno_res";pie_scale=1
  p=.extra_sconline.visPseudocellAnno_cluster(inputData=res,argList=argList,min_effective_size=5,pd=pd,cell_annoCol="cluster_anno_res",return_obj=return_obj)
  
  
  
  return(p)
}


.sconline.clusteringVis_archive=function(argList,cluster_obj,annoCol=NULL,n_clusters=NULL,pie_scale=1,tol_level=0.9){
  
  #argList=.ArgList;cluster_obj=.tmp;n_clusters=NULL
  
  {
    require(ggraph)
    library(igraph)
    require(tidyverse)
    require(ape)
    require(seriation)
  } 
  
  if(is.null(n_clusters)){
    n_clusters=cluster_obj$n_clusters
  }
  if(is.null(n_clusters)){
    stop("Number of desired clusters should be provided!")
  }
  
  
  diff_clust=cluster_obj$cluster_object
  
  d_conMat=cutree(diff_clust,k=n_clusters)
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  if(!is.null(annoCol)){
    
    
    
    if(is.numeric(pd[,annoCol])){
      x=matrix(pd[,annoCol],ncol=1)
    } else {
      x=as.matrix(.myOneHotFn(inputVector=pd[,annoCol]))
    }
    
    if(sum(is.na(pd$UMAP_1))>0|sum(colnames(pd)=="UMAP_1")==0){
      warning("UMAP coordinates are missing for some annotations. there might an issue in the data!")
    }
    res=prop_mat %*% x
    
    res=as.data.frame(res)
    
    res$pseudocell=row.names(prop_mat)
    pd$cluster_anno_res=pd[,annoCol]
    
    #tmp=d_conMat[match(row.names(res),names(d_conMat))]
    #res$cluster_anno_res=tmp
  } else {
    prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
    prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
    prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
    #prop anno
    
    colMax_vals_m=qlcMatrix::colMax(prop_merged)
    colMax_vals_m=.extra_matrix_colNorm(prop_merged,colValues=1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
    prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
    pd$cluster_anno_res=paste0("C",as.character(prop_m_hardCluster$i))
    res=data.frame(pseudocell=names(d_conMat),cluster=paste0("C",d_conMat),stringsAsFactors=F)
  }
  
  #inputData=res;argList=argList;min_effective_size=5;pd=pd;cell_annoCol="cluster_anno_res";pie_scale=1
  p=.extra_sconline.visPseudocellAnno_cluster(inputData=res,argList=argList,min_effective_size=5,pd=pd,cell_annoCol="cluster_anno_res",pie_scale=pie_scale)
  
  
  
  return(p)
}


.sconline.MASCfn<- function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL,
                          verbose = FALSE, jackknife=F,statistical.test="Wald") {
  
  #Adapted from Fonseka et al. PMID: 30333237
  
  # Check inputs
  require(lme4)
  if (is.factor(dataset[[contrast]]) == FALSE&is.numeric(dataset[[contrast]])==FALSE) {
    stop("Specified contrast term should be coded as a factor or numeric in the dataset")
  }
  
  match.arg(statistical.test,c("LRT","Wald"))
  
  
  # Convert cluster assignments to string
  cluster <- as.character(cluster)
  # Prepend design matrix generated from cluster assignments
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  # Create output list to hold results
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]
  
  # Create model formulas
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                          paste0("(1|", random_effects, ")", collapse = " + ")),
                        collapse = " + ")
    if (verbose == TRUE&statistical.test=="LRT") {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else if (!is.null(fixed_effects) && is.null(random_effects)) {
    model_rhs <- paste0(fixed_effects, collapse = " + ")
    if (verbose == TRUE&statistical.test=="LRT") {
      message(paste("Using null model:", "cluster ~", model_rhs))
      # For now, do not allow models without mixed effects terms
      
    }
    stop("No random effects specified")
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
    if (verbose == TRUE&statistical.test=="LRT") {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else {
    model_rhs <- "1" # only includes intercept
    if (verbose == TRUE&statistical.test=="LRT") {
      message(paste("Using null model:", "cluster ~", model_rhs))
      
    }
    stop("No random or fixed effects specified")
  }
  
  # Initialize list to store model objects for each cluster
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]
  
  # Run nested mixed-effects models for each cluster
  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    if (verbose == TRUE) {
      message(paste("Creating logistic mixed models for", test_cluster))
    }
    null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                   model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                   model_rhs), collapse = ""))
    # Run null and full mixed-effects models
    
    full_model <- lme4::glmer(formula = full_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    
    
    # calculate confidence intervals for contrast term beta
    if(is.factor(dataset[[contrast]])){
      contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
      
    } else {
      contrast_lvl2 <- contrast
    }
    
    contrast_ci <- confint.merMod(full_model, method = "Wald",
                                  parm = contrast_lvl2)
    
    
    if(statistical.test=="Wald"){
      pval=summary(full_model)
      
      pval=pval$coefficients[contrast_lvl2,4]
    } else {
      null_model <- lme4::glmer(formula = null_fm, data = dataset,
                                family = binomial, nAGQ = 1, verbose = 0,
                                control = glmerControl(optimizer = "bobyqa"))
      model_lrt <- anova(null_model, full_model)
      pval=model_lrt[["Pr(>Chisq)"]][2]
    }
    
    # Save model objects to list
    cluster_models[[i]]$confint <- contrast_ci
    cluster_models[[i]]$pval <- pval
    cluster_models[[i]]$full_model <- full_model
    #jackknifing
    jk_pvalvec=c()
    jk_coefvec=c()
    jk_stable=1
    if(jackknife){
      for(ibatch in unique(dataset[,random_effects])){
        tmp_dataset=dataset[which(dataset[,random_effects]!=ibatch),]
        
        jk_full_model <- tryCatch({lme4::glmer(formula = full_fm, data = tmp_dataset,
                                               family = binomial, nAGQ = 1, verbose = 0,
                                               control = glmerControl(optimizer = "bobyqa"))},error=function(e){return(F)})
        
        if(class(jk_full_model)!=class(T)){
          jk_coefvec=c(jk_coefvec,fixef(jk_full_model)[[contrast_lvl2]])
          if(statistical.test=="Wald"){
            tmp_pval=summary(jk_full_model)
            tmp_pval=tmp_pval$coefficients[contrast_lvl2,4]
            jk_pvalvec=c(jk_pvalvec,tmp_pval)
          } else {
            jk_null_model <- tryCatch({lme4::glmer(formula = null_fm, data = tmp_dataset,
                                                   family = binomial, nAGQ = 1, verbose = 0,
                                                   control = glmerControl(optimizer = "bobyqa"))},error=function(e) {return(F)})
            
            if(class(jk_null_model)!=class(T)){
              jk_model_lrt <- anova(jk_null_model, jk_full_model)
              # calculate confidence intervals for contrast term beta
              jk_pvalvec=c(jk_pvalvec,jk_model_lrt[["Pr(>Chisq)"]][2])
            } else {
              jk_stable=0
            }
          }
          
        } else {
          jk_stable=0
        }
        
      }
    } else {
      jk_pvalvec=(-1)
      jk_coefvec=(-1)
    }
    
    cluster_models[[i]]$jk_pval_median <- median(jk_pvalvec)
    cluster_models[[i]]$jk_pval_mean <- mean(jk_pvalvec)
    cluster_models[[i]]$jk_pval_max <- max(jk_pvalvec)
    
    cluster_models[[i]]$jk_coef_median <- median(jk_coefvec)
    cluster_models[[i]]$jk_coef_mean <- mean(jk_coefvec)
    cluster_models[[i]]$jk_stable <- jk_stable
    cluster_models[[i]]$jk_coef_min <- jk_coefvec[which(abs(jk_coefvec)==min(abs(jk_coefvec)))[1]]
  }
  
  # Organize results into output dataframe
  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
                       size = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$pval)
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))
  output[[paste(contrast_lvl2,"JK","Min", "OR", sep = ".")]] <- sapply(cluster_models, function(x) {if(x$jk_coef_min==(-1)){-1} else{exp(x$jk_coef_min)}})
  output[[paste(contrast_lvl2,"JK","Mean", "OR", sep = ".")]] <- sapply(cluster_models, function(x) {if(x$jk_coef_mean==(-1)){-1} else {exp(x$jk_coef_mean)}})
  output[[paste(contrast_lvl2,"JK","Median", "OR", sep = ".")]] <- sapply(cluster_models, function(x) {if(x$jk_coef_median==(-1)){-1} else {exp(x$jk_coef_median)}})
  output[[paste(contrast_lvl2,"JK","Max", "pvalue", sep = ".")]] <- sapply(cluster_models, function(x) {if(x$jk_pval_max==(-1)){-1} else {x$jk_pval_max}})
  output[[paste(contrast_lvl2,"JK","Mean", "pvalue", sep = ".")]] <- sapply(cluster_models, function(x) {if(x$jk_pval_mean==(-1)){-1} else {x$jk_pval_mean}})
  output[[paste(contrast_lvl2,"JK","Median", "pvalue", sep = ".")]] <- sapply(cluster_models, function(x) {if(x$jk_pval_median==(-1)){-1} else {x$jk_pval_median}})
  output[[paste(contrast_lvl2,"JK","Stable", sep = ".")]] <- sapply(cluster_models, function(x) x$jk_stable)
  
  return(output)
}

.sconline.MASC=function(argList,contrast,random_effects,fixed_effects=NULL,test_condition=NULL,ctrl_condition=NULL,n_clusters=50,tol_level=0.9){
  #argList=.ArgList;contrast="status";random_effects="anno_batch";fixed_effects="anno_sex";test_condition="Abeta";ctrl_condition="Ctrl";n_clusters=50;tol_level=0.9
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  cluster_obj=.sconline.cluster(argList,n_clusters=10,clustering_method="average")
  diff_clust=cluster_obj$cluster_object
  
  d_conMat=cutree(diff_clust,k=n_clusters)
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  {
    prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
    if(sum(colnames(pd)=="sample")==1&length(setdiff(pd$sample,colnames(prop_mat)))==0){
      prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
    } else {
      prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),row.names(pd)]
    }
    
    prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
    #prop anno
    
    colMax_vals_m=qlcMatrix::colMax(prop_merged)
    colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_merged,colValues=1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
    
  }
  
  prop_mat=prop_m_hardCluster
  
  if(is.null(test_condition)){
    if(is.null(ctrl_condition)){
      test_condition=setdiff(unique(pd[,contrast]),pd[1,contrast])
    } else {
      test_condition=setdiff(unique(pd[,contrast]),ctrl_condition)
    }
  }
  
  if(is.null(ctrl_condition)){
    ctrl_condition=setdiff(unique(pd[,contrast]),test_condition)
    warning(paste(ctrl_condition,"was chosen as control condition!"))
  }
  
  res_list=list()
  for(icontrast in test_condition){
    tmp_res=NULL
    tmp_pd=pd
    tmp_pd[,contrast]=as.character(tmp_pd[,contrast])
    tmp_pd=tmp_pd[which(tmp_pd[,contrast] %in% c(icontrast,ctrl_condition)),]
    tmp_pd[tmp_pd[,contrast]!=icontrast,contrast]="other"
    tmp_pd[,contrast]=factor(tmp_pd[,contrast],levels=c("other",icontrast))
    tmp_prop_mat=prop_mat[,row.names(tmp_pd)]
    for(i in 1:nrow(tmp_prop_mat)){
      #print(i)
      clust_anno=rep("other",ncol(tmp_prop_mat))
      clust_anno[tmp_prop_mat[i,]>0]="cluster"
      tmp=.sconline.MASCfn(dataset=tmp_pd, cluster=clust_anno, contrast=contrast, random_effects = random_effects, fixed_effects = fixed_effects,
                           verbose = FALSE, jackknife=F,statistical.test="Wald")
      tmp=tmp[,!grepl("\\.JK",colnames(tmp))]
      tmp=tmp[!grepl("other",tmp$cluster),]
      tmp$contrast=icontrast
      tmp$cluster=row.names(tmp_prop_mat)[i]
      tmp_res=rbind(tmp_res,tmp)
    }
    
    prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
    prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x = 1 / (rowSums(prop_merged)+0.000000000001)) %*% prop_merged
    centroid_umap=.sconline.fetch_data("umap_pseudocells",argList = argList)
    centroid_umap=prop_merged %*% as.matrix(centroid_umap[match(colnames(prop_merged),centroid_umap$centroid),c("UMAP_1","UMAP_2")])
    row.names(centroid_umap)=row.names(prop_merged)
    centroid_umap=centroid_umap[match(tmp_res$cluster,row.names(centroid_umap)),]
    tmp_res=cbind(tmp_res,centroid_umap)
    or_col=colnames(tmp_res)
    or_col=or_col[grepl("\\.OR",or_col)&(!grepl("OR\\.",or_col))]
    tmp_res$zscore=qnorm(tmp_res$model.pvalue/2,lower.tail = F)*sign(log2(tmp_res[,or_col]))
    
    tmp_p=p=ggplot(tmp_res,aes(UMAP_1,UMAP_2))+geom_density2d(data=tmp_pd,color="gray10")+
      geom_point(aes(size=abs(zscore)^2,fill=zscore),color="black",shape=21)+scale_fill_gradient2(low="steelblue",high="red",mid="white")+
      cowplot::theme_cowplot()+scale_size(range = c(1,3))
    tmp_p2=.sconline.cluster.Vis(cluster_obj = cluster_obj,argList = argList,n_clusters = n_clusters)
    tmp_p=tmp_p+tmp_p2
    tmp_res=list(data=tmp_res,plot=tmp_p)
    
    res_list=c(res_list,list(tmp_res))
    names(res_list)[length(res_list)]=icontrast
  }
  
  return(res_list)
}


.sconline.convertSeuratToExpressionSet=function(object){
  
  if (!.hasSlot(object, "version")){
    stop("We are unable to convert Seurat objects less than version 2.X", 
         "Please use devtools::install_version to install Seurat v2.3.4 and update your object to a 2.X object", 
         call. = FALSE)
  }
  
  if (slot(object = object, name = "version") >= package_version(x = "2.0.0") && 
      slot(object = object, name = "version") < package_version(x = "3.0.0")) {
    
    counts= object@raw.data
    expData=.extraExport2ExpressionSetFn(counts=counts,pd=as.data.frame(object@meta.data),fd=NULL)
  } else {
    expData=.extraExport2ExpressionSetFn(counts=object@assays$RNA@counts,pd=as.data.frame(object@meta.data),fd=as.data.frame(object@assays$RNA@meta.features))
  }
  
  return(expData)
}


.sconline.PseudobulkPctfn=function(argList,n_clusters=NULL, parsing.col.names = c("anno_batch"), use.sconline.cluster4parsing=T,cluster_obj=NULL,
                                        inputExpData=NULL,min_size_limit=20,inputPhenoData=NULL,inputEmbedding=NULL,tol_level=0.9,use.sconline.embeddings=F,nPCs=NULL,ncores=5){
  
  
  #parsing.col.names: The columns in the pheonData that will be used to parse the expression data and generate the pseudocell/pseudobulk data
  #use.sconline.cluster4parsing: use the sconline cluster as factor for the parsing
  #inputPhenoData: in case we want to run the function outside sconline space
  #min_size_limit: minimum acceptable size (ie, #cells) for each pseudobulk
  #nPCs: the dimension of the embedding space for the construction of pseudobulk data
  #inputEmbedding: the embedding space to be used for the generation of the pseudobulk. only needed when pseudocell.size is not null
  #use.sconline.embeddings: use the embedding space used for sconline to generate the pseudobulk samples. Only needed when the pseudocell.size is not set to null
  
  #Function adds/modifies three annotation columns: pseuodcell_size, QC_Gene_total_count, QC_Gene_unique_count
  #QC_Gene_total_count: equivalant to nUMI for the pseudobulk samples
  #QC_Gene_unique_count: equivalant to nGene for the pseudobulk samples
  #use scale(tst$QC_Gene_total_count) and scale(tst$pseudocell_size) as additional covariates for the DE analysis
  
  require(Matrix)
  
  if(is.null(inputPhenoData)&is.null(argList)){
    stop("Both argList and inputPhenoData cannot be null!")
  }
  
  
  if(is.null(argList)){
    use.sconline.cluster4parsing=F
    use.sconline.embeddings=F
  }
  
  
  if(is.null(inputPhenoData)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPhenoData
  }
  
  
  if(use.sconline.cluster4parsing){
    parsing.col.names=c(parsing.col.names,"cluster_anno_res")
    prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
    
    if(!is.null(n_clusters)){
      cat(paste0("Analysis based on ",n_clusters," clusters"))
      if(is.null(cluster_obj)){
        stop("Cluster object needs to be provided!")
      }
      diff_clust=cluster_obj$cluster_object
      d_conMat=cutree(diff_clust,k=n_clusters)
      prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
      prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
      prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
      #prop anno
      
      
    } else {
      cat(paste0("Analysis at the pseudocell level"))
      prop_merged=prop_mat
    }
    
    colMax_vals_m=qlcMatrix::colMax(prop_merged)
    colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
    prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
    pd$cluster_anno_res=paste0("C",as.character(prop_m_hardCluster$i))
  }
  
  if(is.null(inputExpData)){
    if(!file.exists(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))){
      stop("Expression data is missing!")
    }
    inputExpData=qread(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
    inputExpData=.extraExport2ExpressionSetFn(counts=inputExpData@assays$RNA@counts,pd=as.data.frame(inputExpData@meta.data),fd=as.data.frame(inputExpData@assays$RNA@meta.features))
  } else if(tolower(class(inputExpData))=="seurat"){
    inputExpData=.sconline.convertSeuratToExpressionSet(object=inputExpData)
  } else if(class(inputExpData)!="SingleCellExperiment") {
    stop("Unrecognized inputExpression data!")
  }
  
  inputExpData=inputExpData[,colnames(inputExpData) %in% row.names(pd)]
  if(ncol(inputExpData)==0){
    stop("Expression data doesn't match with the phenoData!")
  }
  pd=pd[match(colnames(inputExpData), row.names(pd)),]
  
  for(icol in parsing.col.names){
    if(sum(colnames(colData(inputExpData))==icol)==0){
      if(sum(colnames(pd)==icol)==0){
        stop(paste0(icol," column was not identified!"))
      } else {
        if(class(pd[,icol])[1]==class(factor())){
          colData(inputExpData)[,icol]=as.character(pd[,icol])
        } else {
          colData(inputExpData)[,icol]=pd[,icol]
        }
        
      }
    }
  }
  
  
  inputExpData$lib_anno=apply(as.data.frame(colData(inputExpData)[,parsing.col.names]),1,function(x)paste(x,collapse="_"))
  
  
  design_mat=as.matrix(.myOneHotFn(colData(inputExpData)$lib_anno))
  if(!is.null(min_size_limit)){
    design_mat=design_mat[,colSums(design_mat)>=min_size_limit]
  }
  
  design_mat=t(design_mat)
  design_mat=as(design_mat,"dgCMatrix")
  
  counts(inputExpData)@x=rep(1,length(counts(inputExpData)@x))
  agg_mat=design_mat %*% t(counts(inputExpData))
  agg_mat=t(agg_mat)
  
  agg_mat=sweep(agg_mat,2,rowSums(design_mat),"/")  
  
  
  res_pd=as.data.frame(colData(inputExpData))
  res_pd=res_pd[!duplicated(res_pd$lib_anno),]
  row.names(res_pd)=res_pd$lib_anno
  res_pd=res_pd[match(colnames(agg_mat),row.names(res_pd)),]
  res=SingleCellExperiment(assays = list(counts = agg_mat),colData = res_pd,rowData=as.data.frame(rowData(inputExpData)))
  
  res$QC_Gene_total_count=NA
  res$QC_Gene_unique_count=apply(counts(res),2,function(x) sum(x>0))
  lib_sizes=as.data.frame(table(colData(inputExpData)$lib_anno))
  lib_sizes=lib_sizes[match(colData(res)$lib_anno,lib_sizes[,1]),]
  res$pseudocell_size=lib_sizes[,2]
  
  
  
  return(res)
}


.sconline.PseudobulkGeneration=function(argList=NULL,n_clusters=NULL, parsing.col.names = c("anno_batch"), use.sconline.cluster4parsing=T,cluster_obj=NULL,
                                   pseudocell.size=40,inputExpData=NULL,min_size_limit=20,inputPhenoData=NULL,inputEmbedding=NULL,tol_level=0.9,use.sconline.embeddings=F,nPCs=NULL,ncores=5,calculate.outlier.score=F,rand_pseudobulk_mod=T,organism,cols_to_sum=NULL){
  
  
  #parsing.col.names: The columns in the pheonData that will be used to parse the expression data and generate the pseudocell/pseudobulk data
  #use.sconline.cluster4parsing: use the sconline cluster as factor for the parsing
  #inputPhenoData: in case we want to run the function outside sconline space
  #min_size_limit: minimum acceptable size (ie, #cells) for each pseudobulk
  #pseudocell.size: average pseudocell size.
  #if pseudocell.size=null, function turns to a tranditional pseudobulk method, ie, all cells at each parsing level are combined together
  #nPCs: the dimension of the embedding space for the construction of pseudobulk data
  #inputEmbedding: the embedding space to be used for the generation of the pseudobulk. only needed when pseudocell.size is not null
  #use.sconline.embeddings: use the embedding space used for sconline to generate the pseudobulk samples. Only needed when the pseudocell.size is not set to null
  
  #Function adds/modifies three annotation columns: pseuodcell_size, QC_Gene_total_count, QC_Gene_unique_count
  #QC_Gene_total_count: equivalant to nUMI for the pseudobulk samples
  #QC_Gene_unique_count: equivalant to nGene for the pseudobulk samples
  #use scale(tst$QC_Gene_total_count) and scale(tst$pseudocell_size) as additional covariates for the DE analysis
  
  
  if(is.null(inputPhenoData)&is.null(argList)){
    stop("Both argList and inputPhenoData cannot be null!")
  }
  
  if(!is.null(pseudocell.size)){
    if(pseudocell.size<2){
      stop("pseudocell.size cannot be less than 2!")
    }
  }
  
  if(is.null(inputExpData)&!rand_pseudobulk_mod){
    warning("It's advised to provide the inputEmbedding")
  }
  
  if(is.null(nPCs)&!is.null(pseudocell.size)&!rand_pseudobulk_mod){
    if(!is.null(argList)){
      nPCs=argList$nPCs
    } else if(!is.null(inputEmbedding)){
        warning(paste0("Setting the nPCs to ", nrow(inputEmbedding)," based on the inputEmbedding"))
    } else {
      stop("nPCs argument (number of PCs for the generation of embedding) need to be provided")
    }
    
    if(nPCs>pseudocell.size){
      warning(paste0("nPCs larger than pseuodcell.size is not advised. setting nPCs to ",pseudocell.size-5))
      nPCs=pseudocell.size-5
    }
  }
  
  
  
  if(is.null(argList)){
    use.sconline.cluster4parsing=F
    use.sconline.embeddings=F
  }
  
  
  if(is.null(inputPhenoData)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPhenoData
  }
  
  
  if(use.sconline.cluster4parsing){
    parsing.col.names=c(parsing.col.names,"cluster_anno_res")
    prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
    
    if(!is.null(n_clusters)){
      cat(paste0("Analysis based on ",n_clusters," clusters"))
      if(is.null(cluster_obj)){
        stop("Cluster object needs to be provided!")
      }
      diff_clust=cluster_obj$cluster_object
      d_conMat=cutree(diff_clust,k=n_clusters)
      prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
      prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
      prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
      #prop anno
      
      
    } else {
      cat(paste0("Analysis at the pseudocell level"))
      prop_merged=prop_mat
    }
    
    colMax_vals_m=qlcMatrix::colMax(prop_merged)
    colMax_vals_m=.extra_matrix_colNorm(prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
    prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
    pd$cluster_anno_res=paste0("C",as.character(prop_m_hardCluster$i))
  }
  
  if(is.null(inputExpData)){
    if(!file.exists(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))){
      stop("Expression data is missing!")
    }
    inputExpData=qread(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
    inputExpData=.extraExport2ExpressionSetFn(counts=inputExpData@assays$RNA@counts,pd=as.data.frame(inputExpData@meta.data),fd=as.data.frame(inputExpData@assays$RNA@meta.features))
  } else if(tolower(class(inputExpData))=="seurat"){
    inputExpData=.sconline.convertSeuratToExpressionSet(object=inputExpData)
  } else if(class(inputExpData)!="SingleCellExperiment") {
    stop("Unrecognized inputExpression data!")
  }
  
  inputExpData=inputExpData[,colnames(inputExpData) %in% row.names(pd)]
  if(ncol(inputExpData)==0){
    stop("Expression data doesn't match with the phenoData!")
  }
  pd=pd[match(colnames(inputExpData), row.names(pd)),]
  
  for(icol in parsing.col.names){
    if(sum(colnames(colData(inputExpData))==icol)==0){
      if(sum(colnames(pd)==icol)==0){
        stop(paste0(icol," column was not identified!"))
      } else {
        if(class(pd[,icol])[1]==class(factor())){
          colData(inputExpData)[,icol]=as.character(pd[,icol])
        } else {
          colData(inputExpData)[,icol]=pd[,icol]
        }
        
      }
    }
  }
  
  if(length(unique(parsing.col.names))==1){
    inputExpData$lib_anno=colData(inputExpData)[,unique(parsing.col.names)]
  } else {
    inputExpData$lib_anno=apply(as.data.frame(colData(inputExpData)[,unique(parsing.col.names)]),1,function(x)paste(x,collapse="_"))
  }
  
  
  if(is.null(pseudocell.size)){
    #inputData=inputExpData;colName="lib_anno";mode="sum";cols_to_sum=cols_to_sum
    sl_data=.extra_sconline.PseudobulkFn(inputData=inputExpData,colName="lib_anno",mode="sum",min_size_limit=min_size_limit,cols_to_sum=cols_to_sum,ncores=ncores)
  } else {
    inputExpList=.mySplitObject_v2(inputExpData,colName="lib_anno",min_dataset_size=min_size_limit,ncores=ncores)
    #inputExpList2=.mySplitObject_v2(inputExpData,colName="library_anno2",min_dataset_size=min_size_limit,ncores=ncores)
    if(sum(unlist(lapply(inputExpList,class))=="SingleCellExperiment")!=length(inputExpList)){
      cat("Error: consider increasing RAM! re-trying with lower number of cores")
      inputExpList=.mySplitObject_v2(inputExpData,colName="lib_anno",min_dataset_size=min_size_limit,ncores=1)
      if(sum(unlist(lapply(inputExpList,class))=="SingleCellExperiment")!=length(inputExpList)){
        stop("Persistent error: increase RAM!")
      }
    }
    if(use.sconline.embeddings&!rand_pseudobulk_mod){
      load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
      inputEmbedding=harmony_embeddings[,1:nPCs,drop=F]
    }
    
    data_size=unlist(lapply(inputExpList,ncol))
    inputExpList=inputExpList[data_size>=min_size_limit]
    data_size=unlist(lapply(inputExpList,ncol))
    
    if(F){
      for(icheck in inputExpList){
        #inputExp=icheck;include.outlier.score=calculate.outlier.score
        tst=.extra_sconline.FixedSizeFn(inputExp=icheck,inputEmbedding=inputEmbedding,pseudocell_size=pseudocell.size,nPCs=nPCs,include.outlier.score=calculate.outlier.score)
      }
    }
    
    
    sl_data=suppressWarnings(parallel::mclapply(inputExpList,.extra_sconline.FixedSizeFn,inputEmbedding=inputEmbedding,pseudocell_size=pseudocell.size,nPCs=nPCs,include.outlier.score=calculate.outlier.score,rand_pseudobulk_mod=rand_pseudobulk_mod,mc.cores = ncores,cols_to_sum=cols_to_sum))
    if(sum(unlist(lapply(sl_data,class))=="SingleCellExperiment")!=length(inputExpList)){
      cat("Error: consider increasing RAM! re-trying with lower number of cores")
      sl_data=suppressWarnings(parallel::mclapply(inputExpList,.extra_sconline.FixedSizeFn,inputEmbedding=inputEmbedding,pseudocell_size=pseudocell.size,nPCs=nPCs,include.outlier.score=calculate.outlier.score,cols_to_sum=cols_to_sum,mc.cores = 1))
      if(sum(unlist(lapply(sl_data,class))=="SingleCellExperiment")!=length(inputExpList)){
        stop("Persistent error: increase RAM!")
      }
    }
    sl_data_size=lapply(sl_data,function(x) max(x$pseudocell_size))
    sl_data=sl_data[sl_data_size>=min_size_limit]
    sl_data=lapply(sl_data,function(x) x[,which(x$pseudocell_size>=min_size_limit),drop=F])
    
    sl_data=.mycBindFn(sl_data)
    sl_data$QC_Gene_total_count=apply(counts(sl_data),2,sum)
    sl_data$QC_Gene_unique_count=apply(counts(sl_data),2,function(x) sum(x>0))
  }
  
  sl_data$QC_MT.pct=.extraMitoPctFn(inputData = sl_data,organism = organism)
  
  
  return(sl_data)
}


.sconline.Pseudobulk10=function(inputExpData,embeddings,pseudobulk_split_col,min_dataset_size=4){
  #creates pseudobulk samples of median size 10
  
  if(is.null(inputExpData)){
    stop("inputExpData should be provided!")
  } else if(tolower(class(inputExpData))=="seurat"){
    inputExpData=.sconline.convertSeuratToExpressionSet(object=inputExpData)
  } else if(class(inputExpData)!="SingleCellExperiment") {
    stop("Unrecognized inputExpression data!")
  }
  
  #data.list=.mySplitObject(inputExpData,pseudobulk_split_col)
  data.list=.mySplitObject_v2(object=inputExpData,colName=pseudobulk_split_col,min_dataset_size=min_dataset_size)
  
  
  embeddings=embeddings[match(colnames(inputExpData),row.names(embeddings)),]
  if(sum(is.na(embeddings))>0){
    stop("Error!")
  }
  
  pseudo_res=list()
  for(i in 1:length(data.list)){
    iclust=unique(colData(data.list[[i]])[,pseudobulk_split_col])
    {
      tmpExp=data.list[[i]]#[,colData(data.list[[i]])[,'pseudobulk_split_col']==iclust]
      #tmpExp=SingleCellExperiment(assays = list(counts = tmpExp),colData = as.data.frame(colData(data.list[[i]]@meta.data[data.list[[i]]@meta.data$cluster_name==iclust,]),rowData=as.data.frame(data.list[[i]]@assays$RNA@meta.features))
      if(length(setdiff(colnames(tmpExp),row.names(embeddings)))>0){
        stop("Error in input embeddings file. Embedding info was not found for some cells!")
      }
      tmp_embeddings=embeddings[match(colnames(tmpExp),row.names(embeddings)),]
      #inputExpData=tmpExp;inputPCAembeddings=tmp_embeddings; n.adaptiveKernel=5; nPropIter=3;nPCs=30;verbose=T
      icounter=0
      runCheck=T
      
      best_coverage=0
      sl_solution=NULL
      
      while(icounter<100&runCheck){
        icounter=icounter+1
        runCheck=F
        if(ncol(tmpExp)<=25&ncol(tmpExp)>13){
          s1=sample(ncol(tmpExp),ncol(tmpExp)/2)
          s2=as.matrix(counts(tmpExp))[,-s1]
          s1=as.matrix(counts(tmpExp))[,s1]
          tmp=matrix(0,nrow=nrow(s1),ncol=2)
          tmp[,1]=rowSums(s1)
          tmp[,2]=rowSums(s2)
          colnames(tmp)=paste0("group_",colnames(tmpExp)[1:2])
          sl_solution=tmp
        } else if(ncol(tmpExp)<=13){
          s1=as.matrix(counts(tmpExp))
          tmp=matrix(rowSums(s1),nrow=nrow(s1),ncol=1)
          colnames(tmp)=colnames(tmpExp)[1]
          colnames(tmp)=paste0("group_",colnames(tmp))
          sl_solution=tmp
        } else {
          
          
          tmp=tryCatch({.myPseudoCellfn(inputExpData=tmpExp,inputPCAembeddings=tmp_embeddings, n.adaptiveKernel=5, nPropIter=3,nPCs=30,verbose=F)}, error=function(e) {return(T)})
          if(class(tmp)==class(T)){
            runCheck=T
          } else {
            if(tmp$coverage>best_coverage){
              sl_solution=tmp$collapsedExpData
              best_coverage=tmp$coverage
            }
            
            if(best_coverage<0.75){
              runCheck=T
            }
          }
        }
        
      }
      
      
      
      if(runCheck&is.null(sl_solution)){
        stop("Error!")
      } else {
        tmp=sl_solution
        if(class(tmp)[1]=="numeric"){
          tmp=matrix(tmp,ncol=1)
          colnames(tmp)=colnames(tmpExp)[1]
        } else if(sum(grepl("^group",colnames(tmp)))==1){
          
          tmp2=matrix(tmp[,grepl("^group",colnames(tmp))],ncol=1)
          colnames(tmp2)=gsub("^group_","",colnames(tmp)[grepl("^group",colnames(tmp))])
          tmp=tmp2
        } else {
          tmp=tmp[,grepl("^group",colnames(tmp))]
          colnames(tmp)=gsub("^group_","",colnames(tmp))
        }
      }
      
      
      
      tmp=SingleCellExperiment(assays = list(counts = tmp),colData = as.data.frame(colData(data.list[[i]])[match(colnames(tmp),row.names(colData(data.list[[i]]))),]),rowData=as.data.frame(rowData(data.list[[i]])))
      
    }
    
    pseudo_res=c(pseudo_res,list(tmp))
    names(pseudo_res)[length(pseudo_res)]=iclust
  }
  
  pseudo_res=.mycBindFn(pseudo_res,batchNames = NULL)
  return(pseudo_res)
}


.sconline.MeanExp=function(argList,inputExpdata,cluster_obj=NULL,n_clusters=NULL,normalization.method = NULL,inputPd=NULL,cluster_col=NULL,return.only.unique.cellTypes=F,cluster_selection_regularExp=NULL){
  #normalization.method: NULL, "LogNormalize", "RC"
  #return.only.unique.cellTypes: filter the clustering results to include only one representative from each cell type.
  #argList=.ArgList;inputExpdata=data;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;affinity_param=5;scaling_factor=10;n_clusters=NULL;normalization.method = "RC"
  
  if(is.null(inputPd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPd
  }
  
  anno=NULL
  x_exp=NULL
  x_exp_binary=NULL
  anno_m=NULL
  
  if(is.null(n_clusters)&is.null(cluster_col)){
    stop("Either n_clusters or cluster_col argument should be provided!")
  } else if(!is.null(n_clusters)&!is.null(cluster_col)){
    stop("Both cluster_col and n_clusters can't be provided at the same time!")
  } else if(!is.null(cluster_col)){
    if(sum(colnames(pd)==cluster_col)==0){
      stop("Provided cluster_col column was not found in the meta data!")
    }
    n_clusters=length(unique(pd[,cluster_col]))
    anno_m=data.frame(cluster=pd[,cluster_col])
    row.names(anno_m)=row.names(pd)
    
    if(!is.null(normalization.method)){
      logNormData=t(Seurat:::NormalizeData.default(counts(inputExpdata),normalization.method = normalization.method,verbose = F))
    } else {
      logNormData=t(counts(inputExpdata))
    }
    
    prop_m_hardCluster=reshape2::dcast(cluster~cell,data=data.frame(cluster=anno_m$cluster,cell=row.names(anno_m),count=1,stringsAsFactors = F),value.var="count",drop=F,fun.aggregate=length)
    row.names(prop_m_hardCluster)=prop_m_hardCluster[,1]
    prop_m_hardCluster=as(as.matrix(prop_m_hardCluster[,-1]),"dgCMatrix")
    
    x_exp_binary=t(prop_m_hardCluster %*% logNormData[colnames(prop_m_hardCluster),])
    
  } else {
    diff_clust=cluster_obj$cluster_object
    d_conMat=cutree(diff_clust,k=n_clusters)
    d_conMat=cutree(diff_clust,k=n_clusters)
    
    prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
    prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
    prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
    prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
    
    #prop anno
    if(sum(colnames(pd)=="anno_orig_cellState")==1){
      anno=prop_merged %*% as.matrix(.myOneHotFn(inputVector=pd$anno_orig_cellState))
      
      if(return.only.unique.cellTypes){
        tmp=apply(anno,1,function(x){
          y=which(x==max(x))[1]
          tmp2=data.frame(cluster=colnames(anno)[y],purity=x[y])
          tmp2
        })
        
        tmp=do.call("rbind",tmp)
        tmp$pseudocell=row.names(tmp)
        tmp=tmp[order(tmp$purity,decreasing = T),]
        tmp=tmp[!duplicated(tmp$cluster),]
        prop_merged=prop_merged[row.names(prop_merged) %in% tmp$pseudocell,]
        anno=prop_merged %*% as.matrix(.myOneHotFn(inputVector=pd$anno_orig_cellState))
      }
      
      
      anno=as.data.frame(anno)
      anno$dominant_type=apply(anno,1,function(x) colnames(anno)[which(x==max(x))])
      anno$ps=row.names(prop_merged)
    }
    
    colMax_vals_m=as.numeric(qlcMatrix::colMax(prop_merged))
    colMax_vals_m[which(colMax_vals_m==0)]=1
    colMax_vals_m=prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=Matrix::drop0(colMax_vals_m,tol=0.95)
    anno_m=colMax_vals_m %*% as.matrix(.myOneHotFn(inputVector=pd$anno_orig_cellState))
    anno_m=.extra_matrix_rowNorm(as.matrix(anno_m))#as.matrix(Matrix::Diagonal(x=1/rowSums(anno_m)) %*% anno_m)
    anno_m=as.data.frame(anno_m)
    anno_m$dominant_type=apply(anno_m,1,function(x) colnames(anno_m)[which(x==max(x))])
    anno_m$ps=row.names(prop_merged)
    
    colMax_vals_m=colMax_vals_m[,colSums(colMax_vals_m)==1]
    colMax_vals_m=as.data.frame(as.matrix(t(colMax_vals_m)))
    colMax_vals_m=apply(colMax_vals_m,1,function(x) colnames(colMax_vals_m)[x==1])
    colMax_vals_m=data.frame(cell=names(colMax_vals_m),cluster=colMax_vals_m,stringsAsFactors = F)
    tmp_pd=pd[match(colMax_vals_m$cell,pd$sample),]
    tmp_pd$anno_cluster_res=colMax_vals_m$cluster
    
    if(!is.null(normalization.method)){
      logNormData=t(Seurat:::NormalizeData.default(counts(inputExpdata),normalization.method = normalization.method,verbose = F))
    } else {
      logNormData=t(counts(inputExpdata))
    }
    
    x_exp=t(prop_merged %*% logNormData[colnames(prop_merged),])
    x_exp_binary=t(prop_m_hardCluster %*% logNormData[colnames(prop_merged),])
    
    
  }
  
  return(list(anno_hardCluster=anno_m,anno_softCluster=anno,meanExp_hardCluster=x_exp_binary,meanExp_softCluster=x_exp))
}


.sconline.markerPlot=function(argList,inputExpData,geneName,pseudocell_size_offset=0.2,pseudocell_size_pwr=1,include_seurat_featurePlot=F){
  #metric: zscore, pct.1
  library(dplyr)
  library(purrr)
  library(cowplot)
  library(patchwork)
  
  if(length(which(toupper(geneName)==toupper(row.names(inputExpData))))==0){
    stop("input gene was not found!")
  }
  
  geneName=row.names(inputExpData)[which(toupper(geneName)==toupper(row.names(inputExpData)))]
  
  meta_data=.sconline.fetch_data("meta_z",argList = argList)
  load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  inputExpData=inputExpData[,colnames(inputExpData) %in% row.names(pd)]
  pd=pd[match(colnames(inputExpData),row.names(pd)),]
  input_cell_bkg_umap=.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"))
  
  p1=NULL
  if(include_seurat_featurePlot){
    if(class(inputExpData)!="Seurat"){
      warning("converting the inputExpData to seurat object speeds up the process!")
      inputExpData=.extraExport2SeuratFn(inputData = inputExpData)
      inputExpData=Seurat::NormalizeData(inputExpData,verbose=F)
      
    }
    
    umap.reduction <- Seurat::CreateDimReducObject(embeddings = as.matrix(pd[,c("UMAP_1","UMAP_2")]), 
                                                   key = "UMAP_", assay = "RNA", global = TRUE)
    
    inputExpData[["umap"]]=umap.reduction
    
    p1=Seurat::FeaturePlot(inputExpData, features=geneName,cols=c("lightgrey","#0D0887FF","#9C179EFF","#ED7953FF","#F0F921FF"))
    #p2=.my2dPlot_counts(inputPCA=pd,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=geneName,geneNameCol=NULL,expData=inputExpData,combine_figs=F)
    
    
  }
  
  zscore_df=data.frame(pseudocell=row.names(meta_data$meta_z),zscore=meta_data$meta_z[,geneName],stringsAsFactors = F)
  pct_df=data.frame(pseudocell=row.names(meta_data$med_pct.1),pct.1=meta_data$med_pct.1[,geneName],stringsAsFactors = F)
  sc_df=merge(zscore_df,pct_df,by="pseudocell")
  sc_df=merge(sc_df,UMAP_centroid,by.x="pseudocell",by.y="centroid",all.x=T)
  sc_df=sc_df[sc_df$zscore!=0,]
  sc_df$scaled_pct.1=(sc_df$pct.1+pseudocell_size_offset)^pseudocell_size_pwr
  
  p3=ggplot(sc_df,aes(UMAP_1,UMAP_2,color=zscore,size=scaled_pct.1))+geom_density2d(data=input_cell_bkg_umap,aes(UMAP_1,UMAP_2),color="lightgrey",size=0.4)+geom_point()+theme_cowplot()+ theme(plot.title = element_text(hjust = 0.5))+ggtitle(geneName)+scale_color_gradientn(colors = c("lightblue","darkblue","black","yellow","orange","red"),breaks=c(min(sc_df$zscore),-2,0,3,5,max(sc_df$zscore)))+scale_size_identity()
  
  if(!is.null(p1)){
    p=p1+p3+ plot_layout(nrow=1,ncol=2)
  } else {
    p=p3
  }
  
  
  
  return(p)
}


.sconline.cluster.DotPlot=function (argList,inputExpData, features, cluster_obj, pseudocell_assignments=NULL,n_clusters=NULL, cols = c("blue","white",
                                                                                      "red"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                            group.by = NULL, split.by = NULL, 
                            scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) {
  #cols = c("blue","white",
  #         "red"); col.min = -2.5; col.max = 2.5; dot.min = 0; dot.scale = 6; 
  #         group.by = NULL; split.by = NULL; 
  #         scale = TRUE; scale.by = "radius"; scale.min = NA; scale.max = NA
  
  #Adapted from Seurat dotPlot function
  require(ggplot2)
  require(cowplot)
  require(Matrix)
  
  if(!is.null(pseudocell_assignments)){
    d_conMat=pseudocell_assignments
  } else {
    if(is.null(n_clusters)){
      n_clusters=cluster_obj$n_clusters
    }
    if(is.null(n_clusters)){
      stop("Number of desired clusters should be provided!")
    }
    
    diff_clust=cluster_obj$cluster_object
    if(is.null(cluster_obj$cluster_object)){
      stop("A proper cluster object needs to be provided")
    }
    d_conMat=cutree(diff_clust,k=n_clusters)
  }
  
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),colnames(inputExpData)]
  prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x = 1 / (rowSums(prop_merged)+0.000000000001)) %*% prop_merged
  
  logNormData=t(Seurat:::NormalizeData.default(counts(inputExpData)[,colnames(prop_mat)],normalization.method = "RC",verbose = F))
  logNormData=logNormData[,colnames(logNormData) %in% features,drop=F]
  exp_binary=logNormData
  exp_binary@x=rep(1,length(exp_binary@x))
  
  avg_exp=log2(prop_merged %*% logNormData+1)
  
  pct.1=prop_merged %*% exp_binary
  
  markers_obj=NULL
  for(icol in 1:ncol(avg_exp)){
    markers_obj=rbind(markers_obj,data.frame(features.plot=colnames(avg_exp)[icol],id=row.names(prop_merged),pct.exp=pct.1[,colnames(avg_exp)[icol]],avg.exp=avg_exp[,icol],stringsAsFactors = F))
  }
  
  matWeights=.myEffSizePropMat(prop_merged)
  matWeights=matWeights$effective_sample_size
  matWeights=matWeights[order(matWeights)]
  markers_obj$id=factor(as.character(markers_obj$id),levels = names(matWeights))
  
  
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  data.plot=markers_obj[markers_obj$features.plot %in% features,]
  if(!is.factor(data.plot$id)){
    data.plot$id=factor(data.plot$id)
  }
  #cluster annotations
  id.levels <- levels(x = data.plot$id)
  data.plot$id <- as.vector(x = data.plot$id)
  
  #data.features$id: cluster ids
  
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  
  
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  
  data.plot=split(data.plot,data.plot$features.plot)
  
  data.plot = lapply(X = data.plot, 
                     FUN = function(x) {
                       data.use =x$avg.exp
                       if (scale) {
                         data.use <- scale(x = data.use)
                         data.use <- MinMax(data = data.use, min = col.min, 
                                            max = col.max)
                       }
                       else {
                         data.use <- log(x = data.use)
                       }
                       x$avg.exp.scaled=as.numeric(data.use)
                       return(x)
                     })
  
  
  data.plot=do.call("rbind",data.plot)
  if (split.colors) {
    data.plot$avg.exp.scaled <- as.numeric(x = cut(x = data.plot$avg.exp.scaled, 
                                                   breaks = 20))
  }
  
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                           fill = color.by),shape=21) + scale.func(range = c(0, dot.scale), 
                                                                                                                                                   limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                                             axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) + theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", 
                              space = "free_x", switch = "y")
  }
  if (split.colors) {
    plot <- plot + scale_fill_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_fill_distiller(palette = cols)
  } else if(length(cols)==2) {
    plot <- plot + scale_fill_gradient(low = cols[1], high = cols[2])
  } else if(length(cols)>2) {
    plot <- plot + scale_fill_gradient2(low = cols[1],mid = cols[2], high = cols[3],midpoint = 0)
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  
  if(length(features)==1&sum(colnames(data.plot)=="class_name")>0&is.null(feature.groups)){
    {
      plot=plot+facet_grid(~class_name,scales="free",space="free")+
        theme(axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1),axis.text.y=element_text(size=rel(0.7)),panel.grid.major.x = element_blank())+ylab("")+xlab("")
    }
    
  } else {
    plot=plot+theme(panel.spacing = unit(x = 1, 
                                         units = "lines"), strip.background = element_blank())
  }
  
  return(plot)
}


.sconline.cluster_validation=function(n_clusters,argList,prune_clusters=F,inputExpData=NULL,cluster_obj,cell_anno_col="anno_orig_cellState", inputPd=NULL,tol_level=0.9,combinatorial_pct_tol=1,marker_sig1_thr=3,marker_pct2_thr=0.3,marker_pct_diff_thr=0.2,forgiveness_factor=1) {
  
  if(is.null(n_clusters)){
    n_clusters=cluster_obj$n_clusters
  }
  if(is.null(n_clusters)){
    stop("Number of desired clusters should be provided!")
  }
  
  
  diff_clust=cluster_obj$cluster_object
  
  d_conMat=cutree(diff_clust,k=n_clusters)
  
  if(prune_clusters){
    #cluster_assignments=d_conMat;clust_obj=diff_clust;tol_level=0.9;combinatorial_pct_tol=0.1;marker_sig1_thr=3;marker_pct2_thr=0.3;marker_pct_diff_thr=0.2
    if(is.null(inputExpData)){
      stop("inputExpData needs to be provided for pruning")
    }
    d_conMat=.sconline.cluster_pruning(cluster_assignments=d_conMat,clust_obj=diff_clust,
                                       inputExpData=inputExpData,
                                       argList=argList,
                                       combinatorial_pct_tol=combinatorial_pct_tol,marker_sig1_thr=marker_sig1_thr,marker_pct2_thr=marker_pct2_thr,marker_pct_diff_thr=marker_pct_diff_thr,forgiveness_factor=forgiveness_factor)
  }
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(is.null(inputPd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPd
  }
  
  pd=pd[!is.na(pd[,cell_anno_col]),]
  prop_mat=prop_mat[,colnames(prop_mat) %in% row.names(pd)]
  pd=pd[match(row.names(pd),colnames(prop_mat)),]
  
  
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
  prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  #prop anno
  
  colMax_vals_m=qlcMatrix::colMax(prop_merged)
  colMax_vals_m=.extra_matrix_colNorm(prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
  prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
  #prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
  prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
  prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
  pd$anno_cluster_res=paste0("C",as.character(prop_m_hardCluster$i))
  
  
  
  ari = mclust::adjustedRandIndex(pd$anno_cluster_res, pd[,cell_anno_col])
  nmi=aricode::NMI(pd$anno_cluster_res, pd[,cell_anno_col],variant="sum")
  ami=aricode::AMI(pd$anno_cluster_res, pd[,cell_anno_col])
  ariScores=data.frame(cluster_count=n_clusters,ari=ari,nmi=nmi,ami=ami,stringsAsFactors = F)
  
  
  tmp_pd=pd
  
  tmp_clust_size=as.data.frame(table(tmp_pd[,cell_anno_col]))
  tmp_clust_size=tmp_clust_size[scale(tmp_clust_size[,2])<2,]
  sl_ind=tmp_pd[,cell_anno_col] %in% as.character(tmp_clust_size[,1])
  tmp_pd=tmp_pd[sl_ind,]
  
  
  ari = mclust::adjustedRandIndex(tmp_pd$anno_cluster_res, tmp_pd[,cell_anno_col])
  nmi=aricode::NMI(tmp_pd$anno_cluster_res, tmp_pd[,cell_anno_col],variant="sum")
  ami=aricode::AMI(tmp_pd$anno_cluster_res, tmp_pd[,cell_anno_col])
  ariScores_balanced=data.frame(cluster_count=n_clusters,ari=ari,nmi=nmi,ami=ami,stringsAsFactors = F)
  
  return(list(all_cellTypes=ariScores,rm_LargeCellTypes=ariScores_balanced))
}

.sconline.cluster_validation.purity=function(n_clusters,argList,cluster_obj,cell_anno_col="anno_orig_cellState", inputPd=NULL,tol_level=0.9) {
  
  if(is.null(n_clusters)){
    n_clusters=cluster_obj$n_clusters
  }
  if(is.null(n_clusters)){
    stop("Number of desired clusters should be provided!")
  }
  
  
  diff_clust=cluster_obj$cluster_object
  
  d_conMat=cutree(diff_clust,k=n_clusters)
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(is.null(inputPd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPd
  }
  
  pd=pd[!is.na(pd[,cell_anno_col]),]
  prop_mat=prop_mat[,colnames(prop_mat) %in% row.names(pd)]
  pd=pd[match(row.names(pd),colnames(prop_mat)),]
  
  
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
  prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  #prop anno
  
  colMax_vals_m=qlcMatrix::colMax(prop_merged)
  colMax_vals_m=.extra_matrix_colNorm(prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
  prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
  #prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
  prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
  prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
  pd$anno_cluster_res=paste0("C",as.character(prop_m_hardCluster$i))
  
  
  #inputPhenoData = pd;cluster_col = "anno_cluster_res";cellType_col = cell_anno_col
  purity_res=.sconline.recoveredClusters_resiprocal(inputPhenoData = pd,cluster_col = "anno_cluster_res",cellType_col = cell_anno_col)
  
  return(purity_res)
}

.sconline.seuratClusteringBN=function(argList,ncores=NULL){
  load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  pd=pd[match(row.names(pd),row.names(harmony_embeddings)),]
  library(future)
  options(future.globals.maxSize = 1000 * 1024^4)
  
  if(is.null(ncores)){
    plan("multicore", workers = argList$ncores)
  } else {
    plan("multicore", workers = ncores)
  }
  
  
  seurat_knn=Seurat:::FindNeighbors(object = harmony_embeddings[,1:argList$nPCs,drop=F])
  seurat_clusters=Seurat:::FindClusters.default(object = seurat_knn[["snn"]],resolution = c(seq(0.01,0.1,0.01),seq(0.15,4,0.05)))
  
  clust_count=apply(seurat_clusters,2,function(x) length(unique(x)))
  seurat_clusters=seurat_clusters[,!duplicated(clust_count)]
  
  ariScores=lapply(1:ncol(seurat_clusters),function(x) {
    ari = mclust::adjustedRandIndex(seurat_clusters[,x], pd$anno_orig_cellState)
    nmi=aricode::NMI(seurat_clusters[,x], pd$anno_orig_cellState,variant="sum")
    ami=aricode::AMI(seurat_clusters[,x], pd$anno_orig_cellState)
    data.frame(resolution=colnames(seurat_clusters)[x],cluster_count=length(unique(seurat_clusters[,x])),ari=ari,nmi=nmi,ami=ami,stringsAsFactors = F)
  })
  
  ariScores=do.call("rbind",ariScores)
  
  tmp_pd=pd
  
  tmp_clust_size=as.data.frame(table(tmp_pd$anno_orig_cellState))
  tmp_clust_size=tmp_clust_size[scale(tmp_clust_size[,2])<2,]
  sl_ind=tmp_pd$anno_orig_cellState %in% as.character(tmp_clust_size[,1])
  tmp_pd=tmp_pd[sl_ind,]
  tmp_seurat=seurat_clusters[sl_ind,]
  
  ariScores_balanced=lapply(1:ncol(seurat_clusters),function(x) {
    ari = mclust::adjustedRandIndex(tmp_seurat[,x], tmp_pd$anno_orig_cellState)
    nmi=aricode::NMI(tmp_seurat[,x], tmp_pd$anno_orig_cellState,variant="sum")
    ami=aricode::AMI(tmp_seurat[,x], tmp_pd$anno_orig_cellState)
    data.frame(resolution=colnames(tmp_seurat)[x],cluster_count=length(unique(tmp_seurat[,x])),ari=ari,nmi=nmi,ami=ami,stringsAsFactors = F)
  })
  
  ariScores_balanced=do.call("rbind",ariScores_balanced)
  #.res_seurat=list(all_cellTypes=ariScores,rm_LargeCellTypes=ariScores_balanced)
  return(list(all_cellTypes=ariScores,rm_LargeCellTypes=ariScores_balanced))
  
}

.sconline.GSEA.readGMT=function (file,bkg_genes=NULL,min.gs.size=NULL,max.gs.size=NULL) {
  if (!grepl("\\.gmt$", file)[1]&F) {
    stop("Pathway information must be in a .gmt file format")
  }
  geneSetDB = readLines(file)
  geneSetDB = strsplit(geneSetDB, "\t")
  names(geneSetDB) = sapply(geneSetDB, "[", 1)
  geneSetDB = lapply(geneSetDB, "[", -1:-2)
  geneSetDB = lapply(geneSetDB, function(x) {
    x[which(x != "")]
  })
  
  if(!is.null(bkg_genes)){
    for(i in 1:length(geneSetDB)){
      tmp=geneSetDB[[i]]
      tmp=bkg_genes[which(toupper(bkg_genes) %in% toupper(tmp))]
      geneSetDB[[i]]=tmp
    }
  }
  
  if(!is.null(min.gs.size)){
    size.dist=unlist(lapply(geneSetDB,length))
    geneSetDB=geneSetDB[size.dist>=min.gs.size]
  }
  
  if(!is.null(max.gs.size)){
    size.dist=unlist(lapply(geneSetDB,length))
    geneSetDB=geneSetDB[size.dist<=max.gs.size]
  }
  
  return(geneSetDB)
}

.sconline.GSEA.GSVAfn=function(expData,gs_path,min.gs.size=10,max.gs.size=500,bkg_gene_pct=0.01,ncores=10){
  require(GSVA)
  
  row.names(expData)=toupper(row.names(expData))
  
  gs_elsivier=.sconline.GSEA.readGMT(gs_path,bkg_genes = row.names(expData),min.gs.size=min.gs.size,max.gs.size=max.gs.size)
  
  gs.scores=GSVA::gsva(expData,gset.idx.list=gs_elsivier,min.sz=min.gs.size,max.sz=max.gs.size,parallel.sz=ncores,method="ssgsea")
  
  
  return(gs.scores)
}

#expData=qread("~/myBucket/integrative_Tushar/SC_data/meta_analysis/refined_anno/DEanalysis/pseudo_fixedsize/OPC_human_iNPH_size40.qs")
#=gs_path="~/myBucket/integrative_Tushar/SC_data/meta_analysis/refined_anno/GSEA/genesets/Elsevier_Pathway_Collection.gmt"
#min.gs.size=10;aucell_maxrank_fraction=0.1;max.gs.size=500;bkg_gene_pct=0.01;ncores=10

.sconline.GSEA.AUCellfn=function(expData,gs_path,min.gs.size=10,aucell_maxrank_fraction=0.1,max.gs.size=500,ncores=10){
  require(AUCell)
  
  if(is.null(expData)){
    stop("Expression data is missing!")
  } else if(tolower(class(expData))=="seurat"){
    expData=.sconline.convertSeuratToExpressionSet(object=expData)
  } else if(class(expData)!="SingleCellExperiment") {
    stop("Unrecognized inputExpression data!")
  }
  
  tmp_counts=counts(expData)
  tmp_counts=apply(tmp_counts,1,function(x) sum(x>0))
  row.names(expData)=toupper(row.names(expData))
  
  res=NULL
  gs_size=NULL
  
  gs_elsivier=.sconline.GSEA.readGMT(gs_path,bkg_genes = row.names(expData),min.gs.size=min.gs.size,max.gs.size=max.gs.size)
  if(length(gs_elsivier)>0){
    gs_size=data.frame(gs=names(gs_elsivier),size=unlist(lapply(gs_elsivier,length)))
    
    cells_rankings <- AUCell_buildRankings(as.matrix(counts(expData)))
    cells_rankings@colData=colData(expData)
    
    res= AUCell_calcAUC(gs_elsivier, rankings = cells_rankings, aucMaxRank=nrow(expData)*aucell_maxrank_fraction,nCores = 1)
    res=res@assays@data[[1]]
  }
  
  return(list(res=res,gs_size=gs_size,gs_list=gs_elsivier))
}

#expData=tmp_sl_data[row.names(tmp_sl_data) %in% bkg_genes,];gs_path=file.path(gs_dir,sl_gs);min.gs.size=10;max.gs.size=250;ncores=10
.sconline.GSEA.wilcox=function(expData,gs_path,min.gs.size=10,max.gs.size=250,ncores=10){
  
  require(coin)
  
  expData_counts=expData#@assays$RNA@counts
  
  res=NULL
  gs_size=NULL
  
  gs_elsivier=.sconline.GSEA.readGMT(gs_path,bkg_genes = row.names(expData),min.gs.size=min.gs.size,max.gs.size=max.gs.size)
  if(length(gs_elsivier)>0){
  gs_size=data.frame(gs=names(gs_elsivier),size=unlist(lapply(gs_elsivier,length)))
  
  .mywilcoxFn=function(irow,expData_counts,gs_elsivier){
    is_sl=rep("no",nrow(expData_counts))
    is_sl[tolower(row.names(expData_counts)) %in% tolower(gs_elsivier[[irow]])]="yes"
    tmp_scores=NULL
    if(sum(is_sl=="yes")>=min.gs.size&sum(is_sl=="yes")<=max.gs.size){
      tmp_scores=c()
      for(icol in 1:ncol(expData_counts)){
        tmp_scores=c(tmp_scores,(-1)*statistic(wilcox_test(expData_counts[,icol]~as.factor(is_sl)), "standardized")[1])
      }
      tmp_scores=as.data.frame(t(tmp_scores))
      colnames(tmp_scores)=colnames(expData_counts)
      row.names(tmp_scores)=names(gs_elsivier)[irow]
    }
    return(tmp_scores)
  }
  
  res=parallel::mclapply(1:length(gs_elsivier),.mywilcoxFn,expData_counts=expData_counts,gs_elsivier=gs_elsivier,mc.cores=ncores)
  
  res=do.call("rbind",res)
  }
  return(list(res=as.matrix(res),gs_size=gs_size,gs_list=gs_elsivier))
}

.sconline.GSEA.aREA=function(expData,gs_path,min.gs.size=10,max.gs.size=500,ncores=10,signed=F){
  
  require(viper)
  res=NULL
  gs_size=NULL
  gs_elsivier=.sconline.GSEA.readGMT(gs_path,bkg_genes = row.names(expData),min.gs.size=min.gs.size,max.gs.size=max.gs.size)
  if(length(gs_elsivier)>0){
    gs_size=data.frame(gs=names(gs_elsivier),size=unlist(lapply(gs_elsivier,length)))
    
    expData=as.matrix(expData)
    row.names(expData)=toupper(row.names(expData))
    
    if(signed){
      regulon=lapply(names(gs_elsivier), function(x, gset,expData){
        gene <- gset[[x]]
        gene=intersect(gene,row.names(expData))
        tmp_pc=prcomp(scale(t(expData[row.names(expData) %in% gene,])))$x
        tmp_cor=cor(t(expData[gene,]),tmp_pc[,1])
        tmp_cor=tmp_cor[gene,]
        if(sum(tmp_cor)<0){
          tmp_cor=tmp_cor*(-1)
        }
        list(tfmode=structure(rep(sign(tmp_cor), length(gene)), names=gene), likelihood=rep(1, length(gene)))
      }, gset=gs_elsivier,expData=expData)
      
    } else {
      regulon=lapply(names(gs_elsivier), function(x, gset,expData){
        gene <- gset[[x]]
        gene=intersect(gene,row.names(expData))
        list(tfmode=structure(rep(1, length(gene)), names=gene), likelihood=rep(1, length(gene)))
      }, gset=gs_elsivier,expData=expData)
      
    }
    
    names(regulon)=names(gs_elsivier)
    res=viper::aREA(expData, regulon,minsize=min.gs.size,cores=ncores)$nes
    
  }
  
  return(list(res=res,gs_size=gs_size,gs_list=gs_elsivier))
}


.myTFmode1=function (regulon, expset, method = "spearman") 
{
  expset <- viper:::filterCV(expset)
  regulon <- viper:::updateRegulon(regulon)
  
  regulon <- lapply(regulon, function(x, genes) {
    filtro <- names(x$tfmode) %in% genes
    x$tfmode <- x$tfmode[filtro]
    if (length(x$likelihood) == length(filtro)) 
      x$likelihood <- x$likelihood[filtro]
    return(x)
  }, genes = rownames(expset))
  
  for(ireg in 1:length(regulon)){
    tmp_colors=rep("grey",nrow(expset))
    tmp_colors[row.names(expset) %in% names(regulon[[ireg]]$tfmode)]="blue"
    tmp_pc=WGCNA::moduleEigengenes(t(expset),colors = tmp_colors,excludeGrey = T)
    
    
    cmat <- cor(t(expset[rownames(expset) %in% names(regulon[[ireg]]$tfmode), ]), tmp_pc$eigengenes[,"MEblue"], method = method)
    tgnames=names(regulon[[ireg]]$tfmode)
    regulon[[ireg]]$tfmode=cmat[match(names(regulon[[ireg]]$tfmode),row.names(cmat))]
    names(regulon[[ireg]]$tfmode)=tgnames
    regulon[[ireg]]$pc_scores=tmp_pc$eigengenes[,"MEblue",drop=F]
  }
  
  return(regulon)
}



.extra.color_branches=function (dend, k = NULL, h = NULL, col, groupLabels = NULL, 
                                clusters, warn = dendextend_options("warn"), ...) 
{
  #adapted from dendextend package
  all_unique= function (x, ...) 
  {
    anyDuplicated(x) == 0L
  }
  is.dendrogram=function (x) 
  {
    inherits(x, "dendrogram")
  }
  
  is.hclust=function (x) {
    inherits(x, "hclust")
  }
  
  get_col <- function(col, k) {
    if (is.function(col)) {
      col <- col(k)
    } else {
      if (length(col) < k) {
        warning("Length of color vector was shorter than the number of clusters - color vector was recycled")
        col <- rep(col, length.out = k)
      }
      if (length(col) > k) {
        warning("Length of color vector was longer than the number of clusters - first k elements are used")
        col <- col[seq_len(k)]
      }
    }
    return(col)
  }
  if (missing(col)) 
    col <- rainbow_fun
  if (!missing(clusters)) {
    if (!missing(k)) 
      warning("Both the parameters 'cluster' and 'k' are not missing - k is ignored.")
    if (length(clusters) != nleaves(dend)) {
      warning("'clusters' is not of the same length as the number of labels. The dend is returned as is.")
      return(dend)
    }
    u_clusters <- unique(clusters)
    k <- length(u_clusters)
    col <- get_col(col, k)
    return(branches_attr_by_clusters(dend, clusters, values = col, 
                                     attr = "col", branches_changed_have_which_labels = "all"))
  }
  old_labels <- labels(dend)
  labels_arent_unique <- !all_unique(old_labels)
  if (labels_arent_unique) {
    if (warn) 
      warning("Your dend labels are NOT unique!\n This may cause an un expected issue with the color of the branches.\n Hence, your labels were temporarily turned unique (and then fixed as they were before).")
    labels(dend) <- seq_along(old_labels)
  }
  if (is.null(k) & is.null(h)) {
    if (warn) 
      warning("k (number of clusters) is missing, using the dend size as a default")
    k <- nleaves(dend)
  }
  if (!is.dendrogram(dend) && !is.hclust(dend)) 
    stop("dend needs to be either a dendrogram or an hclust object")
  
  g=cutree(dend,k)
  #g <- dendextend::cutree(dend, k = k, h = h, order_clusters_as_data = FALSE)
  col=col[match(names(g),dend$label)]
  #groupLabels=groupLabels[match(names(g),dend$label)]
  #groupLabels=unique(groupLabels)
  col=unique(col)
  if (is.hclust(dend)) 
    dend <- as.dendrogram(dend)#, hang = 0)
  k <- max(g)
  if (k == 0L) {
    if (warn) 
      warning("dend has only one level - returning the dendrogram with no colors.")
    return(dend)
  }
  col <- get_col(col, k)
  if (!is.null(groupLabels)) {
    if (length(groupLabels) == 1) {
      if (is.function(groupLabels)) {
        groupLabels <- groupLabels(seq.int(length.out = k))
      }
      else if (is.logical(groupLabels)) {
        if (groupLabels) {
          groupLabels <- seq.int(length.out = k)
        }
        else {
          groupLabels <- NULL
        }
      }
    }
    if (!is.null(groupLabels) && length(groupLabels) != k) {
      stop("Must give same number of group labels as clusters")
    }
  }
  addcol <- function(dend_node, col) {
    if (is.null(attr(dend_node, "edgePar"))) {
      attr(dend_node, "edgePar") <- list(col = col)
    }
    else {
      attr(dend_node, "edgePar")[["col"]] <- col
    }
    unclass(dend_node)
  }
  descendTree <- function(sd) {
    groupsinsubtree <- unique(g[labels(sd)])
    if (length(groupsinsubtree) > 1) {
      for (i in seq(sd)) {
        sd[[i]] <- descendTree(sd[[i]])
      }
    }
    else {
      sd <- dendrapply(sd, addcol, col[groupsinsubtree])
      if (!is.null(groupLabels)) {
        attr(sd, "edgetext") <- groupLabels[groupsinsubtree]
        attr(sd, "edgePar") <- c(attr(sd, "edgePar"), 
                                 list(p.border = col[groupsinsubtree]))
      }
    }
    unclass(sd)
  }
  if (!is.character(labels(dend))) 
    labels(dend) <- as.character(labels(dend))
  dend <- descendTree(dend)
  class(dend) <- "dendrogram"
  if (labels_arent_unique) 
    labels(dend) <- old_labels
  
  noLabel <- function(x) {
    if (stats::is.leaf(x)) {
      attr(x, "label") <- NULL }
    return(x)
  }
  dend=stats::dendrapply(dend, noLabel)
  
  dend
}



.sconline.cluster.Dendro= function (clust_obj, clust_count_range=seq(5,50,5),sl_clust=NULL, groupLabels = NULL, rowText = NULL, 
                                         rowTextAlignment = c("left", "center", "right"), rowTextIgnore = NULL, 
                                         textPositions = NULL, setLayout = TRUE, autoColorHeight = TRUE, 
                                         colorHeight = 0.2, colorHeightBase = 0.2, colorHeightMax = 0.6, 
                                         rowWidths = NULL, dendroLabels = T, addGuide = FALSE, 
                                         guideAll = FALSE, guideCount = 50, guideHang = 0.2, addTextGuide = FALSE, 
                                         cex.colorLabels = 0.8, cex.dendroLabels = 0.9, cex.rowText = 0.8, 
                                         marAll = c(1, 5, 3, 1), saveMar = TRUE, abHeight = NULL, 
                                         abCol = "red", ...) 
{
  #modified version of plotDendroAndColors() from WGCNA
  
  
  clust_obj=clust_obj$cluster_object
  
  color_obj=NULL
  for(iclust in clust_count_range){
    tmp=cutree(clust_obj,k=iclust)
    tmp=data.frame(pseudo=names(tmp),clust=tmp,stringsAsFactors = F)
    if(iclust==sl_clust){
      tmp$group_label=tmp[,2]
    }
    tmp$clust=hues::iwanthue(length(unique(tmp[,2])))[tmp[,2]]
    
    colnames(tmp)[2]=paste0("clust_count_",iclust)
    if(is.null(color_obj)){
      color_obj=tmp
    } else {
      color_obj=merge(color_obj,tmp,by="pseudo")
    }
  }
  color_obj=color_obj[match(clust_obj$label,color_obj$pseudo),]
  
  color_obj=color_obj[,-1]
  group_labels=color_obj[,"group_label"]
  group_labels=paste0("C",unique(group_labels))
  color_obj=color_obj[,colnames(color_obj)!="group_label"]
  
  dendro=clust_obj
  colors=color_obj
  
  
  oldMar = par("mar")
  if (!is.null(dim(colors))) {
    nRows = dim(colors)[2]
  } else {nRows = 1}
  if (!is.null(rowText)){
    nRows = nRows + if (is.null(textPositions)) 
      nRows
  } else {length(textPositions)}
  if (autoColorHeight) 
    colorHeight = colorHeightBase + (colorHeightMax - colorHeightBase) * 
      (1 - exp(-(nRows - 1)/6))
  if (setLayout) 
    layout(matrix(c(1:2), 2, 1), heights = c(1 - colorHeight, 
                                             colorHeight))
  par(mar = c(0, marAll[2], marAll[3], marAll[4]))
  
  if(!is.null(sl_clust)){
    dendro2=.extra.color_branches(dendro,k=sl_clust,col=colors[,paste0("clust_count_",sl_clust)],groupLabels = T)
  } else {
    dendro2=dendro
  }
  
  plot(dendro2, labels = dendroLabels, cex = cex.dendroLabels,main = paste("Colors based on selected cluster count of",sl_clust), 
       ...)
  
  if (addGuide) 
    WGCNA::addGuideLines(dendro, count = if (guideAll) {
      length(dendro$height) + 1
    } else {guideCount}, hang = guideHang)
  if (!is.null(abHeight)) 
    abline(h = abHeight, col = abCol)
  par(mar = c(marAll[1], marAll[2], 0, marAll[4]))
  WGCNA::plotColorUnderTree(dendro, colors, groupLabels, cex.rowLabels = cex.colorLabels, 
                            rowText = rowText, rowTextAlignment = rowTextAlignment, 
                            rowTextIgnore = rowTextIgnore, textPositions = textPositions, 
                            cex.rowText = cex.rowText, rowWidths = rowWidths, addTextGuide = addTextGuide)
  
  if (saveMar) 
    par(mar = oldMar)
}



.myaracne2regulon=function (afile, eset, format = "3col", 
          verbose = TRUE,uniform_likelihood=T) 
{
  
  if (is(eset, "ExpressionSet")) 
    eset <- exprs(eset)
  
  if(!uniform_likelihood){
    stop("Yet to be implemented!")
  }
  
  if (verbose) 
    message("\nLoading the dataset...")
  if (length(eset) == 1) {
    tmp <- strsplit(readLines(eset), "\t")
    dset <- t(sapply(tmp[-1], function(x) as.numeric(x[-(1:2)])))
    colnames(dset) <- tmp[[1]][-(1:2)]
    rownames(dset) <- sapply(tmp[-1], function(x) x[1])
    annot <- t(sapply(tmp[-1], function(x) x[1:2]))
  } else {
    dset <- eset
    annot <- rownames(eset)
    names(annot) <- rownames(eset)
    rm(eset)
  }
  
  tmp <- afile
  aracne <- data.frame(tf = tmp[, 1], target = tmp[, 2], 
                       mi = as.numeric(tmp[, 3])/max(as.numeric(tmp[, 3])))
  
  if (verbose) 
    message("Generating the regulon objects...")
  tmp <- aracne[!is.na(aracne$mi), ]
  tmp <- tmp[rowSums(matrix(as.matrix(tmp[, 2]) %in% rownames(dset), 
                            nrow(tmp), 2)) == 2, ]
  aracne <- tapply(1:nrow(tmp), as.vector(tmp$tf), function(pos, 
                                                            tmp) {
    tfmode <- rep(0, length(pos))
    names(tfmode) <- tmp$target[pos]
    list(tfmode = tfmode, likelihood = tmp$mi[pos])
  }, tmp = tmp)
  names(aracne) <- levels(factor(as.vector(tmp$tf)))
  aracne <- .myTFmode1(aracne, dset)
  rm(dset)
  aracne <- aracne[names(aracne) != "NA"]
  aracne <- lapply(aracne, function(x) {
    filtro <- !(names(x$tfmode) == "NA" | is.na(x$tfmode) | 
                  is.na(x$likelihood))
    x$tfmode <- x$tfmode[filtro]
    x$likelihood <- x$likelihood[filtro]
    return(x)
  })
  aracne <- aracne[sapply(aracne, function(x) length(names(x$tfmode))) > 0]
  regul <- viper:::TFscore(aracne, verbose = verbose)
  class(regul) <- "regulon"
  return(regul)
}

.sconline.GSEA.ViperOriginal=function(expData,scExpset,gs_path,viper_balancing_col=NULL,min.gs.size=10,max.gs.size=500,bkg_gene_pct=0.01,ncores=10,uniform_likelihood=T,balancing_quantile=0.25){
  
  require(viper)
  
  gs_elsivier=.sconline.GSEA.readGMT(gs_path,bkg_genes = row.names(expData),min.gs.size=min.gs.size,max.gs.size=max.gs.size)
  
  gs_3col=NULL
  for(igs in 1:length(gs_elsivier)){
    tmp=gs_elsivier[[igs]]
    gs_3col=rbind(gs_3col,data.frame(tf = names(gs_elsivier)[igs], target = tmp,mi=1,stringsAsFactors = F))
  }
  
  expData=as.matrix(expData)
  row.names(expData)=toupper(row.names(expData))
  
  if(!is.null(viper_balancing_col)){
    tmp=as.data.frame(table(colData(scExpset)[,viper_balancing_col]))
    sl_thr=quantile(tmp[,2],balancing_quantile)
    sl_samples=c()
    for(ibatch in unique(colData(scExpset)[,viper_balancing_col])){
      tmp_data=colnames(scExpset)[which(colData(scExpset)[,viper_balancing_col]==ibatch)]
      tmp_data=sample(tmp_data,min(sl_thr,length(tmp_data)))
      sl_samples=c(sl_samples,tmp_data)
    }
    
  } else {
    sl_samples=colnames(scExpset)
  }
  
  regul=.myaracne2regulon(afile = gs_3col,eset = expData[,colnames(expData) %in% sl_samples],uniform_likelihood=uniform_likelihood)
  
  
  res=viper::aREA(expData, regul,minsize=min.gs.size,cores=ncores)$nes
  
  return(res)
}


.sconline.GSEAanalysisFn=function(inputExpData,gs.method="GSVA",exp.norm.method="none",batch_col=NULL,gs_gmt_file,min.gs.size=25,quantile.norm=F,max.gs.size=500,bkg_gene_pct=0.05,ncores=10,viper_signed=F,aucell_maxrank_fraction=0.1,viper_balancing_col=NULL){
  
  #batch_col: for within batch normalization
  
  exp.norm.method=match.arg(exp.norm.method,c("none","scale","rank","mad","within.batch.z","within.batch.mad"))
  
  
  if(gs.method=="AUCell"&exp.norm.method!="none"){
    cat("  setting the exp.norm.method to none for AUCell method!")
    exp.norm.method="none"
  }
  
  if(class(inputExpData)!="Seurat"){
    
    eset=counts(inputExpData)
    tmp_counts=eset
    tmp_counts=apply(tmp_counts,1,function(x) sum(x>0))
    sl_genes=row.names(inputExpData)[which(tmp_counts>max(3,ncol(inputExpData)*bkg_gene_pct))]
    pd=as.data.frame(colData(inputExpData))
    
  } else {
    eset=inputExpData@assays$RNA@counts
    tmp_counts=eset
    tmp_counts=apply(tmp_counts,1,function(x) sum(x>0))
    sl_genes=row.names(inputExpData)[which(tmp_counts>max(3,ncol(inputExpData)*bkg_gene_pct))]
    pd=as.data.frame(inputExpData@meta.data)
  }
  
  
  #eset=edgeR::cpm(as.matrix(eset))
  #eset=eset[sl_genes,]
  row.names(eset)=toupper(row.names(eset))
  
  #if(quantile.norm){
  #  eset=limma::normalizeQuantiles(eset)
  #}
  
  switch(exp.norm.method, scale = {
    tmp_mean=apply(eset,1,mean)
    tmp_sd=apply(eset,1,sd)+0.001
    eset <- sweep(eset,1,tmp_mean,"-")
    eset=sweep(eset,1,tmp_sd,"/")
    if(sum(eset>10)>0){
      eset[which(eset>10)]=10
    }
    if(sum(eset<(-10))>0){
      eset[which(eset<(-10))]=(-10)
    }
  }, rank = {
    eset <- t(apply(eset, 1, rank)) * punif(length(eset), -0.1, 
                                          0.1)
  }, mad = {
    eset <- t(apply(eset, 1, function(x) (x - median(x))/(mad(x)+0.001)))
  }, none = {
    eset <- eset
  }, within.batch.z ={
    res_norm=list()
    for(ibatch in unique(pd[,batch_col])){
      tmp=eset[,pd[,batch_col]==ibatch]
      tmp2=t(scale(t(tmp)))
      colnames(tmp2)=colnames(tmp)
      row.names(tmp2)=row.names(tmp)
      res_norm=c(res_norm,list(tmp2))
    }
    res_norm=do.call("cbind",res_norm)
    eset=res_norm[,colnames(eset)]
  }, within.batch.mad={
    res_norm=list()
    for(ibatch in unique(pd[,batch_col])){
      tmp=eset[,pd[,batch_col]==ibatch]
      tmp2=t(apply(tmp, 1, function(x) (x - median(x))/(mad(x)+0.001)))
      colnames(tmp2)=colnames(tmp)
      row.names(tmp2)=row.names(tmp)
      res_norm=c(res_norm,list(tmp2))
    }
    res_norm=do.call("cbind",res_norm)
    eset=res_norm[,colnames(eset)]
  } )
  
  res=NULL
  if(gs.method=="GSVA"){
    #expData=eset;gs_path=gs_gmt_file;min.gs.size=min.gs.size;max.gs.size=max.gs.size;ncores=ncores
    res=.sconline.GSEA.GSVAfn(expData=eset,gs_path=gs_gmt_file,min.gs.size=min.gs.size,max.gs.size=max.gs.size,ncores=ncores)
  } else if(gs.method=="AUCell"){
    res=.sconline.GSEA.AUCellfn(expData=inputExpData,gs_path=gs_gmt_file,min.gs.size=min.gs.size,aucell_maxrank_fraction=aucell_maxrank_fraction,max.gs.size=max.gs.size,ncores=ncores)
  } else if(gs.method=="aREA"){
    res=.sconline.GSEA.aREA(expData=eset,gs_path=gs_gmt_file,min.gs.size=min.gs.size,max.gs.size=max.gs.size,ncores=ncores,signed=viper_signed)
  } else if(gs.method=="aREAsigned"){
    res=.sconline.GSEA.aREA(expData=eset,gs_path=gs_gmt_file,min.gs.size=min.gs.size,max.gs.size=max.gs.size,ncores=ncores,signed=T)
  } else if(gs.method=="viper_original"){
    res=.sconline.GSEA.ViperOriginal(expData=eset,scExpset=inputExpData,gs_path=gs_gmt_file,min.gs.size=min.gs.size,max.gs.size=max.gs.size,ncores=ncores,viper_balancing_col=viper_balancing_col)
  } else if(gs.method=="Wilcox"){
    res=.sconline.GSEA.wilcox(expData=eset,gs_path=gs_gmt_file,min.gs.size=min.gs.size,max.gs.size=max.gs.size,ncores=ncores)
    
  }
  
  return(res)
}


#query = NULL; k.param = 20; prune.SNN = 1/15; nn.method = "annoy"; n.trees = 50; annoy.metric = "euclidean"; 
#nn.eps = 0; verbose = TRUE; force.recalc = FALSE; 
#cache.index = FALSE; index = NULL


.extra_sconline.sl_pseudocell.densityPeakFn_org=function (object,argList, query = NULL, k.param = 20, 
                                                      prune.SNN = 1/15, nn.method = "annoy", n.trees = 50, annoy.metric = "euclidean", 
                                                      nn.eps = 0, verbose = TRUE, force.recalc = FALSE, 
                                                      cache.index = FALSE, index = NULL) {
  plan("multicore", workers = argList$ncores)
  options(future.globals.maxSize = 1000 * 1024^4)
  cat("Performing densityPeak clustering\n")
  if (is.null(x = dim(x = object))) {
    warning("Object should have two dimensions, attempting to coerce to matrix", 
            call. = FALSE)
    object <- as.matrix(x = object)
  }
  if (is.null(rownames(x = object))) {
    stop("Please provide rownames (cell names) with the input object")
  }
  n.cells <- nrow(x = object)
  if (n.cells < k.param) {
    warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.", 
            call. = FALSE)
    k.param <- n.cells - 1
  }
  
  if(is.null(query)){
    query=object
  }
  
  {
    nn.ranked <- Seurat:::NNHelper(data = object, query = query, k = k.param, 
                                   method = nn.method, n.trees = n.trees, searchtype = "standard", 
                                   eps = nn.eps, metric = annoy.metric, cache.index = cache.index, 
                                   index = index)
    nn.ranked <- Indices(object = nn.ranked)
  }
  
  snn.matrix <- Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
  rownames(x = snn.matrix) <- rownames(x = object)
  colnames(x = snn.matrix) <- rownames(x = object)
  
  tmp=(snn.matrix+t(snn.matrix))
  tmp=.extra_matrix_rowNorm(tmp)#Matrix::Diagonal(x = 1 / (rowSums(tmp)+0.000000000001)) %*% tmp
  tmp=tmp %*% t(tmp)
  tmp=.extra_matrix_rowNorm(tmp,rowValues = 1 / (Matrix::diag(tmp)+0.000000000001))#Matrix::Diagonal(x = 1 / (Matrix::diag(tmp)+0.000000000001)) %*% tmp
  tmp@x[which(tmp@x>=1)]=1
  
  tmp_density=Matrix::rowSums(snn.matrix)
  tmp_density=tmp_density[order(tmp_density,decreasing = T)]
  
  x=tmp[names(tmp_density),names(tmp_density)]
  
  x_dim=subset(summary(x), j < i)
  x_dim=sparseMatrix(x_dim[,"i"],x_dim[,"j"],x=1,dims = c(nrow(x), ncol(x)))
  x_dim=x*x_dim
  x_dim2=qlcMatrix::rowMax(x_dim)
  x_dim2=data.frame(density=tmp_density,similarity=as.matrix(x_dim2),cell=names(tmp_density),stringsAsFactors = F)
  x_dim2=x_dim2[which(x_dim2$density>1),]
  x_dim2=x_dim2[order(x_dim2$similarity*(-1),x_dim2$density,decreasing = T),]
  
  sim_thr=x_dim2$similarity[argList$internal_pseudocell_count]
  tst=x_dim2[x_dim2$similarity<=sim_thr,]
  sl_pseudocells=row.names(tst)
  
  
  
  return(sl_pseudocells)
}

.extra_sconline.sl_pseudocell.densityPeakFn_new=function (object,argList, query = NULL, k.param = 20, 
                                                      prune.SNN = 1/15, nn.method = "annoy", n.trees = 50, annoy.metric = "euclidean", 
                                                      nn.eps = 0, verbose = TRUE, force.recalc = FALSE, 
                                                      cache.index = FALSE, index = NULL,priority_list=NULL) {
  require(qs)
  plan("multicore", workers = argList$ncores)
  options(future.globals.maxSize = 1000 * 1024^4)
  if(is.null(priority_list)){
    cat("Performing densityPeak clustering\n")
  }
  
  if (is.null(x = dim(x = object))) {
    warning("Object should have two dimensions, attempting to coerce to matrix", 
            call. = FALSE)
    object <- as.matrix(x = object)
  }
  if (is.null(rownames(x = object))) {
    stop("Please provide rownames (cell names) with the input object")
  }
  n.cells <- nrow(x = object)
  if (n.cells < k.param) {
    warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.", 
            call. = FALSE)
    k.param <- n.cells - 1
  }
  
  if(is.null(query)){
    query=object
  }
  
  
  if(!file.exists(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))){
    idx=Seurat:::AnnoyBuildIndex(data = object, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = query,k=k.param,include.distance = T,search.k = -1)
    qsave(nn.ranked.1,file=.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
  } else {
    nn.ranked.1=qread(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
  }
  
  if(F){
    tst=Seurat::FindNeighbors(object,prune.SNN =0)
    tst_nn=as(tst[["nn"]],"dgCMatrix")
    tst_snn=as(tst[["snn"]],"dgCMatrix")
    x=tst_nn*tst_snn
    length(x@x)
    
    x2=adj*tst_nn
    length(x2@x)
    length(adj@x)
  }
  
  
  
  
  
  affinities=.extra_matrix_rowNorm(input_mat = nn.ranked.1$nn.dists,rowValues = 1/(nn.ranked.1$nn.dists[,2]+0.000001))#Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,2]+0.000001)) %*% nn.ranked.1$nn.dists
  affinities@x=-1*affinities@x^2
  affinities@x=exp(affinities@x)
  affinities[,1]=affinities[,2]
  j <- as.numeric(t(nn.ranked.1$nn.idx))
  i <- ((1:length(j)) - 1)%/%k.param + 1
  x=as.numeric(t(affinities))
  adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(object),nrow(object)))
  rownames(adj) <- row.names(object)
  colnames(adj)=c(row.names(object))
  adj=Matrix::drop0(adj,tol=0.01)
  adj=.extra_matrix_rowNorm(adj)# Matrix::Diagonal(x=1/rowSums(adj)) %*% adj
  #adj=Matrix::drop0(adj,tol=quantile(adj@x,0.8))
  adj_t=t(adj)
  adj_t=.extra_matrix_rowNorm(adj_t)# Matrix::Diagonal(x=1/rowSums(adj_t)) %*% adj_t
  adj_t=adj %*% adj_t
  
  if(F){
    for(iitr in 1:4){
      adj=adj %*% adj_t*0.7+adj_t *0.3
      adj=Matrix::drop0(adj,tol=quantile(adj@x,0.1))
    }
  }
  
  #adj=adj %*% adj_t
  
  tst=Matrix::drop0(adj,tol=quantile(adj_t@x,0.8))
  tst@x=rep(1,length(tst@x))
  
  tmp_density=Matrix::rowSums(tst)
  if(F){
    if(!is.null(priority_list)){
      tmp_density[names(tmp_density) %in% priority_list]=max(tmp_density)+tmp_density[names(tmp_density) %in% priority_list]
    }  
  }
  
  tmp_density=tmp_density[order(tmp_density,decreasing = T)]
  if(F){
    pd=pd[match(names(tmp_density),row.names(pd)),]
    harmony_embeddings=harmony_embeddings[match(row.names(pd),row.names(harmony_embeddings)),drop=F]
    knet=RANN::nn2(data=harmony_embeddings,query = harmony_embeddings,k=100,eps=0)
    res_arranged=apply(as.data.frame((knet$nn.idx)),2,function(x) pd$anno_orig_cellState[x])
    res_arranged=res_arranged[,-1]
    res_counts=apply(res_arranged,1,function(x) {x=as.numeric(table(x)); max(x)/sum(x)})
    summary(res_counts)
    summary(res_counts>0.9)
    summary(res_counts>0.5)
    pd$purity=res_counts
    
    
    pd2=pd[match(names(tmp_density),row.names(pd)),]
    pd2$anno_orig_cellState=as.character(pd2$anno_orig_cellState)
    head(pd2$anno_orig_cellState,20)
    head(pd2$purity,20)
    pd2$purity[which(pd2$anno_orig_cellState=="Group1")[1:3]]
    pd2$purity[which(pd2$anno_orig_cellState=="Group3")[1:3]]
    which(pd2$anno_orig_cellState=="Group1")[1:3]
  }
  
  
  res_dist=1
  adj2=adj[match(names(tmp_density),row.names(adj)),match(names(tmp_density),colnames(adj))]
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  
  while(nrow(adj2)>max(5*argList$internal_pseudocell_count,11000)){
    sl_ps_list=c()
    batch_size=max(10000,5*argList$internal_pseudocell_count)
    if(nrow(adj2) %% batch_size<3){
      batch_size=batch_size+3
    }
    for(i in seq(1,nrow(adj2),batch_size)){
      adj3=adj2[i:min(nrow(adj2),i+batch_size-1),]
      prop_mat2=myL2normFn(inputMat=adj3)
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
      
      
      
      tst=summary(c_c_aff)
      tst=tst[tst[,1]>tst[,2],]
      tst=sparseMatrix(i = tst[,1], j = tst[,2], x =tst[,3], dims = c(nrow(adj3),nrow(adj3)))
      row.names(tst)=row.names(adj3)
      tst=as.numeric(qlcMatrix::rowMax(tst))
      res_dist=tst
      names(res_dist)=row.names(adj3)
      res_dist=res_dist[order(res_dist,decreasing = F)]
      sl_ps_list=c(sl_ps_list,names(res_dist)[1:min(length(res_dist),argList$internal_pseudocell_count)])
    }
    
    sl_ps_list=sl_ps_list[!is.na(sl_ps_list)]
    adj2=adj2[sl_ps_list,sl_ps_list]
  }
  
  
  prop_mat2=myL2normFn(inputMat=adj2)
  c_c_aff=t(prop_mat2)
  c_c_aff=prop_mat2 %*% c_c_aff
  
  tst=summary(c_c_aff)
  tst=tst[tst[,1]>tst[,2],]
  tst=sparseMatrix(i = tst[,1], j = tst[,2], x =tst[,3], dims = c(nrow(adj2),nrow(adj2)))
  row.names(tst)=row.names(adj2)
  tst=as.numeric(qlcMatrix::rowMax(tst))
  res_dist=tst
  names(res_dist)=row.names(adj2)
  res_dist=res_dist[order(res_dist,decreasing = F)]
  
  if(F){
    load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
    length(unique(pd$anno_orig_cellState[pd$sample %in% names(res_dist)[1:200]]))
    setdiff(pd$anno_orig_cellState,unique(pd$anno_orig_cellState[pd$sample %in% names(res_dist)[1:200]]))
    
    table(pd$anno_orig_cellState[pd$sample %in% names(res_dist)[1:400]])
    
  }
  
  sl_pseudocells=names(res_dist)[1:min(argList$internal_pseudocell_count,sum(res_dist<0.2))]
  
  return(sl_pseudocells)
}


#query = NULL; k.param = 20; prune.SNN = 1/15; nn.method = "annoy"; n.trees = 50; annoy.metric = "euclidean"; nn.eps = 0; verbose = TRUE; force.recalc = FALSE; cache.index = FALSE; index = NULL
.extra_sconline.sl_pseudocell.densityPeakFn=function (object,argList, query = NULL, k.param = 20, 
                                                      prune.SNN = 1/15, nn.method = "annoy", n.trees = 50, annoy.metric = "euclidean", 
                                                      nn.eps = 0, verbose = TRUE, force.recalc = FALSE, 
                                                      cache.index = FALSE, index = NULL,priority_list=NULL,saveFiles=T) {
  require(qs)
  plan("multicore", workers = argList$ncores)
  options(future.globals.maxSize = 1000 * 1024^4)
  if(is.null(priority_list)){
    cat("Performing densityPeak clustering\n")
  }
  
  if (is.null(x = dim(x = object))) {
    warning("Object should have two dimensions, attempting to coerce to matrix", 
            call. = FALSE)
    object <- as.matrix(x = object)
  }
  if (is.null(rownames(x = object))) {
    stop("Please provide rownames (cell names) with the input object")
  }
  n.cells <- nrow(x = object)
  if (n.cells < k.param) {
    warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.", 
            call. = FALSE)
    k.param <- n.cells - 1
  }
  
  if(is.null(query)){
    query=object
  }
  
  
  if(!file.exists(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))|(!saveFiles)){
    idx=Seurat:::AnnoyBuildIndex(data = object, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = query,k=k.param,include.distance = T,search.k = -1)
    if(saveFiles){
      qsave(nn.ranked.1,file=.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
    }
    
  } else {
    nn.ranked.1=qread(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
  }
  
  affinities=.extra_matrix_rowNorm(input_mat = nn.ranked.1$nn.dists,rowValues = 1/(nn.ranked.1$nn.dists[,2]+0.000001))#Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,2]+0.000001)) %*% nn.ranked.1$nn.dists
  affinities@x=-1*affinities@x^2
  affinities@x=exp(affinities@x)
  affinities[,1]=affinities[,2]
  j <- as.numeric(t(nn.ranked.1$nn.idx))
  i <- ((1:length(j)) - 1)%/%k.param + 1
  x=as.numeric(t(affinities))
  adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(object),nrow(object)))
  adj= .extra_matrix_rowNorm(adj)#Matrix::Diagonal(x=1/rowSums(adj)) %*% adj
  rownames(adj) <- row.names(object)
  colnames(adj)=c(row.names(object))
  
  
  adj_t=t(adj)
  adj_t= .extra_matrix_rowNorm(adj_t)#Matrix::Diagonal(x=1/rowSums(adj_t)) %*% adj_t
  adj=adj %*% adj_t
  
  tmp_density=Matrix::rowSums(adj)
  if(!is.null(priority_list)){
    tmp_density[names(tmp_density) %in% priority_list]=max(tmp_density)+tmp_density[names(tmp_density) %in% priority_list]
  }
  tmp_density=tmp_density[order(tmp_density,decreasing = T)]
  
  res_dist=1
  adj2=adj[match(names(tmp_density),row.names(adj)),match(names(tmp_density),colnames(adj))]
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  
  while(nrow(adj2)>max(5*argList$internal_pseudocell_count,11000)){
    sl_ps_list=c()
    batch_size=max(10000,5*argList$internal_pseudocell_count)
    if(nrow(adj2) %% batch_size<3){
      batch_size=batch_size+3
    }
    for(i in seq(1,nrow(adj2),batch_size)){
      adj3=adj2[i:min(nrow(adj2),i+batch_size-1),]
      prop_mat2=myL2normFn(inputMat=adj3)
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
      
      
      
      tst=summary(c_c_aff)
      tst=tst[tst[,1]>tst[,2],]
      tst=sparseMatrix(i = tst[,1], j = tst[,2], x =tst[,3], dims = c(nrow(adj3),nrow(adj3)))
      row.names(tst)=row.names(adj3)
      tst=as.numeric(qlcMatrix::rowMax(tst))
      res_dist=tst
      names(res_dist)=row.names(adj3)
      res_dist=res_dist[order(res_dist,decreasing = F)]
      sl_ps_list=c(sl_ps_list,names(res_dist)[1:min(length(res_dist),argList$internal_pseudocell_count)])
    }
    
    sl_ps_list=sl_ps_list[!is.na(sl_ps_list)]
    adj2=adj2[sl_ps_list,sl_ps_list]
  }
  
  
  prop_mat2=myL2normFn(inputMat=adj2)
  c_c_aff=t(prop_mat2)
  c_c_aff=prop_mat2 %*% c_c_aff
  
  tst=summary(c_c_aff)
  tst=tst[tst[,1]>tst[,2],]
  tst=sparseMatrix(i = tst[,1], j = tst[,2], x =tst[,3], dims = c(nrow(adj2),nrow(adj2)))
  row.names(tst)=row.names(adj2)
  tst=as.numeric(qlcMatrix::rowMax(tst))
  res_dist=tst
  names(res_dist)=row.names(adj2)
  res_dist=res_dist[order(res_dist,decreasing = F)]
  
  if(F){
    load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
    length(unique(pd$anno_orig_cellState[pd$sample %in% names(res_dist)[1:200]]))
    setdiff(pd$anno_orig_cellState,unique(pd$anno_orig_cellState[pd$sample %in% names(res_dist)[1:200]]))
    
    table(pd$anno_orig_cellState[pd$sample %in% names(res_dist)[1:400]])
    
  }
  
  sl_pseudocells=names(res_dist)[1:min(argList$internal_pseudocell_count,sum(res_dist<0.2))]
  
  return(sl_pseudocells)
}




.extra_sconline.sl_pseudocell.kmeansFn=function(harmony_embeddings,argList,saveFiles=T,seed=1){
  
  myPseudoAffinityMakerFn=function(harmony_embeddings,k.param=20,prune.SNN=1/15,n.trees = 50,saveFiles){
    #nn.ranked.1 <- RANN::nn2(harmony_embeddings, k = 10, eps = 0)
    
    if(!file.exists(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))|(!saveFiles)){
      idx=Seurat:::AnnoyBuildIndex(data = harmony_embeddings, metric = "euclidean", 
                                   n.trees = n.trees)
      nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = harmony_embeddings,k=k.param,include.distance = T,search.k = -1)
      if(saveFiles){
        qsave(nn.ranked.1,file=.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
      }
      
    } else {
      nn.ranked.1=qread(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
    }
    
    nn.ranked=nn.ranked.1$nn.idx
    graph= Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(graph) <- rownames(x = harmony_embeddings)
    colnames(graph) <- rownames(x = harmony_embeddings)
    
    
    #graph=Seurat::FindNeighbors(harmony_embeddings,compute.SNN=T)
    #graph=as(graph[["snn"]], "dgCMatrix")
    
    #j <- as.numeric(t(nn.ranked.1$nn.idx))
    #i <- ((1:length(j)) - 1)%/%ncol(nn.ranked.1$nn.idx) + 1
    #k=1#as.numeric(t(nn.ranked.1$nn.dists))
    #graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(harmony_embeddings), nrow(harmony_embeddings)))
    
    #graph=graph+t(graph)
    #graph@x=rep(1,length(graph@x))
    
    if(F){
      graph=graph %*% t(graph)
      diag(graph)=0
      
      
      listCols_sparse<-function(X){
        #converts a sparse Matrix into a list of its columns
        #each list item contains only the nonzero elements of the column
        X<-as(X,"CsparseMatrix")
        res<-split(X@x, findInterval(seq_len(nnzero(X)), X@p, left.open=TRUE))
        names(res)<-colnames(X)
        res
      }
      
      colapply_sparse_nonzero<-function(X,FUN,...,mc.cores=1){
        #apply a function FUN to NONZERO elements of each column of sparse Matrix X
        #for an alternative that operates on all values, see colapply_full
        #mc: should parallel processing be used? Only recommended if FUN is slow
        #... additional args passed to mclapply or to FUN
        #this function always returns a list of length ncol(X)
        if(mc.cores>1){
          res=mclapply(listCols_sparse(X),FUN,...,mc.cores=mc.cores)
        } else {
          res=lapply(listCols_sparse(X),FUN,...)
        }
        res=unlist(res)
        X@x=res
        return(X)
      }
      
      affinities=colapply_sparse_nonzero(X=t(graph),FUN=function(x) exp((-3)*((max(x)-x)/(max(x)+1))^2),mc.cores=argList$ncores)
      affinities=sqrt(t(affinities)*affinities)
      diag(affinities)=0
    } else {
      affinities=graph
    }
    
    return(affinities)
  }
  
  set.seed(seed)
  doClustering=T
  itrClustering=0
  while(doClustering&itrClustering<10){
    itrClustering=itrClustering+1
    res_clust=kmeans(harmony_embeddings[,1:argList$nPCs,drop=F],argList$internal_pseudocell_count,iter.max = 1000) #,algorithm = "Lloyd")
    
    if(sum(is.na(res_clust$centers))==0){
      doClustering=F
    }
  }
  
  if(sum(is.na(res_clust$centers))>0){
    stop("Error in identification of the pseudocells")
  }
  
  res_clusters=data.frame(sample=names(res_clust$cluster),cluster_id=res_clust$cluster,stringsAsFactors = F)
  if(saveFiles){
    save(res_clusters,file=.myFilePathMakerFn("kmeans_res_clusters",argList=argList))
  }
  
  
  pseudo_affinity=myPseudoAffinityMakerFn(harmony_embeddings = harmony_embeddings,saveFiles=saveFiles)
  row.names(pseudo_affinity)=colnames(pseudo_affinity)=row.names(harmony_embeddings)
  
  sl_pseudo=NULL
  for(x in unique(res_clusters$cluster_id)){
    if(sum(res_clusters$cluster_id==x)>5){
      tmp=res_clusters$sample[res_clusters$cluster_id==x]
      tmp=rowSums(as.matrix(pseudo_affinity[tmp,tmp]))
      tmp=tmp[order(tmp,decreasing = T)]
      tmp=names(tmp)[1]
      x=data.frame(cluster=x,pseudocell=tmp,stringsAsFactors = F)
      sl_pseudo=rbind(sl_pseudo,x)
    }
  }
  
  #pca_centroid=res_clust$centers
  pca_centroid=harmony_embeddings[sl_pseudo$pseudocell,,drop=F]
  row.names(pca_centroid)=sl_pseudo$cluster
  pca_centroid=pca_centroid[row.names(res_clust$centers)[row.names(res_clust$centers) %in% row.names(pca_centroid)],]
  
  output=sl_pseudo$pseudocell
  return(output)
}

#prepare the embedding space and defines the pseudocells in it using k-means
#embedding can be provided as input or it can be calculated by the method using pca and optionally harmony
#inputEmbeddings: a matrix of cell x embedding
#run_harmony: logical; if true, approach performs harmony to do batch correction
.sconline.embeddingFn=function(argList,inputEmbeddings,run_harmony,inputBatchCol='anno_batch',pd=NULL,saveFiles=T){
  #argList=.ArgList;inputEmbeddings=NULL;run_harmony=F;inputBatchCol='anno_batch';pd=NULL;saveFiles=F
  require(Matrix)
  require(qlcMatrix)
  reRunCheck=tryCatch({load(.myFilePathMakerFn("pca_centroids",argList=argList));F}, error=function(e) {return(T)})
  
  if(saveFiles){
    if(!(reRunCheck|argList$newRun)){
      return("Done")
    }
  }
  
  
  res="Done"
  if(is.null(inputEmbeddings)){
    cat('Creating Embeddings\n')
    
    library(future)
    plan("multicore", workers = min(parallel::detectCores(), 8, na.rm=T))
    # plan()
    options(future.globals.maxSize = 1000 * 1024^4)
    
    if(is.null(pd)){
      load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    }
    
    if(run_harmony){
      
      
      
      # CHANGED iter.max per https://github.com/immunogenomics/harmony/issues/25 and per line below
      # also changed lloyd
      #devtools::load_all("~/anaconda3/envs/psuedocells/harmony/")
      
      harmony_embeddings <- harmony::HarmonyMatrix(pca_res[,1:argList$nPCs], pd, 'anno_batch', do_pca = FALSE, verbose=FALSE)
      
    } else {
      harmony_embeddings=pca_res[,1:argList$nPCs]
    }
    
    if(saveFiles){
      .mySaveFn(harmony_embeddings,file=.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    } else {
      res=harmony_embeddings
    }
    
    
    
  } else {
    inputEmbeddings=as.matrix(inputEmbeddings)
    
    if(sum(is.na(as.matrix(inputEmbeddings)))>0){
      stop("Error in the inputEmbeddings")
    }
    
    if(saveFiles){
      if(!is.null(pd)&!is.null(inputEmbeddings)){
        pca_res=inputEmbeddings
        save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
      } else {
        stop("Run .sconline.pca function first or provide the inputEmbeddings")
      }
    }
    
    
    if(run_harmony){
      
      # CHANGED iter.max per https://github.com/immunogenomics/harmony/issues/25 and per line below
      # also changed lloyd
      #devtools::load_all("~/anaconda3/envs/psuedocells/harmony/")
      if(is.null(pd)){
        stop("phenoData should be provided")
      }
      inputEmbeddings <- harmony::HarmonyMatrix(inputEmbeddings[,1:argList$nPCs], pd, inputBatchCol, do_pca = FALSE, verbose=FALSE)
      
    }
    
    harmony_embeddings <- inputEmbeddings[,1:argList$nPCs]
    if(saveFiles){
      .mySaveFn(harmony_embeddings,file=.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    } else {
      res=harmony_embeddings
    }
    
    
  }
  
  argList$internal_pseudocell_count=round(nrow(harmony_embeddings)/argList$pseudocell_size)
  
  
  if(argList$pseudocell_selection_method=="kmeans"){
    pseudocell_names=.extra_sconline.sl_pseudocell.kmeansFn(argList=argList,harmony_embeddings=harmony_embeddings,saveFiles = saveFiles)#,kmeansMethod=kmeans_method)
    #.pseudocell_names=pseudocell_names
    if(F){
      knet=RANN::nn2(data=harmony_embeddings,query = harmony_embeddings[pseudocell_names,,drop=F],k=100,eps=0)
      res_arranged=apply(as.data.frame((knet$nn.idx)),2,function(x) pd$anno_orig_cellState[x])
      res_arranged=res_arranged[,-1]
      res_counts=apply(res_arranged,1,function(x) {x=as.numeric(table(x)); max(x)/sum(x)})
      summary(res_counts)
      summary(res_counts>0.9)
      summary(res_counts>0.5)
      pd_sl=pd[pseudocell_names,]
      table(pd_sl$anno_orig_cellState[res_counts<.8])
    }
    
    #object=harmony_embeddings;argList= argList;priority_list=pseudocell_names
    pseudocell_names=.extra_sconline.sl_pseudocell.densityPeakFn(object=harmony_embeddings,argList= argList,priority_list=pseudocell_names,saveFiles=saveFiles)
  } else {
    pseudocell_names=.extra_sconline.sl_pseudocell.densityPeakFn(object=harmony_embeddings,argList= argList,saveFiles=saveFiles)
  }
  
  
  #pca_centroid=res_clust$centers
  
  pca_centroid=harmony_embeddings[pseudocell_names,,drop=F]
  row.names(pca_centroid)=paste0("ps_",1:nrow(pca_centroid))
  sl_pseudo=data.frame(cluster=paste0("ps_",1:nrow(pca_centroid)),pseudocell=pseudocell_names,stringsAsFactors = F)
  
  if(saveFiles){
    save(pca_centroid,file=.myFilePathMakerFn("pca_centroids",argList=argList))
    save(sl_pseudo,file=.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
  } else {
    res=list(harmony_embeddings=res,pca_centroid=pca_centroid,sl_pseudo=sl_pseudo)
  }
  
  
  
  
  
  return(res)
}

#Constructs the argument list specifying the algorithm parameters
#Path to which save the results or load the results from
#n.prop: number of propagations to do
#newRun: to redo all analyses
#min_ds_size: datasets smaller than this size (based on the defined batch structure) are excluded from the analysis
.sconline.arg_creator=function(prefix,do.split.prop=F,min_cluster_size=20,saveDir,n.prop=4,nPCs,HVG_count=3,HVG_list=NULL,indScaling=T,ncores=min(parallel::detectCores(), 8, na.rm=T),pseudocell_selection_method="kmeans",pseudocell_size=200,input_highly_var_genes=NULL,newRun=F,min_ds_size=50,Leng200=F,include.singletons=T,colNormalize=T,singleton.method="fast"){
  
  dfSettings=NULL
  #HVG_count: for a gene to be selected as highly variable, how many datasets should support its variability. can't be larger than the number of datasets 
  #indScaling: performe independent scaling before pca analysis
  #nn.method: "fast", "snn"
  dfSettings=rbind(dfSettings,data.frame(commonExpressed=T,
                                         nPCs=nPCs,
                                         HVG_count=HVG_count,covariates="",indScaling=indScaling,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  inputExpData=data
  runSetting=dfSettings
  runIndx=1
  DE_supportingFractionThr=0.1
  DE_n.adaptiveKernel=20;pval_thr=0.001
  
  
  if(sum(duplicated(row.names(inputExpData)))>0){
    print(paste(sum(duplicated(row.names(inputExpData))),"Duplicate gene ids were found! duplicates were randomly removed from the data"))
    inputExpData=inputExpData[!duplicated(row.names(inputExpData)),]
  }
  
  argList=.extra_sconline.argFn(runIndx=runIndx,min_cluster_size=min_cluster_size,saveDir=saveDir,do.split.prop=do.split.prop,exNonMicCells=F,ncores=ncores,sensitiveSearch=1,include.singletons=include.singletons,colNormalize=colNormalize,includeHuman=F,includeMouse=F,FinnishSbjBased=F,uniformZscore=F,Leng200=Leng200,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=n.prop,pseudocell_selection_method=pseudocell_selection_method,prop.n.neighbors=n.prop,dist_zscore_gamma=F,dist_zscore_norm=T,dist_zscore_nbinom=F,regularize=T,geoMean=F,prefix=prefix,newRun = newRun,inputDf=runSetting,internal_pseudocell_count=NULL,pseudocell_size=pseudocell_size,external_DE_path=NULL,external_DE_name=NULL,singleton.method=singleton.method)
  argList$exclude_non_freq_pseudocells=F
  argList$input_highly_var_genes=input_highly_var_genes
  argList$min_ds_size=min_ds_size
  argList$HVG_list=argList$HVG_count
  
  # JONAH ADDED 
  # .ArgList$nHighlyVar = 3500
  
  #saveDir is the address to the location that main results are saved
  argList$saveDir
  #saveDirGlobal includes address to the expression dataset
  argList$saveDirGlobal
  argList$newRun=newRun
  return(argList)
}

#Performs HVG (Highly variable gene selection) and PCA analsyis
#if HVG_list is provided, it performs pca on different sets of HVGs
#batch_variable: the column name in the meta data that specifies the batch structure of the data

.sconline.pca=function(inputExpData,argList,batch_variable="anno_batch",organism="unknown",addAnno=F,hierarchical_mode=F,set_HVG_count_thr=2000){
  #library(rliger)
  library(Seurat)
  library(SingleCellExperiment)
  library(scran)
  
  data=.extra_sconline.exp_creatorFn(argList=argList,inputExpData,batch_variable=batch_variable,organism = organism,addAnno=addAnno)
  rm(inputExpData)
  
  if(is.null(argList$HVG_list)){
    argList$HVG_list=argList$HVG_count
  }
  
  argList$HVG_list=argList$HVG_list[argList$HVG_list<=length(data$data)]
  if(length(argList$HVG_list)==0){
    if(is.null(set_HVG_count_thr)){
      stop("Error in the input HVG_list")
    }
  }
  
  cat('Loading and preprocessing the datasets\n')
  run_check=tryCatch({!file.exists(.myFilePathMakerFn("pca_anno",argList = argList,pseudoImportant = F))|argList$newRun|hierarchical_mode}, error=function(e) {return(T)})
  if(length(run_check)==0){
    run_check=T
  }
  if(run_check){
    
    cat('Selecting the highly variable genes\n')
    data=.myHighVarGeneSlFn(data,dataorganism=organism,argList = argList,batch_variable=batch_variable)
    
    
    if(length(argList$HVG_list)==0&(!is.null(set_HVG_count_thr))){
      tmp=data$varFeatures
      tmp=tmp[order(tmp[,2],decreasing = T),]
      argList$HVG_list=argList$HVG_count=tmp[set_HVG_count_thr,2]+1
    }
    
    tmp=data[c("varFeatures","allGenes" )]
    if(!hierarchical_mode){
      if(!file.exists(.myFilePathMakerFn("varGenes",argList=argList,varGenes = T))|argList$newRun){
        save(tmp,file=.myFilePathMakerFn("varGenes",argList=argList,varGenes =T))
      }
    }
    
    if(!hierarchical_mode){
      tmp=data$data_m
      if(!file.exists(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))|argList$newRun){
        if(!is.null(data$data_m)){
          qs::qsave(tmp,file=.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
        }
      }
    }
    
    
    
    argList$HVG_list=argList$HVG_list[argList$HVG_list<=max(data$varFeatures[,2])]
    if(length(data$varFeatures[data$varFeatures[,2]>=argList$HVG_list,1])<500){
      if(is.null(set_HVG_count_thr)){
        stop("Error in the input HVG_list")
      } else {
        warning("incompatible input HVG_list; Using the top 2000 HVGs")
        tmp=data$varFeatures
        tmp=tmp[order(tmp[,2],decreasing = T),]
        argList$HVG_count=argList$HVG_list=tmp[set_HVG_count_thr,2]+0.999
        
      }
    }
    cat('Scaling the data and PCA\n')
    res=.myPCAfn(data,argList = argList,saveFiles = (!hierarchical_mode))
    if(!hierarchical_mode){
      res="Done!"
    }
  } else {
    if(!argList$newRun){
      stop("pca file already exists. if need reRun, set the reRun variable in the argList to TRUE")
    }
  }
  
  return(res)
}

#running the sconline propagation method
#inputEmbeddings: optional; a matrix of [cells] x [PC] dimension
#inputExpData: optional; gene expression data; if specified, batch_variable should be defined too (specifies the column in the meta data that indicate the batch structure of the dataset)
#inputPhenoData should be provided if we want to do harmony on the inputEmbeddings
#batchVariable: the column name of the batchVariable in the phenoData
.sconline.runPropagation=function(argList,inputEmbeddings=NULL,inputPhenoData=NULL,inputExpData=NULL,input_UMAP_embedding=NULL,organism,batch_variable="anno_batch",run_harmony=F,addAnno=F,addClusteringModule=F,umap.method='umap-learn',extendedMode=F,L2Norm=T,mergePseudocells=T,generateUMAP=T,merging_strength=0.3,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2){
  
  require(Matrix)
  require(qlcMatrix)
  #argList=.ArgList;inputEmbeddings=NULL;inputPhenoData=NULL;inputExpData=NULL;organism="Mouse";batch_variable="anno_batch";run_harmony=F;addAnno=F;addClusteringModule=F;L2Norm=T;mergePseudocells=T;generateUMAP=F;extendedMode=F;merging_strength=0.3;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;input_UMAP_embedding=NULL
  #run_harmony=T
  #umap.method="uwot";generateUMAP=T
  #pd=inputPhenoData;inputBatchCol=batch_variable
  
  if(!is.null(input_UMAP_embedding)){
    if(sum(colnames(input_UMAP_embedding) %in% c("UMAP_1","UMAP_2"))!=2){
      stop("input_UMAP_embedding is missing the UMAP_1 and UMAP_2 columns!")
    }
  }
  
  res_embeddings=.sconline.embeddingFn(argList,inputEmbeddings=inputEmbeddings,run_harmony=run_harmony,pd=inputPhenoData,inputBatchCol=batch_variable)
  
  res_umap=.sconline.umapFn(argList,umap.method=umap.method,generateUMAP = generateUMAP,input_UMAP_embedding=input_UMAP_embedding)
  #pd=.sconline.fetch_data("annotation",argList)
  #p=.myDimPlotFn(object=pd, dimCols = c("UMAP_1", "UMAP_2"), attCol = "anno_orig_cellState")
  #ggsave(plot=p,file="~/myBucket/torm.pdf")
  
  if(!is.null(inputExpData)){
    data=.extra_sconline.exp_creatorFn(argList=argList,inputExpData = inputExpData,batch_variable = batch_variable,organism = organism,addAnno=addAnno)
    data=data$data
    #expData = data
    res_prop=.myConcensusDEFn_step2(argList,expData = data,addClusteringModule=addClusteringModule,extendedMode = extendedMode,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr)
  } else {
    res_prop=.myConcensusDEFn_step2(argList,extendedMode = extendedMode,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr)
  }
  
  return(res_prop)
}

.sconline.runPropagation_fast=function(argList,inputEmbeddings=NULL,inputPhenoData=NULL,inputExpData=NULL,input_UMAP_embedding=NULL,organism,batch_variable="anno_batch",run_harmony=F,addAnno=F,addClusteringModule=F,umap.method='uwot',extendedMode=F,L2Norm=T,mergePseudocells=T,generateUMAP=T,merging_strength=0.3,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2){
  
  require(Matrix)
  require(qlcMatrix)
  #argList=.ArgList;inputEmbeddings=NULL;inputPhenoData=NULL;inputExpData=NULL;organism="Mouse";batch_variable="anno_batch";
  #run_harmony=F;addAnno=F;addClusteringModule=F;L2Norm=T;mergePseudocells=T;generateUMAP=T;extendedMode=F;merging_strength=0.3;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;input_UMAP_embedding=NULL
  #run_harmony=T
  #umap.method="uwot";generateUMAP=T
  #pd=inputPhenoData;inputBatchCol=batch_variable
  
  if(!is.null(input_UMAP_embedding)){
    if(sum(colnames(input_UMAP_embedding) %in% c("UMAP_1","UMAP_2"))!=2){
      stop("input_UMAP_embedding is missing the UMAP_1 and UMAP_2 columns!")
    }
  }
  
  res_embeddings=.sconline.embeddingFn(argList,inputEmbeddings=inputEmbeddings,run_harmony=run_harmony,pd=inputPhenoData,inputBatchCol=batch_variable)
  
  res_umap=.sconline.umapFn(argList,umap.method=umap.method,generateUMAP = generateUMAP,input_UMAP_embedding=input_UMAP_embedding)
  #pd=.sconline.fetch_data("annotation",argList)
  #p=.myDimPlotFn(object=pd, dimCols = c("UMAP_1", "UMAP_2"), attCol = "anno_orig_cellState")
  #ggsave(plot=p,file="~/myBucket/torm.pdf")
  
  if(!is.null(inputExpData)){
    inputExpData$oneDS="oneDS"
    data=.extra_sconline.exp_creatorFn(argList=argList,inputExpData = inputExpData,batch_variable = "oneDS",organism = organism,addAnno=addAnno)
    data=data$data
    #expData = data
    res_prop=.myConcensusDEFn_step2_fast(argList,expData = data,addClusteringModule=addClusteringModule,extendedMode = extendedMode,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr)
  } else {
    res_prop=.myConcensusDEFn_step2(argList,extendedMode = extendedMode,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr)
  }
  
  return(res_prop)
}

#Creates a Seurat object from the cells
.sconline.create_seurat=function(argList,inputExpData=NULL,n_clusters=NULL,inputEmbeddings=NULL){
  if(is.null(inputExpData)){
    inputExpData=qs::qread(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
  } else if(class(inputExpData)=="SingleCellExperiment"){
    inputExpData=.extraExport2SeuratFn(inputExpData)
  } else if(class(inputExpData)=="liger"){
    inputExpData=.extra_sconline_LigerToExpSet(inputExpData)
    inputExpData=.extraExport2SeuratFn(inputExpData)
  } else {
    stop("unrecognized inputExpData format")
  }
  
  if(is.null(inputEmbeddings)){
    inputEmbeddings=.sconline.fetch_data(dataType = "embeddings",argList = argList)
  }
  
  pd=.sconline.fetch_data("annotation",argList = argList)
  inputExpData=inputExpData[,colnames(inputExpData) %in% row.names(pd)]
  if(length(setdiff(colnames(inputExpData),row.names(pd)))>0){
    stop("Error! phenotype info were not found for some of the cells")
  } else {
    pd=pd[match(colnames(inputExpData),row.names(pd)),]
  }
  
  sl_meta_data=setdiff(colnames(pd),colnames(inputExpData@meta.data))
  if(length(sl_meta_data)>0){
    if(length(sl_meta_data)==1){
      tmp_pd=data.frame(new=pd[,sl_meta_data])
      colnames(tmp_pd)=sl_meta_data
      AddMetaData(inputExpData,tmp_pd) 
    } else {
      inputExpData=AddMetaData(inputExpData,pd[,sl_meta_data]) 
    }
  }
  
  umap_data=.sconline.fetch_data(dataType = "umap",argList = argList)
  umap_data=umap_data[match(colnames(inputExpData),row.names(umap_data)),]
  umap.reduction <- Seurat::CreateDimReducObject(embeddings = umap_data, 
                                         key = "UMAP_", assay = "RNA", global = TRUE)
  
  inputExpData[["umap"]]=umap.reduction
  
  #Adding in the embedding data
  embeddings=inputEmbeddings[match(colnames(inputExpData),row.names(inputEmbeddings)),]
  colnames(embeddings)=gsub("PC_","",colnames(embeddings))
  reduction.data <- Seurat::CreateDimReducObject(embeddings = embeddings, 
                                         assay = "RNA",  
                                         key = "PC_")
  inputExpData[["pca"]] <- reduction.data
  
  if(!is.null(n_clusters)){
    clust_obj=.sconline.cluster(argList,n_clusters=n_clusters,clustering_method="average")
    c_assignments=clust_obj$cell_cluster_assignments
    c_assignments=c_assignments[match(colnames(inputExpData),row.names(c_assignments)),]
    Idents(inputExpData)=c_assignments[,1]
  }
  
  return(inputExpData)
}

#Propagates the cell annotations to the pseudocells
#res_prop: output of .sconline.runPropagation(); .sconline.fetch_data("sconline_arrays")
#collapse_datasets: logical; combine the results across datasets or provide the propagation results at the dataset level
#return_plot: return a ggplot object of the results
.sconline.anno2pseudocell=function(argList,annoCol,collapse_datasets=T,return_plot=T,min_effective_size=5){
  
  #argList=.ArgList;annoCol="anno_orig_cellState";collapse_datasets=T;return_plot=T;min_effective_size=5
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  if(sum(is.na(pd[,annoCol]))>0){
    warning("Cells with missing annotations are excluded from the analysis")
    pd=pd[!is.na(pd[,annoCol]),]
  }
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)# Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  if(collapse_datasets){
    
    if(is.numeric(pd[,annoCol])){
      x=matrix(pd[,annoCol],ncol=1)
    } else {
      x=as.matrix(.myOneHotFn(inputVector=pd[,annoCol]))
    }
    
    if(sum(is.na(pd$UMAP_1))>0|sum(colnames(pd)=="UMAP_1")==0){
      warning("UMAP coordinates are missing for some annotations. there might an issue in the data!")
    }
    res=prop_mat %*% x
    
    res=as.data.frame(res)
    
    res$pseudocell=row.names(prop_mat)
    res$effective_size=matEffectiveSize
    
  } else {
    res=list()
    for(ids in unique(pd$anno_batch)){
      tmp_pd=pd[pd$anno_batch==ids,]
      tmp_prop_mat=prop_mat[,match(tmp_pd$sample,colnames(prop_mat))]
      #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
      tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
      
      tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
      tmp_weights[tmp_effsize<4]=0
      
      if(is.numeric(tmp_pd[,annoCol])){
        x=matrix(tmp_pd[,annoCol],ncol=1)
      } else {
        x=as.matrix(.myOneHotFn(inputVector=tmp_pd[,annoCol]))
      }
      
      if(sum(is.na(tmp_pd$UMAP_1))>0|sum(colnames(tmp_pd)=="UMAP_1")==0){
        warning("UMAP coordinates are missing for some annotations. there might an issue in the data!")
      }
      tmp_res=tmp_prop_mat %*% x
      
      if(median(rowSums(tmp_prop_mat))>0.9){
        tmp_res=sweep(tmp_res,1,tmp_weight,"*")
      }
      
      tmp_res=as.data.frame(tmp_res)
      tmp_res$dataset=ids
      
      tmp_res$pseudocell=row.names(tmp_prop_mat)
      tmp_res$effective_size=tmp_effsize
      res=c(res,list(tmp_res))
      
    }
    
    
    res=do.call(eval(parse(text='plyr::rbind.fill')), res)
    res[is.na(res)]=0
    res=res[,c(setdiff(colnames(res),c("dataset","pseudocell","effective_size")),c("dataset","pseudocell","effective_size"))]
    
    if(collapse_datasets){
      res=res[,-which(colnames(res)=="dataset")]
      absent_pseudocells=apply(res[,-which(colnames(res) %in% c("pseudocell","effective_size"))],1,sum)
      absent_pseudocells=which(absent_pseudocells<0.2)
      if(length(absent_pseudocells)>0){
        res[absent_pseudocells,-which(colnames(res)=="pseudocell")]=NA
      }
      
      res=aggregate(.~pseudocell,data=res,function(x) sum(x,na.rm=T))
    }
  }
  
  
  
  p=""
  if(return_plot){
    #
    if(sum(colnames(res)=="cluster")>0){
      p=.extra_sconline.visPseudocellAnno_cluster(inputData=res,argList = argList,min_effective_size=min_effective_size)
    } else {
      p=.extra_sconline.visPseudocellAnno(inputData=res,argList = argList,min_effective_size=min_effective_size)
    }
    
  }
  return(list(results=res,plot=p))
  
}

.sconline.annoEntropy=function(argList,annoCol="anno_batch",distance_measure="entropy"){
  #distance measure can be entropy or cross_entropy. The difference is that cross entropy considers the bkg proportion of cells in the estimate
  
  match.arg(distance_measure,c("entropy","cross_entropy"))
  
  pd=.sconline.fetch_data("annotation",argList = argList)
  p_dist=as.data.frame(table(pd[,annoCol])/nrow(pd))
  
  pseudo_anno=.sconline.anno2pseudocell(argList = argList,annoCol = annoCol,collapse_datasets=T,return_plot=T)
  q_dist=pseudo_anno$results
  sl_cols=setdiff(colnames(q_dist),c("pseudocell","effective_size"))
  p_dist=p_dist[match(sl_cols,p_dist[,1]),]
  
  q_dist[q_dist==0]=0.00001
  q_dist$cross_entropy=apply(q_dist[,sl_cols],1,function(x){
    (-1)*sum(p_dist[,2]*log2(x))
  })
  
  q_dist$entropy=apply(q_dist[,sl_cols],1,function(x) {
    x=x[x>0]
    y=sum(x*log2(x))
    (-1)*y
  })
  
  
  
  pd_summary=.extra_sconline.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"),n=300)
  UMAP_centroid=.sconline.fetch_data("umap_pseudocells",argList=argList)
  UMAP_centroid=UMAP_centroid[match(q_dist$pseudocell,UMAP_centroid$centroid),]
  q_dist$UMAP_1=UMAP_centroid$UMAP_1
  q_dist$UMAP_2=UMAP_centroid$UMAP_2
  
  low_col="red"
  mid_col="white"
  high_col="green"
  if(distance_measure=="cross_entropy"){
    low_col="green"
    high_col="red"
  }
  
  p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),size=0.3,color="gray") + geom_point(data=q_dist,aes_string(x='UMAP_1', y='UMAP_2', fill=distance_measure,size='effective_size'),color="black",shape=21)+
    coord_equal()+theme_classic()+scale_fill_gradient2(low=low_col,high=high_col,mid=mid_col,midpoint = mean(c(min(q_dist[,distance_measure]),max(q_dist[,distance_measure]))))
  
  
  return(list(data=q_dist[,colnames(q_dist) %in% c("pseudocell","effective_size", "cross_entropy", 
                                                   "entropy", "UMAP_1","UMAP_2")],plot=p))
  
}

.sconline.anno2pseudocell_tmp=function(res_prop,argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T,min_effective_size=5,pd=NULL,subset_pd=F){
  if(is.null(pd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  }
  
  
  res=list()
  for(i in 1:length(res_prop)){
    
    effective_size=res_prop[[i]]$matEffectiveSize
    
    tmp_pd=pd
    
    if(is.numeric(tmp_pd[,annoCol])){
      x=matrix(tmp_pd[,annoCol],ncol=1)
    } else {
      x=as.matrix(.myOneHotFn(inputVector=tmp_pd[,annoCol]))
    }
    tmp_pd=tmp_pd[match(colnames(res_prop[[i]]$prop_mat),row.names(pd)),]
    x=x[match(colnames(res_prop[[i]]$prop_mat),row.names(pd)),]
    if(sum(is.na(tmp_pd$UMAP_1))>0){
      warning("UMAP coordinates are missing for some annotations. there might an issue in the data!")
    }
    tmp_res=res_prop[[i]]$prop_mat %*% x
    tmp_effective_size=res_prop[[i]]$matEffectiveSize
    tmp_effective_size=tmp_effective_size[match(row.names(tmp_res),names(tmp_effective_size))]
    if(median(rowSums(res_prop[[i]]$prop_mat))>0.9){
      tmp_weight=res_prop[[i]]$matWeights
      tmp_weight=tmp_weight[match(row.names(tmp_res),names(tmp_weight))]
      
      tmp_res=sweep(tmp_res,1,tmp_weight,"*")
    }
    
    
    tmp_res=as.data.frame(tmp_res)
    if(sum(names(res_prop)=="zscore")==0){
      tmp_res$dataset=res_prop[[i]]$data$dsName
    } else {
      tmp_res$dataset=names(res_prop)[i]
    }
    
    tmp_res$pseudocell=row.names(res_prop[[i]]$prop_mat)
    tmp_res$effective_size=tmp_effective_size
    res=c(res,list(tmp_res))
  }
  
  res=as.data.frame(data.table::rbindlist(res))
  
  if(collapse_datasets){
    res=res[,-which(colnames(res)=="dataset")]
    absent_pseudocells=apply(res[,-which(colnames(res) %in% c("pseudocell","effective_size"))],1,sum)
    absent_pseudocells=which(absent_pseudocells<0.2)
    if(length(absent_pseudocells)>0){
      res[absent_pseudocells,-which(colnames(res)=="pseudocell")]=NA
    }
    
    res=aggregate(.~pseudocell,data=res,function(x) sum(x,na.rm=T))
  }
  
  p=""
  if(return_plot){
    if(subset_pd){
      p=.extra_sconline.visPseudocellAnno_cluster(inputData=res,argList = argList,min_effective_size=min_effective_size,pd=pd)
    } else {
      p=.extra_sconline.visPseudocellAnno_cluster(inputData=res,argList = argList,min_effective_size=min_effective_size,pd=NULL)
    }
    
  }
  return(list(results=res,plot=p))
  
}


#Fetch the data based on the argList
#if dataType=NULL, prints the datatypes that can be return by the object
.sconline.fetch_data=function(dataType=NULL,argList=NULL){
  result=""
  if(is.null(dataType)){
    cat("*********\nannotation\nembeddings\npca\nexpression\npca_centroids\numap_pseudocells\nmarkers\nsconline_arrays\nmeta_z\n*********\n")
  } else {
    if(is.null(argList)){
      stop("argList should be provided!")
    }
    result = switch(  
      dataType,  
      "annotation"= {
        if(file.exists(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))){
          load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
        } else {
          load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
        }
        
        pd
      },  
      "pca_centroids"= {
        load(.myFilePathMakerFn("pca_centroids",argList=argList))
        pca_centroid
      },  
      "embeddings"= {
        load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
        harmony_embeddings
      },  
      "expression"= {
        expData=qs::qread(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
        expData
      },
      "umap_pseudocells"= {
        load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
        UMAP_centroid
      },
      "pca"= {
        load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
        pca_res
      },
      "umap"= {
        load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
        pd=as.matrix(pd[,c("UMAP_1","UMAP_2")])
        pd
      },
      "markers"= {
        load(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T,propImportant = T))
        res_arranged
      },
      "sconline_arrays"= {
        data=qread(.myFilePathMakerFn("res_dataset_array",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
        data
      },
      "meta_z"= {
        data=qread(.myFilePathMakerFn("res_meta",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
        data
      }
    )  
  }
  return(result)
}

#Correlation of pseudocells at the cell level
#res_prop: output of .sconline.runPropagation(); .sconline.fetch_data("sconline_arrays")
#plot_save_path: if provided, saves a heatmap in the path
.sconline.pseudocellCor=function(argList,annoCol=NULL,plot_save_path=NULL){
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  res_cor=as.matrix(qlcMatrix::cosSparse(t(prop_mat)))
  row.names(res_cor)=colnames(res_cor)=row.names(prop_mat)
  diag(res_cor)=0
  axis_col=rep("white",nrow(res_cor))
  if(!is.null(annoCol)){
    anno_data=.sconline.anno2pseudocell(argList=argList,annoCol = annoCol,collapse_datasets = T,return_plot = F)
    anno_data=anno_data$results[match(row.names(res_cor),anno_data$results$pseudocell),]
    anno_data=anno_data[,-which(colnames(anno_data) %in% c("pseudocell","effective_size"))]
    anno_data=colnames(anno_data)[apply(anno_data,1,function(x) {
      if(max(x)>0.5){
        which(x==max(x))[1]
      } else {
        #NA
        which(x==max(x))[1]
      }
    })]
    anno_data[is.na(anno_data)]="white"
    anno_data2=as.numeric(factor(anno_data))
    axis_col=c("gray",hues::iwanthue(length(unique(anno_data2))))[anno_data2]
    axis_col[anno_data=="white"]="white"
  }
  
  colramp = colorRampPalette(c("dodgerblue","black","yellow"))(9)
  
  singleton_pseudocells=which(rowSums(res_cor)==1)
  if(length(singleton_pseudocells)>0&length(singleton_pseudocells)<nrow(res_cor)){
    res_cor=res_cor[-singleton_pseudocells,-singleton_pseudocells]
    axis_col=axis_col[-singleton_pseudocells]
    
  }
  
  if(!is.null(plot_save_path)){
    pdf(file = "~/myBucket/torm.pdf")
    heatmap(res_cor,col=colramp, RowSideColors = axis_col,distfun = function(x) as.dist(1-cor(t(x))),hclustfun = function(x) hclust(x, method="ward.D2"), ColSideColors = axis_col, margins = c(5,10),scale="none")
    dev.off()
  }
  
  return(res_cor)
}

#Plot the pseudocells in the umap space
.sconline.plot_pseudocell_umap=function(argList,selectedPseudocells=NULL){
  require(ggplot2)
  umap_data=.sconline.fetch_data(dataType="umap_pseudocells",argList=argList)
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  umap_data=umap_data[umap_data$centroid %in% row.names(prop_mat),]
  p=ggplot(umap_data,aes(UMAP_1,UMAP_2,label=centroid))
  
  if(!is.null(selectedPseudocells)){
    p=p+geom_point()+geom_label(data=umap_data[which(umap_data$centroid %in% as.character(selectedPseudocells)),],color="red")
  } else {
    p=p+geom_label(data=umap_data)
  }
  
  p=p+theme_classic()
  
  return(p)
  
}


.sconline.HVGselection=function(argList,inputExpData,batch_variable="anno_batch",run_harmony=F,L2Norm=T){
  
  #argList=.ArgList;inputExpData=data;batch_variable="anno_batch";run_harmony=F;L2Norm=T
  #argList$input_highly_var_genes=NULL
  
  if(!is.null(argList$input_highly_var_genes)){
    stop("Highly variable genes are already specified in the argList")
  }
  
  argList$prefix="tmp_HVG_analysis"
  data=.extra_sconline.exp_creatorFn(argList=argList,inputExpData=inputExpData,batch_variable=batch_variable,organism = "unknown",addAnno=F)
  cat('running pca step\n')
  
  cat('-- Selecting the highly variable genes\n')
  data=.myHighVarGeneSlFn(data,dataorganism="unknwon",batch_variable=batch_variable,argList = argList)
  
  cat('-- Scaling the data and PCA\n')
  pca_data=.myPCAfn(data,argList = argList,saveFiles=F)
  
  argList$pseudocell_size=round((sum(unlist(lapply(data$data,ncol))))/300)
  #inputExpData = data;inputEmbeddings=pca_data$pca_res;inputPhenoData=as.data.frame(pca_data$pd);run_harmony=run_harmony;batch_variable=batch_variable;generateUMAP = F;saveFiles = F;mergePseudocells=T;hierarchical_refinement=F;colNormalize=F
  prop_mat=.sconline.runPropagation(argList = argList,inputExpData = data,inputEmbeddings=pca_data$pca_res,inputPhenoData=as.data.frame(pca_data$pd),run_harmony=run_harmony,batch_variable=batch_variable,generateUMAP = F,saveFiles = F,mergePseudocells=T,hierarchical_refinement=F,colNormalize=F)
  
    
    return(res_prop)
}


.sconline.subset_purityAnalysis=function(argList,sl_cells=NULL,inputPCAdata=NULL,batch_variable="anno_batch",collapse_datasets=T,minCellCountThr=4,analysis_seed=1,extendedMode=F,cluster_count=NULL){
  #Running pca
  library(rliger)
  library(Seurat)
  library(scater)
  library(scran)
  
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
    # TODO jonah parallize buyt make sure chunked
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
  
  #organism="unknown";addAnno=F
  
  if(is.null(inputPCAdata)){
    inputExpData=.sconline.fetch_data("expression",argList = argList)
    if(!is.null(sl_cells)){
      print(paste0("Selecting ",sum(colnames(inputExpData) %in% sl_cells)," out of ",ncol(inputExpData)," cells."))
      if(sum(colnames(inputExpData) %in% sl_cells)==0){
        stop("None of the selected cells was identified in the dataset!")
      }
      inputExpData=inputExpData[,colnames(inputExpData) %in% sl_cells]
    } else {
      stop("sl_cells argument was not provided!")
    }
    
    data=.extra_sconline.exp_creatorFn(argList=argList,inputExpData=inputExpData,batch_variable=batch_variable,organism = "unknown",addAnno=F)
    argList2=argList
    argList2$HVG_list=argList2$HVG_count=min(max(round(length(data$data)/2),1),argList$HVG_count)
    argList2$input_highly_var_genes=NULL
    cat('Re-running pca step\n')
    
    cat('-- Selecting the highly variable genes\n')
    data=.myHighVarGeneSlFn(data,dataorganism="unknwon",argList = argList2)
    
    cat('-- Scaling the data and PCA\n')
    pca_data=.myPCAfn(data,argList = argList2,saveFiles=F)
    
  } else {
    argList$nPCs=ncol(inputPCAdata)-5
    pca_res=Seurat:::RunPCA.default(object=myPrepDR(scaledData=t(inputPCAdata),
                                                    features=colnames(inputPCAdata), verbose = FALSE),
                                    npcs = max(50,argList$nPCs+10))
    pca_data=list(pca_res=pca_res@cell.embeddings)
  }
  
  
  
  cat("Running Propagation step on the subset\n")
  
  set.seed(analysis_seed)
  supportingFractionThr=argList$DE_supportingFractionThr
  n.adaptiveKernel=argList$DE_n.adaptiveKernel
  nPropIter=argList$DE_nPropIter
  
  #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  pd=pd[which(row.names(pd) %in% row.names(pca_data$pca_res)),]
  
  load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
  sl_pseudo=sl_pseudo[which(sl_pseudo$pseudocell %in% row.names(pca_data$pca_res)),]
  pca_centroid=pca_data$pca_res[sl_pseudo$pseudocell,]
  row.names(pca_centroid)=sl_pseudo$cluster
  
  
  harmony_embeddings=pca_data$pca_res[row.names(pd),]
  if(sum(is.na(harmony_embeddings))>0){
    stop("Error in matching Names!")
  }
  
  pcaList=split(as.data.frame(harmony_embeddings),pd[,batch_variable])
  
  
  dataArranged = parallel::mclapply(names(pcaList), function(thisExp){
    
    return(list(
      dsName=thisExp,
      pcaData=pcaList[[thisExp]]))
  })
  
  if(argList$do.split.prop&F){
    res=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final_v2_split_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=argList$prop.n.neighbors,mc.cores = argList$ncores)
  } else {
    #centroidPCAdata=pca_centroid;argList=argList;exCentroids=NULL;runIndx=1;n.neighbors=argList$prop.n.neighbors;batchPCAdata=harmony_embeddings
    res=.myConcensusDEFn_step2_detail_newprop3_final_v11subset(dataArranged=dataArranged,centroidPCAdata=pca_centroid,argList=argList,exCentroids=NULL,n.neighbors=argList$prop.n.neighbors,batchPCAdata=harmony_embeddings,returnPropMat=T,cluster_count=cluster_count)
  }
  
  return(res)
  
}

.extra_hierarchical_clustList=function(x){
  res=NULL
  if(sum(names(x$res)=="object")>0){
    for(i in 1:length(x$res)){
      tmp=.extra_hierarchical_clustList(x$res[[i]])
      tmp$cell_index=paste0(x$index,"_",tmp$cell_index)
      res=rbind(res,tmp)
    }
  } else {
    if(sum(names(x$res)=="res")>0){
      res=x$res$res
    } else {
      res=x$res
    }
  }
  
  return(res)
  
}


.extra_hierarchical_clustMerge=function(argList,meta_data,pd,input_embeddings,prop_mat,input_pca_centroids,
                                        sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,cell_count_diff_thr=0,combinatorial_pct_tol=1,forgiveness_factor=1,tol_level=0.9){
  
  
  
  de_pct_res=.extra_sconline.PctScoreFn(pct_mat=meta_data$med_pct.1,argList = argList,meta_z_mat=meta_data$meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F,cell_count_diff_thr=cell_count_diff_thr,hierarchical_mode=T)
  
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  
  #affinity_param=min(affinity_param,quantile(as.numeric(de_pct_res),0.1))
  #affinity_param=max(affinity_param,2)
  input_embeddings=input_embeddings[match(colnames(prop_mat),row.names(input_embeddings)),]
  
  ps_sim_mat=.extra_sconline.pseudosim_archive11(argList=argList,input_prop_mat=prop_mat,hierarchical_mode=T,binarize = ncol(prop_mat)<100000,cos_dist = T,input_pca_centroids=input_pca_centroids,input_embeddings=input_embeddings)
  
  ps_sim_mat=ps_sim_mat[row.names(de_pct_res),colnames(de_pct_res)]
  
  
  
  sl_param=NULL
  sl_score=0
  for(affinity_param in 2:20){
    de_sim_mat=exp(-1*(de_pct_res)/affinity_param)
    
    tst=cor(as.matrix(de_sim_mat),as.matrix(ps_sim_mat))
    diag(tst)=0
    tmp_score=median(apply(tst,1,max),na.rm = T)
    if(tmp_score>sl_score){
      sl_score=tmp_score
      sl_param=affinity_param
    }
  }
  sl_param1=sl_param
  
  tst=de_pct_res
  diag(tst)=100
  sl_param=apply(as.matrix(tst),1,function(x) {#quantile(x,0.01);
    x=x[order(x,decreasing = F)]
    (x[2]+quantile(x,0.01))/2})
  sl_param=median(sl_param)
  sl_param=max(sl_param,2)
  if(!is.null(sl_param1)){
    sl_param=(sl_param+sl_param1)/2
  }
  
  
  #
  de_sim_mat=exp(-1*(de_pct_res)/sl_param)
  
  diag(ps_sim_mat)=0
  diag(de_sim_mat)=0
  ps_max_vals=qlcMatrix::rowMax(ps_sim_mat)
  de_max_vals=qlcMatrix::rowMax(de_sim_mat)
  ps_sim_mat=.extra_matrix_rowNorm(input_mat = ps_sim_mat,rowValues = 1/(ps_max_vals+0.0000001))#Matrix::Diagonal(x=1/(ps_max_vals+0.0000001)) %*% ps_sim_mat
  de_sim_mat=.extra_matrix_rowNorm(input_mat = de_sim_mat,rowValues = 1/(de_max_vals+0.0000001))#Matrix::Diagonal(x=1/(de_max_vals+0.0000001)) %*% de_sim_mat
  ps_sim_mat=sweep(as.matrix(ps_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
  de_sim_mat=sweep(as.matrix(de_sim_mat),1,(de_max_vals+ps_max_vals)/2,"*")
  diag(ps_sim_mat)=1
  diag(de_sim_mat)=1
  diff=((1-ps_sim_mat)+(1-de_sim_mat))
  
  diff=pmax(diff,0)
  diff=diff+t(diff)
  diag(diff)=0
  res_eff_size=.myEffSizePropMat(prop_mat)$effective_sample_size
  names(res_eff_size)=row.names(prop_mat)
  diff_clust=hclust(as.dist(diff),method = "average",members = res_eff_size[colnames(diff)])
  
  
  d_conMat=1:nrow(prop_mat)
  names(d_conMat)=row.names(prop_mat)
  
  pd=pd[match(colnames(prop_mat),row.names(pd)),]
  d_conMat=.sconline.cluster_pruning(cluster_assignments=d_conMat,clust_obj=diff_clust,
                                     inputExpData=inputExpData,
                                     argList=argList,
                                     combinatorial_pct_tol=combinatorial_pct_tol,marker_sig1_thr=sig1_thr,marker_pct2_thr=pct2_thr,marker_pct_diff_thr=pct_diff_thr,forgiveness_factor=forgiveness_factor,input_meta_data=meta_data,inputPhenoData=pd,input_prop_mat=prop_mat)
  
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),row.names(pd)]
  prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  #prop anno
  
  colMax_vals_m=qlcMatrix::colMax(prop_merged)
  colMax_vals_m=.extra_matrix_colNorm(prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=Matrix::drop0(colMax_vals_m,tol=tol_level)
  prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
  prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
  prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
  prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
  pd2=data.frame(cluster=paste0("C",as.character(prop_m_hardCluster$i)))
  row.names(pd2)=row.names(pd)
  pseudocell_cluster_assignments=data.frame(pseudocell=names(d_conMat),cluster=paste0("C",d_conMat),stringsAsFactors=F)
  
  return(list(cluster_object=diff_clust,distance_matrix=diff,cell_cluster_assignments=pd2,pseudocell_cluster_assignments=pseudocell_cluster_assignments))
  
}



.extra_hierarchical_run=function(argList,inputEmbeddings=NULL,inputPhenoData=NULL,organism,batch_variable,inputExpData,L2Norm,mergePseudocells,merging_strength,sig1_thr,pct2_thr,pct_diff_thr,run_harmony=F,index=NULL,dataset_size_thr=20000,min_cluster_size=500,input_res_embeddings=NULL,input_res_prop=NULL){
  
  
  
  if(is.null(index)){
    index=1
  }
  
  if(ncol(inputExpData[[1]])<min_cluster_size){
    res=list(index=index,res=data.frame(cluster="C1",cell=colnames(inputExpData[[1]]),cell_index=index))
    return(res)
  }
  
  if(is.null(input_res_prop)){
    if(is.null(input_res_embeddings)){
      if(is.null(inputEmbeddings)){
        #inputExpData = inputExpData[[1]];argList = argList;batch_variable = batch_variable;organism=organism;hierarchical_mode = T
        
        pca_res=.sconline.pca(inputExpData = inputExpData[[1]],argList = argList,batch_variable = batch_variable,organism=organism,hierarchical_mode = T)
        pca=pca_res$pca_res
        pd=pca_res$pd
      } else {
        if(is.null(inputPhenoData)){
          stop("inputPhenoData matching with the inputEmbeddings should be provided")
        }
        pca=inputEmbeddings
        pd=inputPhenoData
      }
      
      if(nrow(pca)<1000){
        argList$pseudocell_size=50
      } else if(nrow(pca)<25000){
        argList$pseudocell_size=100
      }
      
      if(run_harmony=="auto"){
        run_harmony=dim(inputExpData[[1]]<20000)
      }
      #inputEmbeddings=pca;run_harmony=run_harmony;inputBatchCol=batch_variable;saveFiles = F
      res_embeddings=.sconline.embeddingFn(argList,inputEmbeddings=pca,run_harmony=run_harmony,pd=pd,inputBatchCol=batch_variable,saveFiles = F)
      
    } else {
      res_embeddings=input_res_embeddings
    }
    
    inputExpData[[1]]=inputExpData[[1]][,colnames(inputExpData[[1]]) %in% row.names(res_embeddings$harmony_embeddings)]
    inputPhenoData=inputPhenoData[match(colnames(inputExpData[[1]]),row.names(inputPhenoData)),]
    res_embeddings$harmony_embeddings=res_embeddings$harmony_embeddings[match(colnames(inputExpData[[1]]),row.names(res_embeddings$harmony_embeddings)),]
    
    .tst=list(res_embeddings=res_embeddings,inputPhenoData=inputPhenoData,inputExpData=inputExpData)
    #qsave(list(res_prop=res_prop,tst=.tst,argList=argList),file="~/torm.qs")
    #argList = argList;expData = inputExpData;hierarchical_mode = T;addClusteringModule=F;extendedMode = F;L2Norm=L2Norm;mergePseudocells=mergePseudocells;merging_strength=merging_strength;sig1_thr=sig1_thr;pct2_thr=pct2_thr;pct_diff_thr=pct_diff_thr;input_embedding=res_embeddings$harmony_embeddings;input_pca_centroid=res_embeddings$pca_centroid;input_pca_centroids_assignments=res_embeddings$sl_pseudo;input_pd=inputPhenoData
    #.tst=list(argList = argList,expData = inputExpData,hierarchical_mode = T,addClusteringModule=F,extendedMode = F,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,input_embedding=res_embeddings$harmony_embeddings,input_pca_centroid=res_embeddings$pca_centroid,input_pca_centroids_assignments=res_embeddings$sl_pseudo,input_pd=inputPhenoData)
    #argList=.tst$argList;expData=.tst$expData;hierarchical_mode = T;addClusteringModule=F;extendedMode = F;L2Norm=L2Norm;mergePseudocells=mergePseudocells;merging_strength=merging_strength;input_embedding=.tst$input_embedding;input_pca_centroid=.tst$input_pca_centroid;input_pca_centroids_assignments=.tst$input_pca_centroids_assignments;input_pd=.tst$input_pd
    res_prop=.myConcensusDEFn_step2_fast(argList = argList,expData = inputExpData,hierarchical_mode = T,addClusteringModule=F,extendedMode = F,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,input_embedding=res_embeddings$harmony_embeddings,input_pca_centroid=res_embeddings$pca_centroid,input_pca_centroids_assignments=res_embeddings$sl_pseudo,input_pd=inputPhenoData)
    
  } else {
    res_prop=input_res_prop
  }
  gc()
  #input_res=res_prop;inputPhenoData = inputPhenoData;input_pca_centroids=res_embeddings$pca_centroid;input_embeddings=res_embeddings$harmony_embeddings
  #.tst=list(argList=argList,input_res=res_prop,inputPhenoData = inputPhenoData,input_pca_centroids=res_embeddings$pca_centroid,input_embeddings=res_embeddings$harmony_embeddings)
  
  if(nrow(res_prop$prop_mat)==1){
    clust_res=data.frame(cluster=clust_res$cell_cluster_assignments[,1],cell=colnames(res_prop$prop_mat),stringsAsFactors = F)
    
  } else if(nrow(res_prop$prop_mat)==2){
    #meta_data=res_prop$res_meta;pd=inputPhenoData;input_embeddings=res_embeddings$harmony_embeddings;prop_mat=res_prop$prop_mat;input_pca_centroids=res_embeddings$pca_centroid
    clust_res=.extra_hierarchical_clustMerge(argList=argList,meta_data=res_prop$res_meta,pd=inputPhenoData,input_embeddings=res_embeddings$harmony_embeddings,prop_mat=res_prop$prop_mat,input_pca_centroids=res_embeddings$pca_centroid,
                                             sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,cell_count_diff_thr=0)
    clust_res=data.frame(cluster=clust_res$cell_cluster_assignments[,1],cell=row.names(clust_res$cell_cluster_assignments),stringsAsFactors = F)
    
  } else {
    clust_res=.extra_hierarchical_split(argList=argList,input_res=res_prop,inputPhenoData = inputPhenoData,input_pca_centroids=res_embeddings$pca_centroid,input_embeddings=res_embeddings$harmony_embeddings)
  }
  
  
  if(F){
    #tmp=SingleCellExperiment(assays = list(counts = res),colData = pd,rowData=fd)
    tmp=.extraExport2SeuratFn(tmp_inputExpData[[1]])
    tmp=Seurat::NormalizeData(tmp)
    tmp=Seurat::FindVariableFeatures(tmp)
    tmp=Seurat::ScaleData(tmp)
    tmp=Seurat::RunPCA(tmp)
    tmp=Seurat::RunUMAP(tmp,dims=1:30)
    Idents(tmp)=tmp$clusters
    table(Idents(tmp))
    p=Seurat::DimPlot(tmp,label = T)
    ggsave(plot=p,file="~/myBucket/torm.pdf")
  }
  
  #tst=variancePartition::fitExtractVarPartModel(exprObj= t(pca), form= ~(1|anno_batch), data=pd )
  
  res=list()
  if(length(unique(clust_res[,"cluster"]))>1){
    clust_indx=0
    for(iclust in unique(clust_res[,"cluster"])){
      clust_indx=clust_indx+1
      sl_cells=clust_res[which(clust_res[,"cluster"]==iclust),2]
      
      tmp_phenoData=inputPhenoData[sl_cells,]
      tmp_inputExpData=inputExpData
      for(iexp in 1:length(tmp_inputExpData)){
        tmp_inputExpData[[iexp]]=tmp_inputExpData[[iexp]][,colnames(tmp_inputExpData[[iexp]]) %in% sl_cells]
      }
      
      tmp_res_embeddings=NULL
      if(length(sl_cells)/ncol(res_prop$prop_mat)>=0.9){
        tmp_res_embeddings=res_embeddings
        tmp_res_embeddings$harmony_embeddings=tmp_res_embeddings$harmony_embeddings[match(colnames(tmp_inputExpData[[1]]),row.names(tmp_res_embeddings$harmony_embeddings)),]
        tmp_res_embeddings$sl_pseudo=tmp_res_embeddings$sl_pseudo[which(tmp_res_embeddings$sl_pseudo$pseudocell %in% colnames(tmp_inputExpData[[1]])),]
        tmp_res_embeddings$pca_centroid=tmp_res_embeddings$pca_centroid[which(row.names(tmp_res_embeddings$pca_centroid) %in% tmp_res_embeddings$sl_pseudo$cluster),]
      }
      
      tmp_res_prop=NULL
      if(length(sl_cells)/ncol(res_prop$prop_mat)>=0.95){
        tmp_res_prop=res_prop
        tmp_res_prop$prop_mat=tmp_res_prop$prop_mat[,match(colnames(tmp_inputExpData[[1]]),colnames(tmp_res_prop$prop_mat))]
        tmp_res_prop$prop_mat=tmp_res_prop$prop_mat[row.names(tmp_res_prop$prop_mat) %in% input_res_embeddings$sl_pseudo$cluster,]
        tmp_res_prop$prop_mat=.extra_matrix_rowNorm(tmp_res_prop$prop_mat)#Matrix::Diagonal(x=1/rowSums(tmp_res_prop$prop_mat)) %*% tmp_res_prop$prop_mat
      }
      tmp_clust_res=tryCatch({.extra_hierarchical_run(argList,inputEmbeddings=NULL,inputPhenoData=tmp_phenoData,organism=organism,batch_variable=batch_variable,inputExpData=tmp_inputExpData,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,run_harmony=run_harmony,index=clust_indx,min_cluster_size=min_cluster_size,input_res_embeddings = tmp_res_embeddings,input_res_prop=tmp_res_prop)}, error=function(e) {return(T)})
      #inputEmbeddings=NULL;inputPhenoData=tmp_phenoData;organism=organism;batch_variable=batch_variable;inputExpData=tmp_inputExpData;L2Norm=L2Norm;mergePseudocells=mergePseudocells;merging_strength=merging_strength;sig1_thr=sig1_thr;pct2_thr=pct2_thr;pct_diff_thr=pct_diff_thr;run_harmony=run_harmony;index=clust_indx;min_cluster_size=min_cluster_size;input_res_embeddings = tmp_res_embeddings;input_res_prop=tmp_res_prop
      if(class(tmp_clust_res)==class(T)){
        qsave(list(sl_cells=sl_cells,pheno=tmp_phenoData,exp=tmp_inputExpData,argList=argList,res_prop=res_prop,res_embeddings=res_embeddings),file="~/torm.qs")
        stop("Error!")
        
        #tmp=qread("~/torm.qs")
        #sl_cells=tmp$sl_cells;tmp_phenoData=tmp$pheno;tmp_inputExpData=tmp$exp;res_prop=tmp$res_prop;res_embeddings=tmp$res_embeddings
      }
      
      
      res=c(res,list(object=list(index=index,res=tmp_clust_res)))
    }
  } else{
    clust_res$cell_index=index
    res=list(object=list(index=index,res=clust_res))
  }
  
  
  if(length(res)>1){
    qsave(res,file="~/torm.qs")
    run_res=list(res=res,index=index)
    clustList=.extra_hierarchical_clustList(run_res)
    if(length(unique(clustList$index))>2){
      prop_mat=clustList$index
      names(prop_mat)=clustList$cell
      prop_mat=t(as.matrix(.myOneHotFn(prop_mat)))
      tmp_exp=inputExpData[[1]]
      tmp_exp=tmp_exp[,match(colnames(prop_mat),colnames(tmp_exp))]
      tmp_exp=list(tmp_exp)
      tmp_pheno=inputPhenoData[match(colnames(prop_mat),row.names(inputPhenoData)),]
      tmp_embedding=res_embeddings$harmony_embeddings
      tmp_embedding=tmp_embedding[match(colnames(prop_mat),row.names(tmp_embedding)),]
      #argList = argList;input_prop_mat = prop_mat;expData = tmp_exp;hierarchical_mode = T;addClusteringModule=F;extendedMode = F;L2Norm=L2Norm;mergePseudocells=mergePseudocells;merging_strength=merging_strength;sig1_thr=sig1_thr;pct2_thr=pct2_thr;pct_diff_thr=pct_diff_thr;input_embedding=tmp_embedding;input_pca_centroid=NULL;input_pca_centroids_assignments=NULL;input_pd=tmp_pheno
      res_prop=.myConcensusDEFn_step2_fast(argList = argList,input_prop_mat = prop_mat,expData = tmp_exp,hierarchical_mode = T,addClusteringModule=F,extendedMode = F,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,input_embedding=tmp_embedding,input_pca_centroid=NULL,input_pca_centroids_assignments=NULL,input_pd=tmp_pheno)
      gc()
      #input_res=res_prop;inputPhenoData = inputPhenoData;input_pca_centroids=res_embeddings$pca_centroid;input_embeddings=res_embeddings$harmony_embeddings
      
      #.tst=list(argList=argList,input_res=res_prop,inputPhenoData = inputPhenoData,input_pca_centroids=res_embeddings$pca_centroid,input_embeddings=res_embeddings$harmony_embeddings,prune_clusters=F)
      #argList=.tst$argList;input_res=.tst$input_res;inputPhenoData = .tst$inputPhenoData;input_pca_centroids=.tst$input_pca_centroids;input_embeddings=.tst$input_embeddings;prune_clusters=T
      
      #argList=argList;meta_data=res_prop$res_meta;pd=inputPhenoData;input_embeddings=res_embeddings$harmony_embeddings;prop_mat=res_prop$prop_mat;input_pca_centroids=res_embeddings$pca_centroid;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;cell_count_diff_thr=0
      clust_res=.extra_hierarchical_clustMerge(argList=argList,meta_data=res_prop$res_meta,pd=inputPhenoData,input_embeddings=res_embeddings$harmony_embeddings,prop_mat=res_prop$prop_mat,input_pca_centroids=res_embeddings$pca_centroid,
                                                        sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,cell_count_diff_thr=0)
      cell_index=clustList[match(row.names(clust_res$cell_cluster_assignments),clustList$cell),]
      res=list(object=list(index=index,res=data.frame(cluster=clust_res$cell_cluster_assignments[,1],cell=row.names(clust_res$cell_cluster_assignments),cell_index=cell_index$cell_index,stringsAsFactors = F)))
    }
    
  }
  
  return(res)
}

.extra_hierarchical_split=function(argList,input_res,inputPhenoData,input_pca_centroids,input_embeddings,forgiveness_factor=1,min_cluster_size=500,prune_clusters=F){
  
  
  #n_clusters=5;clustering_method="average";input_meta_z=input_res$res_meta;input_prop_mat=input_res$prop_mat;hierarchical_mode=T;inputPd=inputPhenoData;input_pca_centroids=input_pca_centroids;input_embeddings=input_embeddings
  if(ncol(input_res$prop_mat)<min_cluster_size){
    return(data.frame(cluster="C1",cell=colnames(input_res$prop_mat)))
  }
  #n_clusters=min(nrow(input_res$prop_mat),5);clustering_method="average";input_meta_z=input_res$res_meta;input_prop_mat=input_res$prop_mat;hierarchical_mode=T;inputPd=inputPhenoData
  clust_res=.sconline.cluster(argList,n_clusters=min(nrow(input_res$prop_mat),5),clustering_method="average",input_meta_z=input_res$res_meta,input_prop_mat=input_res$prop_mat,hierarchical_mode=T,inputPd=inputPhenoData,input_pca_centroids=input_pca_centroids,input_embeddings=input_embeddings,prune_clusters = prune_clusters)
  
  d_conMat=cutree(clust_res$cluster_object,k=min(length(clust_res$cluster_object$order),5))
  
  #inputExpData=inputExpData[[1]];combinatorial_pct_tol=1;marker_sig1_thr=3;marker_pct2_thr=0.3;marker_pct_diff_thr=0.2;cluster_assignments=d_conMat;fast_mode=F;input_meta_z=input_res$res_meta$meta_z;inputPhenoData=inputPhenoData;input_prop_mat=input_res$prop_mat;hierarchical_mode=T
  marker_res=.sconline.marker.combinatorial_count(argList=argList,combinatorial_pct_tol=1,marker_sig1_thr=3,marker_pct2_thr=0.3,marker_pct_diff_thr=0.2,cluster_assignments=d_conMat,fast_mode=F,input_meta_data=input_res$res_meta,inputPhenoData=inputPhenoData,input_prop_mat=input_res$prop_mat)
  
  if(sum(marker_res[,2]<10)>0){
    d_conMat=cutree(clust_res$cluster_object,k=min(nrow(input_res$prop_mat),2))
    
    #inputExpData=inputExpData[[1]];combinatorial_pct_tol=1;marker_sig1_thr=3;marker_pct2_thr=0.3;marker_pct_diff_thr=0.2;cluster_assignments=d_conMat;fast_mode=F;input_meta_z=input_res$res_meta$meta_z;inputPhenoData=inputPhenoData;input_prop_mat=input_res$prop_mat;hierarchical_mode=T
    marker_res=.sconline.marker.combinatorial_count(argList=argList,combinatorial_pct_tol=1,marker_sig1_thr=3,marker_pct2_thr=0.3,marker_pct_diff_thr=0.2,cluster_assignments=d_conMat,fast_mode=F,input_meta_data=input_res$res_meta,inputPhenoData=inputPhenoData,input_prop_mat=input_res$prop_mat)
    if(sum(marker_res[,2]==0)>forgiveness_factor){
      tmp=rep(1,length(d_conMat))
      names(tmp)=names(d_conMat)
      d_conMat=tmp
    }
  }
  
  if(length(unique(d_conMat))>1){
    prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
    prop_mat=.extra_matrix_rowNorm(input_mat = input_res$prop_mat,rowValues = 1/as.numeric(qlcMatrix::rowMax(input_res$prop_mat)))# Matrix::Diagonal(x=1/qlcMatrix::rowMax(input_res$prop_mat)) %*% input_res$prop_mat
    if(sum(colnames(inputPhenoData)=="sample")==1&length(setdiff(inputPhenoData$sample,colnames(prop_mat)))==0){
      prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),inputPhenoData$sample,drop=F]
    } else {
      prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),row.names(inputPhenoData)]
    }
    
    prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
    #prop anno
    
    colMax_vals_m=qlcMatrix::colMax(prop_merged)
    colMax_vals_m=.extra_matrix_colNorm(input_mat = prop_merged,colValues = 1/as.numeric(colMax_vals_m))#prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=0.99)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),,drop=F]
    prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(inputPhenoData),prop_m_hardCluster$j),,drop=F]
    prop_m_hardCluster$cell=colnames(prop_mat)[prop_m_hardCluster[,2]]
    prop_m_hardCluster$cluster=paste0("C",prop_m_hardCluster[,"i"])
    prop_m_hardCluster=prop_m_hardCluster[,c("cluster","cell")]
  } else {
    prop_m_hardCluster=data.frame(cluster="C1",cell=colnames(input_res$prop_mat))
  }
  
  return(prop_m_hardCluster)
  
}

.sconline.sconline.runPropagation_hierarchical=function(argList,inputEmbeddings=NULL,inputPhenoData=NULL,inputExpData=NULL,input_UMAP_embedding=NULL,organism,batch_variable="anno_batch",run_harmony=F,addAnno=F,indScaling=F,addClusteringModule=F,umap.method='uwot',extendedMode=F,L2Norm=T,mergePseudocells=T,generateUMAP=T,merging_strength=0.3,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,min_dataset_size=30,min_cluster_size=500){
  require(Matrix)
  require(qlcMatrix)
  
  #argList=.ArgList;inputEmbeddings=NULL;inputPhenoData=NULL;inputExpData=data;organism="Human";batch_variable="anno_batch";run_harmony=F;addAnno=F;addClusteringModule=F;L2Norm=T;mergePseudocells=T;generateUMAP=F;extendedMode=F;merging_strength=0.3;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;input_UMAP_embedding=NULL
  #run_harmony=T
  #umap.method="uwot";generateUMAP=T
  #pd=inputPhenoData;inputBatchCol=batch_variable
  
  argList$indScaling=indScaling
  argList$min_ds_size=min_dataset_size
  
  if(is.null(inputEmbeddings)){
    tryCatch({.sconline.pca(inputExpData = inputExpData,argList = argList,batch_variable = batch_variable,organism=organism);F}, error=function(e) {return(T)})
  }
  
  
  res_embeddings=.sconline.embeddingFn(argList,inputEmbeddings=inputEmbeddings,run_harmony=run_harmony,pd=inputPhenoData,inputBatchCol=batch_variable,saveFiles = F)
  
  #inputExpData=list(inputExpData);index=1
  #argList=argList;inputEmbeddings=res_embeddings$harmony_embeddings;inputPhenoData=inputPhenoData;organism=organism;inputExpData=list(inputExpData);L2Norm=L2Norm;mergePseudocells=mergePseudocells;merging_strength=merging_strength;sig1_thr=sig1_thr;pct2_thr=pct2_thr;pct_diff_thr=pct_diff_thr;run_harmony=run_harmony;index=0;min_cluster_size=500;input_res_embeddings = res_embeddings;input_res_prop=NULL
  run_res=.extra_hierarchical_run(argList=argList,batch_variable = batch_variable,inputEmbeddings=res_embeddings$harmony_embeddings,inputPhenoData=inputPhenoData,organism=organism,inputExpData=list(inputExpData),L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,run_harmony=run_harmony,index=1,min_cluster_size = min_cluster_size,input_res_embeddings = res_embeddings,input_res_prop = NULL)
  run_res=list(res=run_res,index=0)
  
  res=.extra_hierarchical_clustList(x=run_res)
  
  prop_mat=res$index
  names(prop_mat)=res$cell
  prop_mat=t(as.matrix(.myOneHotFn(prop_mat)))
  tmp_exp=inputExpData
  tmp_exp=tmp_exp[,match(colnames(prop_mat),colnames(tmp_exp))]
  tmp_exp=list(tmp_exp)
  tmp_pheno=inputPhenoData[match(colnames(prop_mat),row.names(inputPhenoData)),]
  tmp_embedding=res_embeddings$harmony_embeddings
  tmp_embedding=tmp_embedding[match(colnames(prop_mat),row.names(tmp_embedding)),]
  #argList = argList;input_prop_mat = prop_mat;expData = tmp_exp;hierarchical_mode = T;addClusteringModule=F;extendedMode = F;L2Norm=L2Norm;mergePseudocells=mergePseudocells;merging_strength=merging_strength;sig1_thr=sig1_thr;pct2_thr=pct2_thr;pct_diff_thr=pct_diff_thr;input_embedding=tmp_embedding;input_pca_centroid=NULL;input_pca_centroids_assignments=NULL;input_pd=tmp_pheno
  res_prop=.myConcensusDEFn_step2_fast(argList = argList,input_prop_mat = prop_mat,expData = tmp_exp,hierarchical_mode = T,addClusteringModule=F,extendedMode = F,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,input_embedding=tmp_embedding,input_pca_centroid=NULL,input_pca_centroids_assignments=NULL,input_pd=tmp_pheno)
  gc()
  #input_res=res_prop;inputPhenoData = inputPhenoData;input_pca_centroids=res_embeddings$pca_centroid;input_embeddings=res_embeddings$harmony_embeddings
  
  #.tst=list(argList=argList,input_res=res_prop,inputPhenoData = inputPhenoData,input_pca_centroids=res_embeddings$pca_centroid,input_embeddings=res_embeddings$harmony_embeddings,prune_clusters=F)
  #argList=.tst$argList;input_res=.tst$input_res;inputPhenoData = .tst$inputPhenoData;input_pca_centroids=.tst$input_pca_centroids;input_embeddings=.tst$input_embeddings;prune_clusters=T
  
  #argList=argList;meta_data=res_prop$res_meta;pd=inputPhenoData;input_embeddings=res_embeddings$harmony_embeddings;prop_mat=res_prop$prop_mat;input_pca_centroids=res_embeddings$pca_centroid;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;cell_count_diff_thr=0
  clust_res=.extra_hierarchical_clustMerge(argList=argList,meta_data=res_prop$res_meta,pd=inputPhenoData,input_embeddings=res_embeddings$harmony_embeddings,prop_mat=res_prop$prop_mat,input_pca_centroids=res_embeddings$pca_centroid,
                                           sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,cell_count_diff_thr=0)
  
  
  
}


