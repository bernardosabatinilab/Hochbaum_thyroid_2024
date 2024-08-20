source("~/code/sconline_code.R")
library(fgsea)
library(qs)

# Function that takes as input any DE result, any gene set, and any ranking variable and outputs fGSEA results
myfGSEAfn=function(input_DE,input_geneSet_filename,input_rankBy,input_scoreType,pval_thresh=0.05,logFC_thresh=c(log2(1.25),log2(1)),robScore_thresh=0.5) {
  DE_anno=.extraMouseGeneAnnoAdderFn(input_DE$gene)
  test=DE_anno[match(input_DE$gene,rownames(DE_anno)),]
  tmp=cbind(input_DE,test)
  DE_prot=tmp[which(tmp$gene_biotype=="protein_coding"),drop=FALSE,]
  
  # Order genes by abs(t-statistic)
  orderedGenes=DE_prot
  orderedGenes=orderedGenes[order(orderedGenes[,input_rankBy],decreasing=TRUE),] 
  rankedVec=orderedGenes[,input_rankBy]
  names(rankedVec)=orderedGenes$gene
  
  # construct gene set 
  minSize=15
  maxSize=500
  gs=.sconline.GSEA.readGMT(file=input_geneSet_filename,bkg_genes=DE_prot$gene,min.gs.size=minSize,max.gs.size=maxSize)
  print(head(gs))
  # run GSEA
  fGSEA_run <- fgsea(pathways = gs, 
                     stats = rankedVec, 
                     scoreType= input_scoreType, 
                     minSize=minSize, 
                     maxSize=maxSize)
  print(head(fGSEA_run))
  fGSEA_run=fGSEA_run[order(fGSEA_run$pval,decreasing = F),]
  fGSEA_run=as.data.frame(fGSEA_run)
  fGSEA_run$leadingEdge=unlist(lapply(fGSEA_run$leadingEdge,function(x) paste(x,collapse = ",")))
  # fGSEA_run=fGSEA_run[fGSEA_run$padj<0.1,]
  # create column(s) for "significant" leadingEdge genes according to specified threshold
  test=apply(fGSEA_run,1,function(x) { # go row by row to filter the leadingEdge genes against those that are significant in input_DE
    genes=as.character(x["leadingEdge"])
    genes=unlist(lapply(genes,function(x) (strsplit(x,","))))
    cols_filteredLeadingEdge=lapply(logFC_thresh,function(i){ # can handle multiple logFC thresholds
      tmp=input_DE
      tmp=tmp[match(genes,tmp$gene),]
      tmp=tmp[(tmp$adj.P.Val<pval_thresh & abs(tmp$cp100k_logFC)>i & abs(tmp$robScore_logFC)>=robScore_thresh),]
      return(tmp$gene)
    })
    return(cols_filteredLeadingEdge)
  })
  
  # after exploring filtering leading edge by logFC thresholds, we ultimately removed this criteria (ie set logFC_thresh=log2(1))
  for (i in seq_along(logFC_thresh)){
    filtColName=paste0("filtLeadingEdge_logFCthresh_",round(logFC_thresh[[i]],2))
    fGSEA_run$newCol=unlist(lapply(test,function(x) paste(x[[i]],collapse=",")))
    colnames(fGSEA_run)[length(fGSEA_run)]=filtColName
    fGSEA_run$newCol=unlist(lapply(test,function(x) length(x[[i]])))
    colnames(fGSEA_run)[length(fGSEA_run)]=paste0("num_",filtColName)
  }
  results=list(fGSEA_run,rankedVec,gs)
  return(results) 
}

