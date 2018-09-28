require(qlcMatrix)
require(Seurat)
RunTruncatedPCA <- function(object, n.genes.pc=40, genescale.method="div.by.max", pc.genes=NULL, ...){
  
  if(is.null(pc.genes)) pc.genes <- rownames(object@data)
  object <- RunPCA(object,pc.genes = pc.genes, ...)
  old.loadings <- GetDimReduction(object, reduction.type = "pca", slot="gene.loadings")
  new.pc.scores <- reCalculatePCScores(object, pc.loadings = old.loadings, n.genes.pc = n.genes.pc, use.binary.weights = F, genescale.method = genescale.method)
  object <- SetDimReduction(object, reduction.type = "pca",slot = "cell.embeddings",new.data = new.pc.scores)
  
  object
}

reCalculatePCScores <- function(object, pc.loadings, n.genes.pc=30, use.binary.weights=F, genescale.method=c("div.by.max","mean.center","scale","seurat.scale")){
  
  data <- object@data[rownames(pc.loadings), ] # retain only genes present in the loadings
  # gene-wise scaling by one of several methods
  if(genescale.method=="div.by.max") data <- data / qlcMatrix::rowMax(data)
  else if(genescale.method=="mean.center") data <- t(scale(t(data), center = T, scale = F))
  else if(genescale.method=="scale") data <- t(scale(t(data)))
  else data <- object@scale.data[rownames(pc.loadings), ] # use Seurat scaled data
  
  # if n.genes.pc == 0, then no cutoff is applied and loading is simply recalculated as input loadings * scaled data
  if(n.genes.pc!=0){
    # apply loading cutoffs so that only genes in top n positive or top n negative loadings for each PC get non-zero loadings
    pc.loadings.reweighted <- apply(pc.loadings, 2, cutoff.by.order, n.cutoff=n.genes.pc)
    
    if(use.binary.weights){
      pc.loadings.reweighted[pc.loadings.reweighted > 0] <- 1
    }
    
    pc.scores <- as.matrix( t(data) %*% pc.loadings.reweighted )
  }
  else{
    pc.scores <- as.matrix( t(data) %*% pc.loadings )
  }
  
  
  # recalculate PC scores (cell-wise scores) based on new gene weights
  return(pc.scores)
}
cutoff.by.order <- function(x, n.cutoff=30){
  x.rank <- rank(x)
  len.x <- length(x)
  x[which( (x.rank > n.cutoff) & (x.rank < (len.x-n.cutoff + 1)) )] <- 0
  
  return(x)
}

# useful plotting functions for eliminating pheno-orthogonal PCs
PlotTopLoadingGenes <- function(object, dim.plot=1, ngenes=4, type=c("both", "pos","neg"),...){

  loadings <- sort(GetGeneLoadings(object=object, reduction.type = "pca",dims.use = dim.plot[1])[,1])
  genes.neg <- names(loadings)[1:ngenes]
  loadings <- sort(loadings, decreasing = T)
  genes.pos <- names(loadings)[1:ngenes]

  FeaturePlot(object=object, features.plot = genes.pos, cols.use = c("grey80", "darkblue"), coord.fixed = T, no.axes = T, ...)#+
    # ggtitle(paste("PC", dim.plot, "pos genes"))
  FeaturePlot(object=object, features.plot = genes.neg, cols.use = c("grey80", "darkblue"), coord.fixed = T, no.axes = T,...)#+
    # ggtitle(paste("PC", dim.plot, "neg genes"))

  # if(type=="both") return(plot_grid(plotlist=c(p1, p2), ncol=2))
  # else if(type=="pos") return(p1)
  # else if(type=="neg") return(p2)

}
