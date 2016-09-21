# Co-expression, collection of scripts

#########################################
# install package LINC from Bioconductor:

## In R-3.3
library(BiocInstaller)
useDevel()
biocValid()              # checks for out of date packages

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("LINC")

# install package LINC from GitHub:
require(devtools)
install_github("ManuelGoepferich/LINC")

########################################
# translate gene ids into entrez genes

# input x has to be a matrix with ENSG gene ids
makeEntrezMatrix <- function(x){
require(mygene)
query <- mygene::getGenes(rownames(x), fields = "entrezgene")
m <- match(rownames(x), query@listData$query)
entrez_query <- query@listData$entrezgene[m]
rownames(x) <- entrez_query
return(x) # return the matrix with entrez gene ids
}

#######################################
# get the gene biotype i.e. protein-coding, lincRNA

# input x has to be a matrix with entrez gene ids
biotypeMatrix <- function(x){
  require(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
l_biotype <- getBM(attributes=c('gene_biotype',
                                'entrezgene'),
                   filters = 'entrezgene',
                   values = rownames(x),
                   mart = ensembl)  
m <- match( rownames(x),l_biotype$entrezgene )
check_id <- identical(as.character(l_biotype$entrezgene[m]), rownames(x))
if(!check_id) warning("Match was not identical, please check this by hand!")
g_biotype <- l_biotype$gene_biotype[m]
return(g_biotype)
}

#############################################
# get the high variance genes

# input x has to be a matrix with genes represented by rows
# input nbest has to be the number of required high variance genes
getHighVar <- function(x, nbest){
varindex <- order(apply(x, 1, var), decreasing = T)
out <- x[varindex[1:nbest], ]
return(out)
}

###############################################
# intersect genes of two matrices

# x and y are matrices, respectively 

function(x, y){
  out <- intersect(rownames(x), rownames(y))
  return(out)
}


##########################################
# example for using package LINC
require(LINC)

LNC_index <- (is.element(GTEX_BREAST_BIOTYPE, c("protein_coding", "lincRNA")))
gtex.breast.linc <- linc(GTEX_BREAST_201SP_ENTREZ[LNC_index, ], codingGenes = GTEX_BREAST_BIOTYPE[LNC_index]   )
gtex.breast.clust <- clusterlinc(gtex.breast.linc, pvalCutOff = 0.00000000000000005, distMethod = "dicedist" )


#########################################
# intersect more than two vectors

# example:
x <- c(1:8)
y <- c(5:10)
z <- c(5, 8, 20)
Reduce(intersect, list(x, y, z) )
# 5 and 8 are elements of vectors x, y and z

########################################
# working on a list of 'LINCcluster' objects

enrichQuery <- function(qlist, query,  size){
  require(ReactomePA)
  out_list <- list()
  for(h in 1: length(qlist)){
    
    index <- is.element(names(qlist[[h]]@results$cluster$neighbours), query)
    if(any(index)){
      set <- qlist[[h]]@results$cluster$neighbours[index]
      set <- unlist(set)[1:size]
      set <- set[!is.na(set)]
    } else {
      set <- NULL
    }
    
    if(!is.null(set)){
      path_promise <- enrichPathway(set) 
      if(!is.null(path_promise) && !is.logical(path_promise)){
        out_list[[h]] <- path_promise@result$Description[1:2] 
      } else {
        out_list[[h]] <- NA
      }
    } else {
      out_list[[h]] <- NA
    }
    
  }
  names(out_list) <- unlist(lapply(qlist, function(x){
    organ <- x@history$customID
    col <- x@history$customCol
    if(col == "blue"){
      type <- "HEALTHY"
    }
    if(col == "red"){
      type <- "CANCER"
    }
    return(paste(organ, type, sep = "_"))
  }))
  
  return(out_list)
}

#path_647979 <- enrichQuery(cluster_list4, query = '647979', size = 100)
path_60674 <- enrichQuery(cluster_list5, query = '60674', size = 100)



########################################
# normalize gene expression
require(preprocessCore)
COLON_normal <- normalize.quantiles( cbind(cluster_list5[[3]]@expression,cluster_list5[[4]]@expression)     )
rownames(COLON_normal) <- rownames(cluster_list5[[3]]@expression)

normal_expression_list5[[3]] <- COLON_normal[,1:ncol(cluster_list5[[3]]@expression)] 
normal_expression_list5[[4]] <- COLON_normal[,(ncol(cluster_list5[[4]]@expression) + 1):ncol(COLON_normal)] 


#######################################
# make the final plots for co-expression patterns

# input cluster_list is a list of 'LINCcluster' objects
# input query is the gene name
# input tissue is a vector of names of objects in cluster_list

dotplot27groups <- function(cluster_list, query, tissue){
  require(ReactomePA)
  comp <- enrich27(cluster_list, query, 100)
  names(comp) <- tissue
  #comp <- list(GAS5 = healthy_gas5,NEAT1 = healthy_neat1, BOTH =  intersection)
  res <- compareCluster(comp, fun="enrichPathway")
  plot(res, showCategory = 4) + theme(axis.text.x = element_text(angle = 90, hjust = 0, 
                                                                 colour = c(rep("red", 13), rep("blue", 14))  ) )
}

dotplot27groups(cluster_list5, "283131", tissue)





