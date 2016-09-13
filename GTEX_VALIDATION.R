# GTEX validation

# variance and expression at all
# empirical p-value
# validation with protein-coding genes
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/GTEX/")
tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2, nrows = 5)

GTEX_test <- (as.matrix(tb_promise))
#genes <- (GTEX_test[,1][2:56319])
samples <- as.vector(GTEX_test[1,])
genes <- 

samples <- as.vector(tb_promise[1,])
m <- match(GTEX_BREAST_IDS, samples)


tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2,
                         colClasses = c(rep("character", 4), rep("NULL", 2680), rep("character", 212), rep("NULL", 5661)))

GTEX_test <- tb_promise


m <- match(gtex_ids_LIVER[,1][12:67], samples)
m[!is.na(m)]


# breast
GTEX_test_mat <- as.matrix(GTEX_test[2:56319, 3:216])
GTEX_test_mat <- matrix(as.numeric(GTEX_test_mat), ncol = 214, nrow = 56318)
colnames(GTEX_test_mat) <- samples[3:216]
genes <- as.vector(GTEX_test[,1])[2:56319]
clear_genes <- vapply(genes, function(x){strsplit(x, '.', fixed = T)[[1]][[1]] }, "ENSG" )
rownames(GTEX_test_mat) <- clear_genes
GTEX_BREAST_CRUDE <- GTEX_test_mat

varindex <- order(apply(GTEX_BREAST_CRUDE, 1, var), decreasing = T)
GTEX_BREAST_CRUDE_15000 <- GTEX_BREAST_CRUDE[varindex[1:9100], ]

minindex <- apply(GTEX_BREAST_CRUDE_15000, 1, function(x){ (max(x) > 10) }  )
GTEX_BREAST_CRUDE_EXPRESSED <- GTEX_BREAST_CRUDE_15000[minindex,]

#GTEX_BREAST_CRUDE_EXPRESSED <- GTEX_BREAST_CRUDE_EXPRESSED[!duplicated(rownames(GTEX_BREAST_CRUDE_EXPRESSED)),]

# make entrez genes
query <- mygene::getGenes(rownames(GTEX_BREAST_CRUDE_EXPRESSED), fields = "entrezgene")
m <- match(rownames(GTEX_BREAST_CRUDE_EXPRESSED), query@listData$query)
entrez_query <- query@listData$entrezgene[m]
rownames(GTEX_BREAST_CRUDE_EXPRESSED) <- entrez_query
GTEX_BREAST_CRUDE_ENTREZ <- GTEX_BREAST_CRUDE_EXPRESSED[(!is.na(rownames(GTEX_BREAST_CRUDE_EXPRESSED))), ]

m <- match(GTEX_BREAST_IDS, colnames(GTEX_BREAST_CRUDE_ENTREZ))
m[!is.na(m)]

GTEX_BREAST_201SP_ENTREZ <- GTEX_BREAST_CRUDE_ENTREZ[, (m[!is.na(m)]) ]
#str(intersect(colnames(GTEX_BREAST_208SP_ENTREZ), GTEX_BREAST_IDS))


ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
l_biotype <- getBM(attributes=c('gene_biotype',
                                'entrezgene'),
                   filters = 'entrezgene',
                   values = rownames(GTEX_BREAST_201SP_ENTREZ),
                   mart = ensembl)  
m <- match( rownames(GTEX_BREAST_201SP_ENTREZ),l_biotype$entrezgene )
#identical(as.character(l_biotype$entrezgene[m]), rownames(GTEX_LIVER_CRUDE_ENTREZ))
GTEX_BREAST_BIOTYPE <- l_biotype$gene_biotype[m]


#str(as.matrix(GTEX_test[2:99, 3:8557]))

# numeric conversion
#GTEX_test_mat <- (data.matrix(tb_promise[2:56319, 3:140]))
#GTEX_test_mat <- (data.matrix(GTEX_test[2:56319, 3:140]))
#GTEX_test_mat <- matrix(as.numeric(GTEX_test[2:56319, 3:140]), ncol = 138, nrow = 56318)
GTEX_test_mat <- as.matrix(GTEX_test[2:56319, 3:140])
GTEX_test_mat <- matrix(as.numeric(GTEX_test_mat), ncol = 138, nrow = 56318)

colnames(GTEX_test_mat) <- samples
# clear genes
clear_genes <- vapply(genes, function(x){strsplit(x, '.', fixed = T)[[1]][[1]] }, "ENSG" )
rownames(GTEX_test_mat) <- clear_genes

GTEX_LIVER_CRUDE <- GTEX_test_mat
# select the genes

varindex <- order(apply(GTEX_LIVER_CRUDE, 1, var), decreasing = T)
GTEX_LIVER_CRUDE_15000 <- GTEX_LIVER_CRUDE[varindex[1:15000], ]

minindex <- apply(GTEX_LIVER_CRUDE_15000, 1, function(x){ (max(x) > 10) }  )
GTEX_LIVER_CRUDE_EXPRESSED <- GTEX_LIVER_CRUDE_15000[minindex,]

# make entrez genes
query <- mygene::getGenes(rownames(GTEX_LIVER_CRUDE_EXPRESSED), fields = "entrezgene")
m <- match(rownames(GTEX_LIVER_CRUDE_EXPRESSED), query@listData$query)
entrez_query <- query@listData$entrezgene[m]
rownames(GTEX_LIVER_CRUDE_EXPRESSED) <- entrez_query
GTEX_LIVER_CRUDE_ENTREZ <- GTEX_LIVER_CRUDE_EXPRESSED[(!is.na(rownames(GTEX_LIVER_CRUDE_EXPRESSED))), ]

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
l_biotype <- getBM(attributes=c('gene_biotype',
                                  'entrezgene'),
                     filters = 'entrezgene',
                     values = rownames(GTEX_LIVER_111SP_ENTREZ),
                     mart = ensembl)  

m <- match( rownames(GTEX_LIVER_111SP_ENTREZ),l_biotype$entrezgene )
#identical(as.character(l_biotype$entrezgene[m]), rownames(GTEX_LIVER_CRUDE_ENTREZ))
GTEX_LIVER_BIOTYPE <- l_biotype$gene_biotype[m]

str(gtex_ids_LIVER[,1])
str(colnames(GTEX_LIVER_CRUDE_ENTREZ))
m <- intersect(colnames(GTEX_LIVER_CRUDE_ENTREZ), gtex_ids_LIVER[,1])

GTEX_LIVER_111SP_ENTREZ <- GTEX_LIVER_CRUDE_ENTREZ[,m]

GTEX_LIVER_111SP_ENTREZ <- GTEX_LIVER_111SP_ENTREZ[!duplicated(rownames(GTEX_LIVER_111SP_ENTREZ)), ]
GTEX_LIVER_BIOTYPE <- GTEX_LIVER_BIOTYPE[!duplicated(rownames(GTEX_LIVER_111SP_ENTREZ))]

biotype_100 <- c(rep(F, 412), rep(T, 9000) )
biotype_20 <- c(rep(F, 412), rep(T, 9000) )
#GTEX_LIVER_111SP_ENTREZ[, (GTEX_LIVER_BIOTYPE == "protein_coding")]
gtex.liver.linc <- linc(GTEX_LIVER_111SP_ENTREZ[(GTEX_LIVER_BIOTYPE == "protein_coding"), ], codingGenes = !is.element(protein_genes, potential_100))
gtex.liver.clust <- clusterlinc(gtex.liver.linc, pvalCutOff = 0.00000000000000005, distMethod = "dicedist" )
gtex.liver.bio <- getbio2(gtex.liver.clust, translate = 'none')                        # 5e-17

# select the best genes

potential_targets <- rownames(GTEX_LIVER_111SP_ENTREZ[(GTEX_LIVER_BIOTYPE == "protein_coding"), ])[1:412]
protein_genes <- rownames(GTEX_LIVER_111SP_ENTREZ[(GTEX_LIVER_BIOTYPE == "protein_coding"), ]

gotrue <- vector()
for(h in 1:412){
  go_promise <- try(getGOInfo(potential_targets[h]), silent = T)
  if(class(go_promise) == "try-error"){
    gotrue[h] <- FALSE
  } else {
    gotrue[h] <- TRUE
  }
}

potential_go <- potential_targets[gotrue]
potential_100 <- potential_go[1:100]

cortested_100 <- gtex.liver.clust@correlation$cortest

# get the best 10 go terms
target_list_100 <- list()
for( f in 1:100){
   entrez20 <- sort(cortested_100[,f], decreasing = F)[1:20]
   
   en_promise <- enrichGO(names(entrez20), ont = "BP")
   if(class(en_promise) == "enrichResult"){
     target_list_100[[f]] <- en_promise@result$ID[1:10]
   } else {
     target_list_100[[f]] <- NA
   }
}
names(target_list_100) <- colnames(cortested_100)
 


str(target_list_100[c(1:42)])

best_list_100 <- vector()
for( f in 1:100){
  best_list_100[f] <- doGOsim(names(target_list_100)[f], target_list_100[[f]] )
}
 

for( f in c(50,57,61,62,68,81,88,98)){
  best_list_100[f] <- doGOsim(names(target_list_100)[f], target_list_100[[f]] )
}


str(best_list_100)

responsible1 <- best_list_100[c(50,57,61,62,68,81,88,98)]
names(responsible1) <- names(target_list_100)[c(50,57,61,62,68,81,88,98)]

# get the empirical p-values
# entrez genes assumed


responsible_PVALUE <- vector()
for(u in 1:length(responsible1)){

genes_in_set <- rownames(GTEX_LIVER_111SP_ENTREZ)
co_size <- 20
score_storage <- vector()
ptarget <- responsible1[u]
best <- names(responsible1)[u]

for(k in 1:1000){
    mok <- sample(genes_in_set, co_size)
    en_promise <- enrichGO(mok, ont = "BP")@result$ID
    if(is.character(en_promise) && length(en_promise) > 0){
      score_storage[k] <- doGOsim(ptarget, en_promise)
    } else {
      score_storage[k] <- 0 
    }
}
ppos <- length(which(best < score_storage))
if(ppos == 0) { pvalue <- 1e-16
} else {
pvalue <- (ppos/ length(score_storage))
}
  
responsible_PVALUE[u] <- pvalue
}


save(responsible_PVALUE, file = "responsible_PVALUE1.RData")
  

#query <- mygene::getGenes(clear_genes, fields = "entrezgene")
#test_genes <- unique(query@listData$entrezgene)

str(gtex.test.bio@results[[1]])


id_promise <- read.table(file ="GTEx_Data_V6_Annotations_SampleAttributesDS.txt", skip = 1, fill = TRUE) # , skip =0, nrows = 10)

apply(id_promise, 1, function(x) { if(length(grep("^GTEX-", x[1])) = 1) { T
          }else{
          F
  }}  )

gtex_mat <- as.matrix(id_promise)
gtex_index <- vapply(gtex_mat[,1], function(x) { if(length(grep("^GTEX-.*", x)) == 1){
  T
}else{
  F }
} , TRUE  )
gtex_mat <- gtex_mat[gtex_index, ]

pancreasgr <- apply(gtex_mat, 1, function(y){
      allcol <- vapply(y, function(x) { if(length(grep("Pancreas", x)) == 1){
        T
      }else{
        F }
      } , TRUE  )
      
      if(any(allcol)){
        T
      } else {
        F
      }
   })

breastgr <- apply(gtex_mat, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Mammary", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(any(allcol)){
    T
  } else {
    F
  }
})


gtex_ids_PANCREAS <- gtex_mat[pancreasgr,] 

str(gtex_ids_LIVER[,7])
table(gtex_ids_LIVER[,7])

gtex_ids_BREAST <- gtex_mat[breastgr,] 


pancreaspos <- apply(gtex_ids_PANCREAS, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Pancreas", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(length(which(allcol)) > 1){
    T
  } else {
    F
  }
})

breastpos <- apply(gtex_ids_BREAST, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Breast", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(length(which(allcol)) > 1){
    T
  } else {
    F
  }
})


gtex_ids_PANCREAS <- gtex_ids_PANCREAS[pancreaspos, ]
GTEX_PANCREAS_IDS <- gtex_ids_PANCREAS[,1]

gtex_ids_BREAST <- gtex_ids_BREAST[breastpos, ]
GTEX_BREAST_IDS <- gtex_ids_BREAST[,1]


GTEX_50LIVER_IDS <- gtex_ids_LIVER[1:50,1]
save(gtex_ids_LIVER,
     GTEX_50LIVER_IDS,
     GTEX_LIVER_CRUDE,
     GTEX_LIVER_BIOTYPE,
     GTEX_LIVER_111SP_ENTREZ,
     responsible1,
     file = "GTEX_LIVER.RData")

save(gtex_ids_BREAST,
     GTEX_BREAST_IDS,
     GTEX_BREAST_CRUDE,
     GTEX_BREAST_201SP_ENTREZ,
     GTEX_BREAST_BIOTYPE,
     file = "GTEX_BREAST.RData")

save(gtex_ids_PANCREAS,
     GTEX_PANCREAS_IDS,
     GTEX_PANCREAS_CRUDE,
     GTEX_PANCREAS_156SP_ENTREZ,
     GTEX_PANCREAS_BIOTYPE,
     file = "GTEX_PANCREAS.RData")


# get all lnc RNAs in the data
#

LNC_index <- (is.element(GTEX_LIVER_BIOTYPE, c("protein_coding", "lincRNA")))
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/LINC/data/")
gtex.liver.linc <- linc(GTEX_LIVER_111SP_ENTREZ[LNC_index, ], codingGenes = GTEX_LIVER_BIOTYPE[LNC_index]   )
gtex.liver.clust <- clusterlinc(gtex.liver.linc, pvalCutOff = 0.00000000000000005, distMethod = "dicedist" )
gtex.liver.bio <- getbio(gtex.liver.clust, translate = 'none')                        # 5e-17
  
  n <- 66
  tmp <- enrichGO(gtex.liver.clust@results[[1]][[5]][[n]], ont = "BP" )
  names(gtex.liver.clust@results[[1]][[5]][n])
  #tmp@result$ID
  tmp@result$Description
  getGOInfo(names(gtex.liver.clust@results[[1]][[5]][n]))
  
  doGOsim("26238", tmp@result$ID[1:10])
  
  
  gtex.liver.linc.pc3 <- linc(GTEX_LIVER_111SP_ENTREZ[LNC_index, ], codingGenes = GTEX_LIVER_BIOTYPE[LNC_index], rmPC = c(4:111)   )
  gtex.liver.clust.pc3 <- clusterlinc(gtex.liver.linc.pc3, pvalCutOff = 0.00000000000000005, distMethod = "dicedist" )  # first 3 over 60 %
  
  
  n <- 65
  tmp <- enrichGO(gtex.liver.clust.pc3@results[[1]][[5]][[n]], ont = "BP" )
  names(gtex.liver.clust.pc3@results[[1]][[5]][n])
  #tmp@result$ID
  tmp@result$Description
 
  
  
   getGOInfo(names(gtex.liver.clust@results.pc3[[1]][[5]][n]))
  
  
  
target <- getGOInfo("84129")  #221178



pca_analysis <- prcomp(GTEX_LIVER_111SP_ENTREZ[LNC_index, ], center = FALSE,
                     scale. = FALSE) 



require(GOSim)
doGOsim <- function(target, predicted){
  
  go_promise <- try(getGOInfo(target), silent = T)
 if(class(go_promise) != "try-error" && length(go_promise) != 0){
  
  t_promise <- grepl("^GO:[0-9]+", go_promise, perl = TRUE)
  go_ref <- unlist(go_promise[t_promise])
  
  if(length(go_ref) > 0){
    
  out <- c()
  for(n in 1:length(go_ref)){
    sim_mat <- getTermSim(c(go_ref[n], predicted))
    out[n] <- max(sim_mat[1,2:(length(predicted) +1 )])
  }
  return(mean(out))
  } else {
  return(NA)  
  }
 } else {
 return(NA)
 }
}


# in work
for(n in 1:length(objective)) {
  
  if(is.na(names(objective)[n])){
    predicted_asvector[n] <- NA
  } else {
    if(is.character(objective[n]) &&
       length(objective[n]) > 0){
      predicted_asvector[n] <- doGOsim(nobjective[n], objective)
    } else {
      predicted_asvector[n] <- NA
    } 
    
  }
  
}


# for list
for(n in 1:length(objective)) {
  
  if(is.na(names(objective)[n])){
    predicted_asvector[n] <- NA
  } else {
    if(is.character(objective[[n]][[1]]) &&
       length(objective[[n]][[1]]) > 0){
      predicted_asvector[n] <- doGOsim(names(objective)[n], objective[[n]][[1]])
    } else {
      predicted_asvector[n] <- NA
    } 
    
  }
  
}


# BREAST TISSUE VALIDATION

LNC_index <- (is.element(GTEX_BREAST_BIOTYPE, c("protein_coding", "lincRNA")))

gtex.breast.linc <- linc(GTEX_BREAST_201SP_ENTREZ[LNC_index, ], codingGenes = GTEX_BREAST_BIOTYPE[LNC_index]   )
gtex.breast.clust <- clusterlinc(gtex.breast.linc, pvalCutOff = 0.00000000000000005, distMethod = "dicedist" )
gtex.breast.bio <- getbio(gtex.breast.clust, translate = 'none')         # 5e-17

j <- 44
getGOInfo("558")

doGOsim("221178", res@result$ID[1:10])
  


pca_analysis <- prcomp(GTEX_BREAST_201SP_ENTREZ[1:20, 1:10], center = FALSE, scale. = FALSE) 

pca_analysis <- prcomp(matrix(rnorm(50), ncol = 5), center = FALSE, scale. = FALSE) 



gtex.breast.linc.pc <- linc(GTEX_BREAST_201SP_ENTREZ[LNC_index, ], codingGenes = GTEX_BREAST_BIOTYPE[LNC_index], rmPC = c(4:201)   )
gtex.breast.clust.pc <- clusterlinc(gtex.breast.linc.pc, pvalCutOff = 0.00000000000000005, distMethod = "dicedist" )
#gtex.breast.bio <- getbio(gtex.breast.clust, translate = 'none')         # 5e-17


n <- 43
set <- gtex.breast.clust.pc@results[[1]]$neighbours[[n]][1:20]
res <- enrichGO(set, ont = "BP")
res@result$Description
names(gtex.breast.clust.pc@results[[1]]$neighbours)[n]


set <- gtex.breast.clust.pc@results[[1]]$neighbours[["283673"]][1:1318]
res <- enrichGO(set, ont = "BP")

x <- GTEX_BREAST_201SP_ENTREZ[LNC_index, ]["3190",]
y <- GTEX_BREAST_201SP_ENTREZ[LNC_index, ]["283673",]
  
# Pancreas
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/GTEX/")
tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2, nrows = 5)

GTEX_test <- (as.matrix(tb_promise))
#genes <- (GTEX_test[,1][2:56319])
samples <- as.vector(GTEX_test[1,])
genes <- 
  
  samples <- as.vector(tb_promise[1,])
m <- match(GTEX_BREAST_IDS, samples)


tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2,
                         colClasses = c(rep("character", 4), rep("NULL", 6115), rep("character", 184), rep("NULL", 2254)))

GTEX_test <- tb_promise


m <- match(gtex_ids_PANCREAS[,1], samples)

match(samples, colnames(GTEX_PANCREAS_CRUDE_ENTREZ) )

m[!is.na(m)]


# Pancreas
GTEX_test_mat <- as.matrix(GTEX_test[2:56319, 3:188])
GTEX_test_mat <- matrix(as.numeric(GTEX_test_mat), ncol = 186, nrow = 56318)
colnames(GTEX_test_mat) <- samples[3:188]
genes <- as.vector(GTEX_test[,1])[2:56319]
clear_genes <- vapply(genes, function(x){strsplit(x, '.', fixed = T)[[1]][[1]] }, "ENSG" )
rownames(GTEX_test_mat) <- clear_genes
GTEX_PANCREAS_CRUDE <- GTEX_test_mat

varindex <- order(apply(GTEX_PANCREAS_CRUDE, 1, var), decreasing = T)
GTEX_PANCREAS_CRUDE_15000 <- GTEX_PANCREAS_CRUDE[varindex[1:9100], ]

minindex <- apply(GTEX_PANCREAS_CRUDE_15000, 1, function(x){ (max(x) > 10) }  )
GTEX_PANCREAS_CRUDE_EXPRESSED <- GTEX_PANCREAS_CRUDE_15000[minindex,]

# make entrez genes
query <- mygene::getGenes(rownames(GTEX_PANCREAS_CRUDE_EXPRESSED), fields = "entrezgene")
m <- match(rownames(GTEX_PANCREAS_CRUDE_EXPRESSED), query@listData$query)
entrez_query <- query@listData$entrezgene[m]
rownames(GTEX_PANCREAS_CRUDE_EXPRESSED) <- entrez_query
GTEX_PANCREAS_CRUDE_ENTREZ <- GTEX_PANCREAS_CRUDE_EXPRESSED[(!is.na(rownames(GTEX_PANCREAS_CRUDE_EXPRESSED))), ]

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
l_biotype <- getBM(attributes=c('gene_biotype',
                                'entrezgene'),
                   filters = 'entrezgene',
                   values = rownames(GTEX_PANCREAS_156SP_ENTREZ),
                   mart = ensembl)  

m <- match( rownames(GTEX_PANCREAS_156SP_ENTREZ),l_biotype$entrezgene )
#identical(as.character(l_biotype$entrezgene[m]), rownames(GTEX_LIVER_CRUDE_ENTREZ))
GTEX_PANCREAS_BIOTYPE <- l_biotype$gene_biotype[m]

#str(gtex_ids_LIVER[,1])
str(colnames(GTEX_PANCREAS_CRUDE_ENTREZ))
colnames(GTEX_PANCREAS_CRUDE_ENTREZ) <- samples[3:188]

m <- intersect(colnames(GTEX_PANCREAS_CRUDE_ENTREZ), gtex_ids_PANCREAS[,1])

GTEX_PANCREAS_156SP_ENTREZ <- GTEX_PANCREAS_CRUDE_ENTREZ[,m]

GTEX_PANCREAS_156SP_ENTREZ <- GTEX_PANCREAS_156SP_ENTREZ[!duplicated(rownames(GTEX_PANCREAS_156SP_ENTREZ)), ]
GTEX_PANCREAS_BIOTYPE <- GTEX_PANCREAS_BIOTYPE[!duplicated(rownames(GTEX_PANCREAS_156SP_ENTREZ))]

#now, validate



# Skin
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/GTEX/")
id_promise <- read.table(file ="GTEx_Data_V6_Annotations_SampleAttributesDS.txt", skip = 1, fill = TRUE) # , skip =0, nrows = 10)

gtex_mat <- as.matrix(id_promise)
gtex_index <- vapply(gtex_mat[,1], function(x) { if(length(grep("^GTEX-.*", x)) == 1){
  T
}else{
  F }
} , TRUE  )
gtex_mat <- gtex_mat[gtex_index, ]

testisgr <- apply(gtex_mat, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Skin", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(any(allcol)){
    T
  } else {
    F
  }
})



skingr <- apply(gtex_mat, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Skin", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  
  allcol2 <- vapply(y, function(x) { if(length(grep('Lower', x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  allcol3 <- vapply(y, function(x) { if(length(grep('Exposed', x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  #if(any(allcol) && any(allcol2)){
  if(any(allcol) && any(allcol2)  && any(allcol3) ){
    T
  } else {
    F
  }
})

gtex_ids_SKIN <- gtex_mat[skingr,] 

testispos <- apply(gtex_ids_TESTIS, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Testis", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(length(which(allcol)) > 1){
    T
  } else {
    F
  }
})

gtex_ids_SKIN <- gtex_ids_SKIN[testispos, ]
GTEX_SKIN_IDS <- gtex_ids_SKIN[,1]

save(gtex_ids_SKIN,
     GTEX_SKIN_IDS,
     GTEX_SKIN_CRUDE,
     #GTEX_SKIN_BIOTYPE,
     #GTEX_SKIN_301SP_ENTREZ,
     
     file = "GTEX_SKIN.RData")

tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2, nrows = 5)

GTEX_test <- (as.matrix(tb_promise))
#genes <- (GTEX_test[,1][2:56319])
samples <- as.vector(GTEX_test[1,])
genes <- 
  
  samples <- as.vector(tb_promise[1,])
m <- match(GTEX_SKIN_IDS, samples)


tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2,
                         colClasses = c(rep("character", 4), rep("NULL", 6745), rep("character", 449), rep("NULL", 1359)))

GTEX_test <- tb_promise


# Ovary
GTEX_test_mat <- as.matrix(GTEX_test[2:56319, 3:449])
GTEX_test_mat <- matrix(as.numeric(GTEX_test_mat), ncol = 447, nrow = 56318)


samples <- as.vector(GTEX_test[1,])

colnames(GTEX_test_mat) <- samples[3:449]
genes <- as.vector(GTEX_test[,1])[2:56319]
clear_genes <- vapply(genes, function(x){strsplit(x, '.', fixed = T)[[1]][[1]] }, "ENSG" )
rownames(GTEX_test_mat) <- clear_genes
GTEX_SKIN_CRUDE <- GTEX_test_mat

varindex <- order(apply(GTEX_SKIN_CRUDE, 1, var), decreasing = T)
GTEX_SKIN_CRUDE_15000 <- GTEX_SKIN_CRUDE[varindex[1:9100], ]

minindex <- apply(GTEX_SKIN_CRUDE_15000, 1, function(x){ (max(x) > 10) }  )
GTEX_SKIN_CRUDE_EXPRESSED <- GTEX_SKIN_CRUDE_15000[minindex,]

# make entrez genes
query <- mygene::getGenes(rownames(GTEX_SKIN_CRUDE_EXPRESSED), fields = "entrezgene")
m <- match(rownames(GTEX_SKIN_CRUDE_EXPRESSED), query@listData$query)
entrez_query <- query@listData$entrezgene[m]
rownames(GTEX_SKIN_CRUDE_EXPRESSED) <- entrez_query
GTEX_SKIN_CRUDE_ENTREZ <- GTEX_SKIN_CRUDE_EXPRESSED[(!is.na(rownames(GTEX_SKIN_CRUDE_EXPRESSED))), ]


ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
l_biotype <- getBM(attributes=c('gene_biotype',
                                'entrezgene'),
                   filters = 'entrezgene',
                   values = rownames(GTEX_SKIN_CRUDE ),
                   mart = ensembl)  

m <- match( rownames(GTEX_SKIN_CRUDE),l_biotype$entrezgene )
#identical(as.character(l_biotype$entrezgene[m]), rownames(GTEX_LIVER_CRUDE_ENTREZ))
GTEX_SKIN_BIOTYPE <- l_biotype$gene_biotype[m]

#str(gtex_ids_LIVER[,1])
str(colnames(GTEX_SKIN_CRUDE_ENTREZ))
#colnames(GTEX_SKIN_CRUDE_ENTREZ) <- samples[3:188]

m <- intersect(colnames(GTEX_SKIN_CRUDE), gtex_ids_SKIN[,1])





# make entrez genes
query <- mygene::getGenes(rownames(GTEX_SKIN_CRUDE_EXPRESSED), fields = "entrezgene")
m <- match(rownames(GTEX_SKIN_CRUDE_EXPRESSED), query@listData$query)
entrez_query <- query@listData$entrezgene[m]
rownames(GTEX_SKIN_CRUDE_EXPRESSED) <- entrez_query
GTEX_SKIN_CRUDE_ENTREZ <- GTEX_SKIN_CRUDE_EXPRESSED[(!is.na(rownames(GTEX_SKIN_CRUDE_EXPRESSED))), ]




GTEX_SKIN_LowerLeg_ENSG <- GTEX_SKIN_CRUDE[,m]

save(GTEX_SKIN_Suprapubic_ENSG,
     GTEX_SKIN_LowerLeg_ENSG,
     file = "GTEX_SKIN_C.RData", compress = "gzip")



GTEX_SKIN_301SP_ENTREZ <- GTEX_SKIN_301SP_ENTREZ[!duplicated(rownames(GTEX_SKIN_301SP_ENTREZ)), ]

Suprapubic

# blood
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/GTEX/")
id_promise <- read.table(file ="GTEx_Data_V6_Annotations_SampleAttributesDS.txt", skip = 1, fill = TRUE) # , skip =0, nrows = 10)

gtex_mat <- as.matrix(id_promise)
gtex_index <- vapply(gtex_mat[,1], function(x) { if(length(grep("^GTEX-.*", x)) == 1){
  T
}else{
  F }
} , TRUE  )
gtex_mat <- gtex_mat[gtex_index, ]

bloodgr <- apply(gtex_mat, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Whole", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  
  allcol2 <- vapply(y, function(x) { if(length(grep("Blood", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  #if(any(allcol) && any(allcol2)){
  if(any(allcol) && length(which(allcol2)) > 1){
    T
  } else {
    F
  }
})

gtex_ids_BLOOD <- gtex_mat[bloodgr,] 

testispos <- apply(gtex_ids_BLOOD, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Testis", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(length(which(allcol)) > 1){
    T
  } else {
    F
  }
})

gtex_ids_BLOOD <- gtex_ids_BLOOD[bloodgr, ]
GTEX_BLOOD_IDS <- gtex_ids_BLOOD[,1]

save(gtex_ids_BLOOD,
     GTEX_BLOOD_IDS,
     GTEX_BLOOD_CRUDE,
     #GTEX_BLOOD_BIOTYPE,
     #GTEX_BLOOD_301SP_ENTREZ,
     
     file = "GTEX_BLOOD.RData")

tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2, nrows = 5)

GTEX_test <- (as.matrix(tb_promise))
#genes <- (GTEX_test[,1][2:56319])
samples <- as.vector(GTEX_test[1,])
genes <- 
  
  samples <- as.vector(tb_promise[1,])
m <- match(GTEX_BLOOD_IDS, samples)


tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2,
                         colClasses = c(rep("character", 4), rep("NULL", 8161), rep("character", 392)))

GTEX_test <- tb_promise



GTEX_test_mat <- as.matrix(GTEX_test[2:56319, 3:392])
GTEX_test_mat <- matrix(as.numeric(GTEX_test_mat), ncol = 390, nrow = 56318)


samples <- as.vector(GTEX_test[1,])

colnames(GTEX_test_mat) <- samples[3:392]
genes <- as.vector(GTEX_test[,1])[2:56319]
clear_genes <- vapply(genes, function(x){strsplit(x, '.', fixed = T)[[1]][[1]] }, "ENSG" )
rownames(GTEX_test_mat) <- clear_genes
GTEX_BLOOD_CRUDE <- GTEX_test_mat



m <- intersect(colnames(GTEX_BLOOD_CRUDE), gtex_ids_BLOOD[,1])

GTEX_BLOOD_CRUDE <- GTEX_BLOOD_CRUDE[,m]



GTEX_BLOOD_CRUDE <- colnames(GTEX_BLOOD_CRUDE)

varindex <- order(apply(GTEX_BLOOD_CRUDE, 1, var), decreasing = T)
GTEX_BLOOD_CRUDE_15000 <- GTEX_BLOOD_CRUDE[varindex[1:9100], ]

minindex <- apply(GTEX_BLOOD_CRUDE_15000, 1, function(x){ (max(x) > 10) }  )
GTEX_BLOOD_CRUDE_EXPRESSED <- GTEX_BLOOD_CRUDE_15000[minindex,]

# make entrez genes
query <- mygene::getGenes(rownames(GTEX_BLOOD_CRUDE_EXPRESSED), fields = "entrezgene")
m <- match(rownames(GTEX_BLOOD_CRUDE_EXPRESSED), query@listData$query)
entrez_query <- query@listData$entrezgene[m]
rownames(GTEX_BLOOD_CRUDE_EXPRESSED) <- entrez_query
GTEX_BLOOD_CRUDE_ENTREZ <- GTEX_BLOOD_CRUDE_EXPRESSED[(!is.na(rownames(GTEX_BLOOD_CRUDE_EXPRESSED))), ]





# thyroid
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/GTEX/")
id_promise <- read.table(file ="GTEx_Data_V6_Annotations_SampleAttributesDS.txt", skip = 1, fill = TRUE) # , skip =0, nrows = 10)

gtex_mat <- as.matrix(id_promise)
gtex_index <- vapply(gtex_mat[,1], function(x) { if(length(grep("^GTEX-.*", x)) == 1){
  T
}else{
  F }
} , TRUE  )
gtex_mat <- gtex_mat[gtex_index, ]

thyroidgr <- apply(gtex_mat, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Thyroid", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  
  
  if(length(which(allcol)) > 1){
    T
  } else {
    F
  }
})

crblgr <- testisgr

gtex_ids_THYROID <- gtex_mat[thyroidgr,] 

testispos <- apply(gtex_ids_THYROID, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Testis", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(length(which(allcol)) > 1){
    T
  } else {
    F
  }
})

gtex_ids_THYROID <- gtex_ids_THYROID[thyroidgr, ]
GTEX_THYROID_IDS <- gtex_ids_THYROID[,1]

save(gtex_ids_THYROID,
     GTEX_THYROID_IDS,
     GTEX_THYROID_CRUDE,
     #GTEX_THYROID_BIOTYPE,
     #GTEX_THYROID_301SP_ENTREZ,
     
     file = "GTEX_THYROID.RData")

tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2, nrows = 5)

GTEX_test <- (as.matrix(tb_promise))
#genes <- (GTEX_test[,1][2:56319])
samples <- as.vector(GTEX_test[1,])
genes <- 
  
  samples <- as.vector(tb_promise[1,])
m <- match(GTEX_THYROID_IDS, samples)


tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2,
                         colClasses = c(rep("character", 4), rep("NULL", 7659), rep("character", 413), rep("NULL", 181)))

GTEX_test <- tb_promise



GTEX_test_mat <- as.matrix(GTEX_test[2:56319, 3:413])
GTEX_test_mat <- matrix(as.numeric(GTEX_test_mat), ncol = 411, nrow = 56318)


samples <- as.vector(GTEX_test[1,])

colnames(GTEX_test_mat) <- samples[3:413]
genes <- as.vector(GTEX_test[,1])[2:56319]
clear_genes <- vapply(genes, function(x){strsplit(x, '.', fixed = T)[[1]][[1]] }, "ENSG" )
rownames(GTEX_test_mat) <- clear_genes
GTEX_THYROID_CRUDE <- GTEX_test_mat



m <- intersect(colnames(GTEX_THYROID_CRUDE), gtex_ids_THYROID[,1])

GTEX_THYROID_CRUDE <- GTEX_THYROID_CRUDE[,m]



# stomach
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/GTEX/")
id_promise <- read.table(file ="GTEx_Data_V6_Annotations_SampleAttributesDS.txt", skip = 1, fill = TRUE) # , skip =0, nrows = 10)

gtex_mat <- as.matrix(id_promise)
gtex_index <- vapply(gtex_mat[,1], function(x) { if(length(grep("^GTEX-.*", x)) == 1){
  T
}else{
  F }
} , TRUE  )
gtex_mat <- gtex_mat[gtex_index, ]

stomachgr <- apply(gtex_mat, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Stomach", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  
  
  if(length(which(allcol)) > 1){
    T
  } else {
    F
  }
})

crblgr <- testisgr

gtex_ids_STOMACH <- gtex_mat[stomachgr,] 

testispos <- apply(gtex_ids_STOMACH, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Testis", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(length(which(allcol)) > 1){
    T
  } else {
    F
  }
})

gtex_ids_STOMACH <- gtex_ids_STOMACH[stomachgr, ]
GTEX_STOMACH_IDS <- gtex_ids_STOMACH[,1]

save(gtex_ids_STOMACH,
     GTEX_STOMACH_IDS,
     GTEX_STOMACH_CRUDE,
     #GTEX_STOMACH_BIOTYPE,
     #GTEX_STOMACH_301SP_ENTREZ,
     
     file = "GTEX_STOMACH.RData")

tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2, nrows = 5)

GTEX_test <- (as.matrix(tb_promise))
#genes <- (GTEX_test[,1][2:56319])
samples <- as.vector(GTEX_test[1,])
genes <- 
  
  samples <- as.vector(tb_promise[1,])
m <- match(GTEX_STOMACH_IDS, samples)


tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2,
                         colClasses = c(rep("character", 4), rep("NULL", 7294), rep("character", 191), rep("NULL", 1068)))

GTEX_test <- tb_promise



GTEX_test_mat <- as.matrix(GTEX_test[2:56319, 3:191])
GTEX_test_mat <- matrix(as.numeric(GTEX_test_mat), ncol = 189, nrow = 56318)


samples <- as.vector(GTEX_test[1,])

colnames(GTEX_test_mat) <- samples[3:191]
genes <- as.vector(GTEX_test[,1])[2:56319]
clear_genes <- vapply(genes, function(x){strsplit(x, '.', fixed = T)[[1]][[1]] }, "ENSG" )
rownames(GTEX_test_mat) <- clear_genes
GTEX_STOMACH_CRUDE <- GTEX_test_mat



m <- intersect(colnames(GTEX_STOMACH_CRUDE), gtex_ids_STOMACH[,1])

GTEX_STOMACH_CRUDE <- GTEX_STOMACH_CRUDE[,m]




# Ovary
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/GTEX/")
id_promise <- read.table(file ="GTEx_Data_V6_Annotations_SampleAttributesDS.txt", skip = 1, fill = TRUE) # , skip =0, nrows = 10)

gtex_mat <- as.matrix(id_promise)
gtex_index <- vapply(gtex_mat[,1], function(x) { if(length(grep("^GTEX-.*", x)) == 1){
  T
}else{
  F }
} , TRUE  )
gtex_mat <- gtex_mat[gtex_index, ]

ovarygr <- apply(gtex_mat, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Ovary", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(any(allcol)){
    T
  } else {
    F
  }
})

gtex_ids_TESTIS <- gtex_mat[ovarygr,] 

ovarypos <- apply(gtex_ids_TESTIS, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Ovary", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(length(which(allcol)) > 1){
    T
  } else {
    F
  }
})

gtex_ids_TESTIS <- gtex_ids_TESTIS[ovarypos, ]
GTEX_TESTIS_IDS <- gtex_ids_TESTIS[,1]

save(gtex_ids_TESTIS,
     GTEX_TESTIS_IDS,
     GTEX_TESTIS_CRUDE,
     #GTEX_TESTIS_BIOTYPE,
     #GTEX_TESTIS_301SP_ENTREZ,
     
     file = "GTEX_TESTIS.RData")

tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2, nrows = 5)

GTEX_test <- (as.matrix(tb_promise))
#genes <- (GTEX_test[,1][2:56319])
samples <- as.vector(GTEX_test[1,])
genes <- 
  
  samples <- as.vector(tb_promise[1,])
m <- match(GTEX_TESTIS_IDS, samples)


tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2,
                         colClasses = c(rep("character", 4), rep("NULL", 6018), rep("character", 106), rep("NULL", 2429)))

GTEX_test <- tb_promise


# Ovary
GTEX_test_mat <- as.matrix(GTEX_test[2:56319, 3:110])
GTEX_test_mat <- matrix(as.numeric(GTEX_test_mat), ncol = 108, nrow = 56318)


samples <- as.vector(GTEX_test[1,])

colnames(GTEX_test_mat) <- samples[3:110]
genes <- as.vector(GTEX_test[,1])[2:56319]
clear_genes <- vapply(genes, function(x){strsplit(x, '.', fixed = T)[[1]][[1]] }, "ENSG" )
rownames(GTEX_test_mat) <- clear_genes
GTEX_TESTIS_CRUDE <- GTEX_test_mat

varindex <- order(apply(GTEX_TESTIS_CRUDE, 1, var), decreasing = T)
GTEX_TESTIS_CRUDE_15000 <- GTEX_TESTIS_CRUDE[varindex[1:9100], ]

minindex <- apply(GTEX_TESTIS_CRUDE_15000, 1, function(x){ (max(x) > 10) }  )
GTEX_TESTIS_CRUDE_EXPRESSED <- GTEX_TESTIS_CRUDE_15000[minindex,]

# make entrez genes
query <- mygene::getGenes(rownames(GTEX_TESTIS_CRUDE_EXPRESSED), fields = "entrezgene")
m <- match(rownames(GTEX_TESTIS_CRUDE_EXPRESSED), query@listData$query)
entrez_query <- query@listData$entrezgene[m]
rownames(GTEX_TESTIS_CRUDE_EXPRESSED) <- entrez_query
GTEX_TESTIS_CRUDE_ENTREZ <- GTEX_TESTIS_CRUDE_EXPRESSED[(!is.na(rownames(GTEX_TESTIS_CRUDE_EXPRESSED))), ]


ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
l_biotype <- getBM(attributes=c('gene_biotype',
                                'entrezgene'),
                   filters = 'entrezgene',
                   values = rownames(GTEX_TESTIS_CRUDE ),
                   mart = ensembl)  

m <- match( rownames(GTEX_TESTIS_CRUDE),l_biotype$entrezgene )
#identical(as.character(l_biotype$entrezgene[m]), rownames(GTEX_LIVER_CRUDE_ENTREZ))
GTEX_TESTIS_BIOTYPE <- l_biotype$gene_biotype[m]

#str(gtex_ids_LIVER[,1])
str(colnames(GTEX_TESTIS_CRUDE_ENTREZ))
#colnames(GTEX_TESTIS_CRUDE_ENTREZ) <- samples[3:188]

m <- intersect(colnames(GTEX_TESTIS_CRUDE), gtex_ids_TESTIS[,1])

GTEX_TESTIS_301SP_ENTREZ <- GTEX_TESTIS_CRUDE_ENTREZ[,m]

GTEX_TESTIS_301SP_ENTREZ <- GTEX_TESTIS_301SP_ENTREZ[!duplicated(rownames(GTEX_TESTIS_301SP_ENTREZ)), ]






LNC_index <- (is.element(GTEX_PANCREAS_BIOTYPE, c("protein_coding", "lincRNA")))
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/LINC/data/")
gtex.pancreas.linc <- linc(GTEX_PANCREAS_156SP_ENTREZ[LNC_index, ], codingGenes = GTEX_PANCREAS_BIOTYPE[LNC_index]   )
gtex.pancreas.clust <- clusterlinc(gtex.pancreas.linc, pvalCutOff = 0.00000000000000005, distMethod = "dicedist" )
#gtex.liver.bio <- getbio(gtex.liver.clust, translate = 'none')                        # 5e-17


n <- 51
set <- gtex.pancreas.clust@results[[1]]$neighbours[[n]][1:20]
res <- enrichGO(set, ont = "CC")
res@result$Description
names(gtex.pancreas.clust@results[[1]]$neighbours)[n]
getGOInfo(names(gtex.pancreas.clust@results[[1]]$neighbours)[n])


set <- gtex.pancreas.clust@results[[1]]$neighbours[[n]][1:100]
res <- enrichPathway(set)

doGOsim("90632", res@result$ID[1:10])

pca_analysis <- prcomp(GTEX_PANCREAS_156SP_ENTREZ[LNC_index, ], center = FALSE,
                       scale. = FALSE) 

gtex.pancreas.linc.pc <- linc(GTEX_PANCREAS_156SP_ENTREZ[LNC_index, ], codingGenes = GTEX_PANCREAS_BIOTYPE[LNC_index], rmPC = c(4:156)   )
gtex.pancreas.clust.pc <- clusterlinc(gtex.pancreas.linc.pc, pvalCutOff = 0.00000000000000005, distMethod = "dicedist" )
#gtex.liver.bio <- getbio(gtex.liver.clust, translate = 'none')                        # 5e-17

n <- 51
set <- gtex.pancreas.clust.pc@results[[1]]$neighbours[[n]][1:20]
res <- enrichGO(set, ont = "BP")
res@result$Description
names(gtex.pancreas.clust.pc@results[[1]]$neighbours)[n]
#getGOInfo(names(gtex.pancreas.clust.pc@results[[1]]$neighbours)[n])

# LUNG
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/GTEX/")
id_promise <- read.table(file ="GTEx_Data_V6_Annotations_SampleAttributesDS.txt", skip = 1, fill = TRUE) # , skip =0, nrows = 10)

gtex_mat <- as.matrix(id_promise)
gtex_index <- vapply(gtex_mat[,1], function(x) { if(length(grep("^GTEX-.*", x)) == 1){
  T
}else{
  F }
} , TRUE  )
gtex_mat <- gtex_mat[gtex_index, ]

lunggr <- apply(gtex_mat, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Lung", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(any(allcol)){
    T
  } else {
    F
  }
})

gtex_ids_LUNG <- gtex_mat[lunggr,] 

lungpos <- apply(gtex_ids_LUNG, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Lung", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(length(which(allcol)) > 1){
    T
  } else {
    F
  }
})

gtex_ids_LUNG <- gtex_ids_LUNG[lungpos, ]
GTEX_LUNG_IDS <- gtex_ids_LUNG[,1]

save(gtex_ids_LUNG,
     GTEX_LUNG_IDS,
     GTEX_LUNG_CRUDE,
     GTEX_LUNG_BIOTYPE,
     GTEX_LUNG_301SP_ENTREZ,
     
     file = "GTEX_LUNG.RData")

tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2, nrows = 5)

GTEX_test <- (as.matrix(tb_promise))
#genes <- (GTEX_test[,1][2:56319])
samples <- as.vector(GTEX_test[1,])
genes <- 
  
  samples <- as.vector(tb_promise[1,])
m <- match(GTEX_LUNG_IDS, samples)


tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2,
                         colClasses = c(rep("character", 4), rep("NULL", 4908), rep("character", 473), rep("NULL", 3172)))

GTEX_test <- tb_promise


# Lung
GTEX_test_mat <- as.matrix(GTEX_test[2:56319, 3:477])
GTEX_test_mat <- matrix(as.numeric(GTEX_test_mat), ncol = 475, nrow = 56318)


samples <- as.vector(GTEX_test[1,])

colnames(GTEX_test_mat) <- samples[3:477]
genes <- as.vector(GTEX_test[,1])[2:56319]
clear_genes <- vapply(genes, function(x){strsplit(x, '.', fixed = T)[[1]][[1]] }, "ENSG" )
rownames(GTEX_test_mat) <- clear_genes
GTEX_LUNG_CRUDE <- GTEX_test_mat

varindex <- order(apply(GTEX_LUNG_CRUDE, 1, var), decreasing = T)
GTEX_LUNG_CRUDE_15000 <- GTEX_LUNG_CRUDE[varindex[1:9100], ]

minindex <- apply(GTEX_LUNG_CRUDE_15000, 1, function(x){ (max(x) > 10) }  )
GTEX_LUNG_CRUDE_EXPRESSED <- GTEX_LUNG_CRUDE_15000[minindex,]

# make entrez genes
query <- mygene::getGenes(rownames(GTEX_LUNG_CRUDE_EXPRESSED), fields = "entrezgene")
m <- match(rownames(GTEX_LUNG_CRUDE_EXPRESSED), query@listData$query)
entrez_query <- query@listData$entrezgene[m]
rownames(GTEX_LUNG_CRUDE_EXPRESSED) <- entrez_query
GTEX_LUNG_CRUDE_ENTREZ <- GTEX_LUNG_CRUDE_EXPRESSED[(!is.na(rownames(GTEX_LUNG_CRUDE_EXPRESSED))), ]


ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
l_biotype <- getBM(attributes=c('gene_biotype',
                                'entrezgene'),
                   filters = 'entrezgene',
                   values = rownames(GTEX_LUNG_301SP_ENTREZ),
                   mart = ensembl)  

m <- match( rownames(GTEX_LUNG_301SP_ENTREZ),l_biotype$entrezgene )
#identical(as.character(l_biotype$entrezgene[m]), rownames(GTEX_LIVER_CRUDE_ENTREZ))
GTEX_LUNG_BIOTYPE <- l_biotype$gene_biotype[m]

#str(gtex_ids_LIVER[,1])
str(colnames(GTEX_LUNG_CRUDE_ENTREZ))
#colnames(GTEX_LUNG_CRUDE_ENTREZ) <- samples[3:188]

m <- intersect(colnames(GTEX_LUNG_CRUDE_ENTREZ), gtex_ids_LUNG[,1])

GTEX_LUNG_301SP_ENTREZ <- GTEX_LUNG_CRUDE_ENTREZ[,m]

GTEX_LUNG_301SP_ENTREZ <- GTEX_LUNG_301SP_ENTREZ[!duplicated(rownames(GTEX_LUNG_301SP_ENTREZ)), ]
#GTEX_PANCREAS_BIOTYPE <- GTEX_PANCREAS_BIOTYPE[!duplicated(rownames(GTEX_PANCREAS_156SP_ENTREZ))]

#now, validate

LNC_index <- (is.element(GTEX_LUNG_BIOTYPE, c("protein_coding", "lincRNA")))
#setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/LINC/data/")
gtex.lung.linc <- linc(GTEX_LUNG_301SP_ENTREZ[LNC_index, ], codingGenes = GTEX_LUNG_BIOTYPE[LNC_index]   )
gtex.lung.clust <- clusterlinc(gtex.lung.linc, pvalCutOff = 0.00000000000000005, distMethod = "dicedist" )
#gtex.

n <- 34
set <- gtex.lung.clust@results[[1]]$neighbours[[n]][1:20]
res <- enrichGO(set, ont = "BP")
res@result$Description
names(gtex.lung.clust@results[[1]]$neighbours)[n]
getGOInfo(names(gtex.lung.clust@results[[1]]$neighbours)[n])

doGOsim("64788", res@result$ID[1:10])

pca_analysis <- prcomp(GTEX_LUNG_301SP_ENTREZ[LNC_index, ], center = FALSE, scale. = FALSE) 


gtex.lung.linc.pc <- linc(GTEX_LUNG_301SP_ENTREZ[LNC_index, ], codingGenes = GTEX_LUNG_BIOTYPE[LNC_index], rmPC = c(8:301)   )
gtex.lung.clust.pc <- clusterlinc(gtex.lung.linc.pc, pvalCutOff = 0.00000000000000005, distMethod = "dicedist" )

n <- "400550"
set <- gtex.lung.clust.pc@results[[1]]$neighbours[[n]][1:100]
res <- enrichGO(set, ont = "CC")

res <- enrichPathway(set)

res@result$Description
names(gtex.lung.clust.pc@results[[1]]$neighbours)[n]

set <- gtex.lung.clust.pc@results[[1]]$neighbours[["400550"]][1:20]
res <- enrichGO(set, ont = "BP")


# Colon
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/GTEX/")
id_promise <- read.table(file ="GTEx_Data_V6_Annotations_SampleAttributesDS.txt", skip = 1, fill = TRUE) # , skip =0, nrows = 10)

gtex_mat <- as.matrix(id_promise)
gtex_index <- vapply(gtex_mat[,1], function(x) { if(length(grep("^GTEX-.*", x)) == 1){
  T
}else{
  F }
} , TRUE  )
gtex_mat <- gtex_mat[gtex_index, ]

colongr <- apply(gtex_mat, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Sigmoid", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(any(allcol)){
    T
  } else {
    F
  }
})

gtex_ids_COLON <- gtex_mat[colongr,] 

colonpos <- apply(gtex_ids_COLON, 1, function(y){
  allcol <- vapply(y, function(x) { if(length(grep("Colon", x)) == 1){
    T
  }else{
    F }
  } , TRUE  )
  
  if(any(allcol)){
    T
  } else {
    F
  }
})


GTEX_COLON_IDS <- gtex_ids_COLON[,1]

save(gtex_ids_COLON,
     GTEX_COLON_IDS,
     GTEX_COLON_CRUDE,
     GTEX_COLON_BIOTYPE,
     GTEX_COLON_138SP_ENTREZ,
     file = "GTEX_COLON.RData")

tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2, nrows = 5)

GTEX_test <- (as.matrix(tb_promise))
#genes <- (GTEX_test[,1][2:56319])
samples <- as.vector(GTEX_test[1,])
genes <- 
  
  samples <- as.vector(tb_promise[1,])
m <- match(GTEX_COLON_IDS, samples)


tb_promise <- read.table(file = "All_Tissue_Site_Details_Analysis.combined.rpkm.gct", skip =2,
                         colClasses = c(rep("character", 4), rep("NULL", 3307), rep("character", 166), rep("NULL", 5080)))

GTEX_test <- tb_promise

GTEX_test_mat <- as.matrix(GTEX_test[2:56319, 3:166])
GTEX_test_mat <- matrix(as.numeric(GTEX_test_mat), ncol = 164, nrow = 56318)


samples <- as.vector(GTEX_test[1,])

colnames(GTEX_test_mat) <- samples[3:166]
genes <- as.vector(GTEX_test[,1])[2:56319]
clear_genes <- vapply(genes, function(x){strsplit(x, '.', fixed = T)[[1]][[1]] }, "ENSG" )
rownames(GTEX_test_mat) <- clear_genes
GTEX_COLON_CRUDE <- GTEX_test_mat

varindex <- order(apply(GTEX_COLON_CRUDE, 1, var), decreasing = T)
GTEX_COLON_CRUDE_15000 <- GTEX_COLON_CRUDE[varindex[1:9100], ]

minindex <- apply(GTEX_COLON_CRUDE_15000, 1, function(x){ (max(x) > 10) }  )
GTEX_COLON_CRUDE_EXPRESSED <- GTEX_COLON_CRUDE_15000[minindex,]

# make entrez genes
query <- mygene::getGenes(rownames(GTEX_COLON_CRUDE_EXPRESSED), fields = "entrezgene")
m <- match(rownames(GTEX_COLON_CRUDE_EXPRESSED), query@listData$query)
entrez_query <- query@listData$entrezgene[m]
rownames(GTEX_COLON_CRUDE_EXPRESSED) <- entrez_query
GTEX_COLON_CRUDE_ENTREZ <- GTEX_COLON_CRUDE_EXPRESSED[(!is.na(rownames(GTEX_COLON_CRUDE_EXPRESSED))), ]


ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
l_biotype <- getBM(attributes=c('gene_biotype',
                                'entrezgene'),
                   filters = 'entrezgene',
                   values = rownames(GTEX_COLON_138SP_ENTREZ),
                   mart = ensembl)  

m <- match( rownames(GTEX_COLON_138SP_ENTREZ),l_biotype$entrezgene )
#identical(as.character(l_biotype$entrezgene[m]), rownames(GTEX_LIVER_CRUDE_ENTREZ))
GTEX_COLON_BIOTYPE <- l_biotype$gene_biotype[m]

#str(gtex_ids_LIVER[,1])
str(colnames(GTEX_COLON_CRUDE_ENTREZ))
#colnames(GTEX_LUNG_CRUDE_ENTREZ) <- samples[3:188]

m <- intersect(colnames(GTEX_COLON_CRUDE_ENTREZ), gtex_ids_COLON[,1])

GTEX_COLON_138SP_ENTREZ <- GTEX_COLON_CRUDE_ENTREZ[,m]

GTEX_COLON_138SP_ENTREZ  <- GTEX_COLON_138SP_ENTREZ[!duplicated(rownames(GTEX_COLON_138SP_ENTREZ )), ]
#GTEX_PANCREAS_BIOTYPE <- GTEX_PANCREAS_BIOTYPE[!duplicated(rownames(GTEX_PANCREAS_156SP_ENTREZ))]

#now, validate

LNC_index <- (is.element(GTEX_COLON_BIOTYPE, c("protein_coding", "lincRNA")))
#setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/LINC/data/")
gtex.colon.linc <- linc(GTEX_COLON_138SP_ENTREZ[LNC_index, ], codingGenes = GTEX_COLON_BIOTYPE[LNC_index]   )
gtex.colon.clust <- clusterlinc(gtex.colon.linc, pvalCutOff = 0.00000000000000005, distMethod = "dicedist" )
#gtex.

n <- 36
set <- gtex.colon.clust@results[[1]]$neighbours[["171391"]][1:20]
res <- enrichGO(set, ont = "BP")
res@result$Description
names(gtex.colon.clust@results[[1]]$neighbours)[n]
getGOInfo(names(gtex.colon.clust@results[[1]]$neighbours)[n])


doGOsim("84976", res@result$ID[1:4])

pca_analysis <- prcomp(GTEX_COLON_138SP_ENTREZ[LNC_index, ], center = FALSE, scale. = FALSE) 

gtex.colon.linc.pc <- linc(GTEX_COLON_138SP_ENTREZ[LNC_index, ], codingGenes = GTEX_COLON_BIOTYPE[LNC_index], rmPC = c(3:138)   )
gtex.colon.clust.pc <- clusterlinc(gtex.colon.linc.pc, pvalCutOff = 0.00000000000000005, distMethod = "dicedist" )

n <- 36
set <- gtex.colon.clust.pc@results[[1]]$neighbours[[n]][1:20]
res <- enrichGO(set, ont = "BP")
res@result$Description
names(gtex.colon.clust.pc@results[[1]]$neighbours)[n]

# intersect of lnc RNAs in all 6 tissues
load("GTEX_LIVER.RData")
load("GTEX_LUNG.RData")
load("GTEX_PANCREAS.RData")
load("GTEX_COLON.RData")
load("GTEX_BREAST.RData")
load("GTEX_PROSTATE.RData")


a <- is.element(GTEX_LIVER_BIOTYPE, "lincRNA")
b <- is.element(GTEX_LUNG_BIOTYPE, "lincRNA")
c <-  is.element(GTEX_PANCREAS_BIOTYPE, "lincRNA")
d <- is.element(GTEX_COLON_BIOTYPE, "lincRNA")
e <- is.element(GTEX_BREAST_BIOTYPE, "lincRNA")
f <- is.element(GTEX_PROSTATE_BIOTYPE, "lincRNA")

liver <- rownames(GTEX_LIVER_111SP_ENTREZ)[a]

lung <- rownames(GTEX_LUNG_301SP_ENTREZ)[b]

pancreas <-  rownames(GTEX_PANCREAS_156SP_ENTREZ)[c]

colon <- rownames(GTEX_COLON_138SP_ENTREZ)[d]

breast <- rownames(GTEX_BREAST_201SP_ENTREZ)[e]

prostate <- rownames(GTEX_PROSTATE_96SP_ENTREZ)[f]

save(liver, lung, pancreas, colon, breast, prostate,
     common_in_all,
     file = "GTEX_SUMMARY.RData")

a <- intersect(liver, colon)
b <- intersect(breast, a)
c <- intersect(prostate, b)
d <- intersect(lung, c)
e <- intersect(pancreas, d)
common_in_all <- e

# GBA Test
# p-value test 

biotype_100 <- c(rep(F, 412), rep(T, 9000) )
biotype_20 <- c(rep(F, 412), rep(T, 9000) )

protein <- sample(rownames(GTEX_LIVER_111SP_ENTREZ[(GTEX_LIVER_BIOTYPE == "protein_coding"), ]), 8413)
rows <- rownames(GTEX_LIVER_111SP_ENTREZ[(GTEX_LIVER_BIOTYPE == "protein_coding"), ])

gtex.liver.linc <- linc(GTEX_LIVER_111SP_ENTREZ[(GTEX_LIVER_BIOTYPE == "protein_coding"), ], codingGenes = is.element(rows, protein) )

objec <- GTEX_LIVER_111SP_ENTREZ[(GTEX_LIVER_BIOTYPE == "protein_coding"), ]

a <- singlelinc(gtex.liver.linc, query = " 140828", testFun = cor.test, coExprCut = 20  )


query = NULL,
onlycor = FALSE,
testFun = NULL,
threshold = 0.05,
underth = TRUE,
coExprCut = NULL,
handleGeneIds = TRUE,
annotateFrom  = 'enrichGO',
verbose = TRUE, ...){

# rows duplicated

pc_promise <- is.element(rows, protein)


if(any(duplicated(rownames(objec)))){
  objec <- objec[(!duplicated(rownames(objec))),]
  pc_promise <- pc_promise[!duplicated(rownames(objec))]
  message("22")
}



gtex.liver.clust <- clusterlinc(gtex.liver.linc, pvalCutOff = 0.05, distMethod = "dicedist" )

cortested <- gtex.liver.clust@correlation$cortest

# get the best 10 go terms
target_list <- list()
real_storage <- vector()

for( f in 1:1000){
  entrez20 <- sort(cortested[,f], decreasing = F)[1:20]
  
  en_promise <- enrichGO(names(entrez20), ont = "BP")
  if(class(en_promise) == "enrichResult"){
    target_list[[f]] <- en_promise@result$ID[1:10]
  } else {
    target_list[[f]] <- NA
  }
  
  if(is.character(target_list[[f]]) && length(target_list[[f]]) > 0){
   real_storage[f] <- doGOsim(colnames(cortested)[f], target_list[[f]])
 } else {
  real_storage[f] <- 0 
 }
  
}
names(real_storage) <- colnames(cortested)
liver_real_storage <- real_storage
save(real_storage, "LIVER_PVALUE_REAL.RData")

# simplify

# get the best 10 go terms
real_storage <- vector()

for( f in 1:1000){
  entrez20 <- sort(cortested[,f], decreasing = F)[1:20]
  
  en_promise <- enrichGO(names(entrez20), ont = "BP")
  if(class(en_promise) == "enrichResult"){
    target <- en_promise@result$ID[1:10]
    target <- target[!is.na(target)]
  } else {
    target <- 0
  }
  
  if(is.character(target) && target != 0 && length(target) > 0 ){
    real_storage[f] <- doGOsim(colnames(cortested)[f], target)
  } else {
    real_storage[f] <- 0 
  }
  
}
#names(real_storage) <- colnames(cortested)
#liver_real_storage <- real_storage
save(real_storage, file = "liver_real_storage.RData")


mok_storage <- vector()

for( f in 1:1000){
  entrez20 <- sample(rownames(cortested), 20)
  
  en_promise <- enrichGO(entrez20, ont = "BP")
  if(class(en_promise) == "enrichResult"){
    target <- en_promise@result$ID[1:10]
    target <- target[!is.na(target)]
  } else {
    target <- 0
  }
  
  if(is.character(target) && target != 0 && length(target) > 0 ){
    mok_storage[f] <- doGOsim(colnames(cortested)[f], target)
  } else {
    mok_storage[f] <- 0 
  }
  
}
save(mok_storage, file = "LIVER_PVALUE_MOK.RData")



### now, moking

target_list <- list()
mok_storage <- vector()

for( f in 1:1000){
  entrez20 <- sample(rownames(cortested), 20)
  
  en_promise <- enrichGO(names(entrez20), ont = "BP")
  if(class(en_promise) == "enrichResult"){
    target_list[[f]] <- en_promise@result$ID[1:10]
  } else {
    target_list[[f]] <- NA
  }
  
  if(is.character(target_list[[f]]) && length(target_list[[f]]) > 0){
    mok_storage[f] <- doGOsim(colnames(cortested)[1], target_list[[f]])
  } else {
    mok_storage[f] <- 0 
  }
  
}
names(mok_storage) <- colnames(cortested)
liver_mok_storage <- mok_storage
save(liver_mok_storage, file = "LIVER_PVALUE_MOK.RData")

t.test(real_storage, mok_storage, paired = T, alternative = "greater")

wilcox.test(real_storage, mok_storage, paired = T, alternative = "greater")


# now plot
real_storage[is.na(real_storage)] <- 0
mok_storage[is.na(mok_storage)] <- 0
t.test(real_storage, mok_storage, pairwise = TRUE, alternative = "greater" )

df_boxplot <- data.frame( SIMILARITY = mok_storage, SIMILARITY2 = real_storage)
gg_freq <- ggplot(df_boxplot, aes(SIMILARITY)) +
  geom_freqpoly(colour = "blue", size = 0.8) +
  theme(
    panel.border = element_rect(color = "grey", fill = NA),   
    panel.background = element_blank()) +
  ggtitle("DISTRIBUTION OF GO TERM SIMILARITY") +
  theme(plot.title = element_text(face = "bold",
                                  color = "black")) +
  geom_freqpoly(colour =  "firebrick", size = 0.8, aes(SIMILARITY2) )



df_boxplot <- data.frame( SIMILARITY = mok_storage, SIMILARITY2 = real_storage)
gg_freq <- ggplot(df_boxplot, aes(SIMILARITY)) +
  geom_histogram(colour = "blue", size = 0.8, fill = "blue", alpha = 0.25, bins = 20) +
  theme(
    panel.border = element_rect(color = "grey", fill = NA),   
    panel.background = element_blank()) +
  ggtitle("DISTRIBUTION OF GO TERM SIMILARITY") +
  theme(plot.title = element_text(face = "bold",
                                  color = "black")) +
  geom_histogram(colour =  "firebrick", size = 0.8, aes(SIMILARITY2), fill = "firebrick", alpha = 0.25, bins = 20 )

gg_freq 




#####################
#### AFTER MASTER

# BREAST
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/GTEX/")
require(biomaRt)

varindex <- order(apply(GTEX_BREAST_CRUDE, 1, var), decreasing = T)
GTEX_BREAST_CRUDE_VAR <- GTEX_BREAST_CRUDE[varindex[1:2500], ]

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
l_biotype <- getBM(attributes=c('gene_biotype',
                                'ensembl_gene_id'),
                   filters = 'ensembl_gene_id',
                   values = rownames(GTEX_BREAST_CRUDE),
                   mart = ensembl)  
m <- match( rownames(GTEX_BREAST_CRUDE_VAR),l_biotype$ensembl_gene_id )
#identical(as.character(l_biotype$entrezgene[m]), rownames(GTEX_LIVER_CRUDE_ENTREZ))
GTEX_BREAST_BIOTYPE <- l_biotype$gene_biotype[m]


l_index <- (l_biotype$gene_biotype == "lincRNA")

ENSG_LINCRNAS <- l_biotype$ensembl_gene_id[l_index]
save(ENSG_LINCRNAS, file = "GTEX_ANNOTATION.RData")

load("GTEX_ANNOTATION.RData")


GTEX_RAW <- list()

GTEX_RAW[[1]] <- GTEX_BREAST_CRUDE[,colnames(GTEX_BREAST_201SP_ENTREZ)]
GTEX_RAW[[2]] <- GTEX_LUNG_CRUDE[,colnames(GTEX_LUNG_301SP_ENTREZ)]
GTEX_RAW[[3]] <- GTEX_COLON_CRUDE[, colnames(GTEX_COLON_138SP_ENTREZ)]
GTEX_RAW[[4]] <- GTEX_PANCREAS_CRUDE
GTEX_RAW[[5]] <- GTEX_STOMACH_CRUDE
GTEX_RAW[[6]] <- GTEX_THYROID_CRUDE
GTEX_RAW[[7]] <- GTEX_LIVER_CRUDE[, colnames(GTEX_LIVER_111SP_ENTREZ)]
GTEX_RAW[[8]] <- GTEX_OVARY_CRUDE
GTEX_RAW[[9]] <- GTEX_TESTIS_CRUDE
GTEX_RAW[[10]] <- GTEX_PROSTATE_CRUDE
GTEX_RAW[[11]] <- GTEX_BLOOD_CRUDE
GTEX_RAW[[12]] <- GTEX_SKIN_LowerLeg_ENSG
GTEX_RAW[[13]] <- GTEX_HPC_CRUDE
GTEX_RAW[[14]] <- GTEX_CRBL_CRUDE
GTEX_RAW[[15]] <- GTEX_CTX_CRUDE

names(GTEX_RAW) <- c("BREAST", "LUNG", "COLON", "PANCREAS", "STOMACH", "THYROID", "LIVER", "OVARY", "TESTIS", "PROSTATE", "BLOOD", "SKIN", "HPC", "CRBL", "CTX")

save(ENSG_LINCRNAS, LINC_ATT, LINC_T_CRITERIA, CANDIDATES, file = "GTEX_ANNOTATION.RData", compress = TRUE)
save(GTEX_RAW, file = "GTEX_RAW.RData", compress = TRUE)

# get attributes of LINCRNAs, median expression, variance
LINC_ATT <- lapply(GTEX_RAW, function(x, ENSG_LINCRNAS){
      all_lincs <- x[ENSG_LINCRNAS,]
      NAMES <- ENSG_LINCRNAS
      MEDIANS <- Biobase::rowMedians(all_lincs)
      VAR <- apply(all_lincs, 1, var)
      entry <- data.frame(as.character(NAMES), MEDIANS, VAR)
      rownames(entry) <- NAMES
      return(entry)
      }, ENSG_LINCRNAS)

#
LINC_ATT_TF <- lapply(LINC_ATT, function(x){
       var_lim <- 0.2*mean(x$VAR)
       EXPR_LEVEL <- x$MEDIANS > 1
       VAR_LIM <- x$VAR > var_lim
       BOTH <- (EXPR_LEVEL + VAR_LIM) > 1
      # entry <- data.frame(EXPR_LEVEL, VAR_LIM, BOTH)
       return(BOTH)
    })

# reduce to candidates
reduceLINCS <- function(LINC_ATT_TF, ENSG_LINCRNAS){
 df  <- as.data.frame(LINC_ATT_TF)
 cand <- apply(df, 1, function(x){ 
   if(length(x[x == TRUE]) > 2){
   tf  <- TRUE
   } else {
   tf  <- FALSE 
   }
   return(tf)
   })
 return(ENSG_LINCRNAS[cand])
}

CANDIDATES <- reduceLINCS(LINC_ATT_TF, ENSG_LINCRNAS)

LINC_RED_TF <- lapply(LINC_ATT_TF, function(x){
  return(x[is.element(ENSG_LINCRNAS, CANDIDATES)])
  })

LINC_T_CRITERIA <- as.data.frame(LINC_RED_TF)
rownames(LINC_T_CRITERIA) <- CANDIDATES

### get expression values for histo

texpr <- mapply(function(x, y){
       texpr <- x[CANDIDATES[y],]
       return(texpr)
       }
       , GTEX_RAW, LINC_RED_TF)
# plot histo
require(ggplot2)
histo_vec <- (unlist(texpr))
ggplot(data.frame(FPKM = histo_vec)) + geom_histogram(aes(FPKM), bins = 1000) + xlim(0, 100) + geom_vline(xintercept = median(histo_vec), colour = "red") + ggtitle(paste("median expression:", round(median(histo_vec), 2)))


toccur <- apply(as.data.frame(LINC_RED_TF), 1, function(x){ length(x[x == TRUE]) } )

ggplot(data.frame(TISSUES = toccur)) + geom_histogram(aes(TISSUES), bins = 15) + xlim(0, 16) + geom_vline(xintercept = median(toccur), colour = "red")

# variance plot
tvar <- mapply(function(x, y){
  tvar <- x[CANDIDATES[y],]
  TVAR <- apply(tvar, 1, var)
  return(TVAR)
}
, GTEX_RAW, LINC_RED_TF)


histo_var <- (unlist(tvar))
ggplot(data.frame(VAR = histo_var)) + geom_histogram(aes(VAR), bins = 100) + xlim(0, 250) + geom_vline(xintercept = median(histo_var), colour = "red") + ggtitle(paste("median variance:", round(median(histo_var), 2)))

# plot samples from GTEx as boxplots
require(ggplot2)
require(reshape2)

for(n in 1:15){
em_promise <- GTEX_RAW[[n]]
df_boxplot  <- suppressMessages(melt(em_promise[1:5000, 1:25 ]))
names(df_boxplot) <- c(NA, "SAMPLES", "VALUE")
gg_box <- ggplot(df_boxplot,
                 aes(SAMPLES, VALUE)) + geom_boxplot(outlier.color =
              "firebrick", colour = "dodgerblue3") +  theme(
          panel.border = element_rect(color = "grey", fill = NA),   
        panel.background = element_blank()) +
        ylim(0, 10) + 
#boxplot(VALUE ~ SAMPLES, data = df_boxplot[1:100000,])

  ggtitle(paste("SAMPLES:", names(GTEX_RAW))) +
  theme(plot.title = element_text(face = "bold",
                                 color = "steelblue4"))
}


