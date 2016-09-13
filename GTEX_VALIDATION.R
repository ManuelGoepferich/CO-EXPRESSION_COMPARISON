# GTEX validation

# set directory (local)
setwd("D:/data/Eigene Dateien/Documents/Dateien/Universität/Master/MASTER/GTEX/")

#require(biomaRt)
#varindex <- order(apply(GTEX_BREAST_CRUDE, 1, var), decreasing = T)
#GTEX_BREAST_CRUDE_VAR <- GTEX_BREAST_CRUDE[varindex[1:2500], ]

#ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#l_biotype <- getBM(attributes=c('gene_biotype',
#                                'ensembl_gene_id'),
#                   filters = 'ensembl_gene_id',
#                   values = rownames(GTEX_BREAST_CRUDE),
#                   mart = ensembl)  
#m <- match( rownames(GTEX_BREAST_CRUDE_VAR),l_biotype$ensembl_gene_id )
#identical(as.character(l_biotype$entrezgene[m]), rownames(GTEX_LIVER_CRUDE_ENTREZ))
#GTEX_BREAST_BIOTYPE <- l_biotype$gene_biotype[m]

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


