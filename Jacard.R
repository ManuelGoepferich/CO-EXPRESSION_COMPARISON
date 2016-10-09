# JACCARD DISTANCE

# measure the jaccard distance for a list of named co-expressed genes

# DEFINITION OF FUNCTION:
measureJaccard <- function(m_list, convert = TRUE){
len <- length(m_list)

measure <- function(x, y){
   dist <- (1 - (length(intersect(x, y))/ length(union(x, y))))
   return(dist)
}

dist_matrix <- matrix(NA, nrow = len, ncol = len)
for(n in 1:len){
  for(m in 1:len){
    dist_matrix[n, m] <- measure(m_list[[n]], m_list[[m]] )
  }
}
colnames(dist_matrix) <- names(m_list)
rownames(dist_matrix) <- names(m_list)

if(convert) dist_matrix <- as.dist(dist_matrix)
return(dist_matrix)
}

# Example 1, most easy
ll <- list(c("A", "F", "H"),
           c("A", "J", "H"),
           c("I", "L", "H"),
           c("P", "F", "L"),
           c("M", "X", "Y"))

names(ll) <- c(1:5)

measureJaccard(ll)
measureJaccard(ll, convert = FALSE)

# Example 2, with the LINC package
# STEP 1: get the co-expressed genes
norad <- lapply(list(gbm_cluster,  
                     ctx_cluster,
                     hpc_cluster,
                     crbl_cluster),
                function(x){
                  getcoexpr(x, '647979')})

# STEP 2: name the list
names(norad) <- c("CANCER_GBM",
                  "HEALTHY_CTX",
                  "HEALTHY_HPC",
                  "HEALTHY_CRBL")

# STEP 3: Jaccard distance
norad_dist <- measureJaccard(norad)
mean(norad_dist) # if this value is high => lincRNA has different interaction partners across tissues
 

