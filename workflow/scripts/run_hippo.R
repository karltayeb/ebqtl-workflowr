#devtools::install_github('https://github.com/karltayeb/HIPPO')
#setwd('../../')

library(zellkonverter)
library(HDF5Array)
library(tictoc)
library(SingleCellExperiment)
library(dplyr)
library(tidyverse)
library(HIPPO)


split <- function(groups, groups2, min.size=50000, max.size=150000){
    groups2 <- as.integer(groups2)
    to_split <- names(which(table(groups) > max.size))
    to_split2 <- names(which(table(groups2) > min.size))
    groups <- dplyr::if_else(
        (groups %in% to_split) & (groups2 %in% to_split2),
        groups2, groups
    )
    return(groups)
}
pre_clusters <- read.table(
    'output/clustering/EB_high_pass_filtered_heirarchecal_kmeans_k3_r5.txt',
    header=TRUE
)
groups <- pre_clusters$r0; table(groups)
groups <- split(groups, pre_clusters$r1); table(groups)
sorted_groups <- sort(table(groups))

tic('loading SingleCellExperiment from H5AD')
path = 'output/proc/EB_high_pass_filtered.h5ad' 
sce <- zellkonverter::readH5AD(
    file=path,
    use_hdf5=TRUE,
    uns=FALSE,
    varm=FALSE, obsm=FALSE,
    varp=FALSE, obsp=FALSE
)
h5 <- H5ADMatrix(path)
toc()

cells <- colnames(sce)
for (cluster_id in names(sorted_groups)){
    save_path = paste0('output/hippo/hippo', cluster_id)
    if (!file.exists(save_path)){
        tic(cluster_id)
        tic('load expression')
        x <- as(h5[,groups==cluster_id], "dgCMatrix")
        sce <- SingleCellExperiment(list(counts=x))
        toc()
        tic('fit HIPPO')
        sce <- hippo(sce, K=1000, max_features=1000, min_cluster_size=100, outlier_proportion=0.0001)
        toc()
        hippo_res <- sce@int_metadata$hippo
        hippo_res$X <- NULL
        hippo_res$clusters <- colData(sce)
        rownames(hippo_res$clusters) <- cells[groups==cluster_id]
        saveRDS(hippo_res, save_path)
        toc()
    }
}

paths <- list.files('../../output/hippo', full.names=T)

path <- paths[1]
hippo_res <- as.data.frame(readRDS(path)$clusters)
hippo_res$hippoGroup <- tail(stringr::str_split(path, '/')[[1]], 1)
pb <- progress::progress_bar$new(total = length(paths)-1)
pb$tick(0)
for (path in tail(paths, -1)){
    pb$tick()
    tmp <- as.data.frame(readRDS(path)$clusters)
    tmp$hippoGroup <- tail(stringr::str_split(path, '/')[[1]], 1)
    hippo_res <- dplyr::bind_rows(hippo_res, tmp)
}

write.csv(hippo_res, '../../output/hippo/clusters.txt')
hippo_res <- as.data.frame(readRDS(path)$clusters)

# Convert to more compact binary cluster matrix
# Also makes tree structure more explicit
makeBinaryClusterMatrix <- function(hippo_res){
    cluster_tables <- apply(hippo_res, 2, table)
    splits <- purrr::map_dbl(
        seq(length(cluster_tables) - 2),
        ~ which.min(head(cluster_tables[[.x + 1]], -1) - cluster_tables[[.x]])[[1]])

    df <- data.frame(matrix(ncol=1, nrow=dim(hippo_res)[1]))
    names(df) <- c('V1')
    rownames(df) <- rownames(hippo_res)
    lev <- rep(1, dim(hippo_res)[2])

    for (k in tail(seq(length(cluster_tables) - 1), -1)){
        kold <- splits[k-1] 
        l <- lev[kold]; lev[kold] <- lev[kold] + 1; lev[k] <- lev[kold]
        idx1 <- (hippo_res[, k] == k) & (hippo_res[, k-1] == kold)
        idx0 <- (hippo_res[, k] == kold) & (hippo_res[, k-1] == kold)
        df[idx0, l] = 0
        df[idx1, l] = 1
    }
    return(df)
}

binary_res <- c()
pb <- progress::progress_bar$new(total = length(paths))
pb$tick(0)
for (path in paths){
    pb$tick()
    tmp <- as.data.frame(readRDS(path)$clusters)
    tmp <- makeBinaryClusterMatrix(tmp)
    tmp$hippoGroup <- tail(stringr::str_split(path, '/')[[1]], 1)
    binary_res <- dplyr::bind_rows(binary_res, tmp)
}
for (path in head(paths, 4)){
    pb$tick()
    tmp <- as.data.frame(readRDS(path)$clusters)
    tmp <- makeBinaryClusterMatrix(tmp)
    tmp$hippoGroup <- tail(stringr::str_split(path, '/')[[1]], 1)
    binary_res <- dplyr::bind_rows(binary_res, tmp)
}
write.csv(binary_res, '../../output/hippo/binary_clusters.txt')

