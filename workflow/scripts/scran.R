library(scran)
library(scRNAseq)
library(zellkonverter)
library(HDF5Array)
library(tictoc)
library(SingleCellExperiment)
library(dplyr)
library(tidyverse)

split <- function(groups, groups2, min.size=20000, max.size=100000){
    groups2 <- as.integer(groups2)
    to_split <- names(which(table(groups) > max.size))
    to_split2 <- names(which(table(groups2) > min.size))
    groups <- dplyr::if_else(
        (groups %in% to_split) & (groups2 %in% to_split2),
        groups2, groups
    )
    return(groups)
}

# TODO: maybe go back and change the pre-clustering
# right now it is just the result of recursive kmeans
# on top 100 PCs of raw count matrix
pre_clusters <- read.table(
    'output/clustering/EB_high_pass_filtered_heirarchecal_kmeans_k3_r5.txt',
    header=TRUE
)

# get coarse clustering from heirarchecal kmeans
# get cluster range from ~25k cells to ~250k cells
# will call quickCluster on each large cluster group
groups <- pre_clusters$r0
sort(table(groups))
groups <- split(groups, pre_clusters$r1)
sort(table(groups))
groups <- split(groups, pre_clusters$r2)
sort(table(groups))
groups <- split(groups, pre_clusters$r3)
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

clusters <- rep("0", length(cells))
names(clusters) <- cells

tic('Making subclusters...')
for (cluster_id in names(sorted_groups)){
    tic(cluster_id)
    x <- as(h5[,groups==cluster_id], "dgCMatrix")
    group_clusters <- quickCluster(x)
    group_clusters  <- paste0(cluster_id, '_', group_clusters)
    clusters[groups==cluster_id] <- group_clusters
    toc()
}
toc()


# save clusters
write.table(clusters, 'output/scran/clusters.txt')
clusters <- read.table('output/scran/clusters.txt')

# run scran
tic('computing size factors with scran...')
size_factors <- calculateSumFactors(h5, clusters=as.factor(clusters[, 1]))
toc()

names(size_factors) <- cells
write.table(size_factors, 'output/scran/sizefactors.txt')

final <- cbind(clusters, size_factors)
colnames(final) <- c('scranCluster', 'sizeFactor')
final %>%
    mutate(scranGroup = stringr::str_split(scranCluster, '_') %>% purrr::map_chr(., 1)) %>%
    relocate(scranCluster, .after=scranGroup) %>%
    rownames_to_column('cell') %>%
    readr::write_csv('output/scran/final.txt')
