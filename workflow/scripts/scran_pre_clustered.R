library(scran)
library(scRNAseq)
library(zellkonverter)
library(tictoc)

# TODO: maybe go back and change the pre-clustering
# right now it is just the result of recursive kmeans
# on top 100 PCs of raw count matrix

tic('loading clustering')
pre_clusters <- read.table(
    'output/clustering/EB_high_pass_filtered_heirarchecal_kmeans_k3_r5.txt',
    header=TRUE
)
toc()

# this was done manually picked clusters with <250k cells
# cluster range from ~25k cells to ~250k cells
pc <- pre_clusters$r0
pc[pc==0] <- pre_clusters[pc==0,]$r1
sorted_cluster_ids <- names(sort(table(pc)))

tic('loading SingleCellExperiment from H5AD')
path = 'output/proc/EB_high_pass_filtered_subsample.h5ad' 
path = 'output/proc/EB_high_pass_filtered.h5ad' 
sce <- zellkonverter::readH5AD(
    file=path,
    use_hdf5=TRUE,
    uns=FALSE,
    varm=FALSE, obsm=FALSE,
    varp=FALSE, obsp=FALSE
)
toc()

scran_subcluster <- function(cluster_id){
    tic(paste0('running scran for pre-cluster: ', cluster_id, '... # of cells: ', sum(pc==cluster_id))) 
    # run scran for each coarse cluster
    scesub <- sce[,pc == cluster_id]

    tic('loading subset into memory');
    x <- as(assay(scesub, 'X'), 'dgCMatrix');
    counts(scesub) <- x
    toc()

    tic('clustering for scran...')
    clusters <- quickCluster(
      x = x,
      use.ranks = FALSE,
      BSPARAM = BiocSingular::IrlbaParam()
    )
    toc()

    tic('running scran')
    scesub <- computeSumFactors(scesub, clusters=as.factor(clusters))
    toc()

    print('saving size factors...')
    sf <- scesub$sizeFactor
    names(sf) <- colnames(scesub)
    out <- rbind(sf, clusters)
    path = paste0('output/scran/sizeFactor_cluster', cluster_id, '.txt')
    write.csv(out, path)
    toc()
}

for (cluster_id in sorted_cluster_ids){
    print(cluster_id)
    scran_subcluster(cluster_id)
}

