setwd('/project2/gilad/ktayeb/ebqtl-workflowr/')
source('.Rprofile')
renv::activate()

library(zellkonverter)
library(tictoc)
library(SingleCellExperiment)
library(CelliD)
library(progress)

adata_path <- snakemake@input[[1]]; print(adata_path)
save_path <- snakemake@output[[1]]; print(save_path)
group_col <- snakemake@params[['group_col']]; print(group_col)
use_hdf5 <- snakemake@params[['use_hdf5']] == 'TRUE'; print(use_hdf5)
nmcs <- snakemake@params[['nmcs']]; print(nmcs)

tic('loading SingleCellExperiment from H5AD')
# read h5ad, assume AnnData.X is already normalized/log transformed
sce <- zellkonverter::readH5AD(
    file=adata_path,
    use_hdf5=use_hdf5,
    uns=FALSE,
    varm=FALSE, obsm=FALSE,
    varp=FALSE, obsp=FALSE
)
toc()

# split by group_col
sce[['all']] <- 1  # default if no real group col specified
groups <- names(sort(table(sce[[group_col]])))
print(sort(table(sce[[group_col]])))

# run for each group
pb <- progress_bar$new(total=length(groups))
base_path <- substr(save_path, 1, nchar(save_path)-4)
for (group in groups){
    pb$tick()
    tmp_path <- paste0(base_path, group, '.rds')
    if(!file.exists(tmp_path)){
        tic('load subset into memory')
        sce_group <- sce[, sce[[group_col]]==group]
        logcounts(sce_group) <- as(assay(sce_group, "X"), 'dgCMatrix')
        toc()
        sce_group <- RunMCA(sce_group, nmcs=nmcs)
        mca <- reducedDim(sce_group, 'MCA')
        saveRDS(list(group=mca), tmp_path)
    } else{
        print(paste0('already fit ', group))
    }
}

results <- list()
for (group in groups){
    tmp_path <- paste0(base_path, group, '.rds')
    tic('load saved fit')
    mca <- readRDS(tmp_path)
    results[[group]] <- mca[[group]]
    toc()
}

tic('saving results')
saveRDS(results, save_path)
toc()

# assay(sce[,sce$donor_id == 'NA19210'], 'X')
# logcounts(sce) <- assay(sce, "X")
# 
# tic('run MCA')
# sce <- RunMCA(sce, nmcs=10)
# mca <- reducedDim(sce, 'MCA')
# saveRDS(mca, save_path)

# marker_mat_path <- 'data/marker_mat.csv'
# attr(mca, 'genesCoordinates')
# marker_mat <- read.csv(marker_mat_path, sep='\t', header=T, row.names=1)
# markers <- marker_mat %>%
#    rownames_to_column() %>%
#    pivot_longer(!rowname) %>%
#    filter(value !=0) %>%
#    group_by(name) %>%
#    summarise(geneset = list(rowname))
#markers <- setNames(markers$geneset, markers$name)
#celltype_hgt <- RunCellHGT(sce, pathways = markers, dims = 1:10, n.features = 200)
#prediction <- rownames(celltype_hgt)[apply(celltype_hgt, 2, which.max)]
#prediction_signif <- ifelse(apply(celltype_hgt, 2, max)>2, yes = prediction, "unassigned")
