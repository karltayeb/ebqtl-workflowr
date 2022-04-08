setwd('/project2/gilad/ktayeb/ebqtl-workflowr/')
source('.Rprofile')

library(zellkonverter)
library(basilisk)
library(anndata)
library(tictoc)
library(SingleCellExperiment)
library(CelliD)
library(progress)
library(tidyverse)

# convert and n x m dgRMatrix
# to an m x n dgCMatrix
dgr2dgc <- function(X){
	require(Matrix)
	X2 <- Matrix(nrow=dim(X)[2], ncol=dim(X)[1], data=0, sparse=T)
	X2@i <- X@j
	X2@p <- X@p
	X2@x <- X@x
	X2@Dimnames <- rev(X@Dimnames) # `dgCMatrix`
	return(X2)
}

# take anndata object
# copy 'X' and adata.obs into 'counts' and 'colData'
lightAnnData2SCE <- function(adata){
	require(SingleCellExperiment)
	sce <- SingleCellExperiment(
		list(X=dgr2dgc(adata$X)),
		colData=DataFrame(adata$obs)
	)
	return(sce)
}

adata_path <- snakemake@input[[1]]; print(adata_path)
cluster_path <- snakemake@input[[2]]; print(cluster_path)

save_dir <- snakemake@params[['prefix']]; print(save_dir)
manifest_path <- snakemake@output[[2]]; print(manifest_path)

use_hdf5 <- snakemake@params[['use_hdf5']] == 'TRUE'; print(use_hdf5)
nmcs <- snakemake@params[['nmcs']]; print(nmcs)

tic('loading H5AD')
proc <- basiliskStart(zellkonverter::zellkonverterAnnDataEnv)
adata <- basiliskRun(proc, function(){read_h5ad(adata_path)})
toc()

# load annotation, make manifest
clusters <- read.table(cluster_path, sep='\t', header=T)
clusters <- dplyr::rename(clusters, cell=1)  
group2cell <- clusters %>%
	filter(cell %in% adata$obs_names) %>%
	group_by(across(c(-cell))) %>%
	chop(cell)
g2c <- group2cell$cell
names(g2c) <- group2cell[,1][[1]]

groups <- clusters %>% select(-cell) %>% table %>% sort
paths <- sapply(names(groups), function(x){paste0(save_dir, '/', x, '.rds')})
paths <- gsub(" ", "_", paths)  # no spaces in file names
manifest <- t(rbind(names(groups),groups, paths))
colnames(manifest) = c('cluster', 'size', 'path')
names(paths) <- names(groups)
print(groups)
print(manifest)

# run for each group
pb <- progress_bar$new(total = length(groups))
for (group in names(groups)) {
    pb$tick()
    tmp_path <- paths[group]
    print(tmp_path)
    if(!file.exists(tmp_path)){
        tic('load subset into SCE')
		cells <- g2c[[group]]
		adata_sub <- adata[cells]
		sce_group <- lightAnnData2SCE(adata_sub)
        logcounts(sce_group) <- as(assay(sce_group, "X"), 'dgCMatrix')
        toc()
        sce_group <- RunMCA(sce_group, nmcs=nmcs)
        mca <- reducedDim(sce_group, 'MCA')
        saveRDS(list(group=mca), tmp_path)
    } else{
        print(paste0('already fit ', group))
    }
}

write.csv(manifest, file=manifest_path)
