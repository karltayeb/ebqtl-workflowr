setwd('/project2/gilad/ktayeb/ebqtl-workflowr/')
source('.Rprofile')
renv::activate()

library(zellkonverter)
library(tictoc)
library(SingleCellExperiment)
library(CelliD)
library(progress)
library(dplyr)
library(tidyverse)
library(Matrix)

loadMarkerGeneSet <- function(){
    require(dplyr)
    require(tidyr)
    require(tibble)
    marker_mat_path <- 'data/marker_mat.csv'
    marker_mat <- read.csv(marker_mat_path, sep='\t', header=T, row.names=1)
    markers <- marker_mat %>%
       rownames_to_column() %>%
       pivot_longer(!rowname) %>%
       filter(value !=0) %>%
       group_by(name) %>%
       summarise(geneset = list(rowname))
    markers <- setNames(markers$geneset, markers$name)
    return(markers)
}

loadCelliD <- function(path, adata_path=NULL){
    cellid_res <- readRDS(path)
    cellCoord <- cellid_res[[1]]
    geneCoord <- attr(cellid_res[[1]], 'genesCoordinates')
    dims = c(dim(geneCoord)[1], dim(cellCoord)[1])
    if (is.null(adata_path)){
        X <- sparseMatrix(i=integer(0), j = integer(0), dims=dims)
        sce <- SingleCellExperiment(X)
        rownames(sce) <- rownames(geneCoord)
        colnames(sce) <- rownames(cellCoord)
    }else{
        sce <- zellkonverter::readH5AD(
            file=adata_path,
            use_hdf5=TRUE,
            uns=FALSE,
            varm=FALSE, obsm=FALSE,
            varp=FALSE, obsp=FALSE
        )
    }
    reducedDim(sce, 'MCA') <- cellCoord
    mainExpName(sce) <- names(cellid_res)[[1]]
    return(sce)
}

predictFromHGT <- function(hgt){
    prediction <- rownames(hgt)[apply(hgt, 2, which.max)]
    pval <- apply(hgt, 2, max)
    prediction_signif <- ifelse(pval>2, yes = prediction, "unassigned")
    return(list(
        prediction=prediction,
        prediction_signif=prediction_signif,
        pval=pval))
}

cellid_path <- 'output/cellid/EB_high_pass_filtered_normalized/EB_high_pass_filtered_normalized.donor_id.cellidNA18489.rds'
cellid_path <- snakemake@input[[1]]
res_path <- snakemake@output[[1]]
res_path2 <- snakemake@output[[2]]


group <- stringr::str_extract_all(cellid_path,"(?<=\\.cellid).+(?=\\.rds)")

print('load cellid results')
cellid_ind <- loadCelliD(cellid_path)


print('load reference cellid results')
cao_cellid_path <- 'output/cellid/fetal_reference.cellid.rds'
cao_adata_path <- 'data/CaoEtAl.Obj.CellsOfAllClusters.ProteinCodingGenes.h5ad'
cellid_cao <- loadCelliD(cao_cellid_path, cao_adata_path)


#' marker gene assignment
print('make marker gene assignments')
marker_gs <- loadMarkerGeneSet()
ind_marker_hgt <- CelliD::RunCellHGT(cellid_ind, pathways=marker_gs)
marker <- c(
    list(group=group, pathways='markers', hgt=ind_marker_hgt),
    predictFromHGT(ind_marker_hgt))

#' fetal reference assignment
print('making fetal reference assignments')
cao_cell_gs <- GetCellGeneSet(cellid_cao, dims=1:50, n.features=200)
cao_group_gs <- GetGroupGeneSet(
  cellid_cao, dims=1:50, n.features=200, group.by='Main_cluster_name')
ind_cao_hgt <- CelliD::RunCellHGT(cellid_ind, pathways=cao_group_gs)
cao <- c(
  list(group=group, pathways='cao', hgt=ind_cao_hgt),
  predictFromHGT(ind_cao_hgt))

print(tail(sort(table(marker$prediction_signif))))
print(tail(sort(table(cao$prediction_signif))))
res <- as_tibble(as.data.frame(rbind(cao, marker))) %>%
    unnest_longer(group) %>% 
    unnest_longer(pathways)

print('save output')
saveRDS(res, res_path)
write.csv(t(as.data.frame(cao$hgt)), res_path2)
