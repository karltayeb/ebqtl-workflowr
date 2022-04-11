# ebqtl 

This repo is a mashup of a workflowr site and a Snakemake workflow. 
All of the code you'll want is under `workflow/`


## Normalization with `scran`

Scran performs normalization by estimating size factors of many overlapping pseudocell samples (sum of counts across multiple cells) **within a cluster of relatively homogeneous cells**, and then deconvolving them into cell level estimates. The pseudocells have lower proportion of zero counts, which is important for the stability of the size factor estimation  procedures typical in RNA-seq.

That said, the main bottleneck for scran on data this large is computing the clustering. `scran` ships with a heirarchecal clustering strategy that is feasible on the order of ~10k-100k cells.

What I did was cluster the raw counts (really the top PCs of the raw counts) to get a rough clustering of the data. Then, I fed each of these rough clusters to scrans clustering (assuming that the default parameters are tailored to the needs of `scran` for size factor estimation). The final clustering is in `output/scran/clusters.txt` which should have ALL cells and clusters like `{kmeans_cluster}_{scran_heirarchecal_cluster}`. 

Then, we run scran on ALL cells with this clutering. The size factors are saved in `output/scran/final.txt`

### Key outputs:
* `output/scran/final.txt` a file with the size factors and cluster ID for each cell. `scran` deconvolution was run within each cluster and then across cluster normalization was applied as described in the scran paper.


## CellId

### Pipeline

If you want to run CellId with a new set of clusters you can do that with the rules in `workflow/rules/celltype_assignment.smk`

There are two rules
* `run_cellid_per_cluster`: takes an anndata object and a tsv assigning cells to clusters. It runs cellid for each cluster seperately and saves the result of each fit to a `.rds`. Since there are a variable number of outputs, we also return a `manifest.txt` which reports all the files generated

* `cellid_assign`: takes the `manifest.txt` from the above rule and extracts celltype assignments mapping to the Cao et. al fetal reference data. This works by extracting the top genes for each cell-type in the fetal reference (informed by a cellid fit to the fetal reference) and then performing a hypergeometric test with the top genes in each cell of EB cells. The cell gets assigned to the cell-type exhibiting the strongest enrichment (smallest p-value).   

### Things you might want to try

#### Cell-to-cell vs cell-to-group

Cellid lets you do either cell-to-cell matching (matching cell in test set to cell in reference set) to do label transfer or cell-to-group matching (matching cell in test set to best annotated cluster in the reference). Right now I have it using 

#### Extending the reference set

We probably want to do celltype annotation not just with the fetal reference data, but with fetal reference + hESCs (or other populations of stem cells). To do this you should:
1. make a combined `AnnData` object with all the reference cells, and run CellId on that.
2. Then go to `workflow/scripts/cellid_assign.R` and change the hard-coded paths `cao_cellid_path` and `cao_adata_path`.
3. We're thinking about running cellid for each germ layer seperately (see below). Katie has a mapping of fetal reference cell-types to germ layer. So it would be good to make a subsetted `AnnData` for the fetal reference in each germ layer.
3. Bonus: The hard coded paths are not the most elegant, since you're going to want to try assignment with multiple reference sources you might want to make these inputs in the snakemake rule.


#### Running with different subpopulations (e.g. run each germ layer seperately)

The genes that are most specific to a cell population (that is, help you identify that cell population) are probably different if you as the question globally (which genes are specific to this cell population out of all cells?) vs (which genes are specific to this cell population, with respect to a subset of all cells). This is one of the principles that motivates clustering methods like HIPPO, and it should inform how we use CellId too, perhaps.

I've already run CellId for each germ layer. To complete the analyis you want to assign them to germ-layer specific references (see above)

### References

The cellid vignette is a useful source for editing these scripts/rules:
https://www.bioconductor.org/packages/release/bioc/vignettes/CelliD/inst/doc/BioconductorVignette.html

## clustering with HIPPO

### Overview and rationale
HIPPO is basically a recursive k-means (or other clustering) where the features selected at each clustering step. Features (genes) are selected by computing the proportion of zeros and selecting the genes with the most significant deviation from expectation (under a poisson model). The rationale for this is that homogeneous cell populations are modelled as having poisson gene counts. When genes are not very highly expressed, changes in the zero proportion are pretty sensitive for detecting inhomogenous rate parameter.

The heirarchecal clustering is appealing because 
1. It gives us a heirarchecal clustering of the cells, which is potentially useful in understanding the relationships/shared lineage of various cell sub-populations.
2. Each split is associated with a set of "most heterogeneous genes", and after each split you can evaluate how well the clustering resolved expression heterogeneity in those genes.
3. Selecting the features at each clustering step means that we are picking features that help differentiate cell-types within a relevant subpopulation, rather than globally. We can hope that looking at the genes that help subdivide cell clusters can give us a hint at their function.

### Implimentation options

The original HIPPO implimentation is quite slow. The main pitfall is that the first few clustering steps are super expensive. This is because (1) you have the most cells and (2) virtually all the genes get selected as features since the EBs are a super heterogeneos cell population. 

I made extensive modification here (capping the number of features per custering step, fast sparse svd rather than coercing to dense, allowing any gene to be selected at each clustering step rather than a strict subset of the genes in the prior step)
https://github.com/karltayeb/HIPPO

and Mengjie has incorporated a lot of my suggestions into a lightweight version of HIPPO here:
https://github.com/ChenMengjie/lightHippo

Moving forward, I suggest working with her version. It seems a lot cleaner. My version I basically forked the repository, removed a bunch of stuff that was not relevant to our particular application, and then made various optimizations/changes (detailed in the README). But it is not super well tested, and quite messy. In contrast, Mengjie's new implimentation is super easy to read, and it looks like she wrote some useful functions for querying the resulting clustering. 

### What's done


