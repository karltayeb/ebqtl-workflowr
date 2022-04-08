# CelliD
rule abstract_run_cellid:
    input:
        'output/proc/EB_high_pass_filtered_normalized_subset.h5ad',
        'output/clustering/hippo_clust.txt'
    output:
        directory('output/cellid/EB_high_pass_filtered_normalized_subset'),
        'output/cellid/EB_high_pass_filtered_normalized_subset_cellid.txt'
    params:
        group_col='all',
        use_hdf5='FALSE',
nmcs=50
    script: '../scripts/cellid_fit.R'

# rule abstract_run_cellid as run_cellid with:
rule run_cellid_by_donor:
    input: 'output/proc/{adata}.h5ad'
    output: 'output/cellid/{adata}/{adata}.{group_col}.cellid.rds'
    params:
        group_col= lambda wildcards: wildcards.group_col,
        use_hdf5='TRUE',
        nmcs=50
    script: '../scripts/cellid_fit.R'

use rule abstract_run_cellid as cao_run_cellid with:
    input:
        'data/CaoEtAl.Obj.CellsOfAllClusters.ProteinCodingGenes.h5ad'
    output:
        'output/cellid/fetal_reference.cellid.rds'


## take 2
# CelliD, run cellid seperately for cells grouped by `cluster`
# since there are varialbe number of outputs, saves
rule run_cellid_per_cluster:
    input:
        adata='data/single_cell/{adata}.h5ad',
        cluster='output/clustering/{cluster}.tsv'
    output:
        manifest='output/cellid/{adata}/{cluster}/manifest.txt'
    params:
        prefix=lambda wildcards, output: '/'.join(output[0].split('/')[:-1]),
        use_hdf5='TRUE',
        nmcs=50
    script: '../scripts/cellid_fit2.R'


rule cellid_assign:
    input:
        'output/cellid/{adata}/fits/{adata}.{group_col}.cellid{sample}.rds',
    output:
        'output/cellid/{adata}/assignments/rds/{adata}.{group_col}.cellid{sample}.assignment.rds',
        'output/cellid/{adata}/assignments/txt/{adata}.{group_col}.cellid{sample}.assignment.txt'
    script:
        '../scripts/cellid_assign.R'


def manifest2cellid(wildcards):
  manifest = f'output/cellid/{wildcards.adata}/{wildcards.cluster}/manifest.txt'
  manifest = pd.read_csv(manifest, sep='\t')
  
  
#rule cellid_assign_term:
#    input:
#      manifest = 'output/cellid/{adata}/{cluster}/manifest.txt'
#    input:
#        expand('output/cellid/EB_high_pass_filtered_normalized/assignments/txt/EB_high_pass_filtered_normalized.donor_id.cellid{sample}.assignment.txt', sample=samples)

