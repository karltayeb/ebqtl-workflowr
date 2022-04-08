import numpy as np
import pandas as pd

def get_prefix(path, c=1):
    return '/'.join(path.split('/')[:-c])

rule process_gtf:
    input:
        gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf",
        tss_loc="data/gencode/gencode.hg38.filtered.tss.tsv",
        bed_loc="data/gencode/gencode.hg38.filtered.tss.bed"
    script:
        "../scripts/staticqtl/gene_locs.R"


rule pseudobulk_adata:
    input:
        adata_path='output/proc/{adata}.h5ad',
        cluster_path='output/clustering/{annotation}.txt'
    output:
        'data/static/pseudobulk/{adata}/{annotation}/manifest.txt'
    params:
        prefix = lambda wildcards, output: get_prefix(output[0]),
        groups=['pb_cluster', 'donor_id']
    conda: '../envs/scvi-scanpy-gpu.yml'
    script:
        '../scripts/pseudobulk.py'

rule genotype_filter:
    input:
	      genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
	      inds="data/static/pseudobulk/{adata}/{annotation}/{type}/individuals.tsv"
    output:
    	  "data/static/pseudobulk/{adata}/{annotation}/{type}/genotypes_filtered.recode.vcf"
    params:
        prefix = lambda wildcards, output: get_prefix(output[0])
    shell:
	      "workflow/scripts/staticqtl/genotype_filter.sh {input.genotypes} {wildcards.annotation} {wildcards.type} {input.inds} {params.prefix}"


rule genotype_012:
    input:
	      "data/static/pseudobulk/{adata}/{annotation}/{type}/genotypes_filtered.recode.vcf"
    output:
    	  expand("data/static/pseudobulk/{{adata}}/{{annotation}}/{{type}}/genotypes_filtered.{out}", out=['012', '012.indv', '012.pos'])
    params:
        prefix = lambda wildcards, output: get_prefix(output[0])
    shell:
	      "workflow/scripts/staticqtl/genotype_012.sh {input} {wildcards.annotation} {wildcards.type} {params.prefix}"


rule genotype_transpose:
    resources:
        mem_mb=50000
    input:
        "data/static/pseudobulk/{adata}/{annotation}/{type}/genotypes_filtered.012"
    output:
        "data/static/pseudobulk/{adata}/{annotation}/{type}/genotypes_filtered.012.transpose"
    shell:
	      "workflow/scripts/staticqtl/genotype_transpose.sh {input} {output}"


rule genotype_reformat:
    resources:
        mem_mb=50000
    input:
        genotypes="data/static/pseudobulk/{adata}/{annotation}/{type}/genotypes_filtered.012.transpose",
        individuals="data/static/pseudobulk/{adata}/{annotation}/{type}/genotypes_filtered.012.indv",
        snp_locs="data/static/pseudobulk/{adata}/{annotation}/{type}/genotypes_filtered.012.pos"
    output:
        snp_locs="data/static/pseudobulk/{adata}/{annotation}/{type}/snp_locs.tsv",
        genotypes="data/static/pseudobulk/{adata}/{annotation}/{type}/genotypes.tsv"
    shell:
        "workflow/scripts/staticqtl/genotype_reformat.sh {input.genotypes} {input.individuals} {input.snp_locs} {output.snp_locs} {output.genotypes}" 

rule matrix_eqtl:
    resources:
        mem_mb=75000,
        time="00:30:00"
    input:
        genotypes="data/static/pseudobulk/{adata}/{annotation}/{type}/genotypes.tsv",
        snp_locs="data/static/pseudobulk/{adata}/{annotation}/{type}/snp_locs.tsv",
        expression="data/static/pseudobulk/{adata}/{annotation}/{type}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv",
        covariates="data/static/pseudobulk/{adata}/{annotation}/{type}/covariates.tsv"
    output:
        eqtls="output/static/pseudobulk/{adata}/{annotation}/{type}/eqtls.tsv",
        df="output/static/pseudobulk/{adata}/{annotation}/{type}/df.tsv"
    script:
        "../scripts/staticqtl/matrixEQTL.R"

def manifest2eqtl(wildcards):
    manifest = (
        f'data/static/pseudobulk/{wildcards.adata}/'
        f'{wildcards.annotation}/manifest.txt'
    )
    df = pd.read_csv(manifest, sep='\t')
    pb_clusters = np.unique(df.pb_cluster)
    return {c:(
            f"output/static/pseudobulk/{wildcards.adata}/"
            f"{wildcards.annotation}/{c}/eqtls.tsv"
            ) for c in pb_clusters}
    
rule matrix_eqtl_all:
    input:
        unpack(manifest2eqtl)
    output:
        manifest = "output/static/pseudobulk/{adata}/{annotation}/manifest.txt"
    shell:
        "{input} > {output.manifest}"
