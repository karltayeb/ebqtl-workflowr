import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()
ZHANG_FILES, = glob_wildcards("data/gwas/zhang/sumstats/{fname}.sumstats")
TRAITS=['TST', 'SHBG']
#CELLTYPES=["IPSC", "MESODERM", "ENDODERM", "EPITHELIAL", "ENDOTHELIAL", "EARLYECTO", "NEUR1", "NEUR2", "NEUR3"]
CELLTYPES=["EarlyEctoderm_neuroepithelium","endoderm", "MesodermEndothelialHematopoetic", "NeuralCrest", "PluripotentUndifferentiated", "PostMitoticNeurons", "tentative_developingEye", "tentative_DevelopingSensoryGanglia", "tentative_neuralTube"]

#rule download_gencode_data:
#    input:
#        HTTP.remote("ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.basic.annotation.gff3.gz"", keep_local=True)
#    output:
#        "data/gencode/gencode.v39.basic.annotation.gff3"
#    run:
#      shell("wget -O data/gencode/gencode.v39.basic.annotation.gff3.gz {input}")
#      shell("gzip -d data/gencode/gencode.v39.basic.annotation.gff3.gz")
#
#rule process_gencode_data:
#    input:
#      "data/gencode/gencode.v39.basic.annotation.gff3"
#    output:
#      "data/gencode/gencode.v39.basic.annotation.filtered.gff3",
#      "data/gencode/gencode.v39.basic.annotation.filtered.bed",
#      "data/gencode/gencode.v39.basic.annotation.filtered.tss.bed"
#    script:
#      "code/mashr/process_gff.R"

rule process_gtf:
    input:
        gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        gtf_loc="data/gencode/gencode.hg38.filtered.gtf",
        tss_loc="data/gencode/gencode.hg38.filtered.tss.tsv",
        bed_loc="data/gencode/gencode.hg38.filtered.tss.bed"
    script:
        "code/mashr/gene_locs.R"

rule list_snps_zhang:
    input:
        expand("data/gwas/zhang/sumstats/{fname}.sumstats", fname=ZHANG_FILES)
    output:
        "data/gwas/{study}/all_snps.txt"
    run:
        shell("cut -f 1 {input}*.sumstats | sort | uniq | sed -e s/SNP//g > {output}")

rule download_fca_sc:
    input:
        counts=HTTP.remote("atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/gene_count_sampled.RDS", keep_local=True),
        cell_metadata=HTTP.remote("atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/df_cell.RDS", keep_local=True),
        gene_metadata=HTTP.remote("atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/df_gene.RDS", keep_local=True)
    output:
        counts="data/fca/counts.sampled.rds",
        cell_metadata="data/fca/cell_metadata.rds",
        gene_metadata="data/fca/gene_metadata.rds"
    run:
        shell("wget -O {output.counts} {input.counts}")
        shell("wget -O {output.cell_metadata} {input.cell_metadata}")
        shell("wget -O {output.gene_metadata} {input.gene_metadata}")

rule download_fca_celltype:
    input:
        expression=HTTP.remote("https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/FCA_RNA_supp_files/gene_expression_celltype.txt"),
        de_genes=HTTP.remote("https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/DE_gene_77_main_cell_type.csv")
    output:
        expression="data/fca/expression.celltype.csv",
        de_genes="data/fca/de_genes.csv"
    run:
        shell("wget -O {output.expression} {input.expression}")
        shell("wget -O {output.de_genes} {input.de_genes}")

rule create_h5ad_fca:
    input:
        "data/fca/counts.sampled.rds",
        "data/fca/cell_metadata.rds",
        "data/fca/gene_metadata.rds"
    output:
        "data/single_cell_objects/fca.sampled.h5ad",
        "data/fca/organ_celltype.tsv"
    script:
        "code/fca/create_h5ad_fca.R"

rule scdrs_score:
    input:
        "data/single_cell_objects/Lowpass.3seqbatches.merged.scvi_processed_and_scaled.h5ad",
        "data/scDRS/gs_files/{geneset}.gs",
        "data/scDRS/gs_files/{geneset}.map"
    output:
        "results/scDRS/{geneset}/{trait}/scores.tsv",
        "results/scDRS/{geneset}/{trait}/scores.ctrl.tsv"
    script:
        "code/scdrs/scDRS.py"

rule mashr_pseudobulk:
    input:
        anndata="data/single_cell_objects/{annotation}.csc.h5ad"
    output:
        expand("data/static/{{aggregation}}/type/{{annotation}}/{type}/{out}.tsv", type=CELLTYPES, out=['expression', 'covariates', 'individuals']),
        sample_summary="data/static/{aggregation}/type/{annotation}/sample_summary.tsv",
        celltype_summary="data/static/{aggregation}/type/{annotation}/samples_per_celltype.tsv"
    conda: "slurmy/r-pseudobulk.yml"
    script:
        "code/mashr/{wildcards.aggregation}_static.R"

rule mashr_pseudobulk2:
    input:
        anndata="output/proc/{anndata}.h5ad",
        annotation="output/clustering/{annotation}.txt"
    output:
        manifest="data/static/pseudobulk/{anndata}/{annotation}/manifest.txt",
        sample_summary="data/static/pseudobulk/{anndata}/{annotation}/sample_summary.tsv",
        celltype_summary="data/static/{pseudobulk}/{annotation}/samples_per_celltype.tsv"
    conda: "slurmy/r-pseudobulk.yml"
    script:
        "code/mashr/pseudobulk_static.R"

rule mashr_genotype_filter:
    input:
	      genotypes="data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz",
	      inds="data/static/pseudobulk/type/{annotation}/{type}/individuals.tsv"
    output:
    	  "data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.recode.vcf"
    shell:
	      "code/mashr/genotype_filter.sh {input.genotypes} {wildcards.annotation} {wildcards.type} {input.inds}"

rule mashr_genotype_012:
    input:
	      "data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.recode.vcf"
    output:
    	  expand("data/static/{{aggregation}}/type/{{annotation}}/{{type}}/genotypes_filtered.{out}", out=['012', '012.indv', '012.pos'])
    shell:
	      "code/mashr/genotype_012.sh {input} {wildcards.annotation} {wildcards.type}"

rule mashr_genotype_transpose:
    resources:
        mem_mb=50000
    input:
	      "data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.012"
    output:
	      "data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.012.transpose"
    shell:
	      "code/mashr/genotype_transpose.sh {input} {output}"

rule mashr_genotype_reformat:
    resources:
        mem_mb=50000
    input:
        genotypes="data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.012.transpose",
        individuals="data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.012.indv",
        snp_locs="data/static/{aggregation}/type/{annotation}/{type}/genotypes_filtered.012.pos"
    output:
        snp_locs="data/static/{aggregation}/type/{annotation}/{type}/snp_locs.tsv",
        genotypes="data/static/{aggregation}/type/{annotation}/{type}/genotypes.tsv"
    shell:
        "code/mashr/genotype_reformat.sh {input.genotypes} {input.individuals} {input.snp_locs} {output.snp_locs} {output.genotypes}" 

rule matrix_eqtl:
    resources:
        mem_mb=75000,
        time="00:30:00"
    input:
        genotypes="data/static/pseudobulk/type/{annotation}/{type}/genotypes.tsv",
        snp_locs="data/static/pseudobulk/type/{annotation}/{type}/snp_locs.tsv",
        expression="data/static/pseudobulk/type/{annotation}/{type}/expression.tsv",
        gene_locs="data/gencode/gencode.hg38.filtered.tss.tsv",
        covariates="data/static/pseudobulk/type/{annotation}/{type}/covariates.tsv"
    output:
        eqtls="results/mashr/{aggregation}/type/{annotation}/eqtls/{type}/eqtls.tsv",
        df="results/mashr/{aggregation}/type/{annotation}/eqtls/{type}/df.tsv"
    conda: "slurmy/r-matrixEQTL.yml"
    script:
        "code/mashr/matrixEQTL.R"

rule mtc_static:
    resources:
        mem_mb=75000,
        time="00:15:00"
    input:
        eqtls="results/mashr/{aggregation}/type/{annotation}/eqtls/{type}/eqtls.tsv",
        df="results/mashr/{aggregation}/type/{annotation}/eqtls/{type}/df.tsv"
    output:
        all_tests="results/mashr/{aggregation}/type/{annotation}/eqtls/{type}/eqtls.mtc.tsv",
        top_tests="results/mashr/{aggregation}/type/{annotation}/eqtls/{type}/eqtls.tophits.tsv"
    conda: "slurmy/r-mashr.yml"
    script:
        "code/mashr/mtc.R"

rule mashr:
    resources:
        mem_mb=75000,
        time="01:00:00"
    input:
        qtls=expand("results/mashr/{{aggregation}}/type/{{annotation}}/eqtls/{type}/eqtls.mtc.tsv", type=CELLTYPES)
    output:
        trained="results/mashr/{aggregation}/type/{annotation}/eqtls/mashr.trained.rds",
        tophits="results/mashr/{aggregation}/type/{annotation}/eqtls/mashr.tophits.rds"
    conda: "slurmy/r-mashr.yml"
    script:
        "code/mashr/mashr.R"
        
rule mashr_no_endothelial:
    resources:
        mem_mb=75000,
        time="01:00:00"
    input:
        qtls=expand("results/mashr/{{aggregation}}/type/{{annotation}}/eqtls/{type}/eqtls.mtc.tsv", type=[t for t in CELLTYPES if t != "ENDOTHELIAL"])
    output:
        trained="results/mashr/{aggregation}/type/{annotation}/eqtls/mashr.trained.noendo.rds",
        tophits="results/mashr/{aggregation}/type/{annotation}/eqtls/mashr.tophits.noendo.rds"
    conda: "slurmy/r-mashr.yml"
    script:
        "code/mashr/mashr.R"
        
###
#rule make_conda: 
#    input:
#        "test.tmp"
#    output:
#        "test.tmp2"
#    conda: "slurmy/r-pseudobulk.yml"
#    shell:
#        "echo booyah > {output}"
#        
#rule download_gs_files:
#    input:
#        HTTP.remote("ndownloader.figshare.com/files/30853708", keep_local=True)
#    output:
#        gsfile="data/scDRS/gs_files/gs_file.zip"
#    run:
#        shell("wget -O data/scDRS/gs_files/gs_file.zip {input}")
###
