rule scanpy_pca_umap_leiden:
    input:
        "output/proc/EB_high_pass_filtered.h5ad"
    output:
        "output/proc/EB_high_pass_filtered_proc.h5ad"
    conda:
        "../envs/scvi-scanpy-gpu.yml"
    script:
        "../scripts/scanpy_process.py"

rule normalized_pca_umap_leiden:
    input:
        "output/proc/EB_high_pass_filtered_normalized.h5ad"
    output:
        "output/proc/EB_high_pass_filtered_normalized_proc.h5ad"
    conda:
        "../envs/scvi-scanpy-gpu.yml"
    script:
        "../scripts/scanpy_process.py"

