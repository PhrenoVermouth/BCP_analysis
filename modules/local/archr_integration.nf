// modules/local/archr_integration.nf
process ARCHR_INTEGRATION {

    tag "$meta.id"
    errorStrategy 'retry'
    maxRetries 1
    // container 'multiome_env.sif' // Left commented for user to uncomment when container is built
    
    // Publish ArchR outputs
    publishDir "$params.outdir/multiome/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(rna_raw_dir), path(atac_frag), path(atac_frag_index), path(whitelist)

    output:
    path "${meta.id}_ArchR_Project", emit: archr_project
    path "*_mqc.png", emit: mqc_plots, optional: true
    path "*_mqc.tsv", emit: mqc_table, optional: true
    path "*.{h5ad,rds}", emit: h5ad_matrices, optional: true
    path "*.csv", emit: peak_csv

    script:
    def wl_arg = whitelist.name != 'NO_FILE' ? "--barcode_translation ${whitelist}" : "--barcode_translation NULL"
    def sam_py = params.sam_python ?: 'python'
    """
    # Cache bust: v10 - collapse duplicate gene symbols for SoupX
    Rscript $projectDir/bin/run_archr_multiome.R \\
        --sample_id ${meta.id} \\
        --rna_matrix_dir ${rna_raw_dir} \\
        --atac_fragments ${atac_frag} \\
        --genome ${params.atac_genome} \\
        --tss_enrichment ${params.atac_tss_enrichment} \\
        --min_frags ${params.atac_min_frags} \\
        --rna_min_umi ${params.rna_min_umi} \\
        ${wl_arg} \\
        --threads ${task.cpus} \\
        --sam_python ${sam_py} \\
        --project_dir ${projectDir} \\
        --outdir .
    """
}
