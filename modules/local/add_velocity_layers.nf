// modules/local/add_velocity_layers.nf

process ADD_VELOCITY_LAYERS {
    tag "$meta.id"
    publishDir "$params.outdir/velocity/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(gzipped_dir), path(counts_h5ad)

    output:
    tuple val(meta), path("${meta.id}_velocity_ready.h5ad"), emit: velocity_h5ad

    script:
    """
    bin/add_velocity_layers.py \
        --counts_h5ad ${counts_h5ad} \
        --velo_dir ${gzipped_dir}/Velocyto/filtered \
        --output ${meta.id}_velocity_ready.h5ad
    """
}