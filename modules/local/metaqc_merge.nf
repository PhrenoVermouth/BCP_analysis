// modules/local/metaqc_merge.nf

process METAQC_MERGE {
    tag "$meta.id"
    publishDir "$params.outdir/metaqc/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(partial_table), path(rho_file)

    output:
    tuple val(meta), path("${meta.id}_total_metaqc.tsv"), emit: metaqc_table

    script:
    """
    python - <<'PY'
    import pandas as pd

    partial = pd.read_csv("${partial_table}", sep=r"\\s+", comment="#")
    rho_df = pd.read_csv("${rho_file}", sep=r"\\s+")

    rho_value = rho_df['Rho'].iloc[0] if not rho_df.empty and 'Rho' in rho_df else ''
    partial['Rho'] = rho_value

    partial.to_csv(f"${meta.id}_total_metaqc.tsv", sep=' ', index=False)
    PY
    """
}
