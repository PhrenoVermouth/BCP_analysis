// modules/local/metaqc_merge.nf

process METAQC_MERGE {
    tag "$meta.id"
    publishDir "$params.outdir/metaqc/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(partial_table), path(rho_file)

    output:
    tuple val(meta), path("${meta.id}_total_metaqc_mqc.tsv"), emit: metaqc_table

    script:
    """
    python - <<'PY'
    import pandas as pd

    partial = pd.read_csv("${partial_table}", sep=r"\\s+", comment="#")
    rho_df = pd.read_csv("${rho_file}", sep=r"\\s+")

    rho_value = rho_df['Rho'].iloc[0] if not rho_df.empty and 'Rho' in rho_df else ''
    partial['Rho'] = rho_value

    output_path = f"${meta.id}_total_metaqc_mqc.tsv"

    with open(output_path, "w") as handle:
        handle.write("# plot_type: 'table'\n")
        handle.write("# section_name: 'Total_metaQC'\n")
        handle.write("# description: 'Integrated QC metrics across Scrublet, filtering, and SoupX'\n")
        handle.write("# pconfig:\n")
        handle.write("#     sortRows: false\n")
        partial.to_csv(handle, sep=' ', index=False)
    PY
    """
}
