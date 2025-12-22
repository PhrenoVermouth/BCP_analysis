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

partial = pd.read_csv("${partial_table}", sep='\t', comment='#')
rho_df = pd.read_csv("${rho_file}", sep='\t')

rho_value = rho_df['Rho'].iloc[0] if not rho_df.empty and 'Rho' in rho_df else ''
partial['Rho'] = rho_value

out_file = f"${meta.id}_total_metaqc.tsv"
header_lines = [
    "# plot_type: 'table'\n",
    "# section_name: 'Total_metaQC'\n",
    "# description: 'Combined QC metrics including ambient RNA correction'\n",
    "# pconfig:\n",
    "#     sortRows: false\n",
]
with open(out_file, 'w') as fh:
    fh.writelines(header_lines)
    partial.to_csv(fh, sep='\t', index=False)
PY
    """
}
