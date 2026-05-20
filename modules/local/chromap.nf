// modules/local/chromap.nf
process CHROMAP {

    tag "$meta.id"
    errorStrategy 'retry'
    maxRetries 2
    maxForks 1

    // We only publish the final fragments to avoid moving huge intermediate files
    publishDir "$params.outdir/chromap/${meta.id}", \
        mode: 'copy', \
        pattern: "*.{tsv.gz,tbi,log}"

    input:
    // r1s, r2s (barcode), r3s (paired end) — lists for multi-lane support
    tuple val(meta), path(r1s), path(r2s), path(r3s), path(ref_fasta), path(index), path(whitelist)

    output:
    tuple val(meta), path("${meta.id}.fragments.tsv.gz"), emit: fragments
    tuple val(meta), path("${meta.id}.fragments.tsv.gz.tbi"), emit: fragments_tbi
    path "${meta.id}.chromap.log", emit: log
    path "${meta.id}_chromap_mqc.json", emit: mqc_json

    script:
    // chromap accepts comma-separated file lists for multi-lane input
    def r1_arg = (r1s instanceof List ? r1s.join(',') : r1s)
    def r2_arg = (r2s instanceof List ? r2s.join(',') : r2s)
    def r3_arg = (r3s instanceof List ? r3s.join(',') : r3s)
    def has_whitelist = whitelist.name != 'NO_FILE'
    """
    # 1. Prepare whitelist and read-format args
    READ_FMT_ARG=""
    WHITELIST_ARG=""

    if [ "${has_whitelist}" = "true" ]; then
        FIRST_R2=\$(echo "${r2_arg}" | cut -d',' -f1)
        WL_LEN=\$(head -1 ${whitelist} | tr -d '\\n' | wc -c)
        BC_LEN=\$(zcat "\$FIRST_R2" | head -2 | tail -1 | tr -d '\\n' | wc -c)
        WHITELIST_FILE="${whitelist}"

        if [ "\$BC_LEN" -gt "\$WL_LEN" ]; then
            BC_START=\$((BC_LEN - WL_LEN))
            BC_END=\$((BC_LEN - 1))
            echo "Chromium X detected: barcode \${BC_LEN}bp > whitelist \${WL_LEN}bp. Using bc:\${BC_START}:\${BC_END} + revcomp whitelist."
            READ_FMT_ARG="--read-format bc:\${BC_START}:\${BC_END},r1:0:-1,r2:0:-1"
            # Generate reverse-complement whitelist
            awk '{s=""; for(i=length(\$1);i>0;i--){c=substr(\$1,i,1); if(c=="A")s=s"T"; else if(c=="T")s=s"A"; else if(c=="C")s=s"G"; else if(c=="G")s=s"C"; else s=s"N"} print s}' ${whitelist} > whitelist_rc.txt
            WHITELIST_FILE="whitelist_rc.txt"
            echo "Generated reverse-complement whitelist (\$(wc -l < whitelist_rc.txt) barcodes)."
        fi
        WHITELIST_ARG="--barcode-whitelist \$WHITELIST_FILE"
    fi

    # 2. Run chromap
    chromap --preset atac \\
            --drop-repetitive-reads 10 \\
            -q 0 \\
            --trim-adapters \\
            -x ${index} \\
            -r ${ref_fasta} \\
            -t ${task.cpus} \\
            -1 ${r1_arg} \\
            -b ${r2_arg} \\
            -2 ${r3_arg} \\
            -o ${meta.id}.fragments.tsv \\
            \$READ_FMT_ARG \\
            \$WHITELIST_ARG 2> ${meta.id}.chromap.log

    # 3. Sort fragments (required by ArchR and tabix)
    sort -k1,1 -k2,2n ${meta.id}.fragments.tsv > ${meta.id}.fragments.sorted.tsv

    # 4. Compress and index
    bgzip -@ ${task.cpus} ${meta.id}.fragments.sorted.tsv
    tabix -p bed ${meta.id}.fragments.sorted.tsv.gz

    # 5. Rename to standard output and cleanup
    mv ${meta.id}.fragments.sorted.tsv.gz ${meta.id}.fragments.tsv.gz
    mv ${meta.id}.fragments.sorted.tsv.gz.tbi ${meta.id}.fragments.tsv.gz.tbi
    rm ${meta.id}.fragments.tsv

    # 6. Generate MultiQC metrics JSON
    python3 $projectDir/bin/chromap_to_mqc.py ${meta.id}.chromap.log
    """
}
