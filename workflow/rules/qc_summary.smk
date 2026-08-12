rule download_qc_summary:
    output:
        tsv=f"{D_WORK}/qc/{{sample}}.qc_summary.tsv"
    params:
        gcs=lambda wc: qc_summary_gcs_for_sample(wc.sample)
    log:
        f"{D_LOGS}/download_qc_summary/{{sample}}.log"
    shell:
        """
        mkdir -p "$(dirname "{output.tsv}")" "$(dirname "{log}")"
        exec &>> "{log}"

        if [ ! -s "{output.tsv}" ]; then
            gcloud storage cp "{params.gcs}" "{output.tsv}"
        else
            echo "Already exists: {output.tsv}"
        fi
        """


rule aggregate_qc_summary:
    input:
        qc=lambda wc: [
            local_qc_summary_for_sample(sample)
            for sample in samples_for_contrast(wc.experiment)
        ]
    output:
        aggregate=f"{D_OUT}/qc/{{experiment}}/qc_aggregate.tsv",
        stats=f"{D_OUT}/qc/{{experiment}}/qc_group_tests.tsv"
    params:
        samples=lambda wc: " ".join(samples_for_contrast(wc.experiment)),
        groups=lambda wc: " ".join(
            SAMPLE_CONTRASTS[sample][wc.experiment]
            for sample in samples_for_contrast(wc.experiment)
        ),
        script=f"{R_EMSEQ}/scripts/aggregate_qc_summary.R"
    conda:
        ENV_METHYLKIT
    log:
        f"{D_LOGS}/aggregate_qc_summary/{{experiment}}.log"
    shell:
        """
        mkdir -p "$(dirname "{output.aggregate}")" "$(dirname "{log}")"
        exec &>> "{log}"

        Rscript "{params.script}" \
          --qc-files "{input.qc}" \
          --samples "{params.samples}" \
          --groups "{params.groups}" \
          --out-aggregate "{output.aggregate}" \
          --out-stats "{output.stats}"
        """



        touch "{output.done}"
        """
