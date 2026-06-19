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
        "../envs/qc_summary.yaml"
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


rule upload_qc_summary:
    input:
        aggregate=f"{D_OUT}/qc/{{experiment}}/qc_aggregate.tsv",
        stats=f"{D_OUT}/qc/{{experiment}}/qc_group_tests.tsv"
    output:
        done=f"{D_OUT}/upload/{{experiment}}.qc_summary.upload.done"
    params:
        bucket=config["meta"]["results_bucket"],
        prefix=config["meta"]["results_prefix"]
    log:
        f"{D_LOGS}/upload_qc_summary/{{experiment}}.log"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.done}")"

        DEST="gs://{params.bucket}/{params.prefix}/qc_summary/{wildcards.experiment}"

        gcloud storage cp "{input.aggregate}" "$DEST/" > "{log}" 2>&1
        gcloud storage cp "{input.stats}" "$DEST/" >> "{log}" 2>&1

        touch "{output.done}"
        """
