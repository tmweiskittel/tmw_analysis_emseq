rule download_methylkit_db:
    message:
        "Download per-sample MethylDackel methylKit tabix database"
    output:
        db=lambda wc: local_methylkit_db_for_sample(wc.sample),
        tbi=lambda wc: local_methylkit_tbi_for_sample(wc.sample)
    params:
        db_gcs=lambda wc: methylkit_db_gcs_for_sample(wc.sample),
        tbi_gcs=lambda wc: methylkit_tbi_gcs_for_sample(wc.sample)
    log:
        cmd=f"{D_LOGS}/download_methylkit_db/{{sample}}.log"
    shell:
        """
        mkdir -p "$(dirname "{output.db}")" "$(dirname "{log.cmd}")"
        exec &>> "{log.cmd}"

        echo "[download-methylkit-db] $(date) sample={wildcards.sample}"

        if [ ! -s "{output.db}" ]; then
            gsutil cp "{params.db_gcs}" "{output.db}"
        else
            echo "Already exists: {output.db}"
        fi

        if [ ! -s "{output.tbi}" ]; then
            gsutil cp "{params.tbi_gcs}" "{output.tbi}"
        else
            echo "Already exists: {output.tbi}"
        fi
        """
