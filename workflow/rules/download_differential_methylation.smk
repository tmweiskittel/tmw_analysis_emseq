rule download_methylkit_db:
    message:
        "Download methylKit tabix database from Google Cloud Storage"
    output:
        db=f"{D_WORK}/methylkit_db/{{sample}}.methyldackel.txt.bgz",
        tbi=f"{D_WORK}/methylkit_db/{{sample}}.methyldackel.txt.bgz.tbi"
    params:
        db_gcs=lambda wc: methylkit_db_for_sample(wc.sample),
        tbi_gcs=lambda wc: methylkit_db_for_sample(wc.sample) + ".tbi"
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
