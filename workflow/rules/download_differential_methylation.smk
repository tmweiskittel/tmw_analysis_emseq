rule download_methylkit_db:
    output:
        db=(
            f"{D_WORK}/methylkit_db/"
            f"{{sample}}.{config['emseq_ref_name']}."
            f"{config['align_method']}.methyldackel.txt.bgz"
        )
    params:
        db_gcs=lambda wc: methylkit_db_gcs_for_sample(wc.sample)
    log:
        cmd=f"{D_LOGS}/download_methylkit_db/{{sample}}.log"
    shell:
        """
        mkdir -p "$(dirname "{output.db}")" "$(dirname "{log.cmd}")"
        exec &>> "{log.cmd}"

        if [ ! -s "{output.db}" ]; then
            gsutil cp "{params.db_gcs}" "{output.db}"
        else
            echo "Already exists: {output.db}"
        fi
        """
