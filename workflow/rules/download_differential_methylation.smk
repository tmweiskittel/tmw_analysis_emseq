rule download_methylkit_raw:
    output:
        db=f"{D_WORK}/methylkit_raw/{{sample}}.CpG.methylKit.gz"
    params:
        db_gcs=lambda wc: (
            f"{gcs_sample_dir(wc.sample)}/"
            f"{wc.sample}.CpG.methylKit.gz"
        )
    log:
        cmd=f"{D_LOGS}/download_methylkit_raw/{{sample}}.log"
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

rule make_single_methylkit_tabix_db:
    input:
        amp=f"{D_WORK}/methylkit_raw/{{sample}}.CpG.methylKit.gz"
    output:
        bgz=(
            f"{D_WORK}/methylkit_db/"
            f"{{sample}}.{config['emseq_ref_name']}."
            f"{config['align_method']}.methyldackel.txt.bgz"
        ),
        tbi=(
            f"{D_WORK}/methylkit_db/"
            f"{{sample}}.{config['emseq_ref_name']}."
            f"{config['align_method']}.methyldackel.txt.bgz.tbi"
        )
     conda:
        ENV_METHYLKIT
    log:
        cmd=f"{D_LOGS}/make_single_methylkit_tabix_db/{{sample}}.log"
    params:
        script=f"{R_EMSEQ}/scripts/make_single_amp_methylkit_obj.R",
        library_id=(
            f"{{sample}}.{config['emseq_ref_name']}."
            f"{config['align_method']}.methyldackel"
        ),
        mincov=lambda wc: config.get("mincov", 5),
        build=config["emseq_ref_name"],
        treatment=1,
        out_dir=f"{D_WORK}/methylkit_db"
    shell:
        """
        mkdir -p "{params.out_dir}" "$(dirname "{log.cmd}")"
        exec &>> "{log.cmd}"

        Rscript "{params.script}" \
          --amp_file "{input.amp}" \
          --library_id "{params.library_id}" \
          --mincov {params.mincov} \
          --out_dir "{params.out_dir}" \
          --treatment {params.treatment} \
          --build "{params.build}"
        """
