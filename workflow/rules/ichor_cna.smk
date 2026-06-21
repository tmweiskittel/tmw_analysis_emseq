rule download_ichor_bam:
    output:
        bam=temp(f"{D_WORK}/ichor/bam/{{sample}}.aligned.sorted.filt.bl.bam"),
        bai=temp(f"{D_WORK}/ichor/bam/{{sample}}.aligned.sorted.filt.bl.bam.bai")
    params:
        bam_gcs=lambda wc: bam_gcs_for_sample(wc.sample),
        bai_gcs=lambda wc: bai_gcs_for_sample(wc.sample)
    log:
        f"{D_LOGS}/download_ichor_bam/{{sample}}.log"
    shell:
        """
        mkdir -p "$(dirname "{output.bam}")" "$(dirname "{log}")"
        exec &>> "{log}"

        if [ ! -s "{output.bam}" ]; then
            gcloud storage cp "{params.bam_gcs}" "{output.bam}"
        fi

        if [ ! -s "{output.bai}" ]; then
            gcloud storage cp "{params.bai_gcs}" "{output.bai}"
        fi
        """


rule ichor_readcounter:
    input:
        bam=f"{D_WORK}/ichor/bam/{{sample}}.aligned.sorted.filt.bl.bam",
        bai=f"{D_WORK}/ichor/bam/{{sample}}.aligned.sorted.filt.bl.bam.bai"
    output:
        wig=temp(f"{D_WORK}/ichor/wig/{{sample}}.wig")
    params:
        window=config.get("ichor", {}).get("window", 1000000),
        quality=config.get("ichor", {}).get("quality", 20)
    conda:
        "../envs/ichor.yaml"
    log:
        f"{D_LOGS}/ichor_readcounter/{{sample}}.log"
    shell:
        """
        mkdir -p "$(dirname "{output.wig}")" "$(dirname "{log}")"
        exec &>> "{log}"

        readCounter \
          --window {params.window} \
          --quality {params.quality} \
          "{input.bam}" > "{output.wig}"
        """


rule ichor_run_sample:
    input:
        wig=f"{D_WORK}/ichor/wig/{{sample}}.wig",
        gc_wig=config["ichor"]["gc_wig"],
        map_wig=config["ichor"]["map_wig"],
        normal_panel=config["ichor"]["normal_panel"]
    output:
        params=temp(f"{D_OUT}/ichor/{{sample}}/{{sample}}.params.txt"),
        seg=temp(f"{D_OUT}/ichor/{{sample}}/{{sample}}.seg"),
        cna=temp(f"{D_OUT}/ichor/{{sample}}/{{sample}}.cna.seg"),
        rds=temp(f"{D_OUT}/ichor/{{sample}}/{{sample}}.rds")
    params:
        out_dir=f"{D_OUT}/ichor/{{sample}}",
        sample="{sample}",
        ploidy=config.get("ichor", {}).get("ploidy", "c(2)"),
        normal=config.get("ichor", {}).get("normal", "c(0.5,0.6,0.7,0.8,0.9,0.95)"),
        genome=config.get("ichor", {}).get("genome", "hg38")
    conda:
        "../envs/ichor.yaml"
    log:
        f"{D_LOGS}/ichor_run_sample/{{sample}}.log"
    shell:
        """
        mkdir -p "{params.out_dir}" "$(dirname "{log}")"
        exec &>> "{log}"

        Rscript "$CONDA_PREFIX/bin/runIchorCNA.R" \
          --id "{params.sample}" \
          --WIG "{input.wig}" \
          --gcWig "{input.gc_wig}" \
          --mapWig "{input.map_wig}" \
          --normalPanel "{input.normal_panel}" \
          --ploidy "{params.ploidy}" \
          --normal "{params.normal}" \
          --genomeBuild "{params.genome}" \
          --outDir "{params.out_dir}"
        """


rule ichor_aggregate_contrast:
    input:
        params=lambda wc: expand(
            f"{D_OUT}/ichor/{{sample}}/{{sample}}.params.txt",
            sample=samples_for_contrast(wc.experiment)
        )
    output:
        summary=temp(f"{D_OUT}/ichor/{{experiment}}/ichor_summary.tsv")
    params:
        script=f"{R_EMSEQ}/scripts/aggregate_ichor_results.R"
    conda:
        "../envs/qc_summary.yaml"
    log:
        f"{D_LOGS}/ichor_aggregate/{{experiment}}.log"
    shell:
        """
        mkdir -p "$(dirname "{output.summary}")" "$(dirname "{log}")"
        exec &>> "{log}"

        Rscript "{params.script}" \
          --params-files "{input.params}" \
          --out "{output.summary}"
        """


rule upload_ichor_results:
    input:
        summary=f"{D_OUT}/ichor/{{experiment}}/ichor_summary.tsv"
    output:
        done=f"{D_OUT}/upload/{{experiment}}.ichor_cna.upload.done"
    params:
        bucket=config["meta"]["results_bucket"],
        prefix=config["meta"]["results_prefix"]
    log:
        f"{D_LOGS}/upload_ichor/{{experiment}}.log"
    shell:
        """
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.done}")"
        exec &>> "{log}"

        DEST="gs://{params.bucket}/{params.prefix}/ichor_cna/{wildcards.experiment}"

        gcloud storage cp "{input.summary}" "$DEST/"

        touch "{output.done}"
        """
