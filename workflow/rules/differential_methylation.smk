# ============================================================================
# differential_methylation.smk
# EM-seq cfDNA downstream differential methylation module
# ============================================================================

rule methylkit_unite:
    message:
        "Unite per-sample methylKit tabix databases into methylBase"
    wildcard_constraints:
        experiment = "[^.]+"
    conda:
        ENV_METHYLKIT
    input:
        mkit_lib_db=lambda wc: [
            local_methylkit_db_for_sample(sample)
            for sample in samples_for_contrast(wc.experiment)
        ]
    log:
        cmd=f"{D_LOGS}/{{experiment}}_methylkit_unite.log"
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_methylkit_unite.tsv"
    params:
        library_id=lambda wc: " ".join(samples_for_contrast(wc.experiment)),
        treatment_list=lambda wc: " ".join(
            map(str, treatment_vector_for_contrast(wc.experiment))
        ),
        script=f"{R_EMSEQ}/scripts/make_methylkit_unite_db.R",
        mincov=lambda wc: contrast_param(wc.experiment, "mincov", 5),
        min_per_group=lambda wc: contrast_param(wc.experiment, "mingroup", 2),
        chunk_size=lambda wc: contrast_param(wc.experiment, "chunksize", 1e6),
        assembly=lambda wc: contrast_param(wc.experiment, "assembly", "hg38"),
    threads:
        32
    output:
        mbase=temp(f"{D_OUT}/dmr/diff/methylBase_{{experiment}}.txt.bgz"),
        mbase_tbi=temp(f"{D_OUT}/dmr/diff/methylBase_{{experiment}}.txt.bgz.tbi"),

    shell:
        """
        exec &>> "{log.cmd}"
        echo "[methylkit-unite] $(date) experiment={wildcards.experiment} threads={threads}"

        rm -f "{output.mbase}"*

        Rscript "{params.script}" \
          --lib_db_list "{input.mkit_lib_db}" \
          --lib_id_list "{params.library_id}" \
          --treatment_list "{params.treatment_list}" \
          --cores {threads} \
          --out_dir "$(dirname "{output.mbase}")" \
          --suffix {wildcards.experiment} \
          --assembly "{params.assembly}" \
          --mincov {params.mincov} \
          --min_per_group {params.min_per_group} \
          --chunk_size {params.chunk_size}
        """


rule methylkit_diff:
    message:
        "Calculate per-base differential methylation using methylKit"
    wildcard_constraints:
        experiment = "[^.]+"
    conda:
        ENV_METHYLKIT
    input:
        mbase=f"{D_OUT}/dmr/diff/methylBase_{{experiment}}.txt.bgz"
    log:
        cmd=f"{D_LOGS}/{{experiment}}_methylkit_diff.log"
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_methylkit_diff.tsv"
    params:
        script=f"{R_EMSEQ}/scripts/make_methylkit_diff_db.R",
        chunk_size=lambda wc: contrast_param(wc.experiment, "chunksize", 1e6),
    threads:
        32
    output:
        mdiff=temp(f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.txt.bgz"),
        mdiff_tbi=temp(f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.txt.bgz.tbi"),

    shell:
        """
        exec &>> "{log.cmd}"
        echo "[methylkit-diff] $(date) experiment={wildcards.experiment} threads={threads}"

        rm -f "{output.mdiff}" "{output.mdiff}.tbi"

        Rscript "{params.script}" \
          --mbase "{input.mbase}" \
          --cores {threads} \
          --out_dir "$(dirname "{output.mdiff}")" \
          --suffix {wildcards.experiment} \
          --chunk_size {params.chunk_size}
        actual_diff=$(ls -t "$(dirname "{output.mdiff}")"/methylDiff_{wildcards.experiment}*.txt.bgz | head -n 1)

        if [ -z "$actual_diff" ]; then
            echo "ERROR: Could not find methylDiff output"
            ls -lh "$(dirname "{output.mdiff}")"
            exit 1
        fi

        mv "$actual_diff" "{output.mdiff}"

        if [ -f "$actual_diff.tbi" ]; then
            mv "$actual_diff.tbi" "{output.mdiff}.tbi"
        fi

        """


rule methylkit_diff_tiled:
    message:
        "Calculate tiled differential methylation using methylKit"
    wildcard_constraints:
        experiment = "[^.]+"
    conda:
        ENV_METHYLKIT
    input:
        mkit_lib_db=lambda wc: [
            local_methylkit_db_for_sample(sample)
            for sample in samples_for_contrast(wc.experiment)
        ]
    log:
        cmd=f"{D_LOGS}/{{experiment}}_methylkit_diff_tiled.log"
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_methylkit_diff_tiled.tsv"
    params:
        library_id=lambda wc: " ".join(samples_for_contrast(wc.experiment)),
        treatment_list=lambda wc: " ".join(
            map(str, treatment_vector_for_contrast(wc.experiment))
        ),
        script=f"{R_EMSEQ}/scripts/make_methylkit_diff_tiled_db.R",
        mincov=lambda wc: contrast_param(wc.experiment, "mincov", 5),
        min_per_group=lambda wc: contrast_param(wc.experiment, "mingroup", 2),
        chunk_size=lambda wc: contrast_param(wc.experiment, "chunksize", 1e6),
        win_size=lambda wc: contrast_param(wc.experiment, "win_size", 1000),
        assembly=lambda wc: contrast_param(wc.experiment, "assembly", "hg38"),
    threads:
        32
    output:
        unite=temp(f"{D_OUT}/dmr/diff/methylBase_{{experiment}}.tiled.txt.bgz"),
        diff=temp(f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.tiled.txt.bgz"),
        unite_tbi=temp(f"{D_OUT}/dmr/diff/methylBase_{{experiment}}.tiled.txt.bgz.tbi"),
        diff_tbi=temp(f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.tiled.txt.bgz.tbi"),
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[methylkit-diff-tiled] $(date) experiment={wildcards.experiment} threads={threads}"

        rm -f "{output.unite}"* "{output.diff}"*

        Rscript "{params.script}" \
          --lib_db_list "{input.mkit_lib_db}" \
          --lib_id_list "{params.library_id}" \
          --treatment_list "{params.treatment_list}" \
          --cores {threads} \
          --out_dir "$(dirname "{output.unite}")" \
          --suffix {wildcards.experiment} \
          --assembly "{params.assembly}" \
          --mincov {params.mincov} \
          --win_size {params.win_size} \
          --min_per_group {params.min_per_group} \
          --chunk_size {params.chunk_size}
        actual_diff=$(ls -t "$(dirname "{output.diff}")"/methylDiff_{wildcards.experiment}.tiled*.txt.bgz | head -n 1)

        if [ -z "$actual_diff" ]; then
            echo "ERROR: Could not find methylDiff tiled output"
            ls -lh "$(dirname "{output.diff}")"
            exit 1
        fi

        mv "$actual_diff" "{output.diff}"

        if [ -f "$actual_diff.tbi" ]; then
            mv "$actual_diff.tbi" "{output.diff}.tbi"
        fi
            
        """


rule methylkit_meth_extract:
    message:
        "Extract percent methylation matrix from methylBase"
    wildcard_constraints:
        experiment = "[^.]+"
    conda:
        ENV_METHYLKIT
    input:
        mbase=f"{D_OUT}/dmr/diff/methylBase_{{experiment}}.txt.bgz"
    log:
        cmd=f"{D_LOGS}/{{experiment}}_methylkit_meth_extract.log"
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_methylkit_meth_extract.tsv"
    params:
        script=f"{R_EMSEQ}/scripts/all_experiment_methylation.R"
    threads:
        1
    resources:
        mem_mb=512000
    output:
        tsv=temp(f"{D_OUT}/dmr/diff/{{experiment}}_pos_meth.tsv")
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[methylkit-meth-extract] $(date) experiment={wildcards.experiment}"

        Rscript "{params.script}" \
          --db_file "{input.mbase}" \
          --out_file "{output.tsv}"
        """


rule methylkit_annotate_cpg:
    message:
        "Annotate methylKit differential methylation results"
    wildcard_constraints:
        experiment = "[^.]+"
    conda:
        ENV_METHYLKIT
    input:
        db=f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.txt.bgz",
        gtf=config["annotation_gtf"]
    log:
        cmd=f"{D_LOGS}/{{experiment}}_methylkit_annotate_cpg.log"
    benchmark:
        f"{D_BENCHMARK}/{{experiment}}_methylkit_annotate_cpg.tsv"
    params:
        script=f"{R_EMSEQ}/scripts/emseq_annotate_methylkit.R"
    threads:
        1
    output:
        tsv=temp(f"{D_OUT}/dmr/annotation/{{experiment}}_annotated.tsv")
    shell:
        """
        exec &>> "{log.cmd}"
        mkdir -p "$(dirname "{output.tsv}")"

        Rscript "{params.script}" \
          --db "{input.db}" \
          --gtf "{input.gtf}" \
          --out "{output.tsv}"
        """

rule visualize_differential_methylation:
    input:
        mdiff=f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.txt.bgz",
        mdiff_tbi=f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.txt.bgz.tbi",
        tiled_mdiff=f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.tiled.txt.bgz",
        tiled_mdiff_tbi=f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.tiled.txt.bgz.tbi",
        matrix=f"{D_OUT}/dmr/diff/{{experiment}}_pos_meth.tsv",
        annotation=f"{D_OUT}/dmr/annotation/{{experiment}}_annotated.tsv"
    output:
        done=f"{D_OUT}/dmr/visualization/{{experiment}}/{{experiment}}.visualization.done",
        summary=f"{D_OUT}/dmr/visualization/{{experiment}}/summary_statistics.csv",
        sig=f"{D_OUT}/dmr/visualization/{{experiment}}/significant_DMCs.tsv",
        hyper=f"{D_OUT}/dmr/visualization/{{experiment}}/hypermethylated_DMCs.tsv",
        hypo=f"{D_OUT}/dmr/visualization/{{experiment}}/hypomethylated_DMCs.tsv",
        meth_diff_hist=f"{D_OUT}/dmr/visualization/{{experiment}}/{{experiment}}_methylation_difference_histogram.png",
        volcano=f"{D_OUT}/dmr/visualization/{{experiment}}/{{experiment}}_volcano_plot.png",
        manhattan=f"{D_OUT}/dmr/visualization/{{experiment}}/{{experiment}}_manhattan_plot.png",
        tiled_meth_diff=f"{D_OUT}/dmr/visualization/{{experiment}}/{{experiment}}_tiled_methylation_difference.png",
    params:
        outdir=lambda wildcards: f"{D_OUT}/dmr/visualization/{wildcards.experiment}",
        benchdir=f"{D_BENCHMARK}/dmr/visualization",
        script=f"{R_EMSEQ}/scripts/figures.R",
        qvalue_cutoff=config.get("dmr_visualization", {}).get("qvalue_cutoff", 0.05),
        meth_diff_cutoff=config.get("dmr_visualization", {}).get("meth_diff_cutoff", 10),
        top_n_heatmap=config.get("dmr_visualization", {}).get("top_n_heatmap", 500)
    log:
        f"{D_LOGS}/dmr/visualization/{{experiment}}.log"
    benchmark:
        f"{D_BENCHMARK}/dmr/visualization/{{experiment}}.tsv"
    conda:
        ENV_METHYLKIT
    threads:
        4
    shell:
        r"""
        mkdir -p "{params.outdir}"
        mkdir -p "$(dirname "{log}")"
        mkdir -p "{params.benchdir}"

        Rscript "{params.script}" \
            --mdiff "{input.mdiff}" \
            --tiled-mdiff "{input.tiled_mdiff}" \
            --matrix "{input.matrix}" \
            --annotation "{input.annotation}" \
            --outdir "{params.outdir}" \
            --contrast "{wildcards.experiment}" \
            --qvalue-cutoff "{params.qvalue_cutoff}" \
            --meth-diff-cutoff "{params.meth_diff_cutoff}" \
            --top-n-heatmap "{params.top_n_heatmap}" \
            > "{log}" 2>&1

        touch "{output.done}"
        """
