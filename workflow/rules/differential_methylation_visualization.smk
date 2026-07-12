rule differential_methylation_visualization:
    input:
        mdiff = f"{D_OUT}/dmr/diff/methylDiff_{{contrast}}.txt.bgz",
        mdiff_tbi = f"{D_OUT}/dmr/diff/methylDiff_{{contrast}}.txt.bgz.tbi",
        tiled_mdiff = f"{D_OUT}/dmr/diff/methylDiff_{{contrast}}.tiled.txt.bgz",
        tiled_mdiff_tbi = f"{D_OUT}/dmr/diff/methylDiff_{{contrast}}.tiled.txt.bgz.tbi",
        matrix = f"{D_OUT}/dmr/diff/{{contrast}}_pos_meth.tsv",
        annotation = f"{D_OUT}/dmr/annotation/{{contrast}}_annotated.tsv"
    output:
        done = f"{D_OUT}/dmr/visualization/{{contrast}}.visualization.done",
        summary = f"{D_OUT}/dmr/visualization/{{contrast}}/summary_statistics.csv",
        sig = f"{D_OUT}/dmr/visualization/{{contrast}}/significant_DMCs.tsv",
        hyper = f"{D_OUT}/dmr/visualization/{{contrast}}/hypermethylated_DMCs.tsv",
        hypo = f"{D_OUT}/dmr/visualization/{{contrast}}/hypomethylated_DMCs.tsv"
    params:
        outdir = lambda wildcards: f"{D_OUT}/dmr/visualization/{wildcards.contrast}",
        script = f"{R_EMSEQ}/scripts/visualize_differential_methylation.R",
        qvalue_cutoff = config.get("dmr_visualization", {}).get("qvalue_cutoff", 0.05),
        meth_diff_cutoff = config.get("dmr_visualization", {}).get("meth_diff_cutoff", 10),
        top_n_heatmap = config.get("dmr_visualization", {}).get("top_n_heatmap", 500)
    log:
        f"{D_LOGS}/dmr_visualization/{{contrast}}.log"
    benchmark:
        f"{D_BENCHMARK}/dmr_visualization/{{contrast}}.tsv"
    conda:
        ENV_METHYLKIT
    threads:
        4
    shell:
        r"""
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {benchmark})

        Rscript {params.script} \
            --mdiff {input.mdiff} \
            --tiled-mdiff {input.tiled_mdiff} \
            --matrix {input.matrix} \
            --annotation {input.annotation} \
            --outdir {params.outdir} \
            --contrast {wildcards.contrast} \
            --qvalue-cutoff {params.qvalue_cutoff} \
            --meth-diff-cutoff {params.meth_diff_cutoff} \
            --top-n-heatmap {params.top_n_heatmap} \
            > {log} 2>&1

        touch {output.done}
        """
