rule upload_differential_methylation_results:
    input:
        mbase=f"{D_OUT}/dmr/diff/methylBase_{{experiment}}.txt.bgz",
        mbase_tbi=f"{D_OUT}/dmr/diff/methylBase_{{experiment}}.txt.bgz.tbi",
        mdiff=f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.txt.bgz",
        mdiff_tbi=f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.txt.bgz.tbi",
        tiled_mbase=f"{D_OUT}/dmr/diff/methylBase_{{experiment}}.tiled.txt.bgz",
        tiled_mbase_tbi=f"{D_OUT}/dmr/diff/methylBase_{{experiment}}.tiled.txt.bgz.tbi",
        tiled_mdiff=f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.tiled.txt.bgz",
        tiled_mdiff_tbi=f"{D_OUT}/dmr/diff/methylDiff_{{experiment}}.tiled.txt.bgz.tbi",
        matrix=f"{D_OUT}/dmr/diff/{{experiment}}_pos_meth.tsv",
        annotation=f"{D_OUT}/dmr/annotation/{{experiment}}_annotated.tsv"
    output:
        done=f"{D_OUT}/upload/{{experiment}}.upload.done"
    params:
        bucket=config["meta"]["results_bucket"],
        prefix=config["meta"]["results_prefix"]
    log:
        f"{D_LOGS}/upload_differential_methylation/{{experiment}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.done}")"

        DEST="gs://{params.bucket}/{params.prefix}/differential_methylation/{wildcards.experiment}"

        gcloud storage cp {input.mbase} "$DEST/" > "{log}" 2>&1
        gcloud storage cp {input.mbase_tbi} "$DEST/" >> "{log}" 2>&1
        gcloud storage cp {input.mdiff} "$DEST/" >> "{log}" 2>&1
        gcloud storage cp {input.mdiff_tbi} "$DEST/" >> "{log}" 2>&1
        gcloud storage cp {input.tiled_mbase} "$DEST/" >> "{log}" 2>&1
        gcloud storage cp {input.tiled_mbase_tbi} "$DEST/" >> "{log}" 2>&1
        gcloud storage cp {input.tiled_mdiff} "$DEST/" >> "{log}" 2>&1
        gcloud storage cp {input.tiled_mdiff_tbi} "$DEST/" >> "{log}" 2>&1
        gcloud storage cp {input.matrix} "$DEST/" >> "{log}" 2>&1
        gcloud storage cp {input.annotation} "$DEST/" >> "{log}" 2>&1

        touch "{output.done}"
        """
