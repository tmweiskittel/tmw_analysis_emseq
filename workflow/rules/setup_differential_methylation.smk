rule download_gencode_gtf:
    output:
        gtf=config["annotation_gtf"]

    params:
        url=config["annotation_gtf_url"]

    shell:
        """
        mkdir -p "$(dirname "{output.gtf}")"
        curl -L "{params.url}" -o "{output.gtf}"
        """
