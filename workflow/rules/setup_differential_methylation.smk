rule download_gencode_gtf:
    output:
        gtf=config["annotation_gtf"]
    shell:
        """
        mkdir -p "$(dirname "{output.gtf}")"
        curl -L "{config[annotation_gtf_url]}" -o "{output.gtf}"
        """
