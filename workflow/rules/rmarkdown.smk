rule update_extrachips:
    input: os.path.join("workflow", "scripts", "update_extrachips.R")
    output: os.path.join("output", rmd_hash + "_extrachips.updated")
    params:
        version = "1.5.5"
    conda: "../envs/rmarkdown.yml"
    threads: 1
    log: log_path + "/rmarkdown/update_extrachips.log"
    shell:
        """
        Rscript --vanilla {input} {output} {params.version} &>> {log}
        """

rule create_site_yaml:
    input:
        config_yaml = "config/config.yml",
        rmd_yaml = "config/rmarkdown.yml",
        r = "workflow/scripts/create_site_yaml.R"
    output: os.path.join(rmd_path, "_site.yml")
    conda: "../envs/rmarkdown.yml"
    retries: 100
    threads: 1
    log: log_path + "/rmarkdown/create_site_yaml.log"
    shell:
        """
        Rscript --vanilla {input.r} {output} &>> {log}
        """

rule create_setup_chunk:
    input:
        config = "config/rmarkdown.yml",
        r = "workflow/scripts/create_setup_chunk.R"
    output:
        rmd = "analysis/setup_chunk.Rmd"
    conda: "../envs/rmarkdown.yml"
    retries: 100
    threads: 1
    log: log_path + "/rmarkdown/create_setup_chunk.log"
    shell:
        """
        Rscript --vanilla {input.r} {output.rmd} &>> {log}
        """

rule compile_annotations_html:
    input:
        blacklist = blacklist,
        extrachips = rules.update_extrachips.output,
        here = here_file,
        rmd = "workflow/modules/annotation_description.Rmd",
        rds = expand(
        os.path.join(annotation_path, "{file}.rds"),
        file = ['all_gr', 'gene_regions', 'seqinfo', 'trans_models', 'tss']
        ),
        scripts = os.path.join("workflow", "scripts", "custom_functions.R"),
        setup = rules.create_setup_chunk.output,
        site_yaml = rules.create_site_yaml.output,
        yaml = expand(
        os.path.join("config", "{file}.yml"),
        file = ['config', 'colours', 'rmarkdown']
        )
    output:
        rmd = "analysis/annotation_description.Rmd",
        rds = os.path.join(annotation_path, "colours.rds"),
        html = "docs/annotation_description.html",
        fig_path = directory(
            os.path.join("docs", "annotation_description_files", "figure-html")
        )
    conda: "../envs/rmarkdown.yml"
    threads: 1
    log: log_path + "/rmarkdown/compile_annotations_html.log"
    resources:
        mem_mb = 4096,
        disk_mb = 4000,
    shell:
        """
        cp {input.rmd} {output.rmd}
        R -e "rmarkdown::render_site('{output.rmd}')" &>> {log}
        """

rule create_index_rmd:
    input:
        os.path.join("workflow", "modules", "index.Rmd")
    output:
        os.path.join(rmd_path, "index.Rmd")
    retries: 100
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

rule compile_index_html:
    input:
        extrachips = rules.update_extrachips.output,
        html = HTML_OUT,
        here = here_file,		
        rmd = os.path.join(rmd_path, "index.Rmd"),
        setup = rules.create_setup_chunk.output,
        site_yaml = rules.create_site_yaml.output,
        rulegraph = 'workflow/rules/rulegraph.dot'
    output:
        html = "docs/index.html"
    conda: "../envs/rmarkdown.yml"
    threads: 1
    log: log_path + "/rmarkdown/compile_index_html.log"
    shell:
        """
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
        """

