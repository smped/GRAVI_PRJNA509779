rule update_extrachips:
    input: os.path.join("workflow", "scripts", "update_extrachips.R")
    output: temp(os.path.join("output", rmd_hash + "_extrachips.updated"))
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

