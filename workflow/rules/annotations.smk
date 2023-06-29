rule download_blacklist:
	output: blacklist
	params:
		url = urllib.parse.urlunparse(
			(
				'https', 'raw.githubusercontent.com',
				'Boyle-Lab/Blacklist/master/lists/' +
				ucsc_build + "-blacklist.v2.bed.gz",
				'', '', ''
			)
		)
	log: log_path + "/downloads/download_blacklist.log"
	shell:
		"""
		curl {params.url} --output {output} 2> {log}
		"""

rule create_annotations:
	input:
		bam = expand(os.path.join(bam_path, "{bam}.bam"), bam = indiv_pre),
		config = ancient(os.path.join("config", "config.yml")),
		extrachips = rules.update_extrachips.output,
		gtf = gtf,
		r = os.path.join("workflow", "scripts", "create_annotations.R"),
		yaml = os.path.join("config", "params.yml")
	output:
		rds = expand(
		  os.path.join(annotation_path, "{file}.rds"),
		  file = ['all_gr', 'gene_regions', 'tss', 'seqinfo', 'trans_models']
		),
		chrom_sizes = chrom_sizes
	conda: "../envs/rmarkdown.yml"
	threads: 1
	resources:
		mem_mb = 16384	
	log: log_path + "/scripts/create_annotations.log"
	shell:
		"""
		Rscript --vanilla {input.r} {input.gtf} {threads} &>> {log}
		"""

