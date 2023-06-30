rule bedgraph_to_bigwig:
	input:
		bedgraph = os.path.join(
			macs2_path, "{target}", "{sample}_treat_pileup.bdg"
		),
		chrom_sizes = chrom_sizes
	output:
		bigwig = os.path.join(
			macs2_path, "{target}", "{sample}_treat_pileup.bw"
		)
	conda: "../envs/bedgraph_to_bigwig.yml"
	log: log_path + "/bedgraph_to_bigwig/{target}/{sample}.log"
	threads: 1
	retries: 1
	resources:
		mem_mb = 16384,
		runtime = "3h"
	shell:
		"""
		echo -e "\nConverting {input.bedgraph} to BigWig\n" >> {log}
		TEMPDIR=$(mktemp -d -t bdgXXXXXXXXXX)
		SORTED_BDG=$TEMPDIR/temp.bdg

		## Sort the file
		echo -e "\nSorting as $SORTED_BDG..." >> {log}
		sort -k1,1 -k2,2n -S {resources.mem_mb}M {input.bedgraph} |\
		  egrep $'^chr[0-9XY]+\t' > $SORTED_BDG

		## Convert the file
		echo -e "Done\nConverting..." >> {log}
		bedGraphToBigWig $SORTED_BDG {input.chrom_sizes} {output.bigwig}
		echo -e "Done" >> {log}

		## Remove the temp sorted file
		rm -rf $TEMPDIR
		"""

rule get_coverage_summary:
	input: 
		bw = rules.bedgraph_to_bigwig.output.bigwig,
		blacklist = blacklist
	output: os.path.join(macs2_path, "{target}", "{sample}_treat_pileup.summary")
	params:
		script = "workflow/scripts/get_bigwig_summary.R"
	conda: "../envs/rmarkdown.yml"
	log: log_path + "/get_coverage_summary/{target}/{sample}.log"
	threads: 1
	resources:
		mem_mb = 16384
	shell:
		"""
		## Create the summary tsv
		Rscript --vanilla \
			{params.script} \
			{input.bw} \
			{output} &>> {log}
		"""