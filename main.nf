

sraids = Channel.from(
       ["SF3B1", "SRR628582"], 
       ["SF3B1", "SRR628583"], 
       ["SF3B1", "SRR628584"], 
       ["WT", "SRR628585"], 
       ["WT", "SRR628586"], 
       ["WT", "SRR628587"], 
       ["WT", "SRR628588"],
       ["WT", "SRR628589"]
       )

process getFile{
	input:
	set val(condition), val(sraid) from sraids
	
	output:
	set val(condition),val(sraid), file("${sraid}_1.fastq.gz"), file("${sraid}_2.fastq.gz")  into fastq

	shell:
	"""
	fastq-dump --gzip --split-files !{sraid}
	"""
}

process align{
	input:
	set val(condition), val(sraid), file(fastq1), file(fastq2) from fastq

	output:
	set val(condition),val(sraid), file("${sraid}.bam*") into bam

	shell:
	"""
	echo "coucou !{condition} : !{sraid}" > !{sraid}.bam
	"""
}

bam.subscribe{
cond, sraid, bamf -> println bamf.name+ "-"+ bamf.name
}
