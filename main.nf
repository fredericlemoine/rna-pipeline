refchrs = Channel.from("chr1", "chr2","chr3","chr4","chr5",
	"chr6","chr7","chr8","chr9","chr10",
	"chr11","chr12","chr13","chr14","chr15",
	"chr16","chr17","chr18","chr19","chr20",
	"chr21","chr22","chrM","chrX","chrY")

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

gtfF   = file([params.datadir, params.gtf].join(File.separator))
tableF = file([params.datadir, params.table].join(File.separator))

process getChromosomes{
	tag "$chr"

	input:
	val(chr) from refchrs

	output:
	file("${chr}.fa.gz") into chrfa

	shell:
	"""
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/!{chr}.fa.gz
	"""
}

process createGenomeIndex{

	input:
	file(chrfas) from chrfa.toList()

	output:
	file("ref.tar") into refindex

	shell:
	"""
	/bin/mkdir ref
	/bin/gunzip -c *.fa.gz > ref.fa
	STAR --runThreadN 4 \
	     --runMode genomeGenerate \
	     --genomeDir ref/ \
	     --genomeFastaFiles ref.fa
	/bin/tar -cvf ref.tar ref/
	/bin/rm -rf ref.fa ref
	"""
}

process getFile{
	tag "$sraid"

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
	tag "$sraid"

	cpus 2

	input:
	set val(condition), val(sraid), file(fastq1), file(fastq2) from fastq
	file(refIndex) from refindex.first()

	output:
	set val(condition),val(sraid), file("${sraid}.bam*") into bam

	shell:
	"""
	/bin/tar -xf !{refIndex}
	STAR --outSAMstrandField intronMotif     \
	    --outFilterMismatchNmax 4            \
	    --outFilterMultimapNmax 10           \
	    --genomeDir  ref/                    \
	    --readFilesIn !{fastq1} !{fastq2}    \
	    --runThreadN 2                       \
	    --outSAMunmapped None                \
	    --outSAMtype BAM SortedByCoordinate  \
	    --outStd BAM_SortedByCoordinate      \
	    --genomeLoad LoadAndKeep             \
	    --outFileNamePrefix .                \
	    --readFilesCommand gunzip -c         \
	    --limitBAMsortRAM 3000000000         \
	    > !{sraid}.bam
	samtools index !{sraid}.bam
	/bin/rm -rf ref/
	"""
}

process prepareAnnotations{
	output:
	file("annotations_DEXSeq.gtf") into annotFile

	shell:
	'''
	/bin/gunzip -c !{gtfF} > gtf
	/bin/gunzip -c !{tableF} > table
	addGeneNameToUcscGFF.pl gtf table > annotations.gtf
	dexseq_prepare_annotation annotations.gtf annotations_DEXSeq.gtf
	/bin/rm -f gtf table annotations.gtf
	'''
}


process countReads {
	input:
	file(annot) from annotFile.first()
	set val(condition),val(sraid), file(bam) from bam

	output:
	set val(condition),val(sraid),file("counts.txt") into counts

	shell:
	'''
	dexseq_count -p yes -r pos -s no -f bam !{annot} !{bam} counts.txt
	'''
}

/*process analyzeSplicing {
	input:
	set val(condition),val(sraid),file("counts.txt") into counts
	
	output:
	file("*.png")

	shell:
	'''
	#!/usr/bin/env Rscript
	
	
	'''
}
*/

counts.subscribe{
	cond, sraid, counts -> println ${cond}+" "+sraid+" "+counts.name
}
