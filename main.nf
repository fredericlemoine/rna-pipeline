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

process getChromosomes{
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
	file("ref.tar.gz") into refindex

	shell:
	"""
	mkdir ref
	gunzip *.fa.gz
	STAR --runThreadN 4 \
	     --runMode genomeGenerate \
	     --genomeDir ref/ \
	     --genomeFastaFiles *.fa \
	tar -cvf ref.tar3B.gz ref/
	rm -rf ref
	"""
}

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
	file(refIndex) from refindex.first()

	output:
	set val(condition),val(sraid), file("${sraid}.bam*") into bam

	shell:
	"""
	tar -xf !{refIndex}
	STAR --outSAMstrandField intronMotif     \
	    --outFilterMismatchNmax 4            \
	    --outFilterMultimapNmax 10           \
	    --genomeDir  ref/                    \
	    --readFilesIn !{fastq1} !{fastq2}    \
	    --runThreadN 4                       \
	    --outSAMunmapped None                \
	    --outSAMtype BAM SortedByCoordinate  \
	    --outStd BAM_SortedByCoordinate      \
	    --genomeLoad LoadAndKeep             \
	    --outFileNamePrefix $ALIGNPATH/${i}/ \
	    --readFilesCommand gunzip -c         \
	    --limitBAMsortRAM 3000000000         \
	    > !{sraid}.bam
	samtools index !{sraid}.bam
	"""
}

bam.subscribe{
cond, sraid, bamf -> println ${cond}+" "+sraid+" "+bamf.name
}
