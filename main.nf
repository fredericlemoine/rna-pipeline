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

params.datadir = ["$baseDir", "data"].join(File.separator)
params.gtf = "refgene_hg19.gtf.gz"
params.table = "refgene_hg19.table.gz"
params.resultDir= 'result'

gtfF   = file([params.datadir, params.gtf].join(File.separator))
tableF = file([params.datadir, params.table].join(File.separator))
resultDir = file(params.resultDir)

resultDir.with {
    /*if( !empty() ) { deleteDir() }*/
    mkdirs()
}

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
	STAR --runThreadN 2 \
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

process prepareAnnotations{
	output:
	file("annotations.gtf") into annotFile

	shell:
	'''
	/bin/gunzip -c !{gtfF} > gtf
	/bin/gunzip -c !{tableF} > table
	addGeneNameToUcscGFF.pl gtf table > annotations.gtf
	/bin/rm -f gtf table
	'''
}

process prepareDexseqAnnotations{
	input:
	file(annotFile)

	output:
	file("annotations_DEXSeq.gtf") into annotdexseqFile1, annotdexseqFile2

	shell:
	'''
	cat !{annotFile} > ann.gtf
	python /usr/local/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py ann.gtf annotations_DEXSeq.gtf
	rm ann.gtf
	'''
}


process align{
	tag "$sraid"

	cpus 10

	input:
	set val(condition), val(sraid), file(fastq1), file(fastq2) from fastq
	file(refIndex) from refindex.first()

	output:
	set val(condition),val(sraid), file("${sraid}.bam*") into bam

	shell:
	"""
	/bin/tar -xf !{refIndex}
	STAR --outSAMstrandField intronMotif --outFilterMismatchNmax 4 --outFilterMultimapNmax 10 --genomeDir ref --readFilesIn <(gunzip -c !{fastq1}) <(gunzip -c !{fastq2}) --runThreadN 2  --outSAMunmapped None   --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate  --genomeLoad NoSharedMemory --limitBAMsortRAM 3000000000  > !{sraid}.bam
	samtools index !{sraid}.bam
	/bin/rm -rf ref/
	"""
}


process countReads {
	input:
	file(annot) from annotdexseqFile1.first()
	set val(condition),val(sraid), file(bam) from bam

	output:
	set val(condition),val(sraid),file("counts.txt") into counts

	shell:
	'''
	python /usr/local/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -s no -f bam !{annot} *.bam counts.txt
	'''
}

process mapCounts{
	input:
        set val(condition),val(sraid),file(counts) from counts

	output:
	file("mapcounts") into mappedcounts

	shell:
	'''
	#!/usr/bin/env bash

	grep -v "^_" !{counts} | awk '{print "!{condition}\\t!{sraid}\\t" $0}' > mapcounts
	'''
}

allcounts = mappedcounts.collectFile(name:'allcounts.txt')

process analyzeSplicing {

	cache false

	input:
	file(acounts) from allcounts
	file(annot) from annotdexseqFile2.first()

	output:
	file("*_out.*") into dexseqout mode flatten

	shell:
	'''
	#!/usr/bin/env Rscript
	library(DEXSeq)
	library(reshape2)
	options(bitmapType='cairo')

	## Count data
	counts<-read.table("!{acounts}")
	colnames(counts)=c("cond","sraid","exon","count")
	widecount=dcast(counts, exon ~ sraid,value.var="count")
	row.names(widecount)=widecount$exon
	widecount=widecount[,-1]

	## Exon and Gene Names
	exons=sapply(strsplit(row.names(widecount), ":"),"[",2)
	genes=sapply(strsplit(row.names(widecount), ":"),"[",1)

	## Sample Annotation
	samples=unique(counts[,c(1,2)])$cond
	sampleTable <- data.frame(lapply(unique(counts[,c(1,2)]), as.character),libType="paired-end",stringsAsFactors=FALSE)
	row.names(sampleTable)=sampleTable$sraid
	sampleTable=sampleTable[,-2]
	colnames(sampleTable)=c("condition","libType")
	# on remet dans l'ordre
	sampleTable=sampleTable[colnames(widecount),]

	# Write into individual files
	countfiles=paste0(colnames(widecount),".txt")
	for(sample in colnames(widecount)){
		write.table(file=paste0(sample,".txt"),data.frame(row.names=row.names(widecount),count=widecount[,sample]),row.names=TRUE,col.names=FALSE,quote=FALSE,sep="\t")
	}

	# Create DEXSeqDataSet
	dxd= DEXSeqDataSetFromHTSeq(countfiles,sampleData=sampleTable,design=~sample+exon+condition:exon,flattenedfile=!{annot})

	# Stat analysis
	dxd=estimateSizeFactors(dxd)
	dxd=estimateDispersions(dxd)

	png("dispersion_out.png")
	plotDispEsts(dxd)
	dev.off()

	dxd=testForDEU(dxd)
	dxd=estimateExonFoldChanges(dxd,fitExpToVar="condition")
	dxr1=DEXSeqResults( dxd )
	dxr1=na.omit(dxr1)
	#table(dxr1$pvalue<0.1)
	write.table(file="diff_exons_out.txt",dxr1[dxr1$padj<0.1,])
	#table(tapply(dxr1$padj<0.1,dxr1$groupID,any))

	png("maplot_out.png")
	plotMA(dxr1,cex=0.8)
	dev.off()

	for(i in unique(dxr1[dxr1$padj<0.1,"groupID"])){
	      png(paste0(i,"_out.png"))
	      plotDEXSeq( dxr1,i,legend=TRUE,cex.axis=1.2,cex=1.3,lwd=2,norCounts=TRUE,splicing=TRUE,displayTranscripts=TRUE)
	      dev.off()
	}
	'''
}

dexseqout.subscribe{
	file->file.copyTo(resultDir.resolve(file.name))
}
