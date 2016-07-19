refgenome = "BX571856.1" /* NCBI genome ACC */

sraids = Channel.from(
       ["NL", "SRR1598811"], /* NL 4*/
       ["NL", "SRR1598810"], /* NL 3*/
       ["BL", "SRR1598809"], /* BL 4*/
       ["BL", "SRR1598808"]  /* BL 3*/
       )

params.datadir = ["$baseDir", "data"].join(File.separator)
params.resultDir= 'result'

gtfF   = file([params.datadir, params.gtf].join(File.separator))
tableF = file([params.datadir, params.table].join(File.separator))
resultDir = file(params.resultDir)

resultDir.with {
    /*if( !empty() ) { deleteDir() }*/
    mkdirs()
}

process getChromosomes{
	tag "$refgenome"

	output:
	file("${refgenome}.fa.gz") into chrfa

	shell:
	'''
	wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=!{refgenome}&rettype=fasta&retmode=text"
	mv efetch* !{refgenome}.fa
	awk '{if($0~/^>/){split($0,a,"|"); print ">" a[4]}else{print $0}}' !{refgenome}.fa q| gzip -c - > !{refgenome}.fa.gz
	rm !{refgenome}.fa
	'''
}

process createGenomeIndex{

	memory '10 GB'
	cpus 1

	input:
	file(chrfas) from chrfa.toList()

	output:
	file("ref/*") into refindex

	shell:
	'''
	/bin/mkdir ref
	/bin/gunzip -c *.fa.gz > ref.fa
	bowtie2-build ref.fa ref/ref
	/bin/rm -rf ref.fa
	'''
}

process getFile{
	tag "$sraid"

	input:
	set val(condition), val(sraid) from sraids
	
	output:
	set val(condition),val(sraid), file("${sraid}_1.fastq.gz"), file("${sraid}_2.fastq.gz")  into fastq

	shell:
	'''
	fastq-dump --gzip --split-files !{sraid}
	'''
}

process prepareAnnotations{
	output:
	file("annotations.gtf") into annotFile

	shell:
	'''
	wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=!{refgenome}&rettype=fasta&retmode=text"
	mv efetch* !{refgenome}.fa
	wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=!{refgenome}&rettype=gb&retmode=text"
	mv efetch* !{refgenome}.gb
	seqret -sequence !{refgenome}.fa -feature -fformat gb -fopenfile !{refgenome}.gb -osformat gff -auto
	mv *.gff annotations.gtf
	'''
}

process align{
	tag "$sraid"

	cpus 10
	memory "10 GB"

	input:
	set val(condition), val(sraid), file(fastq1), file(fastq2) from fastq
	file(refIndex) from refindex.first()

	output:
	set val(condition),val(sraid), file("${sraid}.bam*") into bam

	shell:
	'''
	mkdir ref
	mv !{refIndex} ref/
	bowtie2 -p 10 -x ref/ref -1 !{fastq1} -2 !{fastq2} -S !{sraid}.sam
	samtools view -o !{sraid}_tmp.bam -Shb !{sraid}.sam
	samtools sort -o !{sraid}.bam -m 10G !{sraid}_tmp.bam
	samtools index !{sraid}.bam
	rm -f !{sraid}_tmp.bam !{sraid}.sam
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
	featureCounts -t gene -g gene -s 0 -a !{annot} -o counts.txt *.bam
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

	grep -v "^#" !{counts} | grep -v "^Geneid" | awk '{print "!{condition}\\t!{sraid}\\t" $0}' > mapcounts
	'''
}

allcounts = mappedcounts.collectFile(name:'allcounts.txt')

process diffExpression {

	cache false

	input:
	file(acounts) from allcounts

	output:
	file("*_out.*") into deseqout mode flatten

	shell:
	'''
	#!/usr/bin/env Rscript
	library(DESeq2)
	library(reshape2)
	options(bitmapType='cairo')

	## Count data
	counts<-read.table("!{acounts}")
	colnames(counts)=c("cond","sraid","gene","chr","start","end","strand","length","count")
	widecount=dcast(counts, gene ~ sraid,value.var="count")
	row.names(widecount)=widecount$gene
	widecount=widecount[,-1]

	## Sample Annotation
	samples=unique(counts[,c(1,2)])$cond
	sampleTable <- data.frame(lapply(unique(counts[,c(1,2)]), as.character),libType="paired-end",stringsAsFactors=FALSE)
	row.names(sampleTable)=sampleTable$sraid
	sampleTable=sampleTable[,-2]
	colnames(sampleTable)=c("condition","libType")
	# on remet dans l'ordre
	sampleTable=sampleTable[colnames(widecount),]
	sampleTable$condition=factor(sampleTable$condition, levels=unique(sampleTable$condition))

	dds <- DESeqDataSetFromMatrix(widecount, sampleTable, design=~condition)

	tc <- colSums(counts(dds))
	png("total_counts_out.png")
	barplot(tc, main="Total count", col=as.numeric(sampleTable[names(tc),"condition"])+1)
	legend("topleft", legend = unique(sampleTable[names(tc),"condition"]),col=unique(as.numeric(sampleTable[names(tc),"condition"]))+1,pch=20)
	dev.off()
	
	# Proportion of null counts per sample
	pn <- apply(counts(dds), 2, function(x) mean(x==0))
	png("null_counts_out.png")
	barplot(100*pn, main="Proportion null counts", col=as.numeric(sampleTable[names(pn),"condition"])+1)
	legend("topright", legend = unique(sampleTable[names(tc),"condition"]),col=unique(as.numeric(sampleTable[names(pn),"condition"]))+1,pch=20)
	dev.off()
	# Stat analysis
	dds <- estimateSizeFactors(dds)
	
	sink("size_factors_out.txt")
	sizeFactors(dds)
	sink()  # returns output to the console

	# Difference between before and after normalization
	png("normalization_plots_out.png",width=1000,height=480)
	par(mfrow=c(1,2))
	boxplot(log2(counts(dds)+1), col=as.numeric(sampleTable[colnames(counts(dds)),"condition"])+1,main="Before normalization",cex=0.5)
	legend("topright", legend = unique(sampleTable[names(tc),"condition"]),col=unique(as.numeric(sampleTable[colnames(counts(dds)),"condition"]))+1,pch=20)
	boxplot(log2(counts(dds,normalized=T)+1), col=as.numeric(sampleTable[colnames(counts(dds)),"condition"])+1,main="After normalization",cex=0.5)
	legend("topright", legend = unique(sampleTable[names(tc),"condition"]),col=unique(as.numeric(sampleTable[colnames(counts(dds)),"condition"]))+1,pch=20)
	dev.off()
	
	# estimation of the dispersions
	dds <- estimateDispersions(dds)
	png("dispersion_plots_out.png",width=1000,height=480)
	plotDispEsts(dds)
	dev.off()

	# PCA
	png("pca_out.png")
	plotPCA(varianceStabilizingTransformation(dds))
	dev.off()

	# Hierarchical Clustering
	transCounts <- getVarianceStabilizedData(dds) 
	hc <- hclust(dist(t(transCounts)), method="ward.D")
	png("clustering_out.png")
	plot(hc,main="Sample clustering",xlab="Samples")
	dev.off()
	
	# statistical testing
	dds <- nbinomWaldTest(dds)
	res <- results(dds)

	# MA-plot
	png("maplot_out.png")
	plotMA(res)
	dev.off()

	# Export the results in a easily readable file
	write.table(res, file="results_out.txt")
	sink("results_summary_out.txt")
	summary(res)
	sink()  # returns output to the console
	'''
}

deseqout.subscribe{
	file->file.copyTo(resultDir.resolve(file.name))
}