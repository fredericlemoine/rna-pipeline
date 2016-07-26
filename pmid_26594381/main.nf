refchrs = Channel.from("chr1", "chr2","chr3","chr4","chr5",
	"chr6","chr7","chr8","chr9","chr10",
	"chr11","chr12","chr13","chr14","chr15",
	"chr16","chr17","chr18","chr19","chr20",
	"chr21","chr22","chrM","chrX","chrY")

sraids = Channel.from(
   ["Adult_Colon",   "SRR2012208"], /*  Biochain_Adult_Colon    */
   ["Adult_Colon",   "SRR2012209"], /*  Agilent_Adult_Colon     */
   ["Fetal_Colon",   "SRR2012627"], /*  Agilent_Fetal_Colon     */
   ["Fetal_Colon",   "SRR2014228"], /*  Biochain_Fetal_Colon    */
   ["Fetal_Stomach", "SRR2014229"], /*  Agilent_Fetal_Stomach   */
   ["Fetal_Stomach", "SRR2014231"], /*  Biochain_Fetal_Stomach  */
   ["Adult_Heart",   "SRR2014232"], /*  Agilent_Adult_Heart     */
   ["Adult_Heart",   "SRR2014233"], /*  Biochain_Adult_Heart    */
   ["Adult_Lung",    "SRR2014234"], /*  Agilent_Adult_Lung      */
   ["Adult_Lung",    "SRR2014235"], /*  Biochain_Adult_Lung     */
   ["Adult_Liver",   "SRR2014237"], /*  Agilent_Adult_Liver     */
   ["Adult_Liver",   "SRR2014238"], /*  Biochain_Adult_Liver    */
   ["Adult_Kidney",  "SRR2014239"], /*  Agilent_Adult_Kidney    */
   ["Adult_Kidney",  "SRR2014240"], /*  Biochain_Adult_Kidney   */
   ["Adult_Stomach", "SRR2014241"], /*  Agilent_Adult_Stomach   */
   ["Adult_Stomach", "SRR2014242"], /*  Biochain_Adult_Stomach  */
   ["Adult_Stomach", "SRR2014243"], /*  OriGene_Adult_1_Stomach */
   ["Adult_Stomach", "SRR2014244"], /*  OriGene_Adult_2_Stomach */
   ["Adult_Stomach", "SRR2014245"]  /*  OriGene_Adult_3_Stomach */
)

params.datadir = ["$baseDir", "data"].join(File.separator)
params.gtf = "refgene_hg19.gtf.gz"
params.table = "refgene_hg19.table.gz"
params.resultDir= 'result'

gtfF   = file([params.datadir, params.gtf].join(File.separator))
tableF = file([params.datadir, params.table].join(File.separator))
resultDir = file(params.resultDir)

gtfChan   = Channel.fromPath( gtfF )
tableChan = Channel.fromPath( tableF )

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
	'''
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/!{chr}.fa.gz
	'''
}

process createGenomeIndex{

	memory '30 GB'
	cpus 10

	input:
	file(chrfas) from chrfa.toList()

	output:
	file("ref/*") into refindex

	shell:
	'''
	/bin/mkdir ref
	/bin/gunzip -c *.fa.gz > ref.fa
	STAR --runThreadN 10 \
	     --runMode genomeGenerate \
	     --genomeDir ref/ \
	     --genomeFastaFiles ref.fa
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
	SRAID=!{sraid}
	PREFIX1=${SRAID:0:3}
	PREFIX2=${SRAID:0:6}
	ascp -i $ASCPKEY -k 1 -T -l300m \
	  anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/$PREFIX1/$PREFIX2/!{sraid}/!{sraid}.sra .
	fastq-dump --gzip --split-files ./!{sraid}.sra
	rm ./!{sraid}.sra
	'''
}

process prepareAnnotations{
	input:
	file(gtf) from gtfChan.first()
	file(table) from tableChan.first()

	output:
	file("annotations.gtf") into annotFile

	shell:
	'''
	/bin/gunzip -c !{gtf} > gtf
	/bin/gunzip -c !{table} > table
	addGeneNameToUcscGFF.pl gtf table > annotations.gtf
	/bin/rm -f gtf table
	'''
}

process align{
	tag "$sraid"

	cpus 10
	memory "30 GB"

	input:
	set val(condition), val(sraid), file(fastq1), file(fastq2) from fastq
	file(refIndex) from refindex.first()

	output:
	set val(condition),val(sraid), file("${sraid}.bam") into bam

	shell:
	'''
	mkdir ref
	mv !{refIndex} ref/
	STAR --outSAMstrandField intronMotif \
	     --outFilterMismatchNmax 4 \
	     --outFilterMultimapNmax 10 \
	     --genomeDir ref \
	     --readFilesIn <(gunzip -c !{fastq1}) <(gunzip -c !{fastq2}) \
	     --runThreadN 10 \
	     --outSAMunmapped None \
	     --outSAMtype BAM SortedByCoordinate \
	     --outStd BAM_SortedByCoordinate \
	     --genomeLoad NoSharedMemory \
	     --limitBAMsortRAM 10000000000 \
	     > !{sraid}.bam
	'''
}

process index {
	tag "$sraid"

	cpus 1
	memory "1 GB"

	input:
	set val(condition),val(sraid), file(bam) from bam

	output:
	set val(condition),val(sraid), file(bam), file("*.bai") into bamindexed

	shell:
	'''
	samtools index !{bam}
	'''
}

process countReads {
	input:
	file(annot) from annotFile.first()
	set val(condition),val(sraid), file(bam) from bamindexed

	output:
	set val(condition),val(sraid),file("counts.txt") into counts

	shell:
	'''
	featureCounts  -O --largestOverlap -p -t exon -g gene_id -s 0 -a !{annot} -o counts.txt *.bam
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
	#dds <- nbinomWaldTest(dds)
	#res <- results(dds)

	# MA-plot
	#png("maplot_out.png")
	#plotMA(res)
	#dev.off()

	# Export the results in a easily readable file
	#write.table(res, file="results_out.txt")
	#write.table(res[which(abs(res$log2FoldChange)>log2(2) & res$padj<0.1),], file="deg_results_out.txt")
	#sink("results_summary_out.txt")
	#summary(res)
	#sink()  # returns output to the console
	'''
}

deseqout.subscribe{
	file->file.copyTo(resultDir.resolve(file.name))
}
