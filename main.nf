channel = Channel.from("COUCOU")

process STARHelp{
	input:
	val(message) from channel
	
	output:
	set file("coucou.txt"),file("essai.pdf") into starhelp

	shell:
	"""
	echo !{message} > coucou.txt
	STAR --help >> coucou.txt
	essai.Rscript
	"""
}


starhelp.subscribe{
f1,f2 -> println f1.name+ "-"+ f2.name
}
