rnaseq-util
===========

utilities for rnaseq:

*  counts2df.R : collect [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html) results into a file, ready for [DESeq](http://www.bioconductor.org/packages/release/bioc/html/DESeq.html)/[edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)

* rFastq.R: sample 10e6 paired reads, save to a tmp folder

* paintCuffFlags.R: convert [cufflinks](http://cufflinks.cbcb.umd.edu/) gtf to bed and set itemRgb according to transcript flag

* GetFasta.R: pull sequence in fasta format from GenBank based on genbank id

* CuffCodingPotential.R: predict coding potential using [CPAT](http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/cpat/_build/html/index.html) export as bed with colour coded coding potential.

* unmappedBam2fq.R: Ever wonder where [tophat's](http://tophat.cbcb.umd.edu/) unmapped.bam reads come from? Convert unmapped.bam file to fastq (apply sequence quality and complexity filters). Now you are set to run [Trinity](http://trinityrnaseq.sourceforge.net/) on fastq files.