rnaseq-util
===========

utilities for rnaseq:

*  counts2df.R : collect [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html) results into a file, ready for [DESeq](http://www.bioconductor.org/packages/release/bioc/html/DESeq.html)/[edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)

* rFastq.R: sample 10e6 paired reads, save to a tmp folder

* paintCuffFlags.R: convert [cufflinks](http://cufflinks.cbcb.umd.edu/) gtf to bed and set itemRgb according to transcript flag