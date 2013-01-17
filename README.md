De novo assembly of Illumina reads using ABySS and alignment using BWA
======================================================================

This workshop is designed by Shaun Jackman <sjackman@gmail.com>.

Purpose
=======

In this lab we will use ABySS to assemble a 200 kbp bacterial
artificial chromosome (BAC). The data set is one lane of paired-end
reads from the Illumina platform. The assembled contigs are aligned to
the human reference genome using BWA-SW and variants are called using
bcftools. IGV is used to visualize these alignments and variants.
snpEff is used to determine the effects of these variants.

After this lab, you will have learned how to use ABySS to assemble a
small genome, use BWA and BWA-SW to align reads and contigs to a
reference genome, use IGV to visualize these alignments, and use bcftools
and snpEff to call variants and determine their effect.

Contents
========

* [Getting Started](#getting-started)
* [Exercise 0: Index the reference using BWA](#exercise-0-index-the-reference-using-bwa)
* [Exercise 1: Align the reads to the reference using BWA (optional)](#exercise-1-align-the-reads-to-the-reference-using-bwa-optional)
* [Exercise 2: Inspect the reads](#exercise-2-inspect-the-reads)
* [Exercise 3: Assemble the reads into contigs using ABySS](#exercise-3-assemble-the-reads-into-contigs-using-abyss)
* [Exercise 4: Align the contigs to the reference using web BLAT](#exercise-4-align-the-contigs-to-the-reference-using-web-blat)
* [Exercise 5: Align the contigs to the reference using BWA-SW](#exercise-5-align-the-contigs-to-the-reference-using-bwa-sw)
* [Exercise 6: Browse the contig to reference alignments using samtools tview](#exercise-6-browse-the-contig-to-reference-alignments-using-samtools-tview)
* [Exercise 7: Browse the contig to reference alignments using IGV](#exercise-7-browse-the-contig-to-reference-alignments-using-igv)
* [Exercise 8: View the contig to reference alignments SAM file](#exercise-8-view-the-contig-to-reference-alignments-sam-file)
* [Exercise 9: Call variants of the reads-to-reference alignments using bcftools (optional)](#exercise-9-call-variants-of-the-reads-to-reference-alignments-using-bcftools-optional)
* [Exercise 10: Call variants of the contigs-to-reference alignments using bcftoolss](#exercse-10-call-variants-of-the-contigs-to-reference-alignments-using-bcftools)
* [Exercise 11: Determine the effects of the SNVs](#exercise-11-determine-the-effects-of-the-snvs)
* [Exercise 12: Compare the assembly variants to the read-alignment variants (optional)](#exercise-12-compare-the-assembly-variants-to-the-read-alignment-variants-optional)
* [Exercise 13: Align the reads to the contigs using BWA (optional)](#exercise-13-align-the-reads-to-the-contigs-using-bwa-optional)

Getting Started
===============

System requirements
-------------------

The total run time of the tools alone is approximately 70 minutes on a
2-core 2 GHz system. 4 GB of RAM and 5 GB of disk space is required.
The following software is required:

* ABySS 1.3.4: assemble short reads *de novo*
	http://www.bcgsc.ca/platform/bioinfo/software/abyss
* BWA 0.6.2: align short reads
	http://bio-bwa.sourceforge.net
* IGV 2.2: visualize a genome
	http://www.broadinstitute.org/igv
* Java 6: execute Java programs
	http://www.java.com
* samtools 0.1.18: manipulate SAM/BAM files
	https://github.com/samtools/samtools
* snpEff 3.1: determine the effect of SNPs
	http://snpeff.sourceforge.net
* tabix 0.2.5: index tab-delimited files
	http://sourceforge.net/projects/samtools/files/tabix

Install the software on Mac OS X
--------------------------------

Install [Homebrew](http://mxcl.github.com/homebrew/).

	ruby -e "$(curl -fsSkL raw.github.com/mxcl/homebrew/go)"

Install ABySS, BWA, samtools and tabix using Homebrew.

	brew tap homebrew/science
	brew install abyss bwa samtools tabix

Install IGV. Close IGV once it opens.

	mkdir ~/abyss
	cd ~/abyss
	wget http://www.broadinstitute.org/igv/projects/current/igv.jnlp
	javaws igv.jnlp

Install snpEff.

	wget http://downloads.sourceforge.net/project/snpeff/snpEff_v3_1_core.zip
	unzip snpEff_v3_1_core.zip
	wget http://downloads.sourceforge.net/project/snpeff/databases/v3_1/snpEff_v3_1_GRCh37.68.zip
	unzip -d snpEff_3_1 snpEff_v3_1_GRCh37.68.zip
	echo "data_dir=$PWD/snpEff_3_1/data" >>snpEff_3_1/snpEff.config

Install the software on Ubuntu and Debian
-----------------------------------------

Install ABySS, BWA, samtools and tabix.

	sudo apt-get install abyss bwa samtools tabix

Install IGV and snpEff using the instructions above for Mac OS X.

Set up the environment
----------------------

Create a working directory.

	mkdir ~/abyss
	cd ~/abyss

Download the workshop scripts using git, if you have it installed.

	git clone git://github.com/sjackman/abyss-activity.git
	mv abyss-activity/* abyss-activity/.git .
	rmdir abyss-activity

If you do not have git installed, use wget.

	wget https://github.com/sjackman/abyss-activity/archive/master.tar.gz
	tar --strip 1 -zxf master.tar.gz

The shell script named `environment` will set environment variables
that specify the location of the installed software. Each time that
you open a new terminal, you will need to run the following two
commands to change to the working directory and source the environment
script:

	cd ~/abyss
	source environment

Check that the tools are installed in the PATH.

	source environment
	which abyss-fac abyss-pe bcftools bgzip bwa java samtools tabix vcfutils.pl
	ls $snpeff/snpEff.jar $snpeff/snpEff.config $snpeff/data/hg37

Check that the workshop scripts are in the PATH.

	which run-abyss run-bwa run-bwasw run-snpeff run-bcftools run-bcftools-assembly

Download the FASTQ files.

	wget ftp://ftp.bcgsc.ca/public/sjackman/30CJCAAXX_4_1.fq.gz
	wget ftp://ftp.bcgsc.ca/public/sjackman/30CJCAAXX_4_2.fq.gz

Download the reference file of Human chromosome 3.

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr3.fa.gz
	gunzip chr3.fa.gz

Exercise 0: Index the reference using BWA
=========================================

Index the reference file.

	cd $top
	bwa index $ref

8 min, 1 GB RAM, 260 MB disk space

Exercise 1: Align the reads to the reference using BWA (optional)
=================================================================

Run BWA.

	cd $top
	mkdir bwa
	cd bwa
	run-bwa 2>&1 |tee bwa.log

40 min, 250 MB RAM, 3 GB disk space

While this job is running in the background, open a new terminal and
continue with the next section. Remember that you have to set the
environment. In the new terminal window type:

	cd ~/abyss
	source environment

Exercise 2: Inspect the reads
=============================

There are two FASTQ files for this one lane of paired-end Illumina
data: one file for the forward reads and one file for the reverse
reads.
Look at the first few reads.

	cd $top
	gunzip -c 30CJCAAXX_4_1.fq.gz |head

How long are the reads? (hint: use `wc -L`)

>     gunzip -c 30CJCAAXX_4_1.fq.gz |head |wc -L
> 50 bp

How many lines are there in both files? (hint: use `wc -l`)

>     gunzip -c 30CJCAAXX_4_1.fq.gz 30CJCAAXX_4_2.fq.gz |wc -l
> 40,869,448 lines

How many lines per read?

> 4 lines per read

How many reads are there in both files?

> 40869448 / 4 = 10,217,362 reads

How many bases are sequenced?

> 10217362*50 = 510,868,100 bp

Assuming the BAC is 200 kbp, what is the depth of coverage?

> 510868100/200000 = 2554 fold coverage

Exercise 3: Assemble the reads into contigs using ABySS
=======================================================

Run the assembly.

	cd $top
	k=48 run-abyss contigs 2>&1 |tee abyss.log

5 min, 200 MB RAM, 2 MB disk space

While the assembly is running, view the script in a text editor. You
may substitute whichever text editor you prefer for `gview`.

	gview $top/bin/run-abyss

Look at the option `-n,--dry-run` of abyss-pe. Its output is the
commands that ABySS will run for the assembly.

	k=32 run-abyss -n

The assembly runs in three stages: assemble contigs without paired-end
information, align the paired-end reads to the initial assembly, and
merge contigs joined by paired-end information. You can instruct ABySS
to stop after any of these stages. Use the -n option to see the
commands for each stage.

	k=32 run-abyss unitigs -n
	k=32 run-abyss pe-sam -n
	k=32 run-abyss contigs -n

Once the assembly has completed, view the contigs in a text editor.

	gview k48/HS0674-contigs.fa

Disabling line wrap makes it easier to browse the file. For Emacs,
select the option "Options -> Line Wrapping in this Buffer -> Truncate
Long Lines." For gview, select the option "Edit -> File Settings ->
Toggle Line Wrap", or type `:set nowrap`

How many contigs are longer than 100 bp?

> 6

What is the length of the longest contig (hint: use wc -L)?

>     wc -L k48/HS0674-contigs.fa
> 67530 bp

What is the N50 of the assembly?

	abyss-fac k48/HS0674-contigs.fa

>     n      n:200  n:N50  min    N80    N50    N20    max    sum
>     14     6      2      8044   54373  54746  67484  67484  212392
> 54746 bp

View the assembly log in a text editor.

	gview abyss.log

What portion of the reads align to the assembly?
(hint: search for "Mapped")

> Mapped 7172456 of 10217362 reads (70.2%)

What is the median fragment size and standard deviation of this library?
(hint: search for "median")

> median = 204 bp, sd = 18 bp

Exercise 4: Align the contigs to the reference using web BLAT
=============================================================

Open BLAT in a web browser.

	http://genome.ucsc.edu

View the assembled contigs in a text editor.

	gview k48/HS0674-contigs.fa

Disabling line wrap makes it easier to select the full sequence.

Select the two contigs whose lengths are approximately 8 and 16.5 kbp
and copy-and-paste their sequence into BLAT.

What is the exact length of these two contigs?

> 8,044 bp and 16,561 bp

Click "browser" for the best alignment and then zoom out 10x.
To which chromosome and band do these contigs align?

> chr3q27.3

What are the nearest two genes?

> SST and RTP2

Set the "Common SNPs" track (in Variation and Repeats) to "pack".
A SNV is displayed with a red line. Zoom in on a SNV. Is it in dbSNP?

Zoom in on the gap between the two contigs.
Set the "RepeatMasker" track (in Variation and Repeats) to "full".
What feature overlaps the gap that likely caused the assembly gap?

> Simple repeat (TA)n

Zoom in to see the sequence of the feature.

Zoom out to see the alignment of both contigs.
Unaligned query sequence is shown with a thin purple line at the end
of the alignment. The thin purple line can be difficult to see.
Which contig has unaligned sequence at one end?

> the 16.5 kbp contig

Select that contig, and copy the unaligned sequence to the clipboard.
Aligned sequence is shown in blue upper-case characters, and unaligned
sequence is shown in black lower-case characters.
Open BLAST in a web browser.

	http://blast.ncbi.nlm.nih.gov

Select "nucleotide blast".
Select the database "Nucleotide collection (nr/nt)".
Paste the sequence into the query box. Click BLAST.
To what sequence is the best BLAST hit?

> Cloning vector pTARBAC2.1

What is the cause of this chimeric contig?

> This contig contains the cloning vector as well as the human insert.

Exercise 5: Align the contigs to the reference using BWA-SW
===========================================================

Warning: do not start BWA-SW until BWA has completed unless your
machine has at least 2 GB of RAM.
Run BWA-SW.

	cd $top
	mkdir k48/bwasw
	cd k48/bwasw
	ln -s ../HS0674-contigs.fa .
	run-bwasw

1 min, 800 MB RAM, 1 MB disk space

While the alignment is running, view the script in a text editor.

	gview $top/bin/run-bwasw

Exercise 6: Browse the contig to reference alignments using samtools tview
==========================================================================

Run samtools tview.

	cd $top/k48/bwasw
	samtools tview HS0674-contigs.bam $ref

Go to the region chr3:186,648,940

	g chr3:186,648,940

What variants do you see?

> a T/G SNV and a 24-bp insertion

Go to the region chr3:186,676,730

	g chr3:186,676,730

Notice the Ns in the contig sequence, indicating a scaffold gap.
How many Ns are in the contig?

> 20 Ns

How many Ns should there be for the size of the gap to agree with the reference?

> 16 Ns

Exercise 7: Browse the contig to reference alignments using IGV
===============================================================

Start IGV.

	javaws http://www.broadinstitute.org/igv/projects/current/igv.jnlp

Select "View -> Preferences... -> Alignments" and change
"Visibility range threshold (kb)" to 1000.
Select the "Genomes -> Load Genome From Server... -> Human hg19".
Then select "File->Load from File..." and "k48/bwasw/HS0674-contigs.bam"
Go to the region chr3:186,600,000-187,600,000 by entering it into
the box labeled "Go".
IGV may take up to a minute to load.

What genes overlap the contigs?

> ST6GAL1, SST, RTP2 and BCL6

Add the dbSNP track.
Select "File->Load from Server..."  then expand "Annotations" and
"Variation and Repeats" and select "dbSNP 1.3.1".
Zoom in on a SNV. Is it in dbSNP? Is it coding?

Bonus: Find a coding SNV. What is its dbSNP rs ID?

> rs1973791 (chr3:187,416,634),
> rs11707167 (chr3:187,416,719),
> rs1056932 (chr3:187,447,032)

Exercise 8: View the contig to reference alignments SAM file
============================================================

View the SAM file in a text editor. Disable line wrap.

	cd $top/k48/bwasw
	gview HS0674-contigs.sam

The contig ID is given in the first column, and the position of the
contig on the reference is given in the third and fourth columns.
Which large contig has two alignments, and what are the positions of
these two alignments?

> The contig that has alignments starting at
> chr3:186,698,393 and chr3:187,439,792.

The orientation is given in the second column. The numbers 0 and 16
indicate positive and negative orientation respectively. What is the
position and orientation of these two alignments?

> chr3:186,698,393 (-) and chr3:187,439,792 (+)

In IGV, go to the region chr3:186,600,000-187,600,000
and find these two alignments.
What large-scale structural rearrangement has occurred, and what is
its approximate size?

> a ~800 kbp inversion

Which two genes are fused as a result of this rearrangement?

> ST6GAL1 and BCL6

Exercise 9: Call variants of the reads-to-reference alignments using bcftools (optional)
=============================================================================================

Check that BWA has completed aligning the reads to the reference.
Run bcftools in this terminal.

	cd $top/bwa
	run-bcftools 2>&1 |tee bcftools.log

20 min, 250 MB RAM, 120 MB disk space

While bcftools is running, view the script in a text editor.

	gview $top/bin/run-bcftools

Exercse 10: Call variants of the contigs-to-reference alignments using bcftools
===============================================================================

Run bcftools.

	cd $top/k48/bwasw
	run-bcftools-assembly 2>&1 |tee bcftools.log

Browse the variants using IGV.
Select "File->Load from File... k48/bwasw/HS0674-contigs.var.vcf.gz"
Right-click on the VCF track and select "Color by Allele".

Exercise 11: Determine the effects of the SNVs
==============================================

Run snpEff.

	cd $top/k48/bwasw
	run-snpeff HS0674-contigs.var.vcf.gz >HS0674-contigs.var.snpeff

1 min, 1.0 GB RAM

View the output of snpEff in a text editor.

	gview HS0674-contigs.var.snpeff

Count the number of SNVs in each category of effect.

	sort -uk1,2 HS0674-contigs.var.snpeff |cut -f16 |cut -d: -f1 |sort |uniq -c |sort -rn

Find all the coding SNVs.

	grep CODING HS0674-contigs.var.snpeff |cut -f1-4,16-18 |uniq

Find the non-synonymous SNV. What is its location?

> chr3:187,416,719

What is its dbSNP rs ID? (hint: use IGV or the UCSC genome browser)

> rs11707167

Open dbSNP in a web browser.

	http://www.ncbi.nlm.nih.gov/projects/SNP/

What is the minor allele frequency (MAF) of this SNP?

> T=0.4716

Exercise 12: Compare the assembly variants to the read-alignment variants (optional)
====================================================================================

Browse the variants called by both methods using IGV.

Start IGV.

	javaws http://www.broadinstitute.org/igv/projects/current/igv.jnlp

Select:

	File->Load from File... k48/bwasw/HS0674-contigs.bam
	File->Load from File... k48/bwasw/HS0674-contigs.var.vcf.gz
	File->Load from File... bwa/HS0674.var.vcf.gz

Go to the region chr3:186,648,960

Is the SNV called by both methods?

> Yes.

Is the insertion called by both methods?

> No, the insertion is called only by the assembly.

Hover the mouse cursor over the insertion to see the inserted sequence.
What is the inserted sequence?

> TCTGGGTTCCTTCAAATCCTGCCT

Compare the inserted sequence to the reference sequence at this
location. What sequence do the inserted sequence and the reference
sequence have in common?

> The insertion duplicates this 8-bp sequence: TCTGGGTT

Load the alignments of the reads to the reference.

	File->Load from File... bwa/30CJCAAXX_4.bam

The aligned reads do not show the insertion, but the alignments that
span the insertion have mismatches at the end of the alignment.
Do those mismatched bases agree with the inserted sequence?

> Yes, the mismatched bases show CT at the 5' end and CC
> at the 3' end, which agree with the inserted sequence.

Exercise 13: Align the reads to the contigs using BWA (optional)
================================================================

Run BWA.

	cd $top
	mkdir k48/bwa
	cd k48/bwa
	bwa index ../HS0674-contigs.fa
	ref=../HS0674-contigs.fa run-bwa 2>&1 |tee bwa.log

12 min, 200 MB RAM, 3 GB disk space

Index the assembly FASTA file.

	samtools faidx ../HS0674-contigs.fa

Browse the BAM file using samtools tview.

	samtools tview 30CJCAAXX_4.bam ../HS0674-contigs.fa

Browse the BAM file using IGV.

	File->Import Genome... k48/HS0674-contigs.fa
	File->Load from File... k48/bwa/30CJCAAXX_4.bam

Find a scaffold gap on the largest contig.
Find the two contigs that have consistent mate pairs joining them.

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/f8201d20af390015a0e171ecfe57738f "githalytics.com")](http://githalytics.com/sjackman/abyss-activity)
