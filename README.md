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

Getting Started
===============

System requirements
-------------------

The total run time of the tools alone is approximately 70 minutes on a
2-core 2 GHz system. 4 GB of RAM and 5 GB of disk space is required.
The following software is required:

* ABySS 1.3.0: assemble short reads *de novo*
	http://www.bcgsc.ca/platform/bioinfo/software/abyss
* BWA 0.5.9: align short reads
	http://bio-bwa.sourceforge.net
* IGV 2.0.11: visualize a genome
	http://www.broadinstitute.org/igv
* Java 6: execute Java programs
	http://www.java.com
* samtools 0.1.16: manipulate SAM/BAM files
	http://samtools.sourceforge.net
* snpEff 1.9.5: determine the effect of SNPs
	http://snpeff.sourceforge.net
* tabix 0.2.5: index tab-delimited files
	http://sourceforge.net/projects/samtools/files/tabix

Set up the environment
----------------------

Create a working directory and edit the file `environment'.

	mkdir ~/abyss
	cd ~/abyss
	gvim environment

This shell script will set environment variables to point to the
location of the software. This script is an example. You will need to
set these paths to reflect your directory structure.

	#!/bin/sh
	export top=~/abyss
	PATH=$top/bin:$PATH
	export snpeff=/opt/snpEff-1.9.5
	export ref=$top/chr3.fa

Each time you open a new terminal, you will need to run the following
two commands. Change to the working directory and source the
environment script:

	cd ~/abyss
	source environment

Check that the tools are installed in the PATH.

	source environment
	which abyss-fac abyss-pe bcftools bgzip bwa java samtools tabix vcfutils.pl
	ls $snpeff/snpEff.jar $snpeff/snpEff.config $snpeff/data/hg37

Download the workshop scripts.

	mkdir bin
	wget -P bin ftp://ftp.bcgsc.ca/public/sjackman/bin/*

Check that the workshop scripts are in the PATH.

	which run-abyss run-bwa run-bwasw run-snpeff run-bcftools run-bcftools-assembly

Download the FASTQ files.

	wget ftp://ftp.bcgsc.ca/public/sjackman/30CJCAAXX_4_1.fq.gz
	wget ftp://ftp.bcgsc.ca/public/sjackman/30CJCAAXX_4_2.fq.gz

Exercise 0: Index the reference using BWA
=========================================

Download the reference file of Human chromosome 3.

	cd $top
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr3.fa.gz
	gunzip chr3.fa.gz

Index the reference file.

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
	zcat 30CJCAAXX_4_1.fq.gz |head

How long are the reads? (hint: use `wc -L`)

	A: zcat 30CJCAAXX_4_1.fq.gz |head |wc -L
	A: 50 bp

How many lines are there in both files? (hint: use `wc -l`)

	A: zcat 30CJCAAXX_4_1.fq.gz 30CJCAAXX_4_2.fq.gz |wc -l
	A: 40,869,448 lines

How many lines per read?

	A: 4 lines per read

How many reads are there in both files?

	A: 40869448 / 4 = 10,217,362 reads

How many bases are sequenced?

	A: 10217362*50 = 510,868,100 bp

Assuming the BAC is 200 kbp, what is the depth of coverage?

	A: 510868100/200000 = 2554 fold coverage

Exercise 3: Assemble the reads into contigs using ABySS
=======================================================

Run the assembly.

	cd $top
	k=48 run-abyss 2>&1 |tee abyss.log

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

	k=32 run-abyss se-contigs -n
	k=32 run-abyss pe-sam -n
	k=32 run-abyss pe-contigs -n

Once the assembly has completed, view the contigs in a text editor.

	gview k48/HS0674-contigs.fa

Disabling line wrap makes it easier to browse the file. For Emacs,
select the option "Options -> Line Wrapping in this Buffer -> Truncate
Long Lines." For gview, select the option "Edit -> File Settings ->
Toggle Line Wrap", or type `:set nowrap`

How many contigs are longer than 100 bp?

	A: 5

What is the length of the longest contig (hint: use wc -L)?

	A: wc -L k48/HS0674-contigs.fa
	A: 121927 bp

What is the N50 of the assembly?

	abyss-fac k48/HS0674-contigs.fa
	A: n   n:200  n:N50  min   N80    N50     N20     max     sum
	A: 10  5      1      8044  54747  121861  121861  121861  212397
	A: 122 kbp

View the assembly log in a text editor.

	gview abyss.log

What portion of the reads align to the assembly?
(hint: search for "Aligned")

	A: 7173439 of 10217362 reads (70.2083%)

What is the median fragment size and standard deviation of this library?
(hint: search for "median")

	A: median = 204 bp, sd = 18 bp

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

	A: 8,044 bp and 16,561 bp

Click "browser" for the best alignment and then zoom out 10x.
To which chromosome and band do these contigs align?

	A: chr3q27.3

What are the nearest two genes?

	A: SST and RTP2

Set the "Common SNPs" track (in Variation and Repeats) to "pack".
A SNV is displayed with a red line. Zoom in on a SNV. Is it in dbSNP?

Zoom in on the gap between the two contigs.
Set the "RepeatMasker" track (in Variation and Repeats) to "full".
What feature overlaps the gap that likely caused the assembly gap?

	A: Simple repeat (TA)n

Zoom in to see the sequence of the feature.

Zoom out to see the alignment of both contigs.
Unaligned query sequence is shown with a thin purple line at the end
of the alignment. The thin purple line can be difficult to see.
Which contig has unaligned sequence at one end?

	A: the 16.5 kbp contig

Select that contig, and copy the unaligned sequence to the clipboard.
Aligned sequence is shown in blue upper-case characters, and unaligned
sequence is shown in black lower-case characters.
Open BLAST in a web browser.

	http://blast.ncbi.nlm.nih.gov

Select "nucleotide blast".
Select the database "Nucleotide collection (nr/nt)".
Paste the sequence into the query box. Click BLAST.
To what sequence is the best BLAST hit?

	A: Cloning vector pTARBAC6, complete sequence

What is the cause of this chimeric contig?

	A: This contig contains the cloning vector as well as the human insert.

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

	A: a T/G SNV and a 24-bp insertion

Go to the region chr3:186,676,730

	g chr3:186,676,730

Notice the Ns in the contig sequence, indicating a scaffold gap.
How many Ns are in the contig?

	A: 19 Ns

How many Ns should there be for the size of the gap to agree with the reference?

	A: 15 Ns

Exercise 7: Browse the contig to reference alignments using IGV
===============================================================

Start IGV.

	http://www.broadinstitute.org/igv/projects/current/igv.jnlp

Select "View -> Preferences... -> Alignments" and change
"Visibility range threshold (kb)" to 1000. Restart IGV.
Select the reference "Human hg19".
Then select "File->Load from File..." and "k48/bwasw/HS0674-contigs.bam"
Go to the region chr3:186,600,000-187,600,000 by entering it into
the box labeled "Go".
IGV may take up to a minute to load.

What genes overlap the contigs?

	A: ST6GAL1, SST, RTP2 and BCL6

Add the dbSNP track.
Select "File->Load from Server..."  then expand "hg19" and
"Variation and Repeats" and select "dbSNP 1.3.1".
Zoom in on a SNV. Is it in dbSNP? Is it coding?

Bonus: Find a coding SNV. What is its dbSNP rs ID?

	A: rs1973791 (chr3:187,416,634)
	A: rs11707167 (chr3:187,416,719)
	A: rs1056932 (chr3:187,447,032)

Exercise 8: View the contig to reference alignments SAM file
============================================================

View the SAM file in a text editor. Disable line wrap.

	cd $top/k48/bwasw
	gview HS0674-contigs.sam

The contig ID is given in the first column, and the position of the
contig on the reference is given in the third and fourth columns.
Which contig has two alignments, and what are the positions of these
two alignments?

	A: The contig that has alignments starting at chr3:186,644,066 and chr3:187,439,792.

The orientation is given in the second column. The numbers 0 and 16
indicate positive and negative orientation respectively. What is the
position and orientation of these two alignments?

	A: chr3:186,644,066 (-) and chr3:187,439,792 (+)

In IGV, go to the region chr3:186,600,000-187,600,000
and find these two alignments.
What large-scale structural rearrangement has occurred, and what is
its approximate size?

	A: a ~800 kbp inversion

Which two genes are fused as a result of this rearrangement?

	A: ST6GAL1 and BCL6

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

	A: chr3:187,416,719

What is its dbSNP rs ID? (hint: use IGV or the UCSC genome browser)

	A: rs11707167

Open dbSNP in a web browser.

	http://www.ncbi.nlm.nih.gov/projects/SNP/

What is the minor allele frequency (MAF) of this SNP?

	A: C=0.490

Exercise 12: Compare the assembly variants to the read-alignment variants (optional)
====================================================================================

Browse the variants called by both methods using IGV.

Start IGV.

	http://www.broadinstitute.org/igv/projects/current/igv.jnlp

Select:

	File->Load from File... k48/bwasw/HS0674-contigs.bam
	File->Load from File... k48/bwasw/HS0674-contigs.var.vcf.gz
	File->Load from File... bwa/HS0674.var.vcf.gz

Go to the region chr3:186,648,960

Is the SNV called by both methods?

	A: Yes.

Is the insertion called by both methods?

	A: No, the insertion is called only by the assembly.

Hover the mouse cursor over the insertion to see the inserted sequence.
What is the inserted sequence?

	A: TCTGGGTTCCTTCAAATCCTGCCT

Compare the inserted sequence to the reference sequence at this
location. What sequence do the inserted sequence and the reference
sequence have in common?

	A: The insertion duplicates this 8-bp sequence: TCTGGGTT

Load the alignments of the reads to the reference.

	File->Load from File... bwa/30CJCAAXX_4.bam

The aligned reads do not show the insertion, but the alignments that
span the insertion have mismatches at the end of the alignment.
Do those mismatched bases agree with the inserted sequence?

	A: Yes, the mismatched bases show CT at the 5' end and
	A: CC at the 3' end, which agree with the inserted sequence.

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

This command may fail with the following error:

	[fai_build_core] line length exceeds 65535 in sequence '94'.

How long is the longest line? (hint: use `wc -L`)

	A: wc -L HS0674-contigs.fa
	A: 121,927 bp

Break the long lines into short lines.

	fold ../HS0674-contigs.fa >HS0674-contigs.fa

Index the FASTA file.

	samtools faidx HS0674-contigs.fa

Browse the BAM file using samtools tview.

	samtools tview 30CJCAAXX_4.bam HS0674-contigs.fa

Browse the BAM file using IGV.

	File->Import Genome... k48/HS0674-contigs.fa
	File->Load from File... k48/bwa/30CJCAAXX_4.bam

Find a scaffold gap on the largest contig.
Find the two contigs that have consistent mate pairs joining them.

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/f8201d20af390015a0e171ecfe57738f "githalytics.com")](http://githalytics.com/sjackman/abyss-activity)
