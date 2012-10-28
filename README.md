De novo assembly of Illumina reads using ABySS and alignment using BWA
2011-07-14 Shaun Jackman <sjackman@gmail.com>

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

* System requirements
Total run time of the tools alone is ~70 minutes on a 2-core 2 GHz
system. 4 GB of RAM and 5 GB of disk space is required. The following
software is required:

ABySS 1.2.7
	http://www.bcgsc.ca/platform/bioinfo/software/abyss
BWA 0.5.9
	http://bio-bwa.sourceforge.net
IGV 2.0.3
	http://www.broadinstitute.org/igv
Java 6
	http://www.java.com
samtools 0.1.16
	http://samtools.sourceforge.net
snpEff 1.9.5
	http://snpeff.sourceforge.net
tabix 0.2.5
	http://sourceforge.net/projects/samtools/files/tabix

* Set up the environment
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
	PATH=/opt/abyss-1.2.7/bin:$PATH
	export snpeff=/opt/snpEff-1.9.5
	export ref=/opt/GRCh37/GRCh37.fa
Check that the tools are installed in the PATH.
	source environment
	which abyss-fac abyss-pe bcftools bgzip bwa java samtools tabix vcfutils.pl
	ls $snpeff/snpEff.jar $snpeff/snpEff.config $snpeff/data/hg37
Copy the workshop scripts to your hard drive.
	mkdir bin
	cp /path_to_the_scripts/run-* bin/
Check that the workshop scripts are in the PATH.
	which run-abyss run-bwa run-bwasw run-snpeff run-bcftools run-bcftools-assembly
Copy the FASTQ files to your hard drive.
	cp /path_to_the_reads/30CJCAAXX_4_[12].fq.gz .

* Align the reads to the reference using BWA
Run BWA
	cd $top
	mkdir bwa
	cd bwa
	run-bwa 2>&1 |tee bwa.log
40 min, 2.4 GB RAM, 3 GB disk space
While this job is running in the background, open a new terminal and
continue with the next section.
	gnome-terminal
	source environment

* Inspect the reads
There are two FASTQ files for this one lane of paired-end Illumina
data: one file for the forward reads and one file for the reverse
reads.
Look at the first few reads.
	cd $top
	zcat 30CJCAAXX_4_1.fq.gz |head
How long are the reads? (hint: use wc -L)
	A: zcat 30CJCAAXX_4_1.fq.gz |head |wc -L
	A: 50 bp
How many lines are there in both files? (hint: use wc -l)
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

* Assemble the reads into contigs using ABySS
Run the assembly.
	cd $top
	k=48 run-abyss 2>&1 |tee abyss.log
5 min, 200 MB RAM, 2 MB disk space
While the assembly is running, view the script in a text editor.
	gview $top/bin/run-abyss
Look at the option -n,--dry-run of aybss-pe. Its output is the
commands that ABySS will run for the assembly.
	k=48 run-abyss -n
The assembly runs in three stages: assemble contigs without paired-end
information, align the paired-end reads to the initial assembly, and
merge contigs joined by paired-end information. You can instruct ABySS
to stop after any of these stages. Use the -n option to see the
commands for each stage.
	k=47 run-abyss se-contigs -n
	k=47 run-abyss se-sam -n
	k=47 run-abyss pe-contigs -n
Once the assembly has completed, view the contigs in a text editor.
	gview k48/HS0674-contigs.fa
Disabling line wrap makes it easier to browse the file.
	:set nowrap
How many contigs are longer than 100 bp?
	A: 5
What is the length of the longest contig (hint: use wc -L)?
	A: wc -L k48/HS0674-contigs.fa
	A: 121927 bp
What is the N50 of the assembly?
	abyss-fac k48/HS0674-contigs.fa
	A: n       n:200   n:N50  min    median mean   N50    max     sum
	A: 10      5       1      8044   16561  42479  121861 121861  212397
	A: 122 kbp
View the assembly log in a text editor.
	gview abyss.log
What portion of the reads align to the assembly?
(hint: search for "Aligned")
	A: 7173439 of 10217362 reads (70.2083%)
What is the median fragment size and standard deviation of this library?
(hint: search for "median")
	A: median = 204 bp, sd = 18 bp

* Align the contigs to the reference using web BLAT
Open BLAT in a web browser.
	http://genome.ucsc.edu
View the assembled contigs in a text editor.
	gview k48/HS0674-contigs.fa
Disabling line wrap makes it easier to select the full sequence.
	:set nowrap
Pick the first two contigs longer than 100 bp and shorter than 25 kbp
and copy-and-paste their sequence into BLAT. How long are these two
contigs?
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
	A: contig 93
Select that contig, and copy the unaligned sequence to the clipboard.
Aligned sequence is shown in blue upper-case characters, and unaligned
sequence is shown in black lower-case characters.
Open BLAST in a web browser.
	http://blast.ncbi.nlm.nih.gov
Select "nucleotide blast".
Select the database "Nucleotide collection (nr/nt)"
Paste the sequence into the query box. Click BLAST.
To what sequence is the best BLAST hit?
	A: Cloning vector pTARBAC6, complete sequence
What is the cause of this chimeric contig?
	A: This contigs contains the cloning vector as well as the human insert.

* Align the contigs to the reference using BWA-SW
Warning: do not start BWA-SW until BWA has completed unless your
machine has at least 8 GB of RAM.
Run BWA-SW.
	cd $top
	mkdir k48/bwasw
	cd k48/bwasw
	ln -s ../HS0674-contigs.fa .
	run-bwasw
1 min, 3.6 GB RAM, 1 MB disk space
While the alignment is running, view the script in a text editor.
	gview $top/bin/run-bwasw

* Browse the contig to reference alignments using samtools tview
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

* Browse the contig to reference alignments using IGV
Start IGV.
Select the reference "Human hg19".
File->Load from File... k48/bwasw/HS0674-contigs.bam
Go to the region chr3:186,600,000-187,600,000
IGV may take up to a minute to load.
What genes overlap the contigs?
	A: ST6GAL1, SST, RTP2 and BCL6
Add the dbSNP track.
File->Load from Server...  then expand "hg19", "Variation and Repeats"
and select "dbSNP 1.3.1".
Zoom in on a SNV. Is it in dbSNP? Is it coding?
Bonus: Find a coding SNV. What is its dbSNP rs ID?
	A: chr3:187,416,634 chr3:187,416,719 chr3:187,447,032

* View the contig to reference alignments SAM file
View the SAM file in a text editor.
	cd $top/k48/bwasw
	gview HS0674-contigs.sam
Disable line wrap.
	:set nowrap
The contig ID is given in the first column.
Which large contig has two alignments?
	A: The largest contig, 94.
The location on the reference is given in the third and fourth columns
and the orientation in the second column. 0 and 16 indicate positive
and negative orientation respectively.
What is the location and orientation of these two alignments?
	A: chr3:186,644,066 (-) and chr3:187,439,792 (+)
In IGV, go to the region chr3:186,600,000-187,600,000
and find these two alignments.
What large-scale structural rearrangement has occurred, and what is
its approximate size?
	A: a ~800 kbp inversion
Which two genes are fused as a result of this rearrangement?
	A: ST6GAL1 and BCL6

* Call variants of the reads-to-reference alignments using bcftools
Check that BWA has completed aligning the reads to the reference.
Run bcftools in this terminal.
	cd $top/bwa
	run-bcftools 2>&1 |tee bcftools.log
20 min, 250 MB RAM, 120 MB disk space
While bcftools is running, view the script in a text editor.
	gview $top/bin/run-bcftools

* Call variants of the contigs-to-reference alignments using bcftools
Run bcftools.
	cd $top/k48/bwasw
	run-bcftools-assembly 2>&1 |tee bcftools.log
Browse the variants using IGV.
File->Load from File... k48/bwasw/HS0674-contigs.var.vcf.gz
Right-click on the VCF track and select "Color by Allele".

* Call the effects of the SNVs
Run snpEff.
	cd $top/k48/bwasw
	run-snpeff HS0674-contigs.var.vcf.gz >HS0674-contigs.var.snpeff
1 min, 3.3 GB RAM
View the output of snpEff in a text editor.
	gview HS0674-contigs.var.snpeff
	set :nowrap
Count the number of SNVs in each category of effect.
	cut -f 1,2,16 HS0674-contigs.var.snpeff |sort -uk1,2 |cut -f3 \
		|sed 's/:.*//' |sort |uniq -c |sort -rn
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

* Compare the assembly variants to the read-alignment variants
Browse both variants using IGV.
Start IGV.
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

Supplementary material

* Align the reads to the contigs using BWA
Run BWA.
	cd $top
	mkdir k48/bwa
	cd k48/bwa
	bwa index HS0674-contigs.fa
	ref=../HS0674-contigs.fa run-bwa 2>&1 |tee bwa.log
12 min, 200 MB RAM, 3 GB disk space

Index the assembly FASTA file.
	samtools faidx HS0674-contigs.fa
This command fails with the following error:
	[fai_build_core] line length exceeds 65535 in sequence '94'.
How long is the longest line? (hint: use wc -L)
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
