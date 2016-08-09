# De novo assembly of Illumina reads using ABySS and alignment using BWA

This workshop is designed by [Shaun Jackman](http://sjackman.ca) [\@sjackman](https://twitter.com/sjackman).

# Purpose

We will use ABySS to assemble a 200 kbp bacterial artificial chromosome (BAC) using one lane of paired-end reads from the Illumina platform. BWA-MEM is used to align the assembled contigs to the human reference genome, and bcftools is used to call variants. IGV is used to visualize these alignments and variants. snpEff is used to determine the effects of these variants.

After this lab, you will have learned how to use ABySS to assemble a small genome, use BWA-MEM to align reads and contigs to a reference genome, use IGV to visualize these alignments, and use bcftools and snpEff to call variants and determine their effect.

# Contents

* [Getting Started](#getting-started)
* [Exercise 0: Index the reference using BWA](#exercise-0-index-the-reference-using-bwa)
* [Exercise 1: Align the reads to the reference using BWA (optional)](#exercise-1-align-the-reads-to-the-reference-using-bwa-optional)
* [Exercise 2: Inspect the reads](#exercise-2-inspect-the-reads)
* [Exercise 3: Assemble the reads into contigs using ABySS](#exercise-3-assemble-the-reads-into-contigs-using-abyss)
* [Exercise 4: Align the contigs to the reference using web BLAT](#exercise-4-align-the-contigs-to-the-reference-using-web-blat)
* [Exercise 5: Align the contigs to the reference using BWA-MEM](#exercise-5-align-the-contigs-to-the-reference-using-bwa-mem)
* [Exercise 6: Browse the contig to reference alignments using samtools tview](#exercise-6-browse-the-contig-to-reference-alignments-using-samtools-tview)
* [Exercise 7: Browse the contig to reference alignments using IGV](#exercise-7-browse-the-contig-to-reference-alignments-using-igv)
* [Exercise 8: View the contig to reference alignments SAM file](#exercise-8-view-the-contig-to-reference-alignments-sam-file)
* [Exercise 9: Call variants of the reads-to-reference alignments using bcftools (optional)](#exercise-9-call-variants-of-the-reads-to-reference-alignments-using-bcftools-optional)
* [Exercise 10: Call variants of the contigs-to-reference alignments using bcftoolss](#exercise-10-call-variants-of-the-contigs-to-reference-alignments-using-bcftools)
* [Exercise 11: Determine the effects of the SNVs](#exercise-11-determine-the-effects-of-the-snvs)
* [Exercise 12: Compare the assembly variants to the read-alignment variants (optional)](#exercise-12-compare-the-assembly-variants-to-the-read-alignment-variants-optional)
* [Exercise 13: Align the reads to the contigs using BWA (optional)](#exercise-13-align-the-reads-to-the-contigs-using-bwa-optional)

# Getting Started

## System requirements

The total run time of the tools alone is approximately 70 minutes on a 2-core 2 GHz system. 4 GB of RAM and 5 GB of disk space is required. The following software is required:

+ ABySS 1.9.0: genome sequence assembler for short reads
  http://www.bcgsc.ca/platform/bioinfo/software/abyss
+ bcftools 1.3.1: manipulate VCF/BCF variant call files
  http://www.htslib.org
+ BWA 0.7.15: align short reads and long contigs
  http://bio-bwa.sourceforge.net
+ IGV 2.3.80: visualize a genome
  http://www.broadinstitute.org/igv
+ Java 7 or later: execute Java programs
  http://www.java.com
+ samtools 1.3.1: manipulate SAM/BAM alignment files
  http://www.htslib.org
+ snpEff 4.2: determine the effect of SNVs
  http://snpeff.sourceforge.net

## Install the software on Mac OS

Install [Xcode](macappstores://itunes.apple.com/us/app/xcode/id497799835).

Open Xcode, select "Xcode -> Preferences -> Downloads -> Command Line Tools -> Install".

Install [Homebrew](http://brew.sh).

```sh
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Install GNU Coreutils and the text editor MacVim using Homebrew.

```sh
brew install coreutils
brew cask install macvim
```

Install ABySS, BWA, samtools and snpEff using Homebrew.

```sh
brew tap homebrew/science
brew install abyss bcftools bwa igv samtools snpeff
```

IGV may alternatively be installed using Java Web Start:

```sh
javaws http://www.broadinstitute.org/igv/projects/current/igv.jnlp
```

snpEff may alternatively be installed manually.

```sh
curl -LO http://downloads.sourceforge.net/project/snpeff/snpEff_v4_2_core.zip
unzip snpEff_v4_2_core.zip
curl -LO http://downloads.sourceforge.net/project/snpeff/databases/v4_2/snpEff_v4_2_GRCh38.82.zip
unzip -d snpEff snpEff_v4_2_GRCh38.82.zip
```

## Install the software on Ubuntu and Debian using `apt-get`

Install ABySS, bcftools, BWA, igv, and samtools.

```sh
sudo apt-get install abyss bcftools bwa igv samtools
sudo ln -s /usr/lib/abyss/abyss-fac /usr/local/bin/
```

Install snpEff using the instructions above for Mac OS.

## Set up the environment

Download the workshop scripts using git, if you have it installed.

```sh
git clone git://github.com/sjackman/abyss-activity.git ~/abyss
```

If you do not have git installed, use curl.

```sh
mkdir ~/abyss
cd ~/abyss
curl -L https://github.com/sjackman/abyss-activity/archive/master.tar.gz | tar -zx --strip 1
```

The shell script named `environment` will set environment variables that specify the location of the installed software. Each time that you open a new terminal, you will need to run the following two commands to change to the working directory and source the environment script:

```sh
cd ~/abyss
source environment
```

Check that the tools are installed in the PATH.

```sh
source environment
which abyss-fac abyss-pe bcftools bgzip bwa curl gview java samtools snpeff vcfutils.pl
```

Check that the workshop scripts are in the PATH.

```sh
which run-abyss run-bwa run-bwamem run-snpeff run-bcftools run-bcftools-assembly
```

Download the FASTQ files.

```sh
curl -LO ftp://ftp.bcgsc.ca/public/sjackman/30CJCAAXX_4_1.fq.gz
curl -LO ftp://ftp.bcgsc.ca/public/sjackman/30CJCAAXX_4_2.fq.gz
```

Download the reference file of Human chromosome 3.

```sh
curl -LO http://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr3.fa.gz
gunzip chr3.fa.gz
```

# Exercise 0: Index the reference using BWA

Index the reference file.

```sh
cd $top
bwa index chr3.fa
```

8 min, 1 GB RAM, 260 MB disk space

# Exercise 1: Align the reads to the reference using BWA (optional)

Run BWA.

```sh
cd $top
bwa mem -t2 chr3.fa 30CJCAAXX_4_1.fq.gz 30CJCAAXX_4_2.fq.gz >bwa.sam
```

40 min, 250 MB RAM, 3 GB disk space

While this job is running in the background, open a new terminal and continue with the next section. Remember that you have to set the environment. In the new terminal window type:

```sh
cd ~/abyss
source environment
```

# Exercise 2: Inspect the reads

There are two FASTQ files for this one lane of paired-end Illumina data: one file for the forward reads and one file for the reverse reads. Look at the first few reads.

```sh
cd $top
gunzip -c 30CJCAAXX_4_1.fq.gz | head
```

How long are the reads? (hint: use `awk '{print length}'`)

> ```sh
> gunzip -c 30CJCAAXX_4_1.fq.gz | awk '{print length}' | head -n2
> ```
> 50 bp

How many lines are there in both files? (hint: use `wc -l`)

> ```sh
> gunzip -c 30CJCAAXX_4_1.fq.gz 30CJCAAXX_4_2.fq.gz | wc -l
> ```
> 40,869,448 lines

How many lines per read?

> 4 lines per read

How many reads are there in both files?

> 40869448 / 4 = 10,217,362 reads

How many bases are sequenced?

> 10217362 * 50 = 510,868,100 bp

Assuming the BAC is 200 kbp, what is the depth of coverage?

> 510868100 / 200000 = 2554 fold coverage

# Exercise 3: Assemble the reads into contigs using ABySS

Run the assembly.

```sh
cd $top
mkdir k48
ln -s ../30CJCAAXX_4_1.fq.gz ../30CJCAAXX_4_2.fq.gz k48/
abyss-pe -C k48 name=HS0674 k=48 v=-v in="30CJCAAXX_4_1.fq.gz 30CJCAAXX_4_2.fq.gz" contigs 2>&1 | tee abyss.log
```

You may need to use the uncompressed FASTQ files on a Mac.

```sh
gunzip -k 30CJCAAXX_4_1.fq.gz 30CJCAAXX_4_2.fq.gz
mkdir k48
ln -s ../30CJCAAXX_4_1.fq ../30CJCAAXX_4_2.fq k48/
abyss-pe -C k48 name=HS0674 k=48 v=-v in="30CJCAAXX_4_1.fq 30CJCAAXX_4_2.fq" contigs 2>&1 | tee abyss.log
```

5 min, 200 MB RAM, 2 MB disk space

Look at the option `-n,--dry-run` of abyss-pe. Its output is the
commands that ABySS will run for the assembly.

```sh
abyss-pe name=HS0674 k=48 in="30CJCAAXX_4_1.fq.gz 30CJCAAXX_4_2.fq.gz" contigs -n
```

The assembly runs in three stages: assemble contigs without paired-end information, align the paired-end reads to the initial assembly, and merge contigs joined by paired-end information. You can instruct ABySS to stop after any of these stages. Use the `-n` option to see the commands for each stage.

```sh
abyss-pe name=HS0674 k=48 in="30CJCAAXX_4_1.fq.gz 30CJCAAXX_4_2.fq.gz" unitigs -n
abyss-pe name=HS0674 k=48 in="30CJCAAXX_4_1.fq.gz 30CJCAAXX_4_2.fq.gz" pe-sam -n
abyss-pe name=HS0674 k=48 in="30CJCAAXX_4_1.fq.gz 30CJCAAXX_4_2.fq.gz" contigs -n
```

Once the assembly has completed, view the contigs in a text editor.

```sh
gview k48/HS0674-contigs.fa
```

Disabling line wrap makes it easier to browse the file. For Emacs, select the option "Options -> Line Wrapping in this Buffer -> Truncate Long Lines." For Vim, select the option "Edit -> File Settings -> Toggle Line Wrap", or type `:set nowrap`

How many contigs are longer than 100 bp?

> 6

What is the length of the longest contig (hint: use `awk '{print length}'`)?

> ```sh
> awk '{print length}' k48/HS0674-contigs.fa | sort -n | tail -n1
> ```
> 67,526 bp

What is the N50 of the assembly?

```sh
abyss-fac k48/HS0674-contigs.fa
```

> ```
> n   n:500  L50  min   N80    N50    N20    E-size  max    sum     name
> 13  6      2    8044  54373  54743  67480  51662   67480  212330  k48/HS0674-contigs.fa
> ```
> 54,743 bp

View the assembly log in a text editor.

```sh
gview abyss.log
```

What portion of the reads align to the assembly? (hint: search for "Mapped")

> Mapped 7167472 of 10217362 reads (70.1%)

What is the median fragment size and standard deviation of this library? (hint: search for "median")

> median = 204 bp, sd = 18 bp

# Exercise 4: Align the contigs to the reference using web BLAT

Open BLAT in a web browser: <http://genome.ucsc.edu>

View the assembled contigs in a text editor.

```sh
gview k48/HS0674-contigs.fa
```

Disabling line wrap makes it easier to select the full sequence.

Select the two contigs whose lengths are approximately 8 and 16.5 kbp and copy-and-paste their sequence into BLAT.

What is the exact length of these two contigs?

> 8,044 bp and 16,561 bp

Click "browser" for the best alignment and then zoom out 10x. To which chromosome and band do these contigs align?

> chr3q27.3

What are the nearest two genes?

> SST and RTP2

Set the "Common SNPs" track in Variation to "pack".

A SNV between the assembled sequence and reference sequence is displayed with a red line. Zoom in on a SNV. Is it in dbSNP?

Zoom in on the gap between the two contigs. Set the "RepeatMasker" track in Repeats to "full". What feature overlaps the gap that likely caused the assembly gap?

> Simple_repeat (TA)n

Zoom in to see the sequence of the feature.

Zoom out to see the alignment of both contigs. Unaligned query sequence is shown with a ridiculously thin purple line at the end of the alignment. The one-pixel-wide purple line is absurdly difficult to see. Which contig has unaligned sequence at one end?

> the 16.5 kbp contig

Select that contig, and copy the unaligned sequence to the clipboard. Aligned sequence is shown in blue upper-case characters, and unaligned sequence is shown in black lower-case characters. Open BLAST in a web browser: <http://blast.ncbi.nlm.nih.gov>

Select "Nucleotide BLAST". Select the database "Nucleotide collection (nr/nt)". Paste the sequence into the query box. Click BLAST. Look at the first few best hits, and ignore hits to chimpanzee (*Pan troglodytes*) clones. What is the best alignment of our non-human contig fragment?

> Cloning vector pTARBAC

What is the cause of this chimeric contig?

> This assembled contig contains both the BAC cloning vector and the human DNA insert.

# Exercise 5: Align the contigs to the reference using BWA-MEM

Warning: do not start BWA-MEM until BWA has completed unless your machine has at least 2 GB of RAM.

Run BWA-MEM.

```sh
bwa mem -t2 chr3.fa k48/HS0674-contigs.fa >k48/HS0674-contigs.sam
samtools sort -o k48/HS0674-contigs.bam k48/HS0674-contigs.sam
samtools index k48/HS0674-contigs.bam
```

1 min, 800 MB RAM, 1 MB disk space

# Exercise 6: Browse the contig to reference alignments using samtools tview

Run samtools tview.

```sh
samtools tview k48/HS0674-contigs.bam chr3.fa
```

Go to the region `chr3:186,931,150`

```
g chr3:186,931,150
```

What variants do you see?

> a T/G SNV and a 24-bp insertion

Go to the region `chr3:186,958,940`

```
g chr3:186,958,940
```

Notice the Ns in the contig sequence, indicating a scaffold gap. How many Ns are in the contig?

> 19 Ns

How many Ns should there be for the size of the gap to agree with the reference?

> 16 Ns

# Exercise 7: Browse the contig to reference alignments using IGV

Start IGV. Open the "IGV 2.3" icon on your desktop, or run either `igv` or `javaws http://www.broadinstitute.org/igv/projects/current/igv.jnlp` at the command line.

Select "View -> Preferences -> Alignments" and change "Visibility range threshold (kb)" to 1000. Select the "Genomes -> Load Genome From Server -> Human hg38". Then select "File -> Load from File" and "k48/HS0674-contigs.bam" Go to the region `chr3:186,500,000-188,000,000` by entering it into the box labeled "Go". IGV may take up to a minute to load.

Which four genes overlap the contigs?

> ST6GAL1, SST, RTP2 and BCL6

Add the dbSNP track. Select "File -> Load from Server" then expand "Annotations", "Variation and Repeats" and select "dbSNP 142". Alternatively, select "File -> Load from URL" and enter `ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz`, or download this file and the associated `.gz.tbi` file using `curl` and load it using "File -> Load from File".

Zoom in on a SNV. Is it in dbSNP? Is it coding?

Bonus: Find a coding SNV. What is its dbSNP rs ID?

> rs1973791 (chr3:187,698,846) and rs11707167 (chr3:187,698,931) in the last exon of RTP2.
> rs1056932 (chr3:187,729,244) in exon 5 of BCL6.

# Exercise 8: View the contig to reference alignments SAM file

View the SAM file in a text editor. Disable line wrap.

```sh
gview k48/HS0674-contigs.sam
```

The contig ID is given in the first column, and the position of the contig on the reference is given in the third and fourth columns. Which large contig has two alignments, and what are the positions of these two alignments?

> The contig that has alignments starting at
> chr3:186,698,393 and chr3:187,439,792.

The orientation is given in the second column. The numbers 0 and 2048 both indicate positive orientation, and 16 indicates negative orientation. What is the position and orientation of these two alignments?

> chr3:186,980,605 (-) and chr3:187,722,004 (+)

In IGV, go to the region `chr3:186,500,000-188,000,000`. The alignment highlighted in blue indicates that the alignment of this contig to the reference genome is split in two, suggesting a misassembly or a structural rearrangement. To find the other half of the alignment, hover over the blue alignment and note its "Read name". Right click on that alignment, select "Select by name", and enter the "Read name" that you noted previously. Both alignments of this one contig are now highlighted in red.

What large-scale structural rearrangement has occurred, and what is its approximate size?

> a ~700 kbp inversion

Which two genes are fused as a result of this rearrangement?

> ST6GAL1 and BCL6

# Exercise 9: Call variants of the reads-to-reference alignments using bcftools (optional)

Check that BWA has completed aligning the reads to the reference, or run BWA now.

```sh
bwa mem -t2 chr3.fa 30CJCAAXX_4_1.fq.gz 30CJCAAXX_4_2.fq.gz >30CJCAAXX_4.sam
```

Call variants using bcftools.

```sh
samtools sort -@2 -o 30CJCAAXX_4.bam 30CJCAAXX_4.sam
samtools index 30CJCAAXX_4.bam
samtools mpileup -u -d 9999 -L 9999 -f chr3.fa 30CJCAAXX_4.bam >30CJCAAXX_4.bcf
bcftools view 30CJCAAXX_4.bcf >30CJCAAXX_4.var.vcf
vcfutils.pl varFilter -D9999 30CJCAAXX_4.var.vcf >30CJCAAXX_4.varFilter.vcf
bgzip 30CJCAAXX_4.varFilter.vcf
bcftools index 30CJCAAXX_4.varFilter.vcf.gz
```

20 min, 250 MB RAM, 120 MB disk space

# Exercise 10: Call variants of the contigs-to-reference alignments using bcftools

Run bcftools.

```sh
cd $top/k48/bwamem
run-bcftools-assembly 2>&1 |tee bcftools.log
```

Browse the variants using IGV. Select "File->Load from File... k48/bwamem/HS0674-contigs.var.vcf.gz" Right-click on the VCF track and select "Color by Allele".

# Exercise 11: Determine the effects of the SNVs

Run snpEff.

```sh
cd $top/k48/bwamem
run-snpeff HS0674-contigs.var.vcf.gz >HS0674-contigs.var.snpeff
```

1 min, 2.8 GB RAM

View the output of snpEff in a text editor.

```sh
gview HS0674-contigs.var.snpeff
```

Count the number of SNVs in each category of effect.

```sh
cut -f1,2,16 HS0674-contigs.var.snpeff |cut -d: -f1 |sort -u |cut -f3 |sort |uniq -c
```

Find all the coding SNVs.

```sh
grep CODING HS0674-contigs.var.snpeff |cut -f1-4,16-18 |uniq
```

Find the non-synonymous SNV. What is its location?

> chr3:187,416,719

What is its dbSNP rs ID? (hint: use IGV or the UCSC genome browser)

> rs11707167

Open dbSNP in a web browser: <http://www.ncbi.nlm.nih.gov/projects/SNP/>

What is the minor allele frequency (MAF) of this SNP?

> T=0.4716

# Exercise 12: Compare the assembly variants to the read-alignment variants (optional)

Browse the variants called by both methods using IGV.

Start IGV.

```sh
javaws http://www.broadinstitute.org/igv/projects/current/igv.jnlp
```

Select:

+ File -> Load from File... `k48/bwamem/HS0674-contigs.bam`
+ File -> Load from File... `k48/bwamem/HS0674-contigs.var.vcf.gz`
+ File -> Load from File... `bwa/HS0674.var.vcf.gz`

Go to the region chr3:186,648,960

Is the SNV called by both methods?

> Yes.

Is the insertion called by both methods?

> No, the insertion is called only by the assembly.

Hover the mouse cursor over the insertion to see the inserted sequence. What is the inserted sequence?

> TCTGGGTTCCTTCAAATCCTGCCT

Compare the inserted sequence to the reference sequence at this location. What sequence do the inserted sequence and the reference sequence have in common?

> The insertion duplicates this 8-bp sequence: TCTGGGTT

Load the alignments of the reads to the reference.

File -> Load from File... `bwa/30CJCAAXX_4.bam`

The aligned reads do not show the insertion, but the alignments that span the insertion have mismatches at the end of the alignment. Do those mismatched bases agree with the inserted sequence?

> Yes, the mismatched bases show CT at the 5' end and CC
> at the 3' end, which agree with the inserted sequence.

# Exercise 13: Align the reads to the contigs using BWA (optional)

Run BWA.

```sh
cd $top
mkdir k48/bwa
cd k48/bwa
bwa index ../HS0674-contigs.fa
ref=../HS0674-contigs.fa run-bwa 2>&1 |tee bwa.log
```

12 min, 200 MB RAM, 3 GB disk space

Index the assembly FASTA file.

```sh
samtools faidx ../HS0674-contigs.fa
```

Browse the BAM file using samtools tview.

```sh
samtools tview 30CJCAAXX_4.bam ../HS0674-contigs.fa
```

Browse the BAM file using IGV.

+ File -> Import Genome... `k48/HS0674-contigs.fa`
+ File -> Load from File... `k48/bwa/30CJCAAXX_4.bam`

Find a scaffold gap on the largest contig. Find the two contigs that have consistent mate pairs joining them.
