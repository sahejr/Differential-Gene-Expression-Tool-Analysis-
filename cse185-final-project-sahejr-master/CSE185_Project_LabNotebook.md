# Lab notebook for CSE185 Final Project

## Session 1 (May 22, 2018)

Databases explicitly considered before searching general paper topics:
http://www.ebi.ac.uk/gwas/ (redirected to from https://www.genome.gov/gwastudies/)
SRA

How-to GWAS stats and analysis papers: https://www.ncbi.nlm.nih.gov/pubmed/16983374 and https://www.ncbi.nlm.nih.gov/pubmed/21962507

Potential Datasets (all accession numbers found and data accessed):
1. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4121804/ (Accession  GSE53081) but different method of differential gene expression; problem is data is probably too big (downloaded 2 files and they are each like 800MB; could drop columns but that is a lot of preprocessing)
2. http://aac.asm.org/content/60/9/5515.full (Accession  PRJNA242614) find antibiotic resistant mutations and do Provean score on them
3. http://www.ebi.ac.uk/gwas/search?query=Jun%20GR (Accession GCST004249) GWAS analysis; add statistical test or see if model development possible

## Session 2 (May 24, 2018)

For future reference, good explanation of differential gene expression analysis pipeline: https://www.ebi.ac.uk/training/online/course/functional-genomics-ii-common-technologies-and-data-analysis-methods/differential-gene

Dataset chosen: Bioproject accession number PRJNA253971, accompanies the paper http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0105004#s3

Question: "how do different tools for differential gene expression analysis differ when run on the same input dataset?"

Secondary question that is answered in the process: what are the differentially expressed genes found via this analysis?

## Session 3+4 (made up both on May 31, 2018)

Made directory called `project` for all files

1. Download reference genome for chromosome 1 (GenBank id CM000937.1): `efetch -db nucleotide -id CM000937.1 -format fasta > CM000937.1.fasta`; Used `head CM000937.1.fasta` to read file header and ensure correct file was downloaded.

2. Constructed index using bowtie:
`bowtie-build ./CM000937.1.fasta AC_chr1`

3. Download reads data

Used reads mapping to chromosome 1 as downsampling method; Data is paired end reads:

A. installed seqtk to be able to download random sampling of reads (https://github.com/lh3/seqtk)
```shell
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
```
Above commands were taken from README at above github link and successfully created executable `seqtk` (with warning that variable `lc` in seqtk.c was set but unused, which seems acceptable given several tutorials I researched showing how to successfully install this tool have output with this warning present).

B. Downloaded datasets one at a time. run seqtk on for downsampling (using seed value of 100), and discard original; code produced below and available at script prepReads.sh

```shell
fastq-dump --split-files -O . SRX642055 --gzip
Read 25287520 spots for SRX642055
Written 25287520 spots for SRX642055
```

Note that because of paired reads, there are 2 reads per spot.

Created files SRX642055_1.fastq and SRX642055_2.fastq; accidently unzipped first one when used gunzip instead of zcat or gunzip -c, which unzipped the file and fed nothing to the command after the pipe used. Adjusted and used the following commands:
```shell
../seqtk/seqtk sample -s 100 SRX642055_1.fastq 0.2 > SRX642055_1.fq
../seqtk/seqtk sample -s 100 <(zcat SRX642055_2.fastq) 0.2 > SRX642055_2.fq   
```

Seed value of 100 used for both so in hopes of all read pair mates both being sampled. 0.2 is the fraction of reads so that file sizes are each ~500MB.

Obtain number of reads and unique read lengths:
```shell
wc -l SRX642055_*
  20237796 SRX642055_1.fq
  20237796 SRX642055_2.fq
  40475592 total
cat SRX642055_1.fq | awk 'NR%4==0 {print length}' | sort -n | uniq -c
  5059449 39
cat SRX642055_2.fq | awk 'NR%4==0 {print length}' | sort -n | uniq -c
  5059449 39
```


Articles found and read:

https://www.nature.com/articles/ncomms10033#methods

https://www.ebi.ac.uk/training/online/course/functional-genomics-ii-common-technologies-and-data-analysis-methods/differential-gene

Reference genome from paper: https://www.nature.com/articles/nature10390

## Saturday June 2, 2018

Dividing by 4, that results in 5,059,449 paired reads, all of length 39.

4. Repeat process for SRA SRX642051 but now choosing to have the same number of reads to control for coverage. Using same seed should also allow for same read pairs being available.
```shell
fastq-dump --split-files -O . SRX642051 --gzip
Read 29725407 spots for SRX642051
Written 29725407 spots for SRX642051

../seqtk/seqtk sample -s 100 <(zcat SRX642051_1.fastq) 5059449 > SRX642051_1.fq   
../seqtk/seqtk sample -s 100 <(zcat SRX642051_2.fastq) 5059449 > SRX642051_2.fq   

wc -l SRX642051_*
```

5. Fastqc Analysis:
made directory reads/fastqcOutput and cd'ed into it before running following commands:
```shell
fastqc -o . ../SRX642051_1.fq ../SRX642051_2.fq
fastqc -o . ../SRX642055_1.fq ../SRX642055_2.fq
```
Resulting html files revealed passing label for main tests (the SRX642055_2.fq had failure for overepresented sequences of all Ns and all Ts). Attempted trimming with sickle pe but no signficant improvement was found


6. Download cufflinks 2.2.1
```shell
curl http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz > cufflinks_2.2.1.tar.gz
tar -xzf cufflinks_2.2.1.tar.gz
rm cufflinks_2.2.1.tar.gz 
mv cufflinks-2.2.1.Linux_x86_64/ cufflinksLib/
```


7. Had difficulty getting tophat to work and thus used HISAT, which the creators of both recommend using instead now for speed and accuracy purposes. Then, I aligned RNA-seq reads to reference genome with HISAT

Build index from chr 1 reference genome fasta file; created 8 files (this makes bowtie reference irrelevant)
```shell
./hisat2-2.1.0/hisat2-build CM000937.1.fasta AC_chr1
```


Align reads from SRX642051 to index (and time this process)

```shell
time ./hisat2-2.1.0/hisat2 --dta-cufflink -x AC_chr1 -1 reads/SRX642051_1.fq -2 reads/SRX642051_2.fq -q -S SRX642051.sam

      78304 (1.68%) aligned discordantly 1 time
    ----
    4574672 pairs aligned 0 times concordantly or discordantly; of these:
      9149344 mates make up the pairs; of these:
        8844673 (96.67%) aligned 0 times
        223854 (2.45%) aligned exactly 1 time
        80817 (0.88%) aligned >1 times
12.59% overall alignment rate

real    4m36.290s
user    4m32.120s
sys     0m3.745s
```

Align reads from SRX642055 to index 
```shell
time ./hisat2-2.1.0/hisat2 --dta-cufflink -x AC_chr1 -1 reads/SRX642055_1.fq -2 reads/SRX642055_2.fq -q -S SRX642055.sam
5059449 reads; of these:
  5059449 (100.00%) were paired; of these:
    4625613 (91.43%) aligned concordantly 0 times
    388229 (7.67%) aligned concordantly exactly 1 time
    45607 (0.90%) aligned concordantly >1 times
    ----
    4625613 pairs aligned concordantly 0 times; of these:
      77277 (1.67%) aligned discordantly 1 time
    ----
    4548336 pairs aligned 0 times concordantly or discordantly; of these:
      9096672 mates make up the pairs; of these:
        8763630 (96.34%) aligned 0 times
        225746 (2.48%) aligned exactly 1 time
        107296 (1.18%) aligned >1 times
13.39% overall alignment rate

real	4m18.409s
user	4m13.382s
sys	0m3.403s
```


Also, I need to get gene annotations
Reannotation by same group as main paper was done in this study (https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-49)
"Accession numbers for non strand-specific RNA-Seq and transcript assemblies:...regenerating epithelial tail tip [SRA: SRX158076, TSA: SUB139331], regenerating tail base [SRA: SRX158077, TSA: SUB139332], tail [SRA: SRX158074, TSA: SUB139330]."

## Tuesday June 5, 2018

8. Prep kallisto arm of analysis

Downloaded transcriptome of regenerating tail epithelium here: https://www.ncbi.nlm.nih.gov/Traces/wgs/GACT01?display=contigs&page=1; though 6 major chromosomes are present, there were only 5 files. I downloaded the first one.

Did not use link https://www.ncbi.nlm.nih.gov/Traces/wgs/GAGC01?display=contigs&page=1, which did possess 6 files. Noting this in case I need to go back.

Indexing transcriptome file, generating transcripts.idx index file:
```shell
time kallisto index -i transcripts.idx GACT01.1.fsa_nt.gz 

[build] loading fasta file GACT01.1.fsa_nt.gz
[build] k-mer length: 31
[build] warning: clipped off poly-A tail (longer than 10)
        from 263 target sequences
[build] counting k-mers ... done.
[build] building target de Bruijn graph ...  done   
[build] creating equivalence classes ...  done
[build] target de Bruijn graph has 359906 contigs and contains 89190386 k-mers 


real	7m48.092s
user	7m15.735s
sys	0m8.346s
```

7. Run kallisto quant to perform pseudoalignment with k-mer based approach and quantify gene expression

First, processed gene annotation from GFF3 file format into GTF (equivalent to GFF2) for use by kallisto. Used the gffread tool that comes with CuffLinks library
```shell
./cufflinksLib/gffread AC_annotations.gff -T -o annot.gtf
```

Modified run_kallisto.sh script from lab4 (changed paths, kept thread and boostrap values).
```shell
#!/bin/bash

# Set variables used in each kallisto run
PROJECT=/home/linux/ieng6/cs185s/sdrandha/project
GTF=${PROJECT}/annot.gtf
KINDEX=${PROJECT}/transcripts.idx

# Do a separate kallisto run for each dataset
for prefix in SRX642051 SRX642055
do
    mkdir -p ${PROJECT}/${prefix}
    kallisto quant -t 3 -b 100 \
	-o ${PROJECT}/$prefix --gtf $GTF -i $KINDEX \
	${PROJECT}/reads/${prefix}_1.fq ${PROJECT}/reads/${prefix}_2.fq
done
```

Executed script as `time ./run_kallisto.sh` to time it. Hit errors like `terminate called after throwing an instance of 'std::length_error'` and research was not clear on solution so I rebuilt index.

Downloaded transcriptome data GCF_000090745.1_AnoCar2.0_rna.fna.gz from official NCBI entry (https://www.ncbi.nlm.nih.gov/genome/?term=txid28377[Organism:noexp]) in fna.gz filetype, which stands for FASTA containing nucleotide data and is zipped. Repeat kallisto indexing command (documentation indicates you can pass in the .gz file directly):
```shell
time kallisto index -i transcriptome.idx AC_RNA.fna.gz

[build] loading fasta file AC_RNA.fna.gz
[build] k-mer length: 31
[build] warning: clipped off poly-A tail (longer than 10)
        from 84 target sequences
[build] warning: replaced 542 non-ACGUT characters in the input sequence
        with pseudorandom nucleotides
[build] counting k-mers ... done.
[build] building target de Bruijn graph ...  done 
[build] creating equivalence classes ...  done
[build] target de Bruijn graph has 159033 contigs and contains 61736661 k-mers 


real	8m47.390s
user	7m59.374s
sys	0m8.423s
```

Altered run_kallisto.sh to point to this new transcriptome.idx file and ran successfully this time
```shell
time ./projectRepo/scripts/run_kallisto.sh

[quant] fragment length distribution will be estimated from the data
[index] k-mer length: 31
[index] number of targets: 38,095
[index] number of k-mers: 61,736,661
[index] number of equivalence classes: 83,563
[quant] running in paired-end mode
[quant] will process pair 1: /home/linux/ieng6/cs185s/sdrandha/project/reads/SRX642051_1.fq
                             /home/linux/ieng6/cs185s/sdrandha/project/reads/SRX642051_2.fq
[quant] finding pseudoalignments for the reads ... done
[quant] processed 5,059,449 reads, 1,871,573 reads pseudoaligned
[quant] estimated average fragment length: 112.839
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 1,008 rounds
[bstrp] number of EM bootstraps complete: 100


[quant] fragment length distribution will be estimated from the data
[index] k-mer length: 31
[index] number of targets: 38,095
[index] number of k-mers: 61,736,661
[index] number of equivalence classes: 83,563
[quant] running in paired-end mode
[quant] will process pair 1: /home/linux/ieng6/cs185s/sdrandha/project/reads/SRX642055_1.fq
                             /home/linux/ieng6/cs185s/sdrandha/project/reads/SRX642055_2.fq
[quant] finding pseudoalignments for the reads ... done
[quant] processed 5,059,449 reads, 1,587,247 reads pseudoaligned
[quant] estimated average fragment length: 123.255
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 958 rounds
[bstrp] number of EM bootstraps complete: 100


real	3m18.432s
user	7m59.642s
sys	0m10.345s
```

8. Use sleuth to perform differential gene expression analysis on kallisto output 

Altered my lab4 sleuth script to update directories and filenames. 

New exp_info.txt file contains:
```
sample condition
SRX642051 S1
SRX642055 S5
```

Ran from project directory:
```shell
time Rscript projectRepo/scripts/sleuthResults.R
reading in kallisto results
dropping unused factor levels
..
normalizing est_counts
17305 targets passed the filter
normalizing tpm
merging in metadata
summarizing bootstraps
..
fitting measurement error models
shrinkage estimation
Error in if (sum(valid) == 0) { : missing value where TRUE/FALSE needed
Calls: sleuth_fit ... do.grouped_df -> overscope_eval_next -> get_quantile
Execution halted
```

Error is explained here (https://github.com/pachterlab/sleuth/issues/125), indicating that a lack of replicates makes sleuth analysis impossible because it requires a measure of variance of samples within condition. Based on reading forum posts troubleshooting this, potential solutions include adding more replicates or "creating" a replicate by duplicating my single replicate. I tried to add more replicates submitted by the paper but the downsampling made the similarlity between samples disappear and made the variance extremely high and resulted in no differentially expressed genes. This also placed a great strain regarding data quota.

Thus I turned to the duplication approach. The variance among these replicates will be extremely low but that will be true, still allowing for a comparison of approaches. This resulted in me duplicating the kallisto outputs per sample, yielding directories (SRX642051_1 and SRX642051_2) and (SRX642055_1 and SRX642055_2)

```shell
time Rscript projectRepo/scripts/sleuthResults.R
reading in kallisto results
dropping unused factor levels
....
normalizing est_counts
20273 targets passed the filter
normalizing tpm
merging in metadata
summarizing bootstraps
....
fitting measurement error models
shrinkage estimation
computing variance of betas
fitting measurement error models
shrinkage estimation
computing variance of betas

real	1m20.238s
user	1m44.794s
sys	0m2.980s
```

Running wc -l on sleuth_results.tab gave a line count of 15,543; subtracting 1 for the header, 15,542 is the number of significantly differentially expressed transcripts (not genes). The first one pertained to myosin heavy chain cardiac muscle beta transcript variant X2. This value seemed high and made me realize that this is for the entire genome and not just chromosome 1, which I only was considering based on my TA's suggestion because of data usage and it would speed up the alignment part of my project.

## Wednesday June 6, 2018

Going back to the alignment approach:

Filter out unmapped reads and built transcriptome with cufflinks

In samBam directory, performed filtering and sorting based on script from cufflinks tutorial:
```shell

time samtools view -F 4 SRX642051.sam | sort -k 3,3 -k 4,4n >> SRX642051_mapped.sam
real	0m23.058s
user	0m17.646s
sys	0m3.877s

head -3 SRX642055.sam > SRX642055_mapped.sam
time samtools view -F 4 SRX642055.sam | sort -k 3,3 -k 4,4n >> SRX642055_mapped.sam
real	0m25.256s
user	0m18.807s
sys	0m4.442s

wc -l *_mapped.sam
  1642750 SRX642051_mapped.sam
  1884047 SRX642055_mapped.sam
  3526797 total
```

Now the mapped.sam files contain only mapped reads and are sorted, with the previous files' headers

```shell
time ./cufflinksLib/cufflinks -g annotEdit.gtf -N -u SRX642051_mapped.sam
```

samtools view -S -b | samtools sort > roommate.bam

## Friday, Saturday, Sunday June 8-10, 2018

10. Use cufflinks, cuffmerge, cuffdiff
Cufflinks to generate transcriptome of S1 reads and S5 reads each:
```shell
time ./cufflinksLib/cufflinks -g annotEdit.gtf -N -u samBam/SRX642051_mapped.sam
[bam_header_read] EOF marker is absent. The input is probably truncated.
[bam_header_read] invalid BAM binary header (this is not a BAM file).
File samBam/SRX642051_mapped.sam doesn't appear to be a valid BAM file, trying SAM...
[17:52:37] Loading reference annotation.
[17:52:38] Inspecting reads and determining fragment length distribution.
> Processed 144088 loci.                       [*************************] 100%
> Map Properties:
>	Normalized Map Mass: 660094.97
>	Raw Map Mass: 660094.97
>	Number of Multi-Reads: 54134 (with 230802 total hits)
>	Fragment Length Distribution: Empirical (learned)
>	              Estimated Mean: 127.53
>	           Estimated Std Dev: 97.17
[17:53:25] Assembling transcripts and initializing abundances for multi-read correction.
> Processed 144088 loci.                       [*************************] 100%
[18:04:33] Loading reference annotation.
[18:04:33] Re-estimating abundances with multi-read correction.
> Processed 3523 loci.                         [*************************] 100%

real	13m13.879s
user	4m13.410s
sys	0m17.855s

time ./cufflinksLib/cufflinks -g annotEdit.gtf -N -u samBam/SRX642055_mapped.sam
[18:10:47] Loading reference annotation.
[18:10:47] Inspecting reads and determining fragment length distribution.
> Processed 119093 loci.                       [*************************] 100%
> Map Properties:
>	Normalized Map Mass: 721012.25
>	Raw Map Mass: 721012.25
>	Number of Multi-Reads: 75035 (with 319639 total hits)
>	Fragment Length Distribution: Empirical (learned)
>	              Estimated Mean: 130.74
>	           Estimated Std Dev: 95.08
[18:11:09] Assembling transcripts and initializing abundances for multi-read correction.
> Processed 119093 loci.                       [*************************] 100%
[18:23:14] Loading reference annotation.
[18:23:14] Re-estimating abundances with multi-read correction.
> Processed 3299 loci.                         [*************************] 100%

real	13m24.856s
user	5m43.736s
sys	0m18.882s
```

Put generated 3 files in sample1_links and sample2_links directories
Combine these with cuffmerge into one transcriptome
```shell
time ./cufflinksLib/cuffmerge -s NC_014776.1.fasta -p 3 assembly_GTF_list.txt

[Sun Jun 10 18:40:26 2018] Beginning transcriptome assembly merge
-------------------------------------------

[Sun Jun 10 18:40:26 2018] Preparing output location ./merged_asm/
Warning: no reference GTF provided!
[Sun Jun 10 18:40:26 2018] Converting GTF files to SAM
[18:40:26] Loading reference annotation.
[18:40:26] Loading reference annotation.
[Sun Jun 10 18:40:26 2018] Assembling transcripts
Warning: Could not connect to update server to verify current version. Please check at the Cufflinks website (http://cufflinks.cbcb.umd.edu).
Command line:
cufflinks -o ./merged_asm/ -F 0.05 -q --overhang-tolerance 200 --library-type=transfrags -A 0.0 --min-frags-per-transfrag 0 --no-5-extend -p 3 ./merged_asm/tmp/mergeSam_file28B3a9 
[bam_header_read] EOF marker is absent. The input is probably truncated.
[bam_header_read] invalid BAM binary header (this is not a BAM file).
File ./merged_asm/tmp/mergeSam_file28B3a9 doesn't appear to be a valid BAM file, trying SAM...
[18:40:27] Inspecting reads and determining fragment length distribution.
Processed 3781 loci.                        
> Map Properties:
>	Normalized Map Mass: 8641.00
>	Raw Map Mass: 8641.00
>	Fragment Length Distribution: Truncated Gaussian (default)
>	              Default Mean: 200
>	           Default Std Dev: 80
[18:40:27] Assembling transcripts and estimating abundances.
Processed 3781 loci.                        
[Sun Jun 10 18:40:33 2018] Comparing against reference file None
Warning: Could not connect to update server to verify current version. Please check at the Cufflinks website (http://cufflinks.cbcb.umd.edu).
[Sun Jun 10 18:40:45 2018] Comparing against reference file None
Warning: Could not connect to update server to verify current version. Please check at the Cufflinks website (http://cufflinks.cbcb.umd.edu).

real	0m30.385s
user	0m12.927s
sys	0m13.059s

```

mergedOld.gtf is with NC reference and it complained with warnings
mergedNew.gtf is with CM reference and it didn't throw warnings
reference gtf did not work

Now with merged transcriptome, to run the diff expr analysis with cuffdiff
```shell
time ./cufflinksLib/cuffdiff -p 3 -o cuffDiffOut -L S1,S5 mergedNew.gtf samBam/SRX642051_mapped.sam samBam/SRX642055_mapped.sam
[18:45:37] Loading reference annotation.
Warning: No conditions are replicated, switching to 'blind' dispersion method
[18:45:37] Inspecting maps and determining fragment length distributions.
[18:46:17] Modeling fragment count overdispersion.
> Map Properties:
>	Normalized Map Mass: 299474.50
>	Raw Map Mass: 307168.75
>	Fragment Length Distribution: Empirical (learned)
>	              Estimated Mean: 104.39
>	           Estimated Std Dev: 79.97
> Map Properties:
>	Normalized Map Mass: 299474.50
>	Raw Map Mass: 292439.72
>	Fragment Length Distribution: Empirical (learned)
>	              Estimated Mean: 110.01
>	           Estimated Std Dev: 78.58
[18:47:52] Calculating preliminary abundance estimates
[18:47:52] Testing for differential expression and regulation in locus.
> Processed 3787 loci.                         [*************************] 100%
Performed 5258 isoform-level transcription difference tests
Performed 5078 tss-level transcription difference tests
Performed 4593 gene-level transcription difference tests
Performed 0 CDS-level transcription difference tests
Performed 0 splicing tests
Performed 0 promoter preference tests
Performing 0 relative CDS output tests
Writing isoform-level FPKM tracking
Writing TSS group-level FPKM tracking
Writing gene-level FPKM tracking
Writing CDS-level FPKM tracking
Writing isoform-level count tracking
Writing TSS group-level count tracking
Writing gene-level count tracking
Writing CDS-level count tracking
Writing isoform-level read group tracking
Writing TSS group-level read group tracking
Writing gene-level read group tracking
Writing CDS-level read group tracking
Writing read group info
Writing run info

real	5m1.061s
user	8m41.740s
sys	0m5.126s

grep "yes" -c isoform_exp.diff 
42

grep "yes" -c gene_exp.diff 
44
```

Now to redo kallisto for only chr1, will do Transcriptome -> alignment to chr 1 ref genome -> turn those alignments into fasta -> kallisto index -> kallisto quant -> sleuth

```shell
time ./hisat2-2.1.0/hisat2 -x AC_chr1 -U (<zcat AC_RNA.fna.gz) -f -p 3 -S txnAlign.sam
head -3 txnAlign.sam > txnAlign_mapped.sam
samtools view -F 4 txnAlign.sam | sort -k 3,3 -k 4,4n >> txnAlign_mapped.sam
samtools view -S -b txnAlign_mapped.sam > txnAlign.bam
samtools fasta --threads 3 txnAlign.bam > chr1Reads.fasta
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 3930 reads
rm txnAlign_mapped.sam
```

Now make new kallisto index comprising of just these transcripts
```shell
time kallisto index -i chr1Txn.idx chr1Reads.fasta

[build] loading fasta file chr1Reads.fasta
[build] k-mer length: 31
[build] warning: clipped off poly-A tail (longer than 10)
        from 1 target sequences
[build] warning: replaced 46 non-ACGUT characters in the input sequence
        with pseudorandom nucleotides
[build] counting k-mers ... done.
[build] building target de Bruijn graph ...  done 
[build] creating equivalence classes ...  done
[build] target de Bruijn graph has 10778 contigs and contains 6600652 k-mers 


real	0m26.887s
user	0m25.473s
sys	0m0.993s
```

Then ran kallisto quant on our reads file with run_kallisto.sh
```shell
[quant] fragment length distribution will be estimated from the data
[index] k-mer length: 31
[index] number of targets: 3,930
[index] number of k-mers: 6,600,652
[index] number of equivalence classes: 7,369
[quant] running in paired-end mode
[quant] will process pair 1: /home/linux/ieng6/cs185s/sdrandha/project/reads/SRX642051_1.fq
                             /home/linux/ieng6/cs185s/sdrandha/project/reads/SRX642051_2.fq
[quant] finding pseudoalignments for the reads ... done
[quant] processed 5,059,449 reads, 246,994 reads pseudoaligned
[quant] estimated average fragment length: 117.331
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 1,036 rounds
[bstrp] number of EM bootstraps complete: 100


[quant] fragment length distribution will be estimated from the data
[index] k-mer length: 31
[index] number of targets: 3,930
[index] number of k-mers: 6,600,652
[index] number of equivalence classes: 7,369
[quant] running in paired-end mode
[quant] will process pair 1: /home/linux/ieng6/cs185s/sdrandha/project/reads/SRX642055_1.fq
                             /home/linux/ieng6/cs185s/sdrandha/project/reads/SRX642055_2.fq
[quant] finding pseudoalignments for the reads ... done
[quant] processed 5,059,449 reads, 194,419 reads pseudoaligned
[quant] estimated average fragment length: 127.55
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 876 rounds
[bstrp] number of EM bootstraps complete: 100


real	0m40.423s
user	1m46.418s
sys	0m2.987s
```

Finally, ran sleuth to do diff expr analysis on these results
```shell

wc -l final_sleuth_results.tab 
1622 final_sleuth_results.tab
```

Subtracting 1 for header, that means 1622 transcripts were found as significantly differentially expressed.

#### Random Notes:
Why the anole is a good model species: https://www.ncbi.nlm.nih.gov/genome/?term=AAWZ00000000.2
WHY I USED HISAT2: https://www.biostars.org/p/289418/

other project articles in case of switch needed:
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0035052#pone.0035052.s002
https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1284-z#Abs1
https://www.sciencedirect.com/science/article/pii/S1879625711000629

kallisto vs alignment based methods:
https://www.biostars.org/p/251666/
https://www.nature.com/articles/s41467-017-00050-4
https://cgatoxford.wordpress.com/2016/08/17/why-you-should-stop-using-featurecounts-htseq-or-cufflinks2-and-start-using-kallisto-salmon-or-sailfish/
