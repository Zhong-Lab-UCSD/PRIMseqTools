# PRIMseqTools

## Overview
Protein-RNA Interaction Sequencing (PRIM-seq) is a high-throughput sequencing technology that efficiently maps cell wide RNA-protein associations in vitro. Here, we distribute PRIMseqTools, a standardized data processing pipeline to identify RNA-protein associations from fastq files of PIM-seq experiments.<br />

## Workflow
![alt text](https://github.com/Zhong-Lab-UCSD/PRIMseqTools/blob/main/workflow.png)
1. Raw read pairs from the PRIM-seq experiment are present in `.fastq` files.
2. Cutadpt is applied to remove 3' linker sequences and 5' adapter sequences from the read pairs (Adapter trimming). 
3. Fastp is then applied to remove low-quality reads whose mean quality is lower than Q20 and too short reads whose length is shorter than 20 bp (Quality filtering).
3. The remaing read pairs are output as pre-processed read pairs in `.fastq` files.
4. The pre-processed read pairs are mapped to transcriptome with BWA separately. ‘-a’ option is enabled to keep all found alignments using default threshold of BWA. This is used in the later filtering of potential homologous read pairs (Mapping). 
5. The mapped read pairs are output in `.bed` file with aligned genes and transcriptome alignment information.
6. The transcriptome alignment information of mapped read pairs is utilized to select read pairs whose two ends’ primary alignments are mapped to two different genes (chimeric read pairs identification).
7. The selected read pairs are further checked to see if both ends have over 50% of their read bases matches the reference transcriptome based on the CIGAR string and if both ends have no shared lesser alignments.
8. The read pairs passing the quality checks are then deduplicated based on the external coordinates of their primary alignments (Deduplication).
9. The deduplicated read pairs with one end sensely mapped to a gene and the other end antisensely mapped to a protein-coding gene  were identified as the valid chimeric read pairs from the library (RNA/protein end assignment).
10. The read ids and alignment infomation of the primary alignment of the valid chimeric read pairs are output in `validChimericReadPairs.csv`. 
11. Chi-square test is applied to the valid chimeric read pairs. Benjamini-Hochberg adjustment is applied to correct all the p-values. Gene pairs with an adjusted p-value less than 0.05 (default) and with an odds ratio larger than 1 (default) are kept. Gene pairs with mapped chimeric read pair count in the library larger than 4 (default) times the average number of mapped chimeric read pairs per gene pair in the library are kept. 
12. The kept gene pairs are output as RNA-protein associations in `RNAProteinAssociations.csv`.


## Software Requirements
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) (4.3 or later)
- [fastp](https://github.com/OpenGene/fastp) (0.22.0 or later)
- [bwa](https://github.com/lh3/bwa) (0.7.17-r1188 or later)
- [samtools](http://www.htslib.org/) (1.6 or later)
- [bedtools](https://bedtools.readthedocs.io/en/latest/) (2.30.0 or later)
- Python 3.4 or later, the following python libraries are required:<br />
    - sys
    - collections
    - cigar
    - glob
    - scipy
    - rpy2
    - datetime
## Additional files required
**BWA Index of the transcriptome to be aligned**<br />
You will need to download or build the bwa index of the target trancriptome for PROPERseqTools to use. Here we provide the compressed bwa index built from [RefSeq GRCh38 transcriptome](https://drive.google.com/file/d/1lAV-dVVwVaPi-qVLXibaAvgjtd1e-QMT/view?usp=sharing)

**Transcript, gene and gene type dictionary file**<br />
You will also need a dictionary file that contains the information of transcript ids to their corresponding gene names/gene ids and corresponding gene types in a csv format with the first column being transcript ids, the second column being gene names and the third column being gene types. Here we provide an example dictionary file for [RefSeq GRCh38 genome](https://github.com/Zhong-Lab-UCSD/PROPERseqTools/blob/master/refSeq_tx_gene_type.csv)



## Usage
**Installation**
1. Clone the current github repository to your local machine. For example<br />
`git clone https://github.com/Zhong-Lab-UCSD/PROPERseqTools`
2. Add the following path of the cloned directory to your `.bashrc` file<br />
`export PATH=$PATH:/home/path/to/PROPERseqTools/bin`

**To excute PROPERseqTools, run**
<pre><code>
properseqTools -a /path/to/read1.fastq
               -b /path/to/read2.fastq
               -i /path/to/bwaIndex/transcriptome.fa
               -o /path/to/outputDir
               -g /path/to/refSeq_tx_gene_type.csv
           
</code></pre>



**Required parameters**
<pre><code>
-a     |String, Path to read1 fastq file, fastq.gz also supported
-b     |String, Path to read2 fastq file, fastq.gz also supported
-o     |String, Path to output directory
-i     |String, Path to bwa index of the target transcriptome
-g     |String, Path to transcirpt, gene and gene type dictionary file
</code></pre>
**Other parameters**
<pre><code>
-d     |Float, odds ratio cutoff used to identify RNA-protein associations, default=1
-p     |Float, false discovery rate cutoff used to identify RNA-protein associations, default=0.05
-c     |Float, read count cutoff coefficient used to identify RNA-protein associations, default=4
-j     |String, Job ID to be prepended to the output files and directories, optional, default=PROPERseq"
-t     |Int, Number of working threads, default=2
-r     |Char, (T or F), removal of intermediate files or not, default=T
-h     |Print usage message" 
           
</code></pre>
