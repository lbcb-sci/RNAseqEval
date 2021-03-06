# Example dataset
This folder contains a small example dataset that can be used to test our evaluation tools. Sample dataset contains only one group, and was simulated using PBSIM on _D. Melanogaster_ chromosome 4, with coverage 5 and with PacBio ROI parameters. Simulation files (generated by PBSIM) are placed in a folder __group1__, inside a root simulation folder __simulation\_X5__. Headers in the generated fastq file were adjusted as described in [RNAseq_benchmark/Data_preparation.md](/RNAseq_benchmark/Data_preparation.md)

Sample dataset contains 2955 reads.

The folder also contains all files needed to test the dataset and the tools.

Contents of the folder:
- example.fastq - small test dataset, 
- example.sam   - mapping of the example dataset to the reference genome using GraphMap (https://github.com/isovic/graphmap) and it option to map RNA reads with the help of annotations.
- dmelanogaster_chr_genome.fa  - genomic sequence for _D. Melanogaster_ chromosome 4, this represents a reference genome for our example dataset
- dmelanogaster_chr.gft - gene annotations for _D. Melanogaster_ chromosome 4
- dmelanogaster_chr_trans_F.fa - transcriptome for for _D. Melanogaster_ chromosome 4, this is generated from geneomic sequence and annotations using [generate_transcriptome.py](/generate_transcriptome.py). Designation '\_F' means that the transcriptome has been filtered, removing sequences shorter than 100 base-pairs, as required by PBSIM.

### Determining match rate
Matchrate for SAM file contining mappings can be determined using errorrates.py script from https://github.com/isovic/samscripts.

    python /home/kkrizanovic/src/samscripts/src/errorrates.py base dmelanogaster_chr4_genome.fa sample.sam

__Results:__

                       	mean    std     median  min     max
    Error rate:     	9.16%	5.77%	8.51%	0.57%	71.89%
    Insertion rate: 	3.90%	2.80%	3.46%	0.00%	13.32%
    Deletion rate:  	0.97%	1.46%	0.81%	0.00%	60.52%
    Mismatch rate:  	4.29%	2.67%	4.12%	0.00%	13.23%
    Match rate:     	91.80%	5.38%	92.41%	75.06%	100.00%
    Read length:    	2310.96	1518.89	2051.50	97.00	18737.00
    Difference ratio:	48:42:10 (mismatch:insertion:deletion)

The script errorrates.py calculates a detailed error statistics. For our benchmark of RNAseq mapping tools, we have used only mean value for math rate.

### Runnning Process_pbsim_data.py
Evaluation of simulated data is performed using Process_pbsim_data.py script. It requires files generated by PBSIM during simulation proces.

    python /home/kkrizanovic/src/RNAseqEval/Process_pbsim_data.py process ./simulation_X5/ sample.sam dmelanogaster_chr4.gtf

__Output:__

    Analysis results:
    Original Samlines: 2955
    Usable whole alignments (with valid CIGAR string): 2852
    Annotations: 334
    Multiexon genes: 316
    Number of exon start hits: 13685
    Number of exon end hits: 13653
    Number of exon start and end hits: 86729
    Number of good whole alignments: 1362
    MAF: Hit both ends: 322
    MAF: Hit all parts: 2017
    MAF: Hit at least one part: 2543
    MAF: Equals at least one part: 2069
    MAF: Number of split reads: 2340
    MAF: Hit both ends, SPLIT read: 235
    MAF: Hit all parts, SPLIT read: 1618
    MAF: Hit at least one part, SPLIT read: 2144
    MAF: Equals at least one part, SPLIT read: 1885

We can see that out of total 2955 reads, 2852 were mapped. Out of those, only 322 correctly match beginning and end of the read origin within 5 base-pairs (__Hit both ends__). However, more than 2000 alignments overlap all exons from read origin (__Hit all parts__) and over 2500 reads overlap at least one exon from read origin (__Hit at least one part__). We cab say that, in general, the alignments are quite correct, i nost cases mappint to the approximate position of read origin, but they are not very precise, very rarely mapping to the beginning and end or read origin within 5 base-pairs.

### Running RNAseqEval.py
Evaluation of real data is performed using RNAseqEval.py script. It requires gene annotations for comparison to mappings.

    python /home/kkrizanovic/src/RNAseqEval/RNAseqEval.py eval-mapping dmelanogaster_chr4_genome.fa sample.sam -a dmelanogaster_chr4.gtf --expression

__Output:__

    Reference format: ANNOTATION
    General information:
    Reference length = 1348131 bp
    Number of chromosomes = 1
    Chromosomes:
            chr4: 5537446bp

    Number of alignments in SAM file (total / unique) = 2955 / 2955
    Alignments with / without CIGAR string = 2852 / 103
    Mapping quality without zeroes (avg / min / max) = 11202.77 / 134 / 78995
    Alignments with mapping quality (>0 / =0) = 27545 / 0
    Number of matches / mismatches / inserts / deletes = 5459186 / 870685 / 259341 / 64612
    Percentage of matches / mismatches / inserts / deletes = 0.82 / 0.13 / 0.04 / 0.01

    Annotation statistics:
    Total gene length = 4189315
    Number of Transcripts / Exons (Multiexon transcripts) = 334 / 3131 (316)
    Maximum number of exons in a gene = 46
    Gene size (Min / Max / Avg) = 69 / 52075 / 12542.86
    Exon size (Min / Max / Avg) = 6 / 9812 / 426.99

    Mapping quality information:
    Bases in reads (aligned / total) (percent) = (6589212 / 6590850) (99.98%)
    Transcripts covered / missed / total = 237 / 97 / 334
    Exons covered / missed / total = 2055 / 1076 / 3131
    Alignments on transcript hit / missed = 2852 / 0
    Alignments on exons hit / missed = 2852 / 0
    Alignments hitting an exon (start / end / both) = 13530 / 13543 / 11456
    Contiguous / non contiguous alignments: 1738 (60.92%) / 1114 (39.05%)

    Transcript/exon expression and coverage information:
    Number of expressed Transcripts = 237
    Expressed transcripts: ...

The above contains several sections. It contains some statistical information on each input it receives, genome reference, gene annotations  and SAM file with alignments. Those sections are straightforward. Last two sections contain statistics calculated by comparing the inputs. The first of them asseses the quality of mapping, while the second one contains gene expression information.

Lookng at __Mapping quality information__ part, we can see that all 2852 alignments that were mapped manage to hit a transcript (__Alignments on transcript__) and an exon (__Alignments on exons__). This is a significantly better result than Process_pbsim_data.py script calculated. However, since RNAseqEval.py script doesn't consider actual read origins, we can not know which alignment is accurate, and which is aligned to an incorrect position while still overlapping a exon.

The metric called __Contiguous alignments__ attempts to determine the correct alignments by comparing them to gene annotations. Those alignments that overlap a contiguous series of exons (not skipping any) belonging to a single transcript, and that manage to align correctly (within 5 base-pairs) to all inner exon boundaries (start and ends) are oconsidered contiguous. In the upper report, 60.92% (or 1738) of the mapped alignments are contiguous. We can say with very high probability that those alignments are correct, but since read origins are not known, we still can not be certain. Some of the other alignments are also aligned to a correct approximate position, but they do not match the annotations well enough.
