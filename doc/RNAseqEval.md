# RNAseqEval.py
Run RNAseqEval.py for general evaluation of mappings in SAM format against reference and optionally annotations. This script is intended to evaluate real dataset mapping. Run RNAseqEval.py without any arguments to print options.

Usage:
     
    RNAseqEval.py eval-mapping <reference FASTA file> <input SAM file> options

## Evaluation method
Eventhough it allows other usage, the main purpose of RNAseqEval.py script is to evaluate the quality of RNA alignments by comparing them to a set of annotations and a reference genome. It's intended use is for real data, for which exact origin of each read is not known. To use the script in this way, it has to be run in eval-mapping mode (see below), with reference genome and mapping in SAM format as required inputs, and with annotations as extra input (-a option).

Example for using RNAseqEval.py to evaluate mappings agains annotations and a reference genome:

    RNAseqEval.py eval-mapping dmelanogaster_genome.fa mappings.sam -a dmelanogaster_annotations.gtf

The script evaluates one read (alignment) at a time, and for each read (alignment) goes through the following steps:
1. Find all candidate annotations (annotations with which the alignment overlaps, looking only at start and end of complete annotation to speed the process up)
2. Compare the alignment to all candidate annotations in more detail and find the one with whom the alignment has the largest overlap. Ths annotations is termed _best_match_annotation._
3. Determine the match between the alignments and the _best_match_annotation_ by calculating four maps:
     - Exon hit map - which exons are overlapped by the alignment
     - Exon complete map - which exons exactly match a part of alignment
     - Exon start map - the start of which exon is covered by the alignments
     - Exon end map - the end of which exon is covered by the alignment
4. Using those four maps, several metrics of similaty between the alignments and the anotation are calculated. The most importan metric is whether the alignment is contiguous or not.

_Contiguous alignment_ represent a read that is correctly aligned to the referece genome or more specifically to the _best_match_annotation_. It covers a contiguous subset of exons from the annotation. Whether an alignment is contiguous is determined using the _hit maps_ calculated in the step 4. Contiguous alignments can skip (not overlap) one or more exons at the start and skip one or more exons at the end. However, if two exons are overlapped by the alignment, all exons between those two must also be overlapped for the alignment to be contiguous. This is determined by applying the following rules:
- Exon _hit map_ must not have _holes_ in the middle. 
- Internal _hit_ (or overlapped) exons must exactly match the alignment. 
- The alignment must match the end of the first exon in the _hit map_ and it must match the start of the last exon in the _hit map_.

_Best_match_annotation_, or the annotation with which the alignment will be compared, is determined based on the number of nucleotide bases from allignment which fall inside and outside the annotation. First, all candidate annotations are determined, i.e. all annotation that overlap with the alignment. For each annotation, a score is calculated by adding the number of bases that fall within the annotation and subtracting the number of bases that fall outside the annotation. The annotation with the highest score is chosen as the _best_match_annotation_. Old scoring algorithm only took into account the number of bases falling inside an annotation and can be used by setting the parameter --old_bma_calc.

The script works in multiple threads and will spawn one thread for each cromosome and strand, and will examine anotations and alignments that strand and chromosome, thus significantly speeding up the analysis. 

__IMPORTANT:__ When making calculations, an error of 5 bases is premitted. Similarly, for an overlap to be valid it has to be at least 5 bases. This can be altered by setting appropriate parameters when running the script.

## Usage modes
RNAseqEval.py script can be used in three differents modes, determined by the first argument. Each mode requires different parameters and allowes different options.

### eval-mapping
Used in eval-mapping mode, RNAseqEval.py script is used to evaluate RNAseq mappings against known FASTA reference and annotations. Annotations can be omitted, but in that case the script will provide only basic output.

Usage:

    RNAseqEval.py eval-mapping <reference FASTA file> <input SAM file> options
    
Allowed options:

     -a <file> : a reference annotation (GFF/GTF/BED) file
     -o (--output) <file> : output file to which the report will be written
     -ex (--expression) : calculate and output gene expression
     --no_check_strand : when matching alignments to annotations, do not take strand 
                         into account, only chromosome and position
     --no_per_base_stats : do not calculate per-base stats, such are percentage of
                           matches, mismatches, inserts and deletes
     -sqn (--save_query_names) : save query names for alignments that managed to hit an exon
                                 and for contiguous alignments. Query names are saved to files
                                 with filenames determined from output file. 
                                 For this option output and annotation files must be specified
     -ai (--alowed_inaccuracy) : A maximumn distance in bases between an annotation and an alignment
                                 where the alignment is still considered correct (default 5)
     -mo (--min_overlap) : A minimum overlap between an annotation and an alignment that is considered valid
                           (default 5)
     --old_bma_calc : Calculate best matching annotation only based on maximizing the number of bases an alignment
                      on an annotation. The number of bases outside an annotation is not take into account in this case
     --leave_chrom_names : do not preprocess chromosome names in SAM, reference and annotation files
                           Chromosome names are by default preprocessed to make them match the format "chr[designation]"
     --calc_new_annotations: calculate potential new annotations, if a sufficient number of alignments (default 3)
                             better fits a combination of exons then any existing annotation, that combination
                             of exons is suggested as a new annotation

### eval-annotations
Used in eval-annotations mode, RNAseqEval.py script will print out basic information on an annotations file.

Usage:

    RNAseqEval.py eval-annotations <annotations file> options

Allowed options:

    -o (--output) <file> : output file to which the report will be written

### eval-maplength
Used in eval-maplength mode, RNAseqEval script will return mapped percentage for each read

Usage:

    RNAseqEval.py eval-maplength <input SAM file> options

Options:

    -o (--output) <file> : output file to which the report will be written

Oposed to first two modes which calculate certain statistical information from input files, in eval-maplength mode the script will print out information on each read in CSV format (on the screen or in a file). The folowinf information is printed out:
- readname name (header "QNAME")
- reference name (header "RNAME")
- read length (header "read length")
- the number of bases aligned for that read (header "bases aligned")

## Output for eval-mapping and eval-annotations modes
Depending on the usage mode, RNAseqEval.py script will display various information about input files and the results of the analysis. If option --calc_new_annotations is set, the script will also try to determine possible new transcripts. More on transcript dicovery can be found at [doc/Transcript_discovery.md](doc/Transcript_discovery.md.).

General information on FASTA reference and mapping SAM file:

    - Reference length - In eval-mapping mode this will be the total lenght of all chromosomes in a FASTA rederence, while in eval-annotations mode this will be the total length of all genes.
    - Number of chromosomes
    - List of chromosomes
    - Number of alignments in SAM file (total / unique) - two alignments are not unique if they represent the same read
    - Alignments with / without CIGAR string
    - Mapping quality without zeroes (avg / min / max)
    - Alignments with mapping quality (>0 / =0)
    - Number of matches / mismatches / inserts / deletes - calculated per base in total for all reads
    - Percentage of matches / mismatches / inserts / deletes

Annotation statistics:

    - Total gene length
    - Total number of transcripts
    - Total number of exons
    - Number of multiexon transcripts
    - Maximum number of exons in a gene
    - Gene size (Min / Max / Avg)
    - Exon size (Min / Max / Avg)

Mapping quality information obtained by comparing alignements in a SAM file to given annotations. Only in eval-mapping mode if annotations are provided.

     - Total number and percentage of bases aligned for all reads
     - The number of transcripts (annotations) "hit" by all reads - an annotation is "hit" by a read if the read overlaps it on at least 5 bases
     - Total number of exons "hit" by all reads - an exon is "hit" by a read if the read overlaps it on at least 5 bases
     - Number of alignments with "hit" on transcripts
     - Number of alignments with "hit" on exons
     - Number of alignments matching a beginning and an end of an exon
     - Number of alignments where more than 50% of bases fall within an annotation
     - Number of alignments with low match count
     - Number of contiguous and non contiguous alignments - as described earlier in the text

If so specified by the option -ex (--expression), the script also calculates gene expression and gene/exon coverage information. This option is available only in eval-mapping mode if annotations are provided. The script will output the number of expressed transcripts. A transcript is considered expressed if at least one read is mapped to its position. For each transcript, the script also prints out the following:

    - transcript name
    - number of exons
    - number of reads that align to it
    - total number of bases aligned to it
    - for each exon in the transcript
         - number of reads aligned to it
         - total number of bases aligned to it
