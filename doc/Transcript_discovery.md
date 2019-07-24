# Transcript discovery

While evaluating RNA alignment on real datasets, [RNAseqEval.py](doc/RNAseqEval.md) script can try and discover new transcripts. Thanscript discovery is enabled by setting the option --calc_new_annotations. New transcripts are discovered in two ways, by conbining existing annotations and by removing or skipping small introns. New transcripts are discovered only for those alignments that do not perfectly match any available annotations. Finally, new transcripts are reported to the user only if they are supported by a minimum number of reads (currently set at three).

During the regular evaluation process, a set of candidate annotations is constructed for each alignment, consisting of all annotations that overlap the alignment. From the set of candidate annotations, a _best_match_annotation_ is chosen based on the number of nucleotides from the alignment that fall inside and outside of each annotation. New transcripts are calculated only for those alignments that do not perfectly match the _best_match_annotation_.

After all new transcripts are determined, they are collected and compared, and only those transcripts that are supported by a minimum number of reads (currently set to 3) are reported to the user.

## Output

General evaluation report only contains information on how many new transcripts were discovered and how many reads were used to construct them:

    Found 213 potential new annotations with 4860 alignments
                    
    Detailed report on annotations can be found in an '_annotations.report' file.

Detailed information on discovered transcripts is given in a separate file. If output file is specified (options -o and --output), text '\_annotations.report' is appened to the output filename and if output filename is not specified, detailed annotation report is writen to '\_annotations.report' file.

An exampe of a detailed annotation report is given below.

    Name: New annotation 18 
    Based on:NM_166774
    Strand: +
    Number of reads: 7
    Type:FUSED ANNOTATION
    Reads:
    m160615_181138_42182_c101000182550000001823232709161603_s1_p0/135263/ccs
    m160615_181138_42182_c101000182550000001823232709161603_s1_p0/138904/ccs
    m160713_175918_42182_c101000162550000001823232709161621_s1_p0/136269/ccs
    m160615_181138_42182_c101000182550000001823232709161603_s1_p0/39346/ccs
    m160615_181138_42182_c101000182550000001823232709161603_s1_p0/137514/ccs
    m160713_175918_42182_c101000162550000001823232709161621_s1_p0/14795/ccs
    m160615_181138_42182_c101000182550000001823232709161603_s1_p0/43113/ccs
    Items: [610684, 610897] [611727, 611807] [613615, 613796] [619848, 620178]

New annotation name is generated sequentially. Field __Based on__ contains the transcript new annotation was based on (initial _best_match_annotation_ for the alignment). Field __Type__ determines if the annotation was constructed by combining existing annotations or by intron skipping. The report also lists all reads used to construct the annotation and finaly exons themselves.

## Combining existing annotations

The first method for new transcript discovery tries to combine existing annotations to achieve the _correct_ alignment. The process is visible on the figure below.

Blue annotation is chosen as the _best_match_annotation_, because, compared to the purple annotation, more nucleotide bases from the alignment fall within it and less nucleotide bases fall outside it. However, the first part of the alignment does not perfectly match the first exon in the annotation. Therefore, a set of candidate annotations is searched for an exon that does perfectly match the first part of the alignment. This is the first exon of the purple _candidate alignment_. To construct the new transcript/annotation the first exon from the blue _best_match_annotation_ is replaced by the fist exon from the purple _candidate annotation_.

<img src="/img/combining annotations.png" width="700" height="300" align="middle">

## Skipping small introns

The second method for new transcript discovery creates new annotations by starting with an existing annotation and removing introns smaller than a preset value (currently set to 10). The new annotation simply combines exons seperated by small introns. The process is visible in the figure below.

<img src="/img/skipping introns.png" width="700" height="300" align="middle">

## Results

The algorithm was tested on datasets from the 

| Dataset | Graphmap | Minimap2 | GMap | Common |
| --- | --- | --- | --- | --- |
| 1 | 247 | 213 | 174 | 33 |
| 2 | 217 | 192 | 178 | 18 |
| 4	| 13| 13 | 10 | 2 |

