#! /usr/bin/python

import sys, os
import re
import setup_RNAseqEval, paramsparser

# For copying SAM lines
import copy

from datetime import datetime

# To enable importing from samscripts submodule
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
import utility_sam
import Annotation_formats

from fastqparser import read_fastq
from report import EvalReport, ReportType

# Multiprocessing stuff
import multiprocessing

# import time

DISTANCE_THRESHOLD = 10000
MIN_OVERLAP_BASES = 5
NUM_PROCESS = 12
NEW_ANNOTATION_MIN = 3      # A Minimum number of reads for new annotation to be proposed
MIN_INTRON_SIZE = 10         # Used for detecting new annotations


# TODO: Osim broja readova koji pokrivaju pojedini gen, izracunati i coverage

annotations2check = ['NM_001170371', 'NM_079876', 'NM_166709', 'NM_001015387', 'NM_001316660', 'NM_001258510'] + ['NM_137066']

# New annotation calculation level
NA_NO_NEW_ANNOTATINS = 0
NA_SKIP_INTRONS = 1
NA_COMBINED_ANNOTATIONS_1 = 2
NA_COMBINED_ANNOTATIONS_2 = 4
NA_DROP_EXONS = 8

# Parameter definitions for paramparser
paramdefs = {'-a' : 1,
             '--version' : 0,
             '-v' : 0,
             '-o' : 1,
             '--output' : 1,
             '-ex' : 0,
             '--expression' : 0,
             '-as' : 0,
             '--alternate_splicing' : 0,
             '--print_erroneous_reads' : 0,
             '--no_check_strand' : 0,
             '--no_per_base_stats' : 0,
             '-sqn' : 0,
             '--save_query_names' : 0,
             '--alowed_inaccurycy' : 1,
             '-ai' : 1,
             '--min_overlap' : 1,
             '-mo' : 1,
             '--graphmap' : 0,
             '--old_bma_calc' : 0,
             '--leave_chrom_names': 0,
             '--calc_new_annotations' : 1,
             '--no_multiprocessing' : 0,
             '--check_reads' : 0}


def cleanup():
    pass


# A function that looks at already calculated exon maps and checks if an alignment is good (correct) and spliced
def isGoodSplitAlignment(exonhitmap, exoncompletemap, exonstartmap, exonendmap):

    isGood = True
    isSpliced = False
    if not (len(exonhitmap) == len(exoncompletemap) and len(exonhitmap) == len(exonstartmap) and len(exonhitmap) == len(exonendmap)):
        raise Exception('ERROR: Exon maps have unequal lengths (%d|%d|%d|%d)!' % (len(exonhitmap), len(exoncompletemap), len(exonstartmap), len(exonendmap)))

    for i in exonhitmap.keys():
        if exonhitmap[i] == 0:
            if exoncompletemap[i] <> 0:
                raise Exception('ERROR: HIT map 0 and COMPLETE map nonzero!')
            if exonstartmap[i] <> 0:
                raise Exception('ERROR: HIT map 0 and START map nonzero!')
            if exonendmap[i] <> 0:
                raise Exception('ERROR: HIT map 0 and END map nonzero!')

    # A list of indices of exons for which a hit map is nonzero
    hitlist = [i for i in exonhitmap.keys() if exonhitmap[i] > 0]

    if len(hitlist) == 0:
        return False, False

    starthit = hitlist[0]
    endhit = hitlist[-1]
    middlelist = hitlist[1:-1]

    # For an alignment to be spliced, the hit list has to skip some exons in the middle
    for x in hitlist[:-1]:
        if exonhitmap[x+1] - exonhitmap[x] > 1:
            isSpliced = True
            break           # No need to look further

    # For an alignment to be strictly good, it has to be uninterrupted
    # It has to end the first hit exon, start the last exon, end complete all in the middle
    middleOK = True
    for x in middlelist:
        if exoncompletemap[x] == 0:
            middleOK = False
            break           # No need to look further

    if (not middleOK) or exonstartmap[endhit] == 0 or exonendmap[starthit] == 0:
        isGood = False

    return isGood, isSpliced


# A function that calculates exon maps for a given alignment and annotation 
def calcExonMaps(samline_list, annotation, allowed_inacc, min_overlap):
    # Initializing the maps
    t_exonhitmap = {(i+1):0 for i in xrange(len(annotation.items))}
    t_exoncompletemap = {(i+1):0 for i in xrange(len(annotation.items))}
    t_exonstartmap = {(i+1):0 for i in xrange(len(annotation.items))}
    t_exonendmap = {(i+1):0 for i in xrange(len(annotation.items))}

    # Calculating hit maps
    for samline in samline_list:
        item_idx = 0
        lstartpos = samline.pos
        reflength = samline.CalcReferenceLengthFromCigar()
        lendpos = lstartpos + reflength
        exonhit = False
        for item in annotation.items:
            item_idx += 1
            if item.overlapsItem(lstartpos, lendpos, allowed_inacc, min_overlap):
                exonhit = True
                t_exonhitmap[item_idx] += 1
                if item.equalsItem(lstartpos, lendpos, allowed_inacc, min_overlap):
                    t_exoncompletemap[item_idx] = 1
                    t_exonstartmap[item_idx] = 1
                    t_exonendmap[item_idx] = 1
                elif item.startsItem(lstartpos, lendpos, allowed_inacc, min_overlap):
                    t_exonstartmap[item_idx] = 1
                elif item.endsItem(lstartpos, lendpos, allowed_inacc, min_overlap):
                    t_exonendmap[item_idx] = 1

    return (t_exonhitmap, t_exoncompletemap, t_exonstartmap, t_exonendmap)



# A helper function that extracts a chromosome name from a fasta header (or other similar strings)
# Annotations and reference can use different chromosome designations, so this is used to 
# correctly compare them
# Chromosome names should be either chromosome [designation] or chr[designation]
# In the case header represents a mitochondrion, 'chrM' is returned!
# If there is no need for chromosome name preprocessing, then the parameter processChromNames should be set to False
def getChromName(header, processChromNames = True):

    # If chromosome names do not need to be processed, simply return unchanged argument (up to the first space)
    if not processChromNames:
        pos = header.find(' ')
        if pos == -1:
            return header
        else:
            return header[:pos]

    chromname = ''
    # regular expressions for searching for long and short chromosome names
    longre = r'(chromosome )(\w*)'
    shortre = r'(chr)(\w*)'

    if header.find('mitochondrion') > -1 or header.find('chrM') > -1:
        chromname = 'chrM'
    else:
        match1 = re.search(longre, header)
        match2 = re.search(shortre, header)
        if match1:
            designation = match1.group(2)
            chromname = 'chr%s' % designation
        elif match2:
            designation = match2.group(2)
            chromname = 'chr%s' % designation
        else:
            chromname = 'genome'        # In case I can't detect chromosome designation or mitochondrion
                                        # I decide its a genome

    return chromname



def load_and_process_reference(ref_file, paramdict, report):
    # Reading FASTA reference
    [headers, seqs, quals] = read_fastq(ref_file)
    processChromNames = True
    if '--leave_chrom_names' in paramdict:
        processChromNames = False

    # Analyzing FASTA file

    # Since annotation file and SAM file (mapper output) do not have to have unified sequence names
    # for individual chromosomes, I'm creating a translation dictionary that will help me access the right
    # sequence more quickly. This will be based on the function which will attempt to analyze a string
    # and determine whether it refers to a genome, a particular chromosome or a mitochndia
    # The chromname2seq dictionary will have a inferred name as key, and index of the corresponding
    # sequence and headeer as value
    chromname2seq = {}

    if len(headers) == 1:
        report.reflength = len(seqs[0])
        chromname = getChromName(headers[0], processChromNames)
        report.chromlengths = {chromname : report.reflength}
        chromname2seq[chromname] = 0
    else:
        for i in xrange(len(headers)):
            header = headers[i]
            chromname = getChromName(header, processChromNames)
            if chromname in report.chromlengths:
                raise Exception('\nERROR: Duplicate chromosome name: %s' % chromname)
                # sys.stderr.write('\nERROR: Duplicate chromosome name: %s' % chromname)
                exit()
            else:
                report.chromlengths[chromname] = len(seqs[i])
                report.reflength += len(seqs[i])
                chromname2seq[chromname] = i

    return [chromname2seq, headers, seqs, quals]



# checks if a given samline could form a split alignment with an existing samline list
# if so, adds the line to the list and returns True, otherwise returns False
def join_split_alignment(samline_list, samline):
    split_possible = True
    for sline in samline_list:
        if not utility_sam.possible_split_alignment(samline, sline, threshold = DISTANCE_THRESHOLD):
            split_possible = False
            break

    if split_possible:
        samline_list.append(samline)

    return split_possible



def load_and_process_SAM(sam_file, paramdict, report, BBMapFormat = False):
    # Loading SAM file into hash
    # Keeping only SAM lines with regular CIGAR string, and sorting them according to position
    qnames_with_multiple_alignments = {}
    [sam_hash, sam_hash_num_lines, sam_hash_num_unique_lines] = utility_sam.HashSAMWithFilter(sam_file, qnames_with_multiple_alignments)

    # If variable BBMapFormat is set to true, all samfiles referring to the same query will be collected together
    # Stil have to decide what to do with query names, currently removing '_part'
    if BBMapFormat:
        new_sam_hash = {}
        for (qname, sam_lines) in sam_hash.iteritems():
            pos = qname.find('_part')
            if pos > -1:
                origQname = qname[:pos]
            else:
                origQname = qname
            if origQname not in new_sam_hash:
                new_sam_hash[origQname] = sam_lines
            else:
                # import pdb
                # pdb.set_trace()
                new_sam_lines = sam_lines + new_sam_hash[origQname]
                new_sam_hash[origQname] = new_sam_lines

        sam_hash = new_sam_hash

    # NOTE: This is a quick and dirty solution
    # Setting this to true so that large deletions are turned into Ns
    # BBMap marks intron RNA alignment gaps with deletions!
    BBMapFormat = True

    # If this option is set in parameters, unmapped queries will be listed in the report
    save_qnames = False
    if '-sqn' in paramdict or '--save_query_names' in paramdict or '--split-qnames' in paramdict:
        save_qnames = True

    # Reorganizing SAM lines, removing unmapped queries, leaving only the first alignment and
    # other alignments that possibly costitute a split alignment together with the first one
    samlines = []
    cnt = 0
    pattern = '(\d+)(.)'
    # for samline_list in sam_hash.itervalues():
    for (samline_key, samline_list) in sam_hash.iteritems():
        cnt += 1
        if samline_list[0].cigar <> '*' and samline_list[0].cigar <> '':            # if the first alignment doesn't have a regular cigar string, skip

            if BBMapFormat:
                # All deletes that are 10 or more bases are replaced with Ns of the same length
                operations = re.findall(pattern, samline_list[0].cigar)
                newcigar = ''
                for op in operations:
                    op1 = op[1]
                    op0 = op[0]
                    if op[1] == 'D' and int(op[0]) >= 10:
                        op1 = 'N'
                    newcigar += op0 + op1
                samline_list[0].cigar = newcigar


            operations = re.findall(pattern, samline_list[0].cigar)
            split = False

            for op in operations[1:-1]:             # Ns cannot appear as the first or the last operation
                if op[1] == 'N':
                    split = True
                    break
            # If the first alignment is split (had Ns in the middle), keep only the first alignment and drop the others
            if split:
                report.num_split_alignments += 1
                # Transform split alignments containing Ns into multiple alignments with clipping
                temp_samline_list = []
                posread = 0
                posref = 0      # NOTE: I don't seem to be using this, probably should remove it
                newcigar = ''
                readlength = samline_list[0].CalcReadLengthFromCigar()
                new_samline = copy.deepcopy(samline_list[0])
                mapping_pos = new_samline.pos
                clipped_bases = new_samline.pos - new_samline.clipped_pos
                hclip_seq = 0        # Used with hard clipping, how big part of sequence should be removed
                clip_type = 'S'     # Soft_clipping by default
                for op in operations:
                    if op[1] == 'N' and int(op[0]) > 1:        # Create a new alignment with clipping
                        newcigar += '%dS' % (readlength - posread)      # Always use soft clipping at the end
                        new_samline.cigar = newcigar
                        # After some deliberation, I concluded that this samline doesn't have to have its position changed
                        # The next samline does, and by the size of N operation in cigar string + any operations before
                        temp_samline_list.append(new_samline)
                        new_samline = copy.deepcopy(samline_list[0])
                        mapping_pos += int(op[0])
                        new_samline.pos = mapping_pos
                        new_samline.clipped_pos = new_samline.pos - clipped_bases
                        posref += int(op[0])
                        if clip_type == 'H':
                            new_samline.seq = new_samline.seq[hclip_seq:]
                        newcigar = '%d%c' % (posread, clip_type)
                    else:                   # Expand a current alignment
                        newcigar += op[0] + op[1]
                        if op[1] in ('D', 'N'):
                            posref += int(op[0])
                            mapping_pos += int(op[0])
                        elif op[1] == 'I':
                            posread += int(op[0])
                            # Everything besides deletes and Ns will be clipped in the next partial alignment
                            # Therefore have to adjust both pos and clipped pos
                            clipped_bases += int(op[0])
                            hclip_seq += int(op[0])
                        elif op[1] in ('S', 'H'):
                            clip_type = op[1]
                            # Clipped bases can not appear in the middle of the original cigar string
                            # And they have already been added to the position,
                            # so I shouldn't adjust my mapping_pos and clipped_bases again
                            # TODO: I should probably diferentiate between hars and soft clipping
                            posread += int(op[0])
                            posref += int(op[0])
                        else:
                            posref += int(op[0])
                            posread += int(op[0])
                            clipped_bases += int(op[0])
                            mapping_pos += int(op[0])
                            hclip_seq += int(op[0])

                new_samline.cigar = newcigar
                temp_samline_list.append(new_samline)

                samlines.append(temp_samline_list)
            else:
                temp_samline_list = [samline_list[0]]        # add the first alignment to the temp list
                multi_alignment = False
                for samline in samline_list[1:]:            # look through other alignments and see if they could form a split alignment with the current temp_samline_list
                    if BBMapFormat:
                        # All deletes that are 10 or more bases are replaced with Ns of the same length
                        operations = re.findall(pattern, samline.cigar)
                        newcigar = ''
                        for op in operations:
                            op0 = op[0]
                            op1 = op[1]
                            if op[1] == 'D' and int(op[0]) >= 10:
                                op1 = 'N'
                            newcigar += op0 + op1
                        samline.cigar = newcigar
                    if not join_split_alignment(temp_samline_list, samline):
                        multi_alignment = True

                if multi_alignment:
                    report.num_multi_alignments += 1
                if len(temp_samline_list) > 1:
                    report.num_possibly_split_alignements += 1
                samlines.append(temp_samline_list)
        else:
            # Samline has invalid cigar and is considered unmapped
            if save_qnames:
                report.unmapped_names.append(samline_list[0].qname)
            pass

    # Sorting SAM lines according to the position of the first alignment
    samlines.sort(key = lambda samline: samline[0].pos)

    # Calculate real split alignments
    num_real_split = 0
    for samline_list in samlines:
        if len(samline_list) > 1:
            num_real_split += 1

    report.num_alignments = sam_hash_num_lines
    report.num_unique_alignments = sam_hash_num_unique_lines
    report.num_real_alignments = len(samlines)
    report.num_real_split_alignments = num_real_split
    report.num_non_alignments = report.num_alignments - len(samlines)       # Not sure if this is correct any more

    return samlines



def load_and_process_annotations(annotations_file, paramdict, report):
    processChromNames = True
    if '--leave_chrom_names' in paramdict:
        processChromNames = False

    # Reading annotation file
    annotations = Annotation_formats.Load_Annotation_From_File(annotations_file)

    # Sorting annotations according to position
    # NOTE: Might not be necessary because they are generally already sorted in a file
    annotations.sort(reverse=False, key=lambda annotation: annotation.start)

    # Analyzing annotations

    # Looking at expressed genes, ones that overlap with at least one read in SAM file
    # Storing them in a dictionary together with a number of hits for each exon in the gene
    expressed_genes = {}

    # Calculating gene coverage, how many bases of each gene and exon are covered by ready
    # Bases covered multiple times are taken into account multiple times
    # The structure of this dictionary is similar to expressed_genes above
    # Each gene has one global counter (index 0), and one counter for each exon
    gene_coverage = {}

    report.totalGeneLength = 0
    report.num_genes = len(annotations)
    report.max_exons_per_gene = 1       # Can not be less than 1
    report.num_exons = 0
    sumGeneLength = 0.0
    sumExonLength = 0.0
    for annotation in annotations:
        # Initializing a list of counters for a gene
        # Each gene has one global counted (index 0), and one counter for each exon
        expressed_genes[annotation.genename] = [0 for i in xrange(len(annotation.items) + 1)]
        gene_coverage[annotation.genename] = [0 for i in xrange(len(annotation.items) + 1)]

        if len(annotation.items) > 1:
            report.num_multiexon_genes += 1

        # Determining a maximum number of exons per gene
        if len(annotation.items) > report.max_exons_per_gene:
            report.max_exons_per_gene = len(annotation.items)

        report.totalGeneLength += annotation.getLength()
        chromname = getChromName(annotation.seqname, processChromNames)
        if chromname in report.chromlengths:
            report.chromlengths[chromname] += annotation.getLength()
        else:
            report.chromlengths[chromname] = annotation.getLength()
        report.num_exons += len(annotation.items)
        glength = annotation.getLength()
        if glength < report.min_gene_length or report.min_gene_length == 0:
            report.min_gene_length = glength
        if glength > report.max_gene_length or report.max_gene_length == 0:
            report.max_gene_length = glength
        sumGeneLength += glength
        for item in annotation.items:
            elength = item.getLength()
            if elength < report.min_exon_length or report.min_exon_length == 0:
                report.min_exon_length = elength
            if elength > report.max_exon_length or report.max_exon_length == 0:
                report.max_exon_length = elength
            sumExonLength += elength
    report.avg_gene_length = sumGeneLength / report.num_genes
    report.avg_exon_length = sumExonLength / report.num_exons

    return annotations, expressed_genes, gene_coverage



# If the best_match_annotation does not produce a "correct" alignment check if it contains overly small
# introns. If by joining exons seperated by the samll exon we can achieve correct alignment, propose
# the new annotation!
def proposeNewFusedAnnotation(samline_list, best_match_annotation, bma_score, allowed_inacc, min_overlap, old_bma_calc = False):
    # Create a new annotation from the best match annotation, that removes short introns.
    new_annotation = Annotation_formats.GeneDescription()
    new_annotation.seqname = samline_list[0].qname
    # new_annotation.genename = "New annotation %d " % (len(new_annotations)+1)
    new_annotation.source = best_match_annotation.genename
    new_annotation.strand = best_match_annotation.strand
    new_annotation.transcriptname = "SMALL INTRON DISCARDING"       # New annotation type
    new_annotation.items = []
    fused = False
    isNAGood = False

    # Fusing exons separated by very short introns
    items = sorted(best_match_annotation.items, key = lambda it: it.start)
    current_aitem = copy.copy(items[0])
    for aitem in items[1:]:
        if aitem.start - current_aitem.end < MIN_INTRON_SIZE:
            current_aitem.end = aitem.end
            fused = True
        else:
            new_annotation.items.append(current_aitem)
            current_aitem = copy.copy(aitem)

    new_annotation.items.append(current_aitem)

    if fused:
        (t_exonhitmap, t_exoncompletemap, t_exonstartmap, t_exonendmap) = calcExonMaps(samline_list, new_annotation, allowed_inacc, min_overlap)

        isNAGood, isSpliced2 = isGoodSplitAlignment(t_exonhitmap, t_exoncompletemap, t_exonstartmap, t_exonendmap)
        new_score = calc_annotation_match_score(samline_list, new_annotation, old_bma_calc)
        if new_score > bma_score:
            return (new_annotation, isNAGood)

    return (None, isNAGood)


# If the best_match_annotation does not produce a "correct" alignment check if better annotation
# can be obtained by combining the best match annotation with any of the other candidate annotations
def proposeNewCombinedAnnotation(samline_list, best_match_annotation, candidate_annotations, bma_score, allowed_inacc, min_overlap, old_bma_calc = False):
    # Initializa new annotation
    new_annotation = Annotation_formats.GeneDescription()
    new_annotation.seqname = samline_list[0].qname
    # new_annotation.genename = "New annotation %d " % (len(new_annotations)+1)
    new_annotation.source = best_match_annotation.genename
    new_annotation.strand = best_match_annotation.strand
    new_annotation.transcriptname = "FUSED ANNOTATION"          # New annotation type
    new_annotation.items = []
    proposeNew = False      # Determines whether we want to propose a new annotation
    isNAGood = False

    # Sort samlines according to position
    sll2 = sorted(samline_list, key = lambda sl: sl.pos)

    # Calculating alignment start and end
    startpos = sll2[0].pos
    endpos = sll2[-1].pos + sll2[-1].CalcReferenceLengthFromCigar()

    numAlignments = len(sll2)       # Number of exons in an alignment
    start = True
    end = False

    # Adding best match annotation exons that are before the alignment, to the new alignment
    for aitem in best_match_annotation.items:
        if aitem.end < startpos:
            new_annotation.items.append(copy.copy(aitem))
   
    counter = 1     # Counting the partial alignments to allow different procession
                    # of the first and the last partial alignment
    # Check each inside partial alignment and compare it to exons in best annotation
    for samline in sll2:    # Find an alignment that overlaps the exon
        lstartpos = samline.pos
        reflength = samline.CalcReferenceLengthFromCigar()
        lendpos = lstartpos + reflength
        good = False
        replacementFound = False
        ovlitem = None
        for aitem in best_match_annotation.items:
            if aitem.overlapsItem(lstartpos, lendpos, allowed_inacc, min_overlap):
                # Check if partial alignment perfectly matches the annotation exon
                # The first partial alignment can begin within the exon
                # The last partial aligment can end within the exon
                # Middle partial alignments must be equal to the exon
                if aitem.equalsItem(lstartpos, lendpos, allowed_inacc, min_overlap) or \
                   (counter == 1 and aitem.endsItem(lstartpos, lendpos, allowed_inacc, min_overlap) \
                        and lstartpos + allowed_inacc >= aitem.start) or \
                   (counter == numAlignments and aitem.startsItem(lstartpos, lendpos, allowed_inacc, min_overlap) \
                        and lendpos - allowed_inacc <= aitem.end):
                    ovlitem = aitem
                    good = True
                break

        # If partial aligment doesn't perfectly match the exon
        # Try to find a better match among other candidate annotations
        newItem = None
        if not good:
            for cannotation in candidate_annotations:
                for aitem in cannotation.items:
                    # The same condition for perfect match as above
                    if aitem.equalsItem(lstartpos, lendpos, allowed_inacc, min_overlap) or \
                       (counter == 1 and aitem.endsItem(lstartpos, lendpos, allowed_inacc, min_overlap) \
                            and lstartpos + allowed_inacc >= aitem.start) or \
                       (counter == numAlignments and aitem.startsItem(lstartpos, lendpos, allowed_inacc, min_overlap) \
                            and lendpos - allowed_inacc <= aitem.end):
                        replacementFound = True
                        newItem = aitem

        # If the replacement exon is found, place it in the new annotation, otherwise place old exon (if it exists)
        if replacementFound:
            new_annotation.items.append(copy.copy(newItem))
            proposeNew = True
        elif ovlitem is not None:
            new_annotation.items.append(copy.copy(ovlitem))

        counter += 1

    # Adding best match annotation exons that are after the alignment,to the new alignment
    for aitem in best_match_annotation.items:
        if aitem.start > endpos:
            new_annotation.items.append(copy.copy(aitem))

    # If better exon matches have been found, propose a new annotation
    if proposeNew:
        (t_exonhitmap, t_exoncompletemap, t_exonstartmap, t_exonendmap) = calcExonMaps(samline_list, new_annotation, allowed_inacc, min_overlap)

        isNAGood, isSpliced2 = isGoodSplitAlignment(t_exonhitmap, t_exoncompletemap, t_exonstartmap, t_exonendmap)
        new_score = calc_annotation_match_score(samline_list, new_annotation, old_bma_calc)
        if new_score > bma_score:
            return (new_annotation, isNAGood)

    return (None, isNAGood)


# If the best_match_annotation does not produce a "correct" alignment check if better annotation
# can be obtained by combining the best match annotation with any of the other candidate annotations
# Second version of the function, does not look only of perfect match new exons, but also for those that
# ovelap more with the alignment
def proposeNewCombinedAnnotation2(samline_list, best_match_annotation, candidate_annotations, bma_score, allowed_inacc, min_overlap, old_bma_calc = False):
    # Initializa new annotation
    new_annotation = Annotation_formats.GeneDescription()
    new_annotation.seqname = samline_list[0].qname
    # new_annotation.genename = "New annotation %d " % (len(new_annotations)+1)
    new_annotation.source = best_match_annotation.genename
    new_annotation.strand = best_match_annotation.strand
    new_annotation.transcriptname = "FUSED ANNOTATION"          # New annotation type
    new_annotation.items = []
    proposeNew = False      # Determines whether we want to propose a new annotation
    isNAGood = False

    # Sort samlines according to position
    sll2 = sorted(samline_list, key = lambda sl: sl.pos)

    # Calculating alignment start and end
    startpos = sll2[0].pos
    endpos = sll2[-1].pos + sll2[-1].CalcReferenceLengthFromCigar()

    numAlignments = len(sll2)       # Number of exons in an alignment
    start = True
    end = False    

    # Adding best match annotation exons that are before the alignment, to the new alignment
    for aitem in best_match_annotation.items:
        if aitem.end < startpos:
            new_annotation.items.append(copy.copy(aitem))

    counter = 1     # Counting the partial alignments to allow different procession
                    # of the first and the last partial alignment
    # Check each inside partial alignment and compare it to exons in best annotation
    for samline in sll2:    # Find an alignment that overlaps the exon
        lstartpos = samline.pos
        reflength = samline.CalcReferenceLengthFromCigar()
        lendpos = lstartpos + reflength
        overlap = -1
        good = False
        replacementFound = False
        ovlitem = None
        bma_item = None
        for aitem in best_match_annotation.items:
            if aitem.overlapsItem(lstartpos, lendpos, allowed_inacc, min_overlap):
                bma_item = aitem
                # Check if partial alignment perfectly matches the annotation exon
                # The first partial alignment can begin within the exon
                # The last partial aligment can end within the exon
                # Middle partial alignments must be equal to the exon
                if aitem.equalsItem(lstartpos, lendpos, allowed_inacc, min_overlap) or \
                   (counter == 1 and aitem.endsItem(lstartpos, lendpos, allowed_inacc, min_overlap) \
                        and lstartpos + allowed_inacc >= aitem.start) or \
                   (counter == numAlignments and aitem.startsItem(lstartpos, lendpos, allowed_inacc, min_overlap) \
                        and lendpos - allowed_inacc <= aitem.end):
                    good = True

                # Calculate the overlap between exon and alignment
                overlap = aitem.basesInside(lstartpos, lendpos)
                break

        # If partial aligment doesn't perfectly match the exon
        # Try to find a better match among other candidate annotations
        newItem = None
        if not good:
            for cannotation in candidate_annotations:
                for aitem in cannotation.items:
                    # Looking for an exon with larger overlap with alignment, which can fit into a current
                    # best match annotation (there is no overlap with other exons in BMA)
                    if aitem.overlapsItem(lstartpos, lendpos, allowed_inacc, min_overlap):
                        newoverlap = aitem.basesInside(lstartpos, lendpos)
                        if newoverlap > overlap:
                            canFit = True
                            for aitem2 in best_match_annotation.items:
                                if (bma_item is not None) and not aitem2.equalsItem(bma_item.start, bma_item.end) and \
                                       aitem2.overlapsItem(aitem.start, aitem.end):
                                    canFit = False
                                    break
                                if bma_item is None and aitem2.overlapsItem(aitem.start, aitem.end):
                                    canFit = False
                                    break

                            if canFit:
                                replacementFound = True
                                overlap = newoverlap
                                newItem = aitem

        # If the replacement exon is found, place it in the new annotation, otherwise place old exon (if it exists)
        if replacementFound:
            new_annotation.items.append(copy.copy(newItem))
            proposeNew = True
        elif bma_item is not None:
            new_annotation.items.append(copy.copy(bma_item))

        counter += 1

    # Adding best match annotation exons that are after the alignment,to the new alignment
    for aitem in best_match_annotation.items:
        if aitem.start > endpos:
            new_annotation.items.append(copy.copy(aitem))

    # If better exon matches have been found, propose a new annotation
    if proposeNew:
        (t_exonhitmap, t_exoncompletemap, t_exonstartmap, t_exonendmap) = calcExonMaps(samline_list, new_annotation, allowed_inacc, min_overlap)

        isNAGood, isSpliced2 = isGoodSplitAlignment(t_exonhitmap, t_exoncompletemap, t_exonstartmap, t_exonendmap)
        new_score = calc_annotation_match_score(samline_list, new_annotation, old_bma_calc)
        if new_score > bma_score:
            return (new_annotation, isNAGood)

    return (None, isNAGood)



# If the best_match_annotation does not produce a "correct" alignment check if better annotation
# can be obtained by dropping an exon (or mseveral exons) from the best match annotation
def proposeNewDropExonAnnotation(samline_list, best_match_annotation, bma_score, allowed_inacc, min_overlap, old_bma_calc = False):
    # Initializa new annotation
    new_annotation = Annotation_formats.GeneDescription()
    new_annotation.seqname = samline_list[0].qname
    # new_annotation.genename = "New annotation %d " % (len(new_annotations)+1)
    new_annotation.source = best_match_annotation.genename
    new_annotation.strand = best_match_annotation.strand
    new_annotation.transcriptname = "DROP EXON ANNOTATION"          # New annotation type
    new_annotation.items = []
    proposeNew = False      # Determines whether we want to propose a new annotation
    isNAGood = False

    # Sort samlines according to position
    sll2 = sorted(samline_list, key = lambda sl: sl.pos)

    # Calculating alignment start and end
    startpos = sll2[0].pos
    endpos = sll2[-1].pos + sll2[-1].CalcReferenceLengthFromCigar()

    numAlignments = len(sll2)       # Number of exons in an alignment
    start = True
    end = False

    # Adding best match annotation exons that are before the alignment, to the new alignment
    for aitem in best_match_annotation.items:
        if aitem.end < startpos:
            new_annotation.items.append(copy.copy(aitem))

    # Checking if any annotation exons that fall within the global alignment limits
    # Do not overlap with any of the partial alignments
    # Such exons are droppped, and a new alignment is proposed
    for bma_item in best_match_annotation.items:
        dropItem = True
        if bma_item.end >= startpos and bma_item.start <= endpos:       # If an annotation item overlaps with the aligment
            for samline in sll2:
                lstartpos = samline.pos
                reflength = samline.CalcReferenceLengthFromCigar()
                lendpos = lstartpos + reflength
                if bma_item.overlaps(lstartpos,lendpos, allowed_inacc, min_overlap):
                    dropItem = False

            if dropItem:
                proposeNew = True
            else:
                new_annotation.items.append(copy.copy(bma_item))

    # Adding best match annotation exons that are after the alignment,to the new alignment
    for aitem in best_match_annotation.items:
        if aitem.start > endpos:
            new_annotation.items.append(copy.copy(aitem))

    # If better exon matches have been found, propose a new annotation
    if proposeNew:
        (t_exonhitmap, t_exoncompletemap, t_exonstartmap, t_exonendmap) = calcExonMaps(samline_list, new_annotation, allowed_inacc, min_overlap)

        isNAGood, isSpliced2 = isGoodSplitAlignment(t_exonhitmap, t_exoncompletemap, t_exonstartmap, t_exonendmap)
        new_score = calc_annotation_match_score(samline_list, new_annotation, old_bma_calc)
        if new_score > bma_score:
            return (new_annotation, isNAGood)

    return (None, isNAGood)



# Calculates a score that determines how well an alignment matches an annotation
# Old score calculation return the number of bases on which alignment and annotation match
# New score calculation also subtracts the number of bases from alignment that fall outside the annotation
# Default is NEW score calculation!
def calc_annotation_match_score(samline_list, cannotation, old_bma_calc = False):
    score = 0

    for samline in samline_list:
        start = samline.pos
        reflength = samline.CalcReferenceLengthFromCigar()
        end = start + reflength
        slBasesInside = 0
        for item in cannotation.items:
            bases = item.basesInside(start, end)
            score += bases
            slBasesInside += bases
            if not old_bma_calc:
                # KK: Punishing bases outside the gene item
                score -= reflength - slBasesInside

    return score


# A wrapper around function eval_mapping_part, for multiprocessing
# This function is called inside a separate process
def eval_mapping_part_MP(proc_id, samlines, annotations, paramdict, chromname2seq, out_q):

    sys.stdout.write('\nStarting process %d...\n' % proc_id)

    (report, expressed_genes, gene_coverage) = eval_mapping_part(samlines, annotations, paramdict, chromname2seq)

    out_q.put([report, expressed_genes, gene_coverage])
    sys.stdout.write('\nEnding process %d...\n' % proc_id)




# A function that takes samlines and annotations and evaluates the mappings
# The function that does most of the evaluation work!!!
# It also determines potential new annotations
def eval_mapping_part(samlines, annotations, paramdict, chromname2seq):

    allowed_inacc = Annotation_formats.DEFAULT_ALLOWED_INACCURACY       # Allowing some shift in positions
    min_overlap = Annotation_formats.DEFAULT_MINIMUM_OVERLAP            # Minimum overlap that is considered

    reads2check = []

    processChromNames = True
    if '--leave_chrom_names' in paramdict:
        processChromNames = False

    calcNewAnnotations = 0
    if '--calc_new_annotations' in paramdict:
        calcNewAnnotations = int(paramdict['--calc_new_annotations'][0])

    # Setting allowed_inaccuracy from parameters
    if '--allowed_inacc' in paramdict:
        allowed_inacc = int(paramdict['--allowed_inacc'][0])
    elif '-ai' in paramdict:
        allowed_inacc = int(paramdict['-ai'][0])

    # Setting minimum overlap from parameters
    if '--min_overlap' in paramdict:
        min_overlap = int(paramdict['--min_overlap'][0])
    elif '-mo' in paramdict:
        min_overlap = int(paramdict['-mo'][0])

    expressed_genes = {}
    gene_coverage = {}
    new_annotations = []
    for annotation in annotations:
        expressed_genes[annotation.genename] = [0 for i in xrange(len(annotation.items) + 1)]
        gene_coverage[annotation.genename] = [0 for i in xrange(len(annotation.items) + 1)]

    report = EvalReport(ReportType.TEMP_REPORT)

    check_strand = True
    if '--no_check_strand' in paramdict:
        check_strand = False

    calculate_expression = False
    if '-ex' in paramdict or '--expression' in paramdict:
        calculate_expression = True

    old_bma_calc = False
    if '--old_bma_calc' in paramdict:
        old_bma_calc = True

    per_base_stats = True
    if '--no_per_base_stats' in paramdict:
        per_base_stats = False

    checkReads = False
    if '--check_reads' in paramdict:
        checkReads = True

    num_hithalfbases = 0

    for samline_list in samlines:
        # Initializing information for a single read
        genescovered = []   # genes covered by an alignment
        badsplit = False
        hit = False
        exonHit = False
        isGood = False
        isSpliced = False
        exon_cnt = 0        # counting exons spanned by an alignement
        gene_cnt = 0        # counting genes spanned by an alignement
        num_alignments = len(samline_list)
        num_misses = 0
        isAlmostGood = False
        max_score = 0
        if num_alignments > 1:
            split = True
        else:
            split = False

        # Assuming that all parts of the split alignment are on the same chromosome
        chromname = getChromName(samline_list[0].rname, processChromNames)
        if chromname not in chromname2seq:
            raise Exception('\nERROR: Unknown chromosome name in SAM file! (chromname:"%s", samline.rname:"%s")' % (chromname, samline_list[0].rname))
        chromidx = chromname2seq[chromname]

        # TODO: Separate code for split and contiguous alignments
        #       Might make the code easier to read

        # PLAN:
        # - calculate reference length for a split read
        # - check for genes that it intersects
        # - then iterate over parts of alignment and exons to evaluate how well the alignment captures the transcript

        # Calculating total alignment reference length for all parts of a split read
        # A distance between the start of the first alignment and the end of the last alignment
        # If all split alignments of the read were sorted according to position, this could be done faster
        readrefstart = -1
        readrefend = -1

        # Calculate readlength from the first alignment (should be the same)
        # and then see how many of those bases were actually aligned
        readlength = samline_list[0].CalcReadLengthFromCigar()
        basesaligned = 0
        
        for samline in samline_list:
            # start = samline.pos
            start = samline.pos
            reflength = samline.CalcReferenceLengthFromCigar()
            end = start + reflength

            if readrefstart == -1 or readrefstart > start:
                readrefstart = start
            if readrefend == -1 or readrefend < end:
                readrefend = end

            if per_base_stats:
                # Using regular expressions to find repeating digit and skipping one character after that
                # Used to separate CIGAR string into individual operations
                # Using regular CIGAR, beacuse extended is not needed here
                cigar = samline.cigar
                pattern = '(\d+)(.)'
                operations = re.findall(pattern, cigar)

                for op in operations:
                    if op[1] in ('M', '='):
                        basesaligned += int(op[0])
                    elif op[1] == 'I':
                        basesaligned += int(op[0])
                    elif op[1] == 'D':
                        pass
                    elif op[1] =='X':
                        basesaligned += int(op[0])
                    elif op[1] in ('N', 'S', 'H', 'P'):
                        pass
                    else:       # Should never happen
                        sys.stderr.write('\nERROR: Invalid CIGAR string operation (%s)' % op[1])

        readreflength = readrefend - readrefstart
        startpos = readrefstart
        endpos = readrefend

        if startpos > endpos:
            sys.stderr.write('\nERROR invalid start/end calculation: %s (%d, %d)' % (samline_list[0].qname, startpos, endpos))

        # Assuming all samlines in samline_list have the same strand
        if samline_list[0].flag & 16 == 0:
            readstrand = Annotation_formats.GFF_STRANDFW
        else:
            readstrand = Annotation_formats.GFF_STRANDRV

        # NOTE: If a samline list overlaps with multiple annotations, find the best match
        # 1 - Find overlapping annotations
        # 2 - Choose annotation with the most bases aligned (should allow other criteria)
        # 3 - Calculate everything only for "the best match" annotation

        # Finding candidate annotations
        candidate_annotations = []
        best_match_annotation = None
        for annotation in annotations:
            # If its the same chromosome, the same strand and the read and the gene overlap, then proceed with analysis
            if chromname == getChromName(annotation.seqname, processChromNames) \
                            and (not check_strand or (readstrand == annotation.strand)) \
                            and annotation.overlapsGene(startpos, endpos):
                candidate_annotations.append(annotation)

        if len(candidate_annotations) > 1:
            # Find the best matching candidate            
            max_score = 0
            for cannotation in candidate_annotations:
                if cannotation.genename not in genescovered:
                    genescovered.append(cannotation.genename)
                    gene_cnt += 1

                score = calc_annotation_match_score(samline_list, cannotation, old_bma_calc)

                if score > max_score:
                    max_score = score
                    best_match_annotation = cannotation

        elif len(candidate_annotations) == 1:
            best_match_annotation = candidate_annotations[0]

        bma_score = max_score       # New variable for best match annotation score, with a more suitable name

        if best_match_annotation is not None:
            # KK: Checing some suspicious annotations!!
            # import pdb
            #pdb.set_trace()
            # if best_match_annotation.genename in annotations2check:
            #     import pdb
            #     pdb.set_trace()
            if checkReads and best_match_annotation.genename in annotations2check:
                reads2check.append(samline_list[0].qname)
                # import pdb
                # pdb.set_trace()

            annotation = best_match_annotation      # So that I dont have to refactor the code

            hit = True

            if num_alignments > len(annotation.items):
                # TODO: BAD split!! Alignment is split, but annotation is not!
                badsplit = True
                # sys.stderr.write('\nWARNING: Bad split alignment with more parts then annotation has exons!\n')

            # Checking if the alignment covering annotations encompasses at least half the read
            # max_score represents the number of bases of the read that are aligned within the candidate annotation
            readlength = samline_list[0].CalcReadLengthFromCigar()
            # sys.stderr.write('\nINFO: Maxscore = %d, readlength = %d' % (max_score, readlength))
            if max_score > (readlength / 2):
                num_hithalfbases += 1

            # Updating gene expression
            # Since all initial values for expression and coverage are zero, this could all probably default to case one
            if annotation.genename in expressed_genes.keys():
                expressed_genes[annotation.genename][0] += 1
                gene_coverage[annotation.genename][0] += annotation.basesInsideGene(startpos, endpos)
            else:
                expressed_genes[annotation.genename][0] = 1
                gene_coverage[annotation.genename][0] = annotation.basesInsideGene(startpos, endpos)

            if annotation.insideGene(startpos, endpos):
                partial = False
            else:
                partial = True

            # Initialize exon hit map and exon complete map (also start and end map)
            # Both have one entry for each exon
            # Hit map collects how many times has each exon been hit by an alignment (it should be one or zero)
            # Complete map collects which exons have been completely covered by an alignement
            # Start map collects which exons are correctly started by an alignment (have the same starting position)
            # End map collects which exons are correctly ended by an alignment (have the same ending position)
            # NOTE: test this to see if it slows the program too much
            exonhitmap = {(i+1):0 for i in xrange(len(annotation.items))}
            exoncompletemap = {(i+1):0 for i in xrange(len(annotation.items))}
            exonstartmap = {(i+1):0 for i in xrange(len(annotation.items))}
            exonendmap = {(i+1):0 for i in xrange(len(annotation.items))}
            num_misses = 0      # The number of partial alignments that do not overlap any exons
                                # Partial alignments that are smaller than allowed_inaccuracy, are not counted
            for samline in samline_list:
                item_idx = 0
                lstartpos = samline.pos
                reflength = samline.CalcReferenceLengthFromCigar()
                lendpos = lstartpos + reflength
                exonhit = False
                for item in annotation.items:
                    item_idx += 1
                    if item.overlapsItem(lstartpos, lendpos, allowed_inacc, min_overlap):
                        exonhit = True
                        exonhitmap[item_idx] += 1
                        if item.equalsItem(lstartpos, lendpos, allowed_inacc, min_overlap):
                            exoncompletemap[item_idx] = 1
                            exonstartmap[item_idx] = 1
                            exonendmap[item_idx] = 1
                        elif item.startsItem(lstartpos, lendpos, allowed_inacc, min_overlap):
                            exonstartmap[item_idx] = 1
                        elif item.endsItem(lstartpos, lendpos, allowed_inacc, min_overlap):
                            exonendmap[item_idx] = 1

                        exon_cnt += 1
                        if calculate_expression:
                            expressed_genes[annotation.genename][item_idx] += 1
                            gene_coverage[annotation.genename][item_idx] += item.basesInside(lstartpos, lendpos)
                        exonHit = True
                        if item.insideItem(lstartpos, lendpos, allowed_inacc, min_overlap):
                            exonPartial = False
                        else:
                            exonPartial = True

                # Checking if a partial alignment hit any exons
                if not exonhit and reflength > allowed_inacc:
                    num_misses += 1
                    # KK: Doesn't work within a separate process
                    # import pdb
                    # pdb.set_trace()

                # TODO: What to do if an exon is partially hit?
                # NOTE: Due to information in hit map and complete map
                #       This information might be unnecessary
                #       It can be deduced from exon maps

            # Analyzing exon maps to extract some statistics
            num_exons = len(annotation.items)
            num_covered_exons = len([x for x in exonhitmap.values() if x > 0])      # Exons are considered covered if they are in the hit map
                                                                                    # This means that they only have to be overlapping with an alignment!
            if num_covered_exons > 0:
                report.num_cover_some_exons += 1    # For alignments covering multiple genes, this will be calculated more than once

            if num_covered_exons == num_exons:
                report.num_cover_all_exons += 1

            num_equal_exons = len([x for x in exoncompletemap.values() if x > 0])
            report.num_equal_exons += num_equal_exons
            report.num_partial_exons += num_covered_exons - num_equal_exons

            # Exons covered by more than one part of a split alignment
            multicover_exons = len([x for x in exonhitmap.values() if x > 1])
            report.num_multicover_exons += multicover_exons

            # Not sure what to do with this
            report.num_undercover_alignments = 0
            report.num_overcover_alignments = 0

            # Exon start and end position
            num_good_starts = len([x for x in exonstartmap.values() if x > 0])
            num_good_ends = len([x for x in exonendmap.values() if x > 0])
            report.num_good_starts += num_good_starts
            report.num_good_ends += num_good_ends

            isGood, isSpliced = isGoodSplitAlignment(exonhitmap, exoncompletemap, exonstartmap, exonendmap)

            # KK: This should probably be included in the isGoodSplitAlignment function
            isAlmostGood = False
            if num_misses > 0:
                if isGood:
                    isAlmostGood = True
                isGood = False

            #KK: If the alignment starts or ends outside annotation, classify it as not good!
            #    This was previously incorrect
            if (startpos < best_match_annotation.start - allowed_inacc) or (endpos > best_match_annotation.end + allowed_inacc):
                if isGood:
                    report.num_extended_good += 1
                    # import pdb
                    # pdb.set_trace()
                isGood = False

        else:
            # No matching annotations were found
            # TODO: Check if anything needs to be done here
            # sys.stderr.write('\nNo matching annotations for %s' % samline_list[0].qname)
            report.num_cover_no_exons += 1

            if isSpliced:
                report.num_possible_spliced_alignment += 1

        new_annotation = None
        isNAGood = False

        # Checking for possible new annotations
        # 1. Fusing exons separated by small introns
        if (calcNewAnnotations & NA_SKIP_INTRONS != 0) and best_match_annotation is not None and not isGood:
            # import pdb
            # pdb.set_trace()
            # Create a new annotation from the best match annotation, that removes short introns.
            (new_annotation, isNAGood) = proposeNewFusedAnnotation(samline_list, best_match_annotation, bma_score, allowed_inacc \
                                                     , min_overlap, old_bma_calc)
            if new_annotation is not None:
                new_annotation.genename = "New annotation %d " % (len(new_annotations)+1)
                new_annotations.append(new_annotation)
                if isNAGood:
                    report.num_correct_new_annotations += 1

        # 2. If the best_match_annotation does not produce a "correct" alignment, maybe it can be improved by combinig it 
        # with other annotation from the candidate annotations set
        # For each exon that is not correctly aligned (first exon can have misaligned beginning and last exon can have 
        # misaligned end, other exons need to have correctly aligned beginning and end)
        # For each exon that is not correctly aligned to the best match annotation, check another annotations to see if 
        # it can be correctly aligned to any of them
        # checking maps with values 0 or 1 for each exon: exonhitmap, exoncompletemap, exonstartmap, exonendmap
        if (calcNewAnnotations & NA_COMBINED_ANNOTATIONS_1 != 0) and best_match_annotation is not None and new_annotation is None and not isGood:
            (new_annotation, isNAGood) = proposeNewCombinedAnnotation(samline_list, best_match_annotation, candidate_annotations \
                                                        , bma_score, allowed_inacc, min_overlap, old_bma_calc)
            if new_annotation is not None:
                new_annotation.genename = "New annotation %d " % (len(new_annotations)+1)
                new_annotations.append(new_annotation)
                if isNAGood:
                    report.num_correct_new_annotations += 1
        elif (calcNewAnnotations & NA_COMBINED_ANNOTATIONS_2 != 0) and best_match_annotation is not None and new_annotation is None and not isGood:
            (new_annotation, isNAGood) = proposeNewCombinedAnnotation2(samline_list, best_match_annotation, candidate_annotations \
                                                        , bma_score, allowed_inacc, min_overlap, old_bma_calc)
            if new_annotation is not None:
                new_annotation.genename = "New annotation %d " % (len(new_annotations)+1)
                new_annotations.append(new_annotation)
                if isNAGood:
                    report.num_correct_new_annotations += 1

        report.num_halfbases_hit = num_hithalfbases

        if num_misses > 0:
            report.num_partial_exon_miss += 1

        if isGood:
            report.num_good_alignment += 1
            report.contig_names.append(samline_list[0].qname)
            if isSpliced:
                report.num_hit_all += 1
        else:
            report.num_bad_alignment += 1

        if isAlmostGood:
            report.num_almost_good += 1

        if exon_cnt > 1:
            report.num_multi_exon_alignments += 1
        elif exon_cnt == 0:
            report.num_cover_no_exons += 1

        if len(genescovered) > 1:
            report.num_multi_gene_alignments += 1

        if badsplit:
            report.num_bad_split_alignments += 1

            # This is obsolete
            # Partial alignment hits are not calculated any more
            # if hit and not partial:
            #     report.num_hit_alignments += 1
            # elif hit and partial:
            #     report.num_partial_alignments += 1
            # else:
            #     report.num_missed_alignments += 1

        if hit:
            report.num_hit_alignments += 1
        else:
            report.num_missed_alignments += 1

            # This is obsolete
            # Partial alignment hits are not calculated any more
            # if exonHit and not exonPartial:
            #     report.num_exon_hit += 1
            # elif exonHit and exonPartial:
            #     report.num_exon_partial += 1
            # else:
            #     report.num_exon_miss += 1

        if exonHit:
            report.num_exon_hit += 1
            report.hitone_names.append(samline_list[0].qname)
        else:
            report.num_exon_miss += 1
            report.incorr_names.append(samline_list[0].qname)

        if hit and not exonHit:
            report.num_inside_miss_alignments += 1

            # IMPORTANT: Double counted !!!!!
            # if len(genescovered) == 1 and not badsplit:
            #     report.num_good_alignment += 1
            # else:
            #    report.num_bad_alignment += 1

        if per_base_stats:
            if isGood:
                report.sum_read_length_correct += readlength
                report.sum_bases_aligned_correct += basesaligned
            if exonHit:
                report.sum_read_length_hitone += readlength
                report.sum_bases_aligned_hitone += basesaligned

    report.pot_new_annotations = new_annotations

    if len(reads2check) > 0:
        report.reads2check = reads2check
    
    return (report, expressed_genes, gene_coverage)



# Calculating chosen quality statistics
# Separataing it from other analysis for clearer code
def calculateQualityStats(report, samlines):
    numq = 0
    sumq = 0.0

    for samline_list in samlines:
        for samline in samline_list:
            quality = samline.chosen_quality
            if quality > 0:
                report.num_good_quality += 1
                if report.max_mapping_quality == 0 or report.max_mapping_quality < quality:
                    report.max_mapping_quality = quality
                if report.min_mapping_quality == 0 or report.min_mapping_quality > quality:
                    report.min_mapping_quality = quality
                numq += 1
                sumq += quality
            else:
                report.num_zero_quality += 1

    if numq > 0:
        report.avg_mapping_quality = sumq / numq


# Calculating general per-base mapping statistics
# Match/Mismatch/Insert/Delete
def calculateGeneralMappingStats(report, samlines, seqs, chromname2seq, processChromNames, correct_gm=False):
    numMatch = 0
    numMisMatch = 0
    numInsert = 0
    numDelete = 0
    numLowMatchCnt = 0

    total_read_length = 0
    total_bases_aligned = 0
    percent_bases_aligned = 0.0

    sys.stderr.write('\n(%s) Analyzing CIGAR strings ...  ' % datetime.now().time().isoformat())
    sys.stderr.write('\nProgress: | 1 2 3 4 5 6 7 8 9 0 |')
    sys.stderr.write('\nProgress: | ')
    numsamlines = len(samlines)
    progress = 0
    currentbar = 0.1
    for samline_list in samlines:
        # Calculating progress
        progress += 1
        if float(progress)/numsamlines >= currentbar:
            sys.stderr.write('* ')
            currentbar += 0.1
        
        # For checking cigar strings
        t_numMatch = 0
        t_numInsert = 0
        t_numDelete = 0
        t_numMisMatch = 0

        # Calculate readlength from the first alignment (should be the same)
        # and then see how many of those bases were actually aligned
        readlength = samline_list[0].CalcReadLengthFromCigar()
        basesaligned = 0
        for samline in samline_list:
            chromname = getChromName(samline.rname, processChromNames)
            if chromname not in chromname2seq:
                # import pdb
                # pdb.set_trace()
                raise Exception('\nERROR: Unknown chromosome name in SAM file! (chromname:"%s", samline.rname:"%s")' % (chromname, samline.rname))
            chromidx = chromname2seq[chromname]

            # Testing code
            try:
                if correct_gm and samline.flag & 16 != 0:
                    samline.pos += 1
                cigar = samline.CalcExtendedCIGAR(seqs[chromidx])
                pos = samline.pos
                quals = samline.qual

                # Using regular expressions to find repeating digit and skipping one character after that
                # Used to separate CIGAR string into individual operations
                pattern = '(\d+)(.)'
                operations = re.findall(pattern, cigar)

                for op in operations:
                    if op[1] in ('M', '='):
                        t_numMatch += int(op[0])
                        basesaligned += int(op[0])
                    elif op[1] == 'I':
                        t_numInsert += int(op[0])
                        basesaligned += int(op[0])
                    elif op[1] == 'D':
                        t_numDelete += int(op[0])
                    elif op[1] =='X':
                        t_numMisMatch += int(op[0])
                        basesaligned += int(op[0])
                    elif op[1] in ('N', 'S', 'H', 'P'):
                        pass
                    else:
                        sys.stderr.write('\nERROR: Invalid CIGAR string operation (%s)' % op[1])
            except Exception, e:
                # import pdb
                # pdb.set_trace()
                sys.stderr.write('ERROR: querry/ref/pos/message = %s/%s/%d/%s \n' % (samline.qname, samline.rname, samline.pos, e.message))
                pass

        numMatch += t_numMatch
        numInsert += t_numInsert
        numDelete += t_numDelete
        numMisMatch += t_numMisMatch

        # Checking CIGAR strings for low match reads
        if (t_numMatch < t_numMisMatch + t_numInsert + t_numDelete):
            strand = '+'
            if samline_list[0].flag & 16 != 0:
                strand = '-'
            numLowMatchCnt += 1
            # sys.stderr.write('\nDEBUG: strand / match / mismatch / insert / delete: %c / %d / %d / %d / %d' % (strand, t_numMatch, t_numMisMatch, t_numInsert, t_numDelete))
            # import pdb
            # pdb.set_trace()
        else:
            pass

        total_read_length += readlength
        total_bases_aligned += basesaligned
        if basesaligned > readlength:
            raise Exception('\nERROR counting aligned and total bases!')
            # TODO: See what happens here
            pass

    # Closing progress bar
    sys.stderr.write('|')
    sys.stderr.write('\nDone!')

    report.num_match = numMatch
    report.num_mismatch = numMisMatch
    report.num_insert = numInsert
    report.num_delete = numDelete

    total = numMatch + numMisMatch + numInsert + numDelete

    report.num_lowmatchcnt = numLowMatchCnt

    if total > 0:
        report.match_percentage = float(report.num_match)/total
        report.mismatch_percentage = float(report.num_mismatch)/total
        report.insert_percentage = float(report.num_insert)/total
        report.delete_percentage = float(report.num_delete)/total

    if total_read_length == 0:
        percent_bases_aligned = -1
    else:
        percent_bases_aligned = 100 * float(total_bases_aligned) / total_read_length

    report.sum_read_length = total_read_length
    report.sum_bases_aligned = total_bases_aligned
    report.percent_bases_aligned = percent_bases_aligned



# TODO: Refactor code, place some code in functions
#       Rewrite analyzing SAM file, detecting multi and split alignments
def eval_mapping_annotations(ref_file, sam_file, annotations_file, paramdict):

    sys.stderr.write('\n')
    sys.stderr.write('\n(%s) START: Evaluating mapping with annotations:' % datetime.now().time().isoformat())

    processChromNames = True
    if '--leave_chrom_names' in paramdict:
        processChromNames = False

    calcNewAnnotations = 0
    if '--calc_new_annotations' in paramdict:
        calcNewAnnotations = paramdict['--calc_new_annotations'][0]

    report = EvalReport(ReportType.MAPPING_REPORT)
    if '-ex' in paramdict or '--expression' in paramdict:
        report.output_gene_expression = True

    check_strand = True
    if '--no_check_strand' in paramdict:
        check_strand = False

    per_base_stats = True
    if '--no_per_base_stats' in paramdict:
        per_base_stats = False

    correct_gm = False
    if '--graphmap' in paramdict:
        sys.stderr.write('\n(%s) Using option --graphmap ... ' % datetime.now().time().isoformat())        
        correct_gm = True

    MP = True
    if '--no_multiprocessing' in paramdict:
        sys.stderr.write('\n(%s) Using option --no_multiprocessing ... ' % datetime.now().time().isoformat())        
        MP = False

    sys.stderr.write('\n(%s) Loading and processing FASTA reference ... ' % datetime.now().time().isoformat())
    [chromname2seq, headers, seqs, quals] = load_and_process_reference(ref_file, paramdict, report)

    sys.stderr.write('\n(%s) Loading and processing SAM file with mappings ... ' % datetime.now().time().isoformat())
    samlines = load_and_process_SAM(sam_file, paramdict, report)

    # The number of alignments (alignment group) after preprocessing
    report.num_evaluated_alignments = len(samlines)

    sys.stderr.write('\n(%s) Loading and processing annotations file ... ' % datetime.now().time().isoformat())
    annotations, expressed_genes, gene_coverage = load_and_process_annotations(annotations_file, paramdict, report)

    sys.stderr.write('\n(%s) Analyzing mappings against annotations ... ' % datetime.now().time().isoformat())
    # Looking at SAM lines to estimate general mapping quality
    # TODO: This used to take a long time, but I managed to speed it up
    #       should be looked at a bit more to see if additional improvements could be made.

    sys.stderr.write('\n(%s) Calculating chosen quality statistics ... ' % datetime.now().time().isoformat())
    # Calculating chosen quality statistics
    calculateQualityStats(report, samlines)    

    # Calculating general per-base mapping statistics
    # Match/Mismatch/Insert/Delete
    # Setting up some sort of a progress bar
    if per_base_stats:
        calculateGeneralMappingStats(report, samlines, seqs, chromname2seq, processChromNames, correct_gm)


    # NOTE: refactoring the code for multiprocessing
    #     - Each process will handle reads and annotations for a single chomosome and strand
    # Workflow:
    # 1. Sort Annotations and mappings (SAM fle) according to chromosome and strand
    # 2. For each dataset spawn a new processing
    # 3. Collect data returned by multiple processes

    if MP:
        # Separating Annotations and Mappings (SAM lines) according to chromosome and strand
        partlist = []       # A list of keys of parts for processing
                            # Each part represents a single chromosome strand
        part_samlines = {}          # A dictionarry containing a list (or deper hierarchy) of samlines for each part
        part_annotations = {}       # A dictionarry containing a list of annotations for each part

        # If separating for the strand
        if check_strand:
            for chromname in chromname2seq.keys():      # Initializing
                partlist.append(chromname + '+')
                partlist.append(chromname + '-')
                part_samlines[chromname + '+'] = []
                part_samlines[chromname + '-'] = []
                part_annotations[chromname + '+'] = []
                part_annotations[chromname + '-'] = []

            # Separating SAM lines
            for samline_list in samlines:
                samline = samline_list[0]           # Looking only at the first samline in the list
                                                    # Due to previous processing, assuming that all
                                                    # others correspond to the same chromosome and strand
                partname = ''
                chromname = getChromName(samline.rname, processChromNames)
                if samline.flag & 16 == 0:
                    readstrand = Annotation_formats.GFF_STRANDFW
                    partname = chromname + '+'
                else:
                    readstrand = Annotation_formats.GFF_STRANDRV
                    partname = chromname + '-'

                part_samlines[partname].append(samline_list)

            # Separating annotations expressed genes and gene coverage
            for annotation in annotations:
                partname = ''
                chromname = getChromName(annotation.seqname, processChromNames)
                if annotation.strand == Annotation_formats.GFF_STRANDFW:
                    partname = chromname + '+'
                else:
                    partname = chromname + '-'
                part_annotations[partname].append(annotation)
                genename = annotation.genename

        # Not separating according to the strand
        else:
            for chromname in chromname2seq.keys():      # Initializing
                partlist.append(chromname)
                part_samlines[chromname] = []
                part_annotations[chromname] = []

            # Separating SAM lines
            for samline_list in samlines:
                samline = samline_list[0]           # Looking only at the first samline in the list
                                                    # Due to previous processing, assuming that all
                                                    # others correspond to the same chromosome and strand
                try:
                    chromname = getChromName(samline.rname, processChromNames)
                    partname = chromname
                    part_samlines[partname].append(samline_list)
                except Exception:
                    import pdb
                    pdb.set_trace()


            # Separating annotations expressed genes and gene coverage
            for annotation in annotations:
                chromname = getChromName(annotation.seqname, processChromNames)
                partname = chromname
                try:
                    # simply skipping annotations that are not on the reference
                    if partname in part_annotations:
                        part_annotations[partname].append(annotation)
                except Exception:
                    import pdb
                    pdb.set_trace()

                genename = annotation.genename

        # import pdb
        # pdb.set_trace()

        # Spawning and calling processes
        out_q = multiprocessing.Queue()
        jobs = []
        proc_id = 0
        for partname in partlist:
            proc_id += 1
            t_samlines = part_samlines[partname]
            t_annotations = part_annotations[partname]
            proc = multiprocessing.Process(name=partname, target=eval_mapping_part_MP, args=(proc_id, t_samlines, t_annotations, paramdict, chromname2seq, out_q,))
            jobs.append(proc)
            proc.start()

        # TODO: Summarize results from different processes
        sys.stderr.write('\n(%s) Collecting results!' % datetime.now().time().isoformat())

        expressed_genes = {}
        gene_coverage = {}
        for i in xrange(len(jobs)):
            [t_report, t_expressed_genes, t_gene_coverage] = out_q.get()
            expressed_genes.update(t_expressed_genes)
            gene_coverage.update(t_gene_coverage)
            report.num_cover_some_exons += t_report.num_cover_some_exons
            report.num_cover_all_exons += t_report.num_cover_all_exons
            report.num_equal_exons += t_report.num_equal_exons
            report.num_partial_exons += t_report.num_partial_exons
            report.num_multicover_exons += t_report.num_multicover_exons
            report.num_undercover_alignments = t_report.num_undercover_alignments
            report.num_overcover_alignments = t_report.num_overcover_alignments
            report.num_good_starts += t_report.num_good_starts
            report.num_good_ends += t_report.num_good_ends
            report.num_possible_spliced_alignment += t_report.num_possible_spliced_alignment
            report.num_good_alignment += t_report.num_good_alignment
            report.num_bad_alignment += t_report.num_bad_alignment
            report.num_multi_exon_alignments += t_report.num_multi_exon_alignments
            report.num_cover_no_exons += t_report.num_cover_no_exons
            report.num_multi_gene_alignments += t_report.num_multi_gene_alignments
            report.num_bad_split_alignments += t_report.num_bad_split_alignments
            report.num_hit_alignments += t_report.num_hit_alignments
            report.num_partial_alignments += t_report.num_partial_alignments
            report.num_missed_alignments += t_report.num_missed_alignments
            report.num_exon_hit += t_report.num_exon_hit
            report.num_exon_partial += t_report.num_exon_partial
            report.num_exon_miss += t_report.num_exon_miss
            report.num_halfbases_hit += t_report.num_halfbases_hit
            report.num_lowmatchcnt = t_report.num_lowmatchcnt
            report.num_inside_miss_alignments += t_report.num_inside_miss_alignments
            report.num_partial_exon_miss += t_report.num_partial_exon_miss
            report.num_almost_good += t_report.num_almost_good
            report.num_extended_good += t_report.num_extended_good
            report.num_hit_all += t_report.num_hit_all
            report.hitone_names += t_report.hitone_names
            report.hithalfbases_names += t_report.hithalfbases_names
            report.contig_names += t_report.contig_names
            report.incorr_names += t_report.incorr_names
            report.unmapped_names += t_report.unmapped_names
            report.pot_new_annotations += t_report.pot_new_annotations
            report.alignments_with_pna = len(report.pot_new_annotations)
            report.num_correct_new_annotations += t_report.num_correct_new_annotations

            report.sum_read_length_correct += t_report.sum_read_length_correct
            report.sum_bases_aligned_correct += t_report.sum_bases_aligned_correct
            report.sum_read_length_hitone += t_report.sum_read_length_hitone
            report.sum_bases_aligned_hitone += t_report.sum_bases_aligned_hitone

            report.reads2check += t_report.reads2check

        # Wait for all processes to end
        for proc in jobs:
            proc.join()

    else:
        # Single process version
        (t_report, expressed_genes, gene_coverage) = eval_mapping_part(samlines, annotations, paramdict, chromname2seq)
        report.num_cover_some_exons += t_report.num_cover_some_exons
        report.num_cover_all_exons += t_report.num_cover_all_exons
        report.num_equal_exons += t_report.num_equal_exons
        report.num_partial_exons += t_report.num_partial_exons
        report.num_multicover_exons += t_report.num_multicover_exons
        report.num_undercover_alignments = t_report.num_undercover_alignments
        report.num_overcover_alignments = t_report.num_overcover_alignments
        report.num_good_starts += t_report.num_good_starts
        report.num_good_ends += t_report.num_good_ends
        report.num_possible_spliced_alignment += t_report.num_possible_spliced_alignment
        report.num_good_alignment += t_report.num_good_alignment
        report.num_bad_alignment += t_report.num_bad_alignment
        report.num_multi_exon_alignments += t_report.num_multi_exon_alignments
        report.num_cover_no_exons += t_report.num_cover_no_exons
        report.num_multi_gene_alignments += t_report.num_multi_gene_alignments
        report.num_bad_split_alignments += t_report.num_bad_split_alignments
        report.num_hit_alignments += t_report.num_hit_alignments
        report.num_partial_alignments += t_report.num_partial_alignments
        report.num_missed_alignments += t_report.num_missed_alignments
        report.num_exon_hit += t_report.num_exon_hit
        report.num_exon_partial += t_report.num_exon_partial
        report.num_exon_miss += t_report.num_exon_miss
        report.num_halfbases_hit += t_report.num_halfbases_hit
        report.num_lowmatchcnt = t_report.num_lowmatchcnt
        report.num_inside_miss_alignments += t_report.num_inside_miss_alignments
        report.num_partial_exon_miss += t_report.num_partial_exon_miss
        report.num_almost_good += t_report.num_almost_good
        report.num_extended_good += t_report.num_extended_good
        report.num_hit_all += t_report.num_hit_all
        report.hitone_names += t_report.hitone_names
        report.hithalfbases_names += t_report.hithalfbases_names
        report.contig_names += t_report.contig_names
        report.incorr_names += t_report.incorr_names
        report.unmapped_names += t_report.unmapped_names
        report.pot_new_annotations += t_report.pot_new_annotations
        report.alignments_with_pna = len(report.pot_new_annotations)
        report.num_correct_new_annotations += t_report.num_correct_new_annotations

        report.sum_read_length_correct += t_report.sum_read_length_correct
        report.sum_bases_aligned_correct += t_report.sum_bases_aligned_correct
        report.sum_read_length_hitone += t_report.sum_read_length_hitone
        report.sum_bases_aligned_hitone += t_report.sum_bases_aligned_hitone

        report.reads2check += t_report.reads2check

    # Annotations to check are stored globally
    report.annotations2check = annotations2check

    if report.sum_read_length_correct > 0:
        report.percent_bases_aligned_correct = 100.0 * float(report.sum_bases_aligned_correct) / report.sum_read_length_correct
    else :
        report.percent_bases_aligned_correc = -1.0
    if report.sum_read_length_hitone > 0:
        report.percent_bases_aligned_hitone = 100.0 * float(report.sum_bases_aligned_hitone) / report.sum_read_length_hitone
    else:
        report.percent_bases_aligned_hitone = -1.0

    # Collecting new annotations (there should be a lot of duplicates)
    # New annotations details are written to a file: annotations.report
    if calcNewAnnotations:
        collected_new_annotations = []
        cna_count = {}
        cna_readlist = {}
        for pna in report.pot_new_annotations:
            alreadyExists = False
            for cna in collected_new_annotations:
                if pna.itemsEqual(cna):
                    alreadyExists = True
                    cna.seqname += ', ' + pna.seqname
                    cna_count[cna.genename] += 1
                    cna_readlist[cna.genename].append(pna.seqname)
                    break

            if not alreadyExists:
                collected_new_annotations.append(pna)
                cna_count[pna.genename] = 0
                cna_readlist[pna.genename] = [pna.seqname]

        # Preserve only annotations that are generated by at least NEW_ANNOTATION_MIN reads
        filtered_new_annotations = []
        for cna in collected_new_annotations:
            name = cna.genename
            if cna_count[name] >= NEW_ANNOTATION_MIN:
                filtered_new_annotations.append(cna)

        report.pot_new_annotations = filtered_new_annotations
        report.cna_count = cna_count
        report.cna_readlist = cna_readlist
        report.detect_new_annotations = True


    # KK: +1s are for testing to avoid division by zero
    report.good_alignment_percent = 100.0 * float(report.num_good_alignment)/(report.num_good_alignment + report.num_bad_alignment + 1)
    report.bad_alignment_percent = 100.0 * float(report.num_bad_alignment)/(report.num_good_alignment + report.num_bad_alignment + 1)
    report.hit_all_percent = 100.0 * float(report.num_hit_all)/(report.num_good_alignment + report.num_bad_alignment + 1)


    # How many genes were covered by alignments
    report.num_genes_covered = 0
    report.num_exons_covered = 0
    for genecnt in expressed_genes.itervalues():
        if genecnt[0] > 0:
            report.num_genes_covered += 1
        for cnt in genecnt[1:]:
            if cnt > 0:
                report.num_exons_covered += 1

    # Pass gene expression and coverage information to report
    report.expressed_genes = expressed_genes
    report.gene_coverage = gene_coverage

    sys.stderr.write('\n(%s) Done!' % datetime.now().time().isoformat())
    sys.stderr.write('\n')

    return report



def eval_mapping_fasta(ref_file, sam_file, paramdict):

    sys.stderr.write('\n')
    sys.stderr.write('\n(%s) START: Evaluating mapping with FASTA reference only:' % datetime.now().time().isoformat())

    report = EvalReport(ReportType.FASTA_REPORT)

    sys.stderr.write('\n(%s) Loading and processing FASTA reference ... ' % datetime.now().time().isoformat())
    [chromname2seq, headers, seqs, quals] = load_and_process_reference(ref_file, paramdict, report)

    sys.stderr.write('\n(%s) Loading and processing SAM file with mappings ... ' % datetime.now().time().isoformat())
    samlines = load_and_process_SAM(sam_file, paramdict, report)

    sys.stderr.write('\n(%s) Calculating chosen quality statistics ... ' % datetime.now().time().isoformat())
    # Calculating chosen quality statistics
    calculateQualityStats(report, samlines)

    # Analyzing mappings
    sys.stderr.write('\n(%s) Analyzing mappings against FASTA reference ... ' % datetime.now().time().isoformat())

    per_base_stats = True
    if '--no_per_base_stats' in paramdict:
        per_base_stats = False

    processChromNames = True
    if '--leave_chrom_names' in paramdict:
        processChromNames = False

    correct_gm = False
    if '--graphmap' in paramdict:
        sys.stderr.write('\n(%s) Using option --graphmap ... ' % datetime.now().time().isoformat())        
        correct_gm = True

    if per_base_stats:
        calculateGeneralMappingStats(report, samlines, seqs, chromname2seq, processChromNames, correct_gm)

    sys.stderr.write('\n(%s) Done!' % datetime.now().time().isoformat())
    sys.stderr.write('\n')

    return report


def eval_mapping(ref_file, sam_file, paramdict):

    out_filename = ''
    hitone_filename = ''
    hithalfbases_filename = ''
    contig_filename = ''
    incorr_filename = ''
    unmapped_filename = ''
    out_file = None

    save_qnames = False
    if '-sqn' in paramdict or '--save_query_names' in paramdict:
        save_qnames = True
        if '-o' not in paramdict and '--output' not in paramdict or '-a' not in paramdict:
            sys.stderr.write('\nInvalid parameters. Paramater --save_query_names must be used with paramters --output and -a')
            exit()

    calcNewAnnotations = 0
    if '--calc_new_annotations' in paramdict:
        calcNewAnnotations = paramdict['--calc_new_annotations'][0]

    if '-o' in paramdict:
        out_filename = paramdict['-o'][0]
    elif '--output' in paramdict:
        out_filename = paramdict['--output'][0]

    if save_qnames:
        hitone_filename = out_filename + '_hit1.names'
        hithalfbases_filename = out_filename + '_hithalfbases.names'
        contig_filename = out_filename + '_ctg.names'
        incorr_filename = out_filename + '_bad.names'
        unmapped_filename = out_filename + '_unmmapped.names'

    if out_filename != '':
        out_file = open(out_filename, 'w+')
    else:
        out_file = sys.stdout


    if '-a' in paramdict:
        annotations_file = paramdict['-a'][0]
        report = eval_mapping_annotations(ref_file, sam_file, annotations_file, paramdict)
    else:
        report = eval_mapping_fasta(ref_file, sam_file, paramdict)

    report.commandline = paramdict['command']

    if calcNewAnnotations:
        annotation_report_filename = out_filename + '_annotation.report'
        with open(annotation_report_filename, 'w+') as annr_file:
            annr_file.write(report.getAnnotationReport())
            annr_file.close()

    # For debugging purposes
    if len(report.reads2check) > 0:
        with open('reads2check.txt', 'w+') as r2c_file:
            for line in report.reads2check:
                r2c_file.write(line + '\n')
            r2c_file.close()

        with open('annotations2check.txt', 'w+') as a2c_file:
            for line in report.annotations2check:
                a2c_file.write(line + '\n')
            a2c_file.close()        

    if save_qnames:
        with open(hitone_filename, 'w+') as hitone_file:
            hitone_file.write(report.get_hitone_names())
            hitone_file.close()

        with open(hithalfbases_filename, 'w+') as hithalf_file:
            hithalf_file.write(report.get_hithalfbases_names())
            hithalf_file.close()

        with open(contig_filename, 'w+') as contig_file:
            contig_file.write(report.get_contig_names())
            contig_file.close()

        with open(incorr_filename, 'w+') as incorr_file:
            incorr_file.write(report.get_incorr_names())
            incorr_file.close()

        with open(unmapped_filename, 'w+') as unmapped_file:
            unmapped_file.write(report.get_unmapped_names())
            unmapped_file.close()

    out_file.write(report.toString())
    out_file.close()


def eval_annotations(annotations_file, paramdict):

    out_filename = ''
    out_file = None

    if '-o' in paramdict:
        out_filename = paramdict['-o'][0]
    elif '--output' in paramdict:
        out_filename = paramdict['--output'][0]

    if out_filename != '':
        out_file = open(out_filename, 'w+')
    else:
        out_file = sys.stdout

    sys.stderr.write('\n')
    sys.stderr.write('\n(%s) START: Evaluating mapping with FASTA reference only:' % datetime.now().time().isoformat())

    report = EvalReport(ReportType.ANNOTATION_REPORT)

    sys.stderr.write('\n(%s) Loading and processing annotations file ... ' % datetime.now().time().isoformat())
    annotations, expressed_genes, gene_coverage = load_and_process_annotations(annotations_file, paramdict, report)

    report.commandline = paramdict['command']

    # Analyzing annotations to discover alternate splicings
    # Groupign annotations which overlap and are on the same strand
    start_new_group = True
    grouped_annotations = []

    for new_annotation in annotations:
        if start_new_group:
            annotation_group = []
            annotation_group.append(new_annotation)
            group_start = new_annotation.start
            group_end = new_annotation.end
            group_strand = new_annotation.strand
            group_chrom = new_annotation.seqname
            start_new_group = False
        else:
            if new_annotation.overlapsGene(group_start, group_end) and group_strand == new_annotation.strand and group_chrom == new_annotation.seqname:
                # Add annotation to current group
                annotation_group.append(new_annotation)
            else:
                # Save the current group and start the new one
                grouped_annotations.append(annotation_group)
                annotation_group = []
                annotation_group.append(new_annotation)
                group_start = new_annotation.start
                group_end = new_annotation.end
                group_strand = new_annotation.strand
                group_chrom = new_annotation.seqname


    # At the end, add last group if it exists
    if len(annotation_group) > 0:
        grouped_annotations.append(annotation_group)

    report.num_annotation_groups = len(grouped_annotations)

    # Analyze groupped annotations
    for ann_group in grouped_annotations:
        num_spliced_alignments = len(ann_group)
        if num_spliced_alignments > 1:
            group_name = ann_group[0].genename + ('(%d)' % len(ann_group))
            group_splicing_info = ''
            report.num_alternate_spliced_genes += 1
            if report.max_spliced_alignments == 0 or report.max_spliced_alignments < num_spliced_alignments:
                report.max_spliced_alignments = num_spliced_alignments
            if report.min_spliced_alignments == 0 or report.min_spliced_alignments > num_spliced_alignments:
                report.min_spliced_alignments = num_spliced_alignments
            for annotation in ann_group:
                num_exons = len(annotation.items)
                transcript = annotation.genename
                group_splicing_info += '%s(%d), ' % (transcript, num_exons)
                if report.max_spliced_exons == 0 or report.max_spliced_exons < num_exons:
                    report.max_spliced_exons = num_exons
                if report.min_spliced_exons == 0 or report.min_spliced_exons > num_exons:
                    report.min_spliced_exons = num_exons

            report.alternate_splicing[group_name] = group_splicing_info

    if '-as' in paramdict or '--alternate_splicing' in paramdict:
        report.output_alternate_splicing = True
    else:
        report.output_alternate_splicing = False

    out_file.write(report.toString())


def eval_maplength(sam_file, paramdict):

    out_filename = ''
    out_file = None

    if '-o' in paramdict:
        out_filename = paramdict['-o'][0]
    elif '--output' in paramdict:
        out_filename = paramdict['--output'][0]

    if out_filename != '':
        out_file = open(out_filename, 'w+')
    else:
        out_file = sys.stdout

    # Writing headers to output file
    out_file.write('QNAME,RNAME,read length,bases aligned\n')

    sys.stderr.write('\n')
    sys.stderr.write('\n(%s) START: Evaluating mapping Lengths:' % datetime.now().time().isoformat())

    # Do not really need a report for this but using one because a function requires it
    report = EvalReport(ReportType.MAPPING_REPORT)

    # I'm still not sure if I need a reference for this, or just SAM file is enough
    # sys.stderr.write('\n(%s) Loading and processing FASTA reference ... ' % datetime.now().time().isoformat())
    # [chromname2seq, headers, seqs, quals] = load_and_process_reference(ref_file, paramdict, report)

    sys.stderr.write('\n(%s) Loading and processing SAM file with mappings ... ' % datetime.now().time().isoformat())
    samlines = load_and_process_SAM(sam_file, paramdict, report)

    # Setting up some sort of a progress bar
    sys.stderr.write('\n(%s) Analyzing mapping lengths ...  ' % datetime.now().time().isoformat())
    sys.stderr.write('\nProgress: | 1 2 3 4 5 6 7 8 9 0 |')
    sys.stderr.write('\nProgress: | ')
    numsamlines = len(samlines)
    progress = 0
    currentbar = 0.1
    for samline_list in samlines:
        # Callculating progress
        progress += 1
        if float(progress)/numsamlines >= currentbar:
            sys.stderr.write('* ')
            currentbar += 0.1
        # Calculate readlength from the first alignment (should be the same)
        # and then see how many of those bases were actually aligned
        readlength = samline_list[0].CalcReadLengthFromCigar()
        basesaligned = 0
        for samline in samline_list:

            # chromname = getChromName(samline.rname)
            # if chromname not in chromname2seq:
            #     import pdb
            #     pdb.set_trace()
            #     raise Exception('\nERROR: Unknown chromosome name in SAM file! (chromname:"%s", samline.rname:"%s")' % (chromname, samline.rname))
            # chromidx = chromname2seq[chromname]

            # Do not need extended CIGAR for this, regular is enough
            # cigar = samline.CalcExtendedCIGAR(seqs[chromidx])
            # pos = samline.pos
            # quals = samline.qual

            cigar = samline.cigar

            # Using regular expressions to find repeating digit and skipping one character after that
            # Used to separate CIGAR string into individual operations
            pattern = '(\d+)(.)'
            operations = re.findall(pattern, cigar)

            for op in operations:
                if op[1] in ('M', '='):
                    basesaligned += int(op[0])
                elif op[1] == 'I':
                    basesaligned += int(op[0])
                elif op[1] =='X':
                    basesaligned += int(op[0])
                elif op[1] in ('N', 'S', 'H', 'P', 'D'):
                    pass
                else:
                    sys.stderr.write('\nERROR: Invalid CIGAR string operation (%s)' % op[1])

        if basesaligned > readlength:
            # import pdb
            # pdb.set_trace()
            raise Exception('\nERROR counting aligned and total bases!')
            # TODO: See what happens here
            pass

        # Writing bases aligned to output
        rname = samline_list[0].rname
        qname = samline_list[0].qname

        out_file.write('%s,%s,%d,%d\n' % (qname, rname, readlength, basesaligned))

    # Closing progress bar
    out_file.close()
    sys.stderr.write('|')
    sys.stderr.write('\nDone!')

    pass


def verbose_usage_and_exit():
    sys.stderr.write('RNAseqEval - A tool for evaulating RNAseq results.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    # sys.stderr.write('\t\tsetup\n')
    # sys.stderr.write('\t\tcleanup\n')
    sys.stderr.write('\t\teval-mapping\n')
    sys.stderr.write('\t\teval-annotations\n')
    sys.stderr.write('\t\teval-maplength\n')
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if (mode == 'setup'):
        if (len(sys.argv) != 2):
            sys.stderr.write('Setup the folder structures and install necessary tools.\n')
            sys.stderr.write('Requires no additional parameters to run.\n')
            sys.stderr.write('\n')
            exit(1)

        setup_RNAseqEval.setup_all()

    elif (mode == 'cleanup'):
        if (len(sys.argv) != 2):
            sys.stderr.write('Cleans up intermediate files.\n')
            sys.stderr.write('Requires no additional parameters to run.\n')
            sys.stderr.write('\n')
            exit(1)

        cleanup()

    elif (mode == 'eval-mapping'):
        if (len(sys.argv) < 4):
            sys.stderr.write('Evaluates RNAseq mapping from a SAM file.\n')
            sys.stderr.write('Can use annotations if provided.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('%s %s <reference FASTA file> <input SAM file> options\n'% (sys.argv[0], sys.argv[1]))
            sys.stderr.write('options:"\n')
            sys.stderr.write('-a <file> : a reference annotation (GFF/GTF/BED) file\n')
            sys.stderr.write('-o (--output) <file> : output file to which the report will be written\n')
            sys.stderr.write('-ex (--expression) : calculate and output gene expression\n')
            sys.stderr.write('--no_check_strand : when matching alignments to annotations, do not take strand \n')
            sys.stderr.write('                    into account, only chromosome and position\n')
            sys.stderr.write('--no_per_base_stats : do not calculate per-base stats, such are percentage of\n')
            sys.stderr.write('                      matches, mismatches, inserts and deletes\n')
            sys.stderr.write('-sqn (--save_query_names) : save query names for alignments that managed to hit an exon\n')
            sys.stderr.write('                            and for contiguous alignments. Query names are saved to files\n')
            sys.stderr.write('                            with filenames determined from output file. \n')
            sys.stderr.write('                            For this option output and annotation files must be specified\n')
            sys.stderr.write('-ai (--alowed_inaccuracy) : A maximumn distance in bases between an annotation and an alignment\n')
            sys.stderr.write('                            where the alignment is still considered correct (default 5)\n')
            sys.stderr.write('-mo (--min_overlap) : A minimum overlap between an annotation and an alignment that is considered valid\n')
            sys.stderr.write('                      (default 5)\n')
            sys.stderr.write('--graphmap : correct for a bug in GraphMap RNA mapping on reverse strand\n')
            sys.stderr.write('              taken into account only when calculating the percentage of matches\n')
            sys.stderr.write('--old_bma_calc : Calculate best matching annotation only based on maximizing the number of bases an alignment\n')
            sys.stderr.write('                 on an annotation. The number of bases outside an annotation is not take into account in this case\n')
            sys.stderr.write('--leave_chrom_names : do not preprocess chromosome names in SAM, reference and annotation files\n')
            sys.stderr.write('                      Chromosome names are by default preprocessed to make them match the format "chr[designation]"\n')
            sys.stderr.write('--calc_new_annotations [LEVEL]: calculate potential new annotations, if a sufficient number of alignments (default 3)\n')
            sys.stderr.write('                                better fits a combination of exons then any existing annotation, that combination\n')
            sys.stderr.write('                                of exons is suggested as a new annotation\n')
            sys.stderr.write('                        LEVEL:  NA_NO_NEW_ANNOTATINS = 0 (default)\n')
            sys.stderr.write('                                NA_SKIP_INTRONS = 1\n')
            sys.stderr.write('                                NA_COMBINED_ANNOTATIONS_1 = 2\n')
            sys.stderr.write('                                NA_COMBINED_ANNOTATIONS_2 = 4\n')
            sys.stderr.write('--no_multiprocessing: disable multiprocessing - slower execution, but usefull for testing\n')
            sys.stderr.write('\n')
            exit(1)

        ref_file = sys.argv[2]
        sam_file = sys.argv[3]

        pparser = paramsparser.Parser(paramdefs)
        paramdict = pparser.parseCmdArgs(sys.argv[4:])

        ref_file = sys.argv[2]
        sam_file = sys.argv[3]

        paramdict['command'] = ' '.join(sys.argv)

        eval_mapping(ref_file, sam_file, paramdict)

    elif (mode == 'eval-annotations'):
        if (len(sys.argv) < 3):
            sys.stderr.write('Evaluates gene annotation from a BED or GTF/GFF file.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('%s %s <annotations file> options\n'% (sys.argv[0], sys.argv[1]))
            sys.stderr.write('options:"\n')
            sys.stderr.write('-o (--output) <file> : output file to which the report will be written\n')
            sys.stderr.write('\n')
            exit(1)


        annotation_file = sys.argv[2]

        pparser = paramsparser.Parser(paramdefs)
        paramdict = pparser.parseCmdArgs(sys.argv[3:])
        paramdict['command'] = ' '.join(sys.argv)

        eval_annotations(annotation_file, paramdict)

    elif (mode == 'eval-maplength'):
        if (len(sys.argv) < 3):
            sys.stderr.write('Evaluates mapping length returning mapped percentage for each read.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('%s %s <input SAM file> options\n'% (sys.argv[0], sys.argv[1]))
            sys.stderr.write('options:"\n')
            sys.stderr.write('-o (--output) <file> : output file to which the report will be written\n')
            sys.stderr.write('\n')
            exit(1)

        sam_file = sys.argv[2]

        pparser = paramsparser.Parser(paramdefs)
        paramdict = pparser.parseCmdArgs(sys.argv[3:])
        paramdict['command'] = ' '.join(sys.argv)

        eval_maplength(sam_file, paramdict)

    else:
        print 'Invalid mode!'
