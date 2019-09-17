#! /usr/bin/python

# This file contains mathods for manipulating fasta and annotation files
# To make them appropriate for testing RNAseq mapping

# All headers will be changed to chr[ID]

import sys, os

# To enable importing from samscripts submodule
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
# import utility_sam
import Annotation_formats

from fastqparser import read_fastq

NEW_ANNOTATION_MIN = 3      # A Minimum number of reads for new annotation to be proposed


def loadNewAnnotations(newAnnotationsFile):
    new_annotations = []
    with open(newAnnotationsFile, 'rU') as newann_file:
        while (True):
            line = newann_file.readline()
            if line == '':      # End of file reached (empty lines will have \n)
                break

            if line.startswith('Name:'):
                # Loading a new annotations
                # import pdb
                # pdb.set_trace()
                new_annotation = Annotation_formats.GeneDescription()
                new_annotation.seqname = line[6:-1]
                line = newann_file.readline()
                new_annotation.source = line[9:-1]
                line = newann_file.readline()
                new_annotation.strand = line[8:-1]
                line = newann_file.readline()
                noreads = int(line[16:-1])
                new_annotation.items = []

                # Reading annotation items
                while (True):
                    line = newann_file.readline()
                    if line == '':
                        sys.stderr.write('Error reading new annotations file %s!' % newAnnotationsFile);
                        import pdb
                        pdb.set_trace()
                        exit()
                    if line.startswith('Items:'):
                        # Extracting annotation items
                        line = line[7:-1]       # Removing starting text and \n at the end
                        line = line[1:-1]       # Removing starting and ending bracket for easier splitting
                        elements = line.split('] [')     # Splitting into separate items
                        for element in elements:
                            pos = element.find(',')
                            start = int(element[:pos])
                            end = int(element[pos+2:])
                            new_item = Annotation_formats.GeneItem()
                            new_item.start = start
                            new_item.end = end
                            new_annotation.items.append(new_item)
                        break

                if noreads >= NEW_ANNOTATION_MIN:
                    new_annotations.append(new_annotation)

    return new_annotations


def compare(filelist):
    annotation_dict = {}
    for filename in filelist:
        tmp_annotations = loadNewAnnotations(filename)
        annotation_dict[filename] = tmp_annotations
    
    sys.stdout.write('\nCOMPARING ANNOTATIONS:')
    for filename, annlist in annotation_dict.items():
        sys.stdout.write('\nFile %s contains %d annotations supported by at least %d reads!' % (filename, len(annlist), NEW_ANNOTATION_MIN))

    num_commonannotations = 0
    filename1 = filelist[0]
    annotations1 = annotation_dict[filename1]
    for annotation1 in annotations1:
        commonannotation = True
        for otherfile in filelist[1:]:
            otherannotations = annotation_dict[otherfile]
            foundinfile = False
            for oannotation in otherannotations:
                if annotation1.itemsEqual(oannotation):
                    foundinfile = True
                    break
            if not foundinfile:
                commonannotation = False
                break

        if commonannotation:
            num_commonannotations += 1

    sys.stdout.write('\n\nFound %d common annotations!\n' % num_commonannotations)


def analyze(annotations_file):

    filename, file_extension = os.path.splitext(annotations_file)

    if file_extension.lower() in ['.gtf', '.gff']:
        filetype = 'GTF'
    elif file_extension.lower() in ['.bed']:
        filetype = 'BED'
    else:
        raise Exception('Invalid annotation file type: %s' % file_extension)

    # Reading annotation file
    # annotations = Annotation_formats.Load_Annotation_From_File(annotations_file, check_duplicates = True)
    annotations = Annotation_formats.Load_Annotation_From_File(annotations_file)

    # for annotation in annotations:
    #     if len(annotation.items) > 1 and annotation.genename[0] == 'Q':
    #         import pdb
    #         pdb.set_trace()

    # Analyzing annotations to discover alternate splicings
    # Grouping annotations which overlap and are on the same strand
    annotation_groups = {}
    group_found = True
    gene_start = gene_end = trcnt = iden = 0

    for annotation in annotations:
        group_found = False
        for idgroup, group in annotation_groups.iteritems():
            gene_start = group[0]
            gene_end = group[1]
            trcnt = group[2]
            iden = idgroup
            if annotation.overlapsGene(gene_start, gene_end):
                group_found = True
                break

        if group_found:
            if annotation.start < gene_start:
                gene_start = annotation.start
            if annotation.end > gene_end:
                gene_end = annotation.end
            trcnt += 1
            annotation_groups[iden] = (gene_start, gene_end, trcnt)
        else:
            iden = annotation.start
            annotation_groups[iden] = (annotation.start, annotation.end, 1)


    new = True
    new_groups = {}
    groupid = group_start = group_end = 0
    for iden, group in sorted(annotation_groups.iteritems(), key=lambda(k,v):v[0]):
        # group = annotation_groups[iden]

        if new:
            groupid = group[0]
            group_start = group[0]
            group_end = group[1]
            trcnt = group[2]
            new = False
        else:
            # If overlaps with the current group join it
            if  not (group[0] > group_end or group[1] < group_start):
                if group[0] < group_start:
                    group_start = group[0]
                if group[1] > group_end:
                    group_end = group[1]
                trcnt += group[2]

            # And if it doesnt overlap, add old group to the new group dictionary
            # And start a new group
            else:
                new_groups[groupid] = (group_start, group_end, trcnt)
                groupid = group[0]
                group_start = group[0]
                group_end = group[1]
                trcnt = group[2]

    # Add last group to thenew  dictionary
    new_groups[groupid] = (group_start, group_end, trcnt)


    sys.stderr.write("\nWritting annotation groups (%d)\n" % len(new_groups))
    sys.stdout.write("ID\tSTART\tEND\tTRCNT\n")
    for idgroup in sorted(new_groups.iterkeys()):
        group = new_groups[idgroup]
        sys.stdout.write("%d\t%d\t%d\t%d\n" % (idgroup, group[0], group[1], group[2]))

        

def verbose_usage_and_exit():
    sys.stderr.write('This script is used to analyze gene annotations.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    sys.stderr.write('\t\tanalyze [filename] - analyze a single annotation file in GFF or BED format\n')
    sys.stderr.write('\t\tcompare [filename] ... - compare multiple files with new annotations (proprietary format)\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) < 3):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if (mode == 'analyze'):
        annotations_file = sys.argv[2]
        analyze(annotations_file)

    elif mode == 'compare':
        filelist = sys.argv[2:]
        compare(filelist)

    else:
        print 'Invalid mode: %s!' % mode
