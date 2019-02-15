#!/usr/bin/env python

"""
Name:    simple_sv_annotation.py

Purpose: Simplify SV annotations from snpEff to highlight exon impact.

Input:   vcf file with SVs annotated with snpEff 4.3 or higher, text file with known fusion pairs, text file with genes of interest

Output:  vcf file with SV annotation field with simplified annotation

Usage:   simple_sv_annotation.py [OPTIONS]

Authors:  David Jenkins (david.jenkins1@astrazeneca.com/dfj@bu.edu), Miika Ahdesmaki (miika.ahdesmaki @ astrazeneca.com / live.fi)

Notes:
    - Simplified with SnpEff 4.3 updates that add better annotation for fusion events resulting from all BND/INV/DUP/DEL
    - No simplified annotation for within-exon deletions (frameshifts, in-frame dels) as not classified as SVs
    - Soft filters intergenic events that are not involved with genes directly

Current scheme with priority 1(high)-3(low)
- exon loss
   - on prioritisation gene list (2)
   - other (3)
- gene_fusion
   - paired (hits two genes)
      - on list of known pairs (1)
      - one gene is a known promiscuous fusion gene (1)
      - other:
         - one or two genes on prioritisation gene list (2)
         - neither gene on prioritisation gene list (3)
   - unpaired (hits one gene)
       - on prioritisation gene list (2)
       - others (3)
- upstream or downstream
   - on prioritisation gene list genes (2)  - e.g. one gene is got into control of another gene's promoter and get overexpressed (oncogene) or underexpressed (tsgene)
- LOF
   - in a tumor suppressor (3)
- other (4)
- missing ANN or SVTYPE: "MissingAnn" in FILTER (4)

Populates:
 - INFO/SIMPLE_ANN that looks like: SIMPLE_ANN=INV|GENE_FUSION|ALK&EML4|NM_...&NM_...|KNOWN_FUSION|1|MODERATE
 - INFO/SV_HIGHEST_TIER (1..4)

See the provided README.md file for more information
"""

from __future__ import print_function
import argparse, os, vcf, sys
from collections import defaultdict
try:
    import pysam
except ImportError:
    pysam = None

def main(vcf_in, outfile, exon_nums, args):
    """
    Adds additional header information, opens the infile for reading, opens
    the outfile for writing, and iterates through the records to identify
    records that are candidates for adding the simple annotation
    """
    vcf_reader = vcf.Reader(filename=vcf_in) if vcf_in != "-" else vcf.Reader(sys.stdin)

    # add a new info field into the header of the output file:
    vcf_reader.infos['SIMPLE_ANN'] = vcf.parser._Info(id="SIMPLE_ANN", num=".", type="String",
          desc="Simplified human readable structural variant annotation: 'SVTYPE | ANNOTATION | GENE(s) | TRANSCRIPT | "
               "DETAIL (exon losses, KNOWN_FUSION, ON_PRIORITY_LIST, NOT_PRIORITISED) | PRIORITY (1-3) '",
          source=None, version=None)

    # Add an info field for highest priority (given multiple annotations per entry):
    vcf_reader.infos['SV_HIGHEST_TIER'] = vcf.parser._Info(id="SV_HIGHEST_TIER", num=1, type="Integer",
           desc="Highest priority tier for the effects of a variant entry", source=None, version=None)

    # Add filters headers:
    # vcf_reader.filters["MissingAnn"] = vcf.parser._Filter(id="MissingAnn", desc="Rejected in SV-prioritize (missing ANN/BND)")
    # vcf_reader.filters["Intergenic"] = vcf.parser._Filter(id="Intergenic", desc="Rejected in SV-prioritize (purely intergenic, or a small event)")
    vcf_writer = vcf.Writer(open(outfile, 'w'), vcf_reader) if outfile != "-" else vcf.Writer(sys.stdout, vcf_reader)

    # Read in gene lists
    known_fusions, tier2_fusion_pairs, known_promiscuous, prioritised_genes = read_gene_lists(
        args.known_fusion_pairs, args.tier2_fusion_pairs, args.known_fusion_promiscuous, args.gene_list)

    for record in vcf_reader:
        if record.FILTER is None:
            record.FILTER = []

        if 'SVTYPE' in record.INFO and 'ANN' in record.INFO:
            vcf_writer.write_record(simplify_ann(record, exon_nums, known_fusions, tier2_fusion_pairs, known_promiscuous, prioritised_genes))
        else: 
            # record.FILTER.append("MissingAnn")
            vcf_writer.write_record(record)
    vcf_writer.close()

def read_gene_lists(known_fusion_pairs, tier2_fusion_pairs, known_fusion_promiscuous, gene_list):
    known_pairs = []
    tier2_pairs = []
    known_promiscuous = []
    gl = []
    if known_fusion_pairs and os.path.isfile(known_fusion_pairs):
        with open(known_fusion_pairs, 'r') as myfhandle:
            for line in myfhandle:
                genes = line.strip().split(",")
                if len(genes) == 2 and len(genes[0]) > 0 and len(genes[1]) > 0:
                    known_pairs.append([genes[0].strip(), genes[1].strip()])
    if tier2_fusion_pairs and os.path.isfile(tier2_fusion_pairs):
        with open(tier2_fusion_pairs, 'r') as myfhandle:
            for line in myfhandle:
                genes = line.strip().split(",")
                if len(genes) == 2 and len(genes[0]) > 0 and len(genes[1]) > 0:
                    tier2_pairs.append([genes[0].strip(), genes[1].strip()])
    if known_fusion_promiscuous and os.path.isfile(known_fusion_promiscuous):
        with open(known_fusion_promiscuous, 'r') as myfhandle:
            for line in myfhandle:
                gene = line.strip()
                known_promiscuous.append(gene)

    if gene_list and os.path.isfile(gene_list):
        with open(gene_list, 'r') as myghandle:
            for line in myghandle:
                gene = line.strip()
                if len(gene) > 0:
                    gl.append(gene)
    return known_pairs, tier2_pairs, known_promiscuous, gl

def simplify_ann(record, exon_nums, known_fusions, tier2_fusion_pairs, known_promiscuous, prioritised_genes):
    """
    Find any annotations that can be simplified and call the method
    to annotate it.
    """
    # marching order is: 'exon_loss_variant', fusions, others
    # to-do: CNV and INS?

    exon_losses_by_tid = defaultdict(list)
    sv_best_tier = 4
    simple_annos = set()
    svtype = record.INFO['SVTYPE']

    lofs = record.INFO.get('LOF', [])  # (GRK7|ENSG00000114124|1|1.00),(DRD3|ENSG00000151577|4|1.00),...
    lof_by_gene = {l.strip('(').strip(')').split('|')[0]: l for l in lofs}
    simple_lofs = set()

    for anno in record.INFO.get('ANN', []):
        anno_fields = anno.split('|')
        # T|splice_acceptor_variant&splice_region_variant&intron_variant|HIGH|NCOA1|ENSG00000084676|transcript|ENST00000348332|
        #   protein_coding|4/20|c.257-52_257-2delTGGAAATAAGCTCTTTTCAGATATGTGATTTTTTTAAGTTTCTTTATTATA||||||INFO_REALIGN_3_PRIME,
        allele, effect, impact, gene, geneid, feature, featureid, biotype, rank, c_change, p_change = anno_fields[:11]
        effects = effect.split('&')
        genes = gene.split('&')
        on_priority_list = any(g in prioritised_genes for g in genes)

        if "exon_loss_variant" in effects:
            # collecting exon losses to process them together further for the variant
            exon_losses_by_tid[featureid].append(anno_fields)

        else:
            # deside of to report the annotation straigh away
            var_detail = ''
            var_priority = 4
            is_fusion = any(e in effects for e in ["gene_fusion", "bidirectional_gene_fusion"])
            is_downstream_upstream = any(e in effects for e in ["downstream_gene_variant", "upstream_gene_variant"])
            if is_fusion:
                # This could be 'gene_fusion', 'bidirectional_gene_fusion' but not 'feature_fusion'
                # 'gene_fusion' could lead to a coding fusion whereas
                # 'bidirectional_gene_fusion' is likely non-coding (opposing frames, _if_ inference correct)
                # result_annos, tier = annotate_fusion(record, genes, effect, featureid, known_fusions, known_promiscuous,
                #                                      prioritised_genes, tier)
                var_priority = 3
                var_detail = "NOT_PRIORITISED"

                # on list of known fusions
                if len(genes) > 1:
                    g1, g2 = genes[:2]
                    if (g1, g2) in known_fusions or (g2, g1) in known_fusions or g1 in known_promiscuous or g2 in known_promiscuous:
                        var_priority = 1 if "bidirectional_gene_fusion" not in effects else 2
                        var_detail = "KNOWN_FUSION"

                    elif (g1, g2) in tier2_fusion_pairs or (g2, g1) in tier2_fusion_pairs:
                        var_priority = 2 if "bidirectional_gene_fusion" not in effects else 3
                        var_detail = "KNOWN_FUSION"

                # one of the genes is of interest
                elif on_priority_list:
                    var_priority = 2 if "bidirectional_gene_fusion" not in effects else 3
                    var_detail = "ON_PRIORITY_LIST"

            elif on_priority_list:
                if gene in lof_by_gene:
                    var_priority = 3
                    var_detail = "ON_PRIORITY_LIST&LOF"
                    simple_lofs.add(gene)

                if is_downstream_upstream:
                    # get SVs affecting prioritised genes or regions up/downstream of them
                    var_priority = 2
                    var_detail = "ON_PRIORITY_LIST"

                # elif impact == "HIGH":
                #     var_priority = 3
                #     var_detail = "ON_PRIORITY_LIST"
                #     if feature == 'interaction':
                #         featureid = featureid.split(':')[-1]

            if var_priority < 4:
                simple_annos.add((svtype, effect, gene, featureid, var_detail, var_priority))
                sv_best_tier = min(var_priority, sv_best_tier)

    if len(exon_losses_by_tid) > 0:
        losses = annotate_exon_loss(record.POS, record.INFO['END'], exon_losses_by_tid, exon_nums, prioritised_genes)
        for (gene, transcriptid, deleted_exons, var_priority) in losses:
            simple_annos.add(('DEL', 'EXON_DEL', gene, transcriptid, deleted_exons, var_priority))
            sv_best_tier = min(var_priority, sv_best_tier)

    if simple_annos:
        simple_annos = sorted(simple_annos)
        record.INFO['SIMPLE_ANN'] = ['|'.join(map(str, a)) for a in simple_annos]

    record.INFO['SV_HIGHEST_TIER'] = sv_best_tier
    if simple_lofs:
        record.INFO['LOF'] = ','.join(sorted(simple_lofs))
    elif 'LOF' in record.INFO:
        del record.INFO['LOF']

    return record

def uniq_list(inlist):
    """Remove unique elements from a list"""
    inset = set(inlist)
    return list(inset)

def find_deleted_exons(annotations):
    """Take the annotations for a particular transcript
    and parse them to find the numbers of the exons that have
    been deleted
    """
    exons = []
    gene = ''
    for anno_fields in annotations:
        allele, effect, impact, g, geneid, feature, featureid, biotype, rank, c_change, p_change = anno_fields[:11]
        gene = gene or g
        try:
            exons.append(int(rank.split('/')[0]))
        except ValueError:
            pass
    return exons, gene

def find_alt_deleted_exons(start, end, exon_nums, annotations):
    """In the case where the user has provided a file of alternate
    transcript numbers, use those numbers to find the exons that have been
    deleted
    """
    exons = []
    gene = ''
    for ann_fields in annotations:
        allele, effect, impact, g, geneid, feature, featureid, biotype, rank, c_change, p_change = ann_fields[:11]
        if gene == '':
            gene = g
        else:
            break
    for i in exon_nums:
        if int(i[0]) > start and int(i[1]) <= end:
            exons.append(int(i[2]))
    return exons, gene

def annotate_exon_loss(start, end, exon_loss_anno_by_tid, exon_nums, prioritised_genes):
    """Create the exon loss simple annotation from the exon dict created
    in simplify_ann

    For each transcript with exon losses, find the numbers for each exon
    and create the annotation
    Example: DEL|EXON_DEL|BLM|NM_001287247.1|Exon2-12del
    """
    annos = set()
    for transcript, annotations in exon_loss_anno_by_tid.items():
        # Remove version number if it exists
        if transcript.split('.')[0] in exon_nums:
            # Use alternate exon numbers for transcript
            exons, gene = find_alt_deleted_exons(start, end, exon_nums[transcript.split('.')[0]], annotations)
        else:
            # se snpEff numbers for transcript
            exons, gene = find_deleted_exons(annotations)
        exons = uniq_list(exons)
        if len(exons) == 0:
            return None
        if max(exons)-min(exons)+1 == len(exons):
            if len(exons) == 1:
                deleted_exons = "Exon"+str(exons[0])+"del"
            else:
                deleted_exons = "Exon"+str(min(exons))+"-"+str(max(exons))+"del"
        else:
            deleted_exons = "Exon"+str(min(exons))+"-"+str(max(exons))+"del"
        var_priority = 2 if gene in prioritised_genes else 3

        annos.add((gene, transcript, deleted_exons, var_priority))

    return annos

# def annotate_fusion(record, genes, effect, featureid, known_fusions, known_promiscuous, prioritised_genes, tier):

# UNUSED FUNCTION FOR NOW
#def annotate_intergenic_var(record, ann_a):
#    """Create a simplified version of the annotation field for an intergenic var
#
#    Regardless of sv type, the simple annotation for an intergenic variant
#    looks like: SIMPLE_ANN=INV|INTERGENIC|LETM2-FGFR1||
#    """
#    simple_ann = "%s|INTERGENIC|%s||" % (record.INFO['SVTYPE'],ann_a[4])
#    try:
#        if simple_ann not in record.INFO['SIMPLE_ANN']: # avoid duplicate entries that bloat the output
#            record.INFO['SIMPLE_ANN'].append(simple_ann)
#    except KeyError:
#        record.INFO['SIMPLE_ANN'] = [simple_ann]

def create_exon_numDict(infile):
    """Create a dictionary of exon numbers based on an input file of alternate
    chromosome numberings.

    example: BRCA1 has no exon 4 (http://www.medscape.com/viewarticle/567639_2),
    making exon numbers from snpEff incorrect.
    """
    exons = {}
    f = open(infile)
    for line in f:
        la = line.rstrip("\n").split("\t")
        name = la[3].split("|")
        transcript = name[0]
        if transcript in exons:
            exons[transcript].append((la[1],la[2],name[1]))
        else:
            exons[transcript] = [(la[1],la[2],name[1])]
    return exons

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Simplify SV annotations from snpEff to highlight exon impact. Requires the pyvcf module.")
    parser.add_argument('vcf', help='VCF file with snpEff annotations')
    parser.add_argument('--gene_list', '-g', help='File with names of genes (one per line) for prioritisation', required=False, default=None)
    parser.add_argument('--known_fusion_pairs', '-k', help='File with known fusion gene pairs, one pair per line delimited by comma', required=False, default=None)
    parser.add_argument('--tier2_fusion_pairs', '--k2', help='File with purative known fusion gene pairs, one pair per line delimited by comma', required=False, default=None)
    parser.add_argument('--known_fusion_promiscuous', '-p', help='File with known promiscuous fusion genes, one gene name per line', required=False, default=None)
    parser.add_argument('--output', '-o', help='Output file name (must not exist). Does not support bgzipped output. Use "-" for stdout. [<invcf>.simpleann.vcf]', required=False)
    parser.add_argument('--exonNums', '-e', help='List of custom exon numbers. A transcript listed in this file will be '
                                                 'annotated with the numbers found in this file, not the numbers found in the snpEff result')
    #parser_excl = parser.add_mutually_exclusive_group(required=False)
    args = parser.parse_args()
    if args.output:
        outfile = args.output
    else:
        outfile = sys.stdout
    exonNumDict = {}
    if args.exonNums:
        exonNumDict = create_exon_numDict(args.exonNums) 
    main(args.vcf, outfile, exonNumDict, args)
