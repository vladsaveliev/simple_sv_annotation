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
    vcf_reader.infos['SV_TOP_TIER'] = vcf.parser._Info(id="SV_TOP_TIER", num=1, type="Integer",
           desc="Highest priority tier for the effects of a variant entry", source=None, version=None)

    # Add filters headers:
    # vcf_reader.filters["MissingAnn"] = vcf.parser._Filter(id="MissingAnn", desc="Rejected in SV-prioritize (missing ANN/BND)")
    # vcf_reader.filters["Intergenic"] = vcf.parser._Filter(id="Intergenic", desc="Rejected in SV-prioritize (purely intergenic, or a small event)")
    vcf_writer = vcf.Writer(open(outfile, 'w'), vcf_reader) if outfile != "-" else vcf.Writer(sys.stdout, vcf_reader)

    # Read in gene lists
    known_pairs, tier2_pairs, known_promiscuous, prio_genes, ts_genes = read_gene_lists(args)

    for record in vcf_reader:
        if record.FILTER is None:
            record.FILTER = []

        if 'SVTYPE' in record.INFO and 'ANN' in record.INFO:
            vcf_writer.write_record(simplify_ann(record, exon_nums, known_pairs, tier2_pairs,
                                                 known_promiscuous, prio_genes, ts_genes))
        else: 
            # record.FILTER.append("MissingAnn")
            vcf_writer.write_record(record)
    vcf_writer.close()

def read_gene_lists(args):
    def _read_list(path):
        out_set = set()
        if path:
            with open(path, 'r') as myfhandle:
                for line in myfhandle:
                    genes = line.strip().split(",")
                    if len(genes) == 1 and len(genes[0].strip()) > 0:
                        out_set.add(genes[0].strip())
                    if len(genes) == 2 and len(genes[0].strip()) > 0 and len(genes[1].strip()) > 0:
                        out_set.add((genes[0].strip(), genes[1].strip()))
        return out_set

    known_pairs = _read_list(args.known_fusion_pairs)
    tier2_pairs = _read_list(args.tier2_fusion_pairs)
    known_promiscuous = _read_list(args.known_fusion_promiscuous)
    gl = _read_list(args.gene_list)
    tsgenes = _read_list(args.tumor_suppressors)
    return known_pairs, tier2_pairs, known_promiscuous, gl, tsgenes

def simplify_ann(record, exon_nums, known_pairs, tier2_pairs, known_promiscuous, prio_genes, ts_genes):
    """
    Find any ANN or LOF annotations that can be prioritized
    """
    # marching order is: 'exon_loss_variant', fusions, others
    # to-do: CNV and INS?

    svtype = record.INFO['SVTYPE']

    # Annotate from ANN
    exon_losses_by_tid = defaultdict(list)
    sv_top_tier = 4
    simple_annos = set()
    for anno in record.INFO.get('ANN', []):
        anno_fields = anno.split('|')
        # T|splice_acceptor_variant&splice_region_variant&intron_variant|HIGH|NCOA1|ENSG00000084676|transcript|ENST00000348332|
        #   protein_coding|4/20|c.257-52_257-2delTGGAAATAAGCTCTTTTCAGATATGTGATTTTTTTAAGTTTCTTTATTATA||||||INFO_REALIGN_3_PRIME,
        allele, effect, impact, gene, geneid, feature, featureid, biotype, rank, c_change, p_change = anno_fields[:11]
        effects = set(effect.split('&'))
        genes = set(gene.split('&'))

        ann_tier = 4
        ann_detail = ''

        if effects & {"exon_loss_variant"}:
            # collecting exon losses to process them together further for the variant
            exon_losses_by_tid[featureid].append(anno_fields)

        elif effects & {"gene_fusion"}:
            # This could be 'gene_fusion', but not 'bidirectional_gene_fusion' or 'feature_fusion'
            # ('gene_fusion' could lead to a coding fusion whereas 'bidirectional_gene_fusion' is
            # likely non-coding (opposing frames, _if_ inference correct))
            if len(genes) == 2:
                g1, g2 = genes
                if {(g1, g2), (g2, g1)} & known_pairs:
                    ann_tier = 1
                    ann_detail = "known_pair"

                elif {g1, g2} & known_promiscuous:
                    ann_tier = 1
                    ann_detail = "known_promiscuous"

                elif {(g1, g2), (g2, g1)} & tier2_pairs:
                    ann_tier = 2
                    ann_detail = "known_pair"

            # One of the genes is of interest
            elif genes & prio_genes:
                ann_tier = 2
                ann_detail = "key_gene"

            else:
                ann_tier = 3
                ann_detail = "unknown"

        # "downstream_gene_variant" and "upstream_gene_variant" can also turn out to be a fusion
        # when gene A falls into control of a promoter of gene B
        elif effects & {"downstream_gene_variant", "upstream_gene_variant"}:
            if len(genes) == 2:
                g1, g2 = genes
                if {(g1, g2), (g2, g1)} & known_pairs:
                    ann_tier = 2
                    ann_detail = "known_pair"

                elif {g1, g2} & known_promiscuous:
                    ann_tier = 2
                    ann_detail = "known_promiscuous"

                elif {(g1, g2), (g2, g1)} & tier2_pairs:
                    ann_tier = 3
                    ann_detail = "known_pair"

            # One of the genes is of interest
            elif genes & prio_genes:
                ann_tier = 3
                ann_detail = "near_key_gene"

        # any other event
        else:
            ann_detail = "unprioritized"
            ann_tier = 4
            featureid = ''
            gene = '&'.join(genes&prio_genes)
            # gene = ''
            # if genes&prio_genes:
            #     if len(genes&prio_genes) > 2:
            #         gene = f'{len(genes&prio_genes)}_key_genes'
            #     else:
            #         gene = '&'.join(genes&prio_genes)
            # if genes-prio_genes:
            #     if gene:
            #         gene += '&'
            #     if len(genes-prio_genes) > 2:
            #         gene += f'{len(genes-prio_genes)}_other_genes'
            #     else:
            #         gene += '&'.join(genes-prio_genes)

        sv_top_tier = min(ann_tier, sv_top_tier)
        simple_annos.add((svtype, effect, gene, featureid, ann_detail, ann_tier))

    if len(exon_losses_by_tid) > 0:
        losses = annotate_exon_loss(record.POS, record.INFO['END'], exon_losses_by_tid, exon_nums, prio_genes)
        for (gene, transcriptid, deleted_exons, ann_tier) in losses:
            simple_annos.add((svtype, 'exon_loss', gene, transcriptid, deleted_exons, ann_tier))
            sv_top_tier = min(ann_tier, sv_top_tier)

    # Annotate from LOF
    lofs = record.INFO.get('LOF', [])  # (GRK7|ENSG00000114124|1|1.00),(DRD3|ENSG00000151577|4|1.00),...
    lof_genes = {l.strip('(').strip(')').split('|')[0] for l in lofs}
    if lof_genes & ts_genes:
        lof_genes = lof_genes & ts_genes
        if lof_genes & prio_genes:
            ann_tier = 2
            ann_detail = 'key_tsgene'
            lof_genes = lof_genes & prio_genes
        else:
            ann_tier = 3
            ann_detail = 'tsgene'

        simple_annos.add((svtype, 'LOF', '&'.join(lof_genes), '', ann_detail, ann_tier))
        sv_top_tier = min(ann_tier, sv_top_tier)

    if simple_annos:
        record.INFO['SIMPLE_ANN'] = ['|'.join(map(str, a)) for a in simple_annos]

    record.INFO['SV_TOP_TIER'] = sv_top_tier

    if 'ANN' in record.INFO:
        del record.INFO['ANN']
    if 'LOF' in record.INFO:
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
    Example: DEL|exon_loss|BLM|NM_001287247.1|exon2-12del
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
                deleted_exons = "exon"+str(exons[0])+"del"
            else:
                deleted_exons = "exon"+str(min(exons))+"-"+str(max(exons))+"del"
        else:
            deleted_exons = "exon"+str(min(exons))+"-"+str(max(exons))+"del"
        var_priority = 1 if gene in prioritised_genes else 2

        annos.add((gene, transcript, deleted_exons, var_priority))

    return annos

def create_exon_num_dict(infile):
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
    parser.add_argument('--tumor_suppressors', '--ts', help='File with names of tumor suppressor genes (one per line) to prioritize LoF events', required=False, default=None)
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
        exonNumDict = create_exon_num_dict(args.exonNums)
    main(args.vcf, outfile, exonNumDict, args)
