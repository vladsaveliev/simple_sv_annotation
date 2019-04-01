simple_sv_annotation
====================

A tool for simplifying SnpEff annotations, and prioritizing fusions, exon losses, and disruption of tumor suppressor genes.

## Requirements

1. python3
2. Python modules (available on conda): cyvcf2, click, [ngs_utils](https://github.com/vladsaveliev/NGS_Utils/) for gene lists
3. a [VCF](https://vcftools.github.io/specs.html) file annotated with [snpEff v4.3+](http://snpeff.sourceforge.net/)

The input VCF should work with the [ANN](http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf) annotation field rather than the older EFF field.

No specific VCF caller requirements, however VCF must have the following fields:

1. SVTYPE
2. MATEID (for SVTYPE=BND)
3. END (for whole exon deletions)

## Usage

```
simple_sv_annotation input.vcf > output.vcf
```

## Licence

This program is distributed under the MIT licence save for the exception below.

## Example Output

Primary output for ```simple_sv_annotation```:

#### 1. Add SIMPLE_ANN field

```simple_sv_annotation.py``` will not alter the ANN field provided by SnpEff. Instead an additional field called SIMPLE_ANN will be added to the SV call. 
A SIMPLE_ANN will only be added to variants that can be simplified, other variants are not altered.

There are six fields in the SIMPLE_ANN tag separated by "`|`".

1. SV type (deletion, duplication, insertion, breakend)
2. Effect (fusion, exon loss, intergenic, intronic)
3. Gene name
4. Feature ID (often transcript name)
5. Prioritization detail. For exon loss variants, deleted exon numbers (Exon5del).
6. Priority of the event (1 highest, 4 lowest)

example:

```
before:

chr17  41258467  del_5  ATATACCTTTTGGTTATATCATTCTTACATAAAGGACACTGTGAAGGCCCTTTCTTCTGGTTGAGAAGTTTCAGCATGCAAAATCTATA  A  .  .  END=41258555;SVTYPE=DEL;SVLEN=-88;UPSTREAM_PAIR_COUNT=0;DOWNSTREAM_PAIR_COUNT=0;PAIR_COUNT=0;ANN=A|exon_loss_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&splice_region_variant&splice_region_variant&splice_region_variant&intron_variant&intron_variant|HIGH|BRCA1|BRCA1|transcript|NM_007294.3|Coding|4/23|c.135-5_212+5delTATAGATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGTATA||||||

after:

chr17  41258467  del_5  ATATACCTTTTGGTTATATCATTCTTACATAAAGGACACTGTGAAGGCCCTTTCTTCTGGTTGAGAAGTTTCAGCATGCAAAATCTATA  A  .  .  END=41258555;SVTYPE=DEL;SVLEN=-88;UPSTREAM_PAIR_COUNT=0;DOWNSTREAM_PAIR_COUNT=0;PAIR_COUNT=0;ANN=A|exon_loss_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&splice_region_variant&splice_region_variant&splice_region_variant&intron_variant&intron_variant|HIGH|BRCA1|BRCA1|transcript|NM_007294.3|Coding|4/23|c.135-5_212+5delTATAGATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGTATA||||||;SIMPLE_ANN=DEL|EXON_DEL|BRCA1|NM_007294.3|Exon5del
```

[![DOI](https://zenodo.org/badge/40991130.svg)](https://zenodo.org/badge/latestdoi/40991130)
