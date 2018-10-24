---
feature_text: |
  ## Precision Medicine
title: Germline Variant Interpretation
categories:
    - Module-07-Clinical
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0007-03-01
---

### Introduction to ClinVar

### Filter VEP annotated VCF further for variants of potential clinical relevance

Use the `filter_vep` tool to filter to clinically interesting variants.

Note - I'm obtaining different results depending on the order of filters supplied if using separate "--filter" options. This should not be the case. Combining into a single expression seems to work though.

Note - there is a difference between the information available for filtering in VEP VCF output vs VEP tabular output. Namely some VCF fields (e.g., FILTER) are not included in the tabular output. Therefore, if you wish to use filter_vep on tabular output (for ease of reading) make sure to complete any vcf-specific filtering first (e.g., using GATK SelectVariants). Alternatively, it may be possible to proceed with filtering on the VCF and then reannotate with VEP specifying tabular output to convert. Yet another option would be to use VCF annotation tools ([vatools.org](http://vatools.org)).

Note - the filter_vep tool immediately and automatically limits to variants with a CSQ entry when filtering VCFs. Keep in this mind if you have annotated your VCF with the coding_only option which only adds CSQ entries for coding region alterations.


```

#Filter hard-filtered exome results for clinical relevance
#Filter VEP VCF
~/bin/ensembl-vep/filter_vep --format vcf -i /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.filtered.PASS.vep.vcf -o /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.filtered.PASS.vep.interesting.vcf --filter "(MAX_AF < 0.001 or not MAX_AF) and ((IMPACT is HIGH) or (IMPACT is MODERATE and (SIFT match deleterious or PolyPhen match damaging)))" --force_overwrite

#Filter tabular VEP
~/bin/ensembl-vep/filter_vep --format tab -i /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.filtered.PASS.vep.tsv -o /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.filtered.PASS.vep.interesting.tsv --filter "(MAX_AF < 0.001 or not MAX_AF) and ((IMPACT is HIGH) or (IMPACT is MODERATE and (SIFT match deleterious or PolyPhen match damaging)))" --force_overwrite

#Filter hard-filtered WGS results for clinical relevance
#Filter VEP VCF
~/bin/ensembl-vep/filter_vep --format vcf -i /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.filtered.PASS.vep.vcf -o /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.filtered.PASS.vep.interesting.vcf --filter "(MAX_AF < 0.001 or not MAX_AF) and ((IMPACT is HIGH) or (IMPACT is MODERATE and (SIFT match deleterious or PolyPhen match damaging)))" --force_overwrite

#Filter tabular VEP
~/bin/ensembl-vep/filter_vep --format tab -i /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.filtered.PASS.vep.tsv -o /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.filtered.PASS.vep.interesting.tsv --filter "(MAX_AF < 0.001 or not MAX_AF) and ((IMPACT is HIGH) or (IMPACT is MODERATE and (SIFT match deleterious or PolyPhen match damaging)))" --force_overwrite


#Filter VQSR-filtered exome results for clinical relevance
#Filter VEP VCF
~/bin/ensembl-vep/filter_vep --format vcf -i /home/ubuntu/data/germline_variants/Exome_Norm_GGVCFs_jointcalls_recalibrated.PASS.vep.vcf -o /home/ubuntu/data/germline_variants/Exome_Norm_GGVCFs_jointcalls_recalibrated.PASS.vep.interesting.vcf --filter "(MAX_AF < 0.001 or not MAX_AF) and ((IMPACT is HIGH) or (IMPACT is MODERATE and (SIFT match deleterious or PolyPhen match damaging)))" --force_overwrite

#Filter tabular VEP
~/bin/ensembl-vep/filter_vep --format tab -i /home/ubuntu/data/germline_variants/Exome_Norm_GGVCFs_jointcalls_recalibrated.PASS.vep.tsv -o /home/ubuntu/data/germline_variants/Exome_Norm_GGVCFs_jointcalls_recalibrated.PASS.vep.interesting.tsv --filter "(MAX_AF < 0.001 or not MAX_AF) and ((IMPACT is HIGH) or (IMPACT is MODERATE and (SIFT match deleterious or PolyPhen match damaging)))" --force_overwrite

#Filter VQSR-filtered WGS results for clinical relevance
#Filter VEP VCF
~/bin/ensembl-vep/filter_vep --format vcf -i /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated.PASS.vep.vcf -o /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated.PASS.vep.interesting.vcf --filter "(MAX_AF < 0.001 or not MAX_AF) and ((IMPACT is HIGH) or (IMPACT is MODERATE and (SIFT match deleterious or PolyPhen match damaging)))" --force_overwrite

#Filter tabular VEP
~/bin/ensembl-vep/filter_vep --format tab -i /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated.PASS.vep.tsv -o /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated.PASS.vep.interesting.tsv --filter "(MAX_AF < 0.001 or not MAX_AF) and ((IMPACT is HIGH) or (IMPACT is MODERATE and (SIFT match deleterious or PolyPhen match damaging)))" --force_overwrite

```

### Explore the filtered results

The above filtering has limited the results to approximately 200 germline variants of potential interest.

- Are there any variants with known clinical significance (See CLIN_SIG column)? Use their HGVS IDs (See HGVSc or HGVSp columns) to search the [ClinGen Allele Registry](http://reg.clinicalgenome.org), and then follow the links (if any) to ClinVar or other resources of interest.
- Are there any variants in known cancer genes? Try intersecting the remaining genes with variants with the [Cancer Gene Census](https://cancer.sanger.ac.uk/census) genes?
- Manually review any variants of interest in IGV

