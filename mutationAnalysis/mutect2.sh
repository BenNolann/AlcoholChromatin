#!/bin/bash
set -euo pipefail

REF=hg38_noaltchromosomes.fa
UT_BAM=UT_merged.rg.bam
E6_BAM=E6_merged.rg.bam
PEAKS=E6vUT.tsv
BLACKLIST=hg38-blacklist.v2.bed
GNOMAD=af-only-gnomad.hg38.vcf.gz

UT_SAMPLE=UT
E6_SAMPLE=E6

# --- 0. Sort peaks, convert to BED ---
sort -k1,1 -k2,2n "$PEAKS" > E6vUT.sorted.tsv
awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' E6vUT.sorted.tsv > E6vUT.sorted.bed

# --- 1. Call somatic variants, E6 as tumor, UT as matched normal ---
gatk Mutect2 \
  -R "$REF" \
  -I "$E6_BAM" -tumor "$E6_SAMPLE" \
  -I "$UT_BAM" -normal "$UT_SAMPLE" \
  -L E6vUT.sorted.bed \
  --germline-resource "$GNOMAD" \
  -O E6_vs_UT_somatic_unfiltered.vcf.gz \
  --f1r2-tar-gz f1r2.tar.gz

# --- 2. Learn read-orientation model ---
gatk LearnReadOrientationModel \
  -I f1r2.tar.gz \
  -O read-orientation-model.tar.gz

# --- 3. Apply statistical filtering ---
gatk FilterMutectCalls \
  -R "$REF" \
  -V E6_vs_UT_somatic_unfiltered.vcf.gz \
  --ob-priors read-orientation-model.tar.gz \
  -O E6_vs_UT_somatic.filtered.vcf.gz

# --- 4. Keep only PASS calls ---
bcftools view -f PASS E6_vs_UT_somatic.filtered.vcf.gz -Oz -o E6_vs_UT_somatic.pass.vcf.gz
bcftools index E6_vs_UT_somatic.pass.vcf.gz

# --- 5. Biallelic SNPs only ---
bcftools view -v snps -m2 -M2 E6_vs_UT_somatic.pass.vcf.gz -Oz -o E6_vs_UT_somatic.snps.vcf.gz
bcftools index E6_vs_UT_somatic.snps.vcf.gz

# --- 6. Blacklist subtraction ---
bedtools intersect -v -a E6_vs_UT_somatic.snps.vcf.gz -b "$BLACKLIST" -header \
  > E6_vs_UT_somatic.snps.noblacklist.vcf
bgzip -f E6_vs_UT_somatic.snps.noblacklist.vcf
bcftools index E6_vs_UT_somatic.snps.noblacklist.vcf.gz

# --- 7. Extract candidate table ---
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%DP\t%AF\t]\n' \
  E6_vs_UT_somatic.snps.noblacklist.vcf.gz \
  > somatic_candidates_E6vUT.tsv

echo "Candidates surviving full pipeline:"
wc -l somatic_candidates_E6vUT.tsv


awk -F'\t' 'BEGIN{OFS="\t"}
{
  e6_af = $8;  ut_af = $12;
  diff = e6_af - ut_af;
  abs_diff = (diff<0 ? -diff : diff);
  e6_dp = $7;  ut_dp = $11;

  if (abs_diff >= 0.10 && e6_dp >= 15 && ut_dp >= 15) {
    print $1, $2, $3, $4, $5, $6, $7, e6_af, $9, $10, $11, ut_af, diff
  }
}' somatic_candidates_E6vUT.tsv > af_diff_candidates.tsv

wc -l af_diff_candidates.tsv