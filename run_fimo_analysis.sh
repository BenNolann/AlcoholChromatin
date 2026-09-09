#!/bin/bash

INPUT_XLS=$1
OUT_PREFIX=$2
MOTIF=$3
GENOME_FA="/Zulu/bnolan/Annotations/human/human/fa/ucsc_PI_noscaffolds/hg38_noaltchromosomes.fa"
#GENOME_FA="/Zulu/bnolan/Annotations/mouse/mouse/mm10.fa"
echo "--- Processing 300bp windows around summits ---"
sed '1d' "$INPUT_XLS" | awk -v OFS="\t" '{print $1,$4-150,$4+150,$0}' | cut -f 1-3,8- > "${OUT_PREFIX}_300bp.summit.tsv"
#sed '1d' "$INPUT_XLS" | awk -v OFS="\t" '{print $1,$2,$3}' > "${OUT_PREFIX}_300bp.summit.tsv"

echo "--- Extracting Fasta Sequences ---"
bedtools getfasta -fi "$GENOME_FA" -bed "${OUT_PREFIX}_300bp.summit.tsv" -fo "${OUT_PREFIX}.fa"

echo "--- Running FIMO ---"
fimo --parse-genomic-coord --oc "${OUT_PREFIX}_fimo_out" "${MOTIF}" "${OUT_PREFIX}.fa"

#echo "--- Filtering for Motif $MOTIF ---"
cut -f 3- "${OUT_PREFIX}_fimo_out/fimo.tsv" | sed '1d' > "${OUT_PREFIX}_fimo.tsv"

echo "--- Motif Count Distribution ---"
# Generate the distribution output you requested
intersectBed -a "${OUT_PREFIX}_300bp.summit.tsv" -b "${OUT_PREFIX}_fimo.tsv" -c > tmp_counts.tsv
cut -f 10 tmp_counts.tsv | sort | uniq -c

echo "--- Filtering for Peaks with exactly 1 Motif ---"
awk '{if($10==1) print $0}' tmp_counts.tsv > "${OUT_PREFIX}_300bp.summit.1motif.tsv"

echo "--- Final Intersection ---"
intersectBed -a "${OUT_PREFIX}_300bp.summit.1motif.tsv" -b "${OUT_PREFIX}_fimo.tsv" -wb > "${OUT_PREFIX}_300bp.summit.1motif.fimo.bed"

# Cleanup
rm tmp_counts.tsv
echo "Done. Results saved with prefix: $OUT_PREFIX"