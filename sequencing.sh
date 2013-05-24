tophat --version
# Tophat v2.0.0

# filone
for LANE in 0 0_5 2 6 18; do
    tophat --library-type fr-unstranded \
           --solexa1.3-quals \
           --segment-length 18 \
           --no-coverage-search \
           --no-novel-juncs \
           -G "gencode.v14.annotation.gtf" \
           -o "$OUT_DIR" \
           hg19 "${LANE}hpi.fastq"
done

# yang_a

for LANE in a_0 a_0_5 a_1 a_2 a_4 ; do
    tophat --library-type fr-secondstrand \
           --segment-length 25 \
           --color \
           --bowtie1 \
           --quals \
           --no-coverage-search \
           --no-novel-juncs \
           -G "gencode.v14.annotation.gtf" \
           -o "$OUT_DIR" \
           hg19 "${LANE}hpi.csfasta" "${LANE}hpi_QV.qual"
done



