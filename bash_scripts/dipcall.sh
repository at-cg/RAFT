TRUTH=$SCRATCH/data/vcf/human/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
QUERY=../output.dip.vcf.gz
REF=$SCRATCH/data/genomes/human/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
TRUTH_BED=$SCRATCH/dip2.bed

shifter --image=pkrusche/hap.py:latest /opt/hap.py/bin/hap.py $TRUTH $QUERY -f $TRUTH_BED -r $REF -o output --engine=vcfeval --pass-only