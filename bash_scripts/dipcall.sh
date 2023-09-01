# Variant calling using hap.py

TRUTH=data/vcf/human/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
QUERY=dipcall_run/output.dip.vcf.gz
REF=data/genomes/human/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
TRUTH_BED=dip2.bed # Restriction to confident call regions of GIAB and dipcall

shifter --image=pkrusche/hap.py:latest /opt/hap.py/bin/hap.py $TRUTH $QUERY -f $TRUTH_BED -r $REF -o output --engine=vcfeval --pass-only