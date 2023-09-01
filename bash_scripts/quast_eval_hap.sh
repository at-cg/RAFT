# Using quast to obtain length of assembly, NG50, etc for both haplotype phased assemblies taken together

EXE=tools/quast/quast.py
CONTIGS=asm/output.bp.hap1hap2.p_ctg.fa
PREFIX=output

rm -r ${PREFIX}

/usr/bin/time python $EXE -t 128 -o ${PREFIX} --est-ref-size 6200000000 --large $CONTIGS