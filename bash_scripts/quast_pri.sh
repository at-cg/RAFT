# Using quast to get assembly length, NG50, etc for primary assembly

EXE=tools/quast/quast.py
CONTIGS=asm/output.bp.p_ctg.fa
PREFIX=output

rm -r ${PREFIX}

/usr/bin/time python $EXE -t 256 -o ${PREFIX} --est-ref-size 3099922541 --large $CONTIGS