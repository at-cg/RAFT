EXE=$SCRATCH/tools/quast/quast.py
CONTIGS=../output.bp.p_ctg.fa
PREFIX=output

rm -r ${PREFIX}

/usr/bin/time python $EXE -t 256 -o ${PREFIX} --est-ref-size 3099922541 --large $CONTIGS