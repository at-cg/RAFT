# QV for primary assembly

EXE=tools/yak/yak
READS=data/input/reads.fasta
ASM=asm/output.bp.p_ctg.fa
PREFIX=output.yak

rm -rf ${PREFIX}*

/usr/bin/time $EXE count -b37 -t 128 -o ${PREFIX} $READS
/usr/bin/time $EXE qv -t 128 ${PREFIX} $ASM