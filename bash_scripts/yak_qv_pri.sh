EXE=$SCRATCH/tools/yak/yak
READS=../../reads.fasta
ASM=../output.bp.p_ctg.fa
PREFIX=output.yak

rm -rf ${PREFIX}*

/usr/bin/time $EXE count -b37 -t 128 -o ${PREFIX} $READS
/usr/bin/time $EXE qv -t 128 ${PREFIX} $ASM