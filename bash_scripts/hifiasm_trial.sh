# Basic hifiasm run

EXE=tools/hifiasm/hifiasm
READS=data/input/reads.fasta
PREFIX=output

rm -f ${PREFIX}*

/usr/bin/time $EXE -o $PREFIX -t 256 -r 1 $READS

awk '/^S/{print ">"$2;print $3}' ${PREFIX}.bp.p_ctg.gfa > ${PREFIX}.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${PREFIX}.bp.hap1.p_ctg.gfa > ${PREFIX}.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${PREFIX}.bp.hap2.p_ctg.gfa > ${PREFIX}.bp.hap2.p_ctg.fa
cat ${PREFIX}.bp.hap1.p_ctg.fa ${PREFIX}.bp.hap2.p_ctg.fa > ${PREFIX}.bp.hap1hap2.p_ctg.fa