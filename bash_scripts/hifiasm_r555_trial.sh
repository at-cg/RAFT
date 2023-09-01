# Computing overlaps (cis + trans)

EXE=$SCRATCH/tools/hifiasm/hifiasm
READS=../reads.fasta
PREFIX=output

rm -f ${PREFIX}*

/usr/bin/time $EXE -o $PREFIX -t 256 --dbg-ovec $READS
cat output.0.ovlp.paf > output.overlaps.paf
cat output.1.ovlp.paf >> output.overlaps.paf

rm output.0.ovlp.paf output.1.ovlp.pa