EXE=$HOME/Desktop/GitHub/RAFT/raft
READS=$HOME/Desktop/GitHub/RAFT/reads.fasta
OVERLAPS=$HOME/Desktop/GitHub/RAFT/overlaps.paf
OUTPUTSEQ=output_reads.fasta
PREFIX=output

rm ${PREFIX}*

$EXE -e 20 -m 1.3 -o ${PREFIX} $READS $OVERLAPS

grep ">" $OUTPUTSEQ > ${OUTPUTSEQ}.headers
cat ${OUTPUTSEQ}.headers | awk -F '[=,-]' '{print $9"\t"$5"\t"$6}' | sort -k 1,1 -k2,2n -k3,3nr > ${OUTPUTSEQ}.headers.bed