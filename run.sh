EXE=$HOME/Desktop/GitHub/chopper/chopper
READS=$HOME/Desktop/GitHub/chopper/reads.fasta
OVERLAPS=$HOME/Desktop/GitHub/chopper/overlaps.paf
KMERS=$HOME/Desktop/GitHub/chopper/repeats_k21.txt
OUTPUTSEQ=output_reads.fasta
PREFIX=output

rm ${PREFIX}*

$EXE -e 20 -m 1.3 -p 0.1 -l 1000 -o ${PREFIX} $READS $OVERLAPS $KMERS

grep ">" $OUTPUTSEQ > ${OUTPUTSEQ}.headers
cat ${OUTPUTSEQ}.headers | awk -F '[=,-]' '{print $9"\t"$5"\t"$6}' | sort -k 1,1 -k2,2n -k3,3nr > ${OUTPUTSEQ}.headers.bed