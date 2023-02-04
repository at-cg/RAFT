EXE=$HOME/Desktop/GitHub/chopper/chopper
READS=$HOME/Desktop/GitHub/chopper/reads.fasta
OVERLAPS=$HOME/Desktop/GitHub/chopper/overlaps.paf
OUTPUTSEQ=output_reads.fasta
PREFIX=output

rm ${PREFIX}*

/usr/bin/time $EXE -e 20 -o ${PREFIX} -m 1.3 $READS $OVERLAPS 

grep ">" $OUTPUTSEQ > ${OUTPUTSEQ}.headers
cat ${OUTPUTSEQ}.headers | awk -F '[=,-]' '{print $9"\t"$5"\t"$6}' | sort -k 1,1 -k2,2n -k3,3nr > ${OUTPUTSEQ}.headers.bed