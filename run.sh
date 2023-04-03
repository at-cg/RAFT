EXE=$HOME/Desktop/GitHub/chopper/chopper
READS=$HOME/Desktop/GitHub/chopper/reads.fasta
KMERS=$HOME/Desktop/GitHub/chopper/repeats_k21.txt
OUTPUTSEQ=output_reads.fasta
PREFIX=output

rm ${PREFIX}*

$EXE -o ${PREFIX} -p 0.05 -l 50 $READS $KMERS 

grep ">" $OUTPUTSEQ > ${OUTPUTSEQ}.headers
cat ${OUTPUTSEQ}.headers | awk -F '[=,-]' '{print $9"\t"$5"\t"$6}' | sort -k 1,1 -k2,2n -k3,3nr > ${OUTPUTSEQ}.headers.bed