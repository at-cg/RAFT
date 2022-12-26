EXE=$HOME/Desktop/GitHub/chopper/chopper
READS=$HOME/Desktop/GitHub/chopper/reads.fasta
OVERLAPS=$HOME/Desktop/GitHub/chopper/overlaps.paf

OUTPUTSEQ=output_reads.fasta

$EXE -e 20 -l 5000 $READS $OVERLAPS 

# $EXE -a 0 $READS $OVERLAPS 

grep ">" $OUTPUTSEQ > ${OUTPUTSEQ}.headers
cat ${OUTPUTSEQ}.headers | awk -F '[=,-]' '{print $9"\t"$5"\t"$6}' | sort -k 1,1 -k2,2n -k3,3nr > ${OUTPUTSEQ}.headers.bed
python countcontainedreads.py ${OUTPUTSEQ}.headers.bed