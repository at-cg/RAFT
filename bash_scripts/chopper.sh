EXE=$SCRATCH/code/chopper3/chopper/chopper
READS=$SCRATCH/data/simulated_reads/human_diploid_30x_hifiasm.trio.0.16.1_hifi/reads.fasta
OVERLAPS=$SCRATCH/data/simulated_reads/human_diploid_30x_hifiasm.trio.0.16.1_hifi/hifiasm_r555_trial/output.overlaps.paf
SEQKIT=$SCRATCH/tools/seqKit-v2.0.0/seqkit
BEDTOOLS=$SCRATCH/tools/bedtools2/bedtools
SEQTK=$SCRATCH/tools/seqtk-1.3/seqtk
GENOMESIZE=$SCRATCH/data/genomes/human/HG002.hifiasm.trio.0.16.1.size
OUTPUTSEQ=output_reads.fasta
PREFIX=chopper

rm ${PREFIX}*

/usr/bin/time $EXE -e 30 -i 5000 -p 5000 -o ${PREFIX} $READS $OVERLAPS 

grep ">" $OUTPUTSEQ > ${OUTPUTSEQ}.headers
cat ${OUTPUTSEQ}.headers | awk -F '[=,-]' '{print $9"\t"$5"\t"$6}' | sort -k 1,1 -k2,2n -k3,3nr > ${OUTPUTSEQ}.headers.bed
$BEDTOOLS genomecov -i ${OUTPUTSEQ}.headers.bed -g $GENOMESIZE -bga | awk '{if ($4==0) print $0}' > ${OUTPUTSEQ}.headers.nocov.bed
$BEDTOOLS genomecov -i ${OUTPUTSEQ}.headers.bed -g $GENOMESIZE > ${OUTPUTSEQ}.headers.cov.bed
cat $OUTPUTSEQ | head -n 1000000 | $SEQKIT seq -g -m 1000 | $SEQKIT watch --fields ReadLen -Q -O readlen_min1k.hist.pdf
cat $OUTPUTSEQ | $SEQKIT stats -a > ${OUTPUTSEQ}.stats