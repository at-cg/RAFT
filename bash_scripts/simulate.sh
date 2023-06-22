EXE=$SCRATCH/tools/seqrequester/seqrequester
HAP1=$SCRATCH/data/genomes/human/HG002.hifiasm.trio.0.16.1.hap1.fa
HAP2=$SCRATCH/data/genomes/human/HG002.hifiasm.trio.0.16.1.hap2.fa
COVERAGE=15
DIST=$SCRATCH/data/dist/HIFI/m64043_200904_190723.ccs.fastq.hist.txt
OUTPUTSEQ=reads.fasta
GENOMESIZE=$SCRATCH/data/genomes/human/HG002.hifiasm.trio.0.16.1.size
SEQKIT=$SCRATCH/tools/seqKit-v2.0.0/seqkit
BEDTOOLS=$SCRATCH/tools/bedtools2/bedtools

SIZE=`grep -v ">" $HAP1 | wc | awk '{print $3-$1}'`
$EXE simulate -truncate -genome $HAP1 -genomesize $SIZE -coverage $COVERAGE -distribution $DIST > $OUTPUTSEQ.hap1
SIZE=`grep -v ">" $HAP2 | wc | awk '{print $3-$1}'`
$EXE simulate -truncate -genome $HAP2 -genomesize $SIZE -coverage $COVERAGE -distribution $DIST > $OUTPUTSEQ.hap2

cat $OUTPUTSEQ.hap1 > $OUTPUTSEQ
cat $OUTPUTSEQ.hap2 >> $OUTPUTSEQ

grep ">" $OUTPUTSEQ > ${OUTPUTSEQ}.headers
cat ${OUTPUTSEQ}.headers | awk -F '[=,-]' '{print $9"\t"$5"\t"$6}' | sort -k 1,1 -k2,2n -k3,3nr > ${OUTPUTSEQ}.headers.bed
$BEDTOOLS genomecov -i ${OUTPUTSEQ}.headers.bed -g $GENOMESIZE -bga | awk '{if ($4==0) print $0}' > ${OUTPUTSEQ}.headers.nocov.bed
$BEDTOOLS genomecov -i ${OUTPUTSEQ}.headers.bed -g $GENOMESIZE > ${OUTPUTSEQ}.headers.cov
cat $OUTPUTSEQ | head -n 1000000 | $SEQKIT seq -g -m 1000 | $SEQKIT watch --fields ReadLen -Q -O readlen_min1k.hist.pdf
cat $OUTPUTSEQ | $SEQKIT stats -a > ${OUTPUTSEQ}.stats

rm $OUTPUTSEQ.hap1 $OUTPUTSEQ.hap2