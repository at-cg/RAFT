# Computing switch error rate using yak trioeval

EXE=tools/yak/yak
READS=data/input/reads.fasta
HAP1=asm/output.bp.hap1.p_ctg.fa
HAP2=asm/output.bp.hap2.p_ctg.fa
ASM=asm/output.bp.hap1hap2.p_ctg.fa # ASM file is the concatenation of HAP1 and HAP2

grep -A 1 "h1tg" $READS > hap1.fasta
grep -A 1 "h2tg" $READS > hap2.fasta

/usr/bin/time $EXE count -b37 -t 128 -o hap1.yak hap1.fasta
/usr/bin/time $EXE count -b37 -t 128 -o hap2.yak hap2.fasta
/usr/bin/time $EXE trioeval -t 128 hap1.yak hap2.yak $ASM > both.out
/usr/bin/time $EXE trioeval -t 128 hap1.yak hap2.yak $HAP1 > hap1.out
/usr/bin/time $EXE trioeval -t 128 hap1.yak hap2.yak $HAP2 > hap2.out
rm hap1.fasta hap2.fasta
