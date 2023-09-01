# Calculate QV for Hap1 and Hap2

EXE=tools/yak/yak
ASM=asm/output.bp.hap1hap2.p_ctg.fa
KMERS=yak_trial/output.yak

/usr/bin/time $EXE qv -t 128 -p -K3.2g -l100k $KMERS $ASM