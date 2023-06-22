HAP1=../output.bp.hap1.p_ctg.fa
HAP2=../output.bp.hap2.p_ctg.fa
EXE=$SCRATCH/tools/minimap2-2.23_x64-linux/minimap2
REF=$SCRATCH/data/genomes/human/Homo_sapiens.GRCh38.cdna.all.fa
ALN=$SCRATCH/data/genomes/human/cdna_asm_alignment/asmgene.grch38.cdna.paf
PAFTOOLS=$SCRATCH/tools/minimap2-2.23_x64-linux/paftools.js
PREFIX=output

rm -f ${PREFIX}*

export PATH="$PATH:$SCRATCH/tools/minimap2-2.23_x64-linux"
$EXE -cxsplice:hq -t 128 $HAP1 $REF > asmgene.hap1.cdna.paf
$EXE -cxsplice:hq -t 128 $HAP2 $REF > asmgene.hap2.cdna.paf

$PAFTOOLS asmgene -a $ALN asmgene.hap1.cdna.paf  > asmgene.hap1.stats
$PAFTOOLS asmgene -a $ALN asmgene.hap2.cdna.paf  > asmgene.hap2.stats