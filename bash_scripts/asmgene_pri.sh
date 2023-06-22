ASM=../output.bp.p_ctg.fa
EXE=$SCRATCH/tools/minimap2-2.23_x64-linux/minimap2
REF=$SCRATCH/data/genomes/human/Homo_sapiens.GRCh38.cdna.all.fa
ALN=$SCRATCH/data/genomes/human/cdna_asm_alignment/asmgene.grch38.cdna.paf
PAFTOOLS=$SCRATCH/tools/minimap2-2.23_x64-linux/paftools.js
PREFIX=output

rm -f ${PREFIX}*

export PATH="$PATH:$SCRATCH/tools/minimap2-2.23_x64-linux"
$EXE -cxsplice:hq -t 128 $ASM $REF > asmgene.asm.cdna.paf

$PAFTOOLS asmgene -a $ALN asmgene.asm.cdna.paf  > asmgene.asm.stats