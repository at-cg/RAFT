# Running asmgene to get multicopy genes retained and gene completeness in primary assembly

ASM=asm/output.bp.p_ctg.fa
EXE=tools/minimap2-2.23_x64-linux/minimap2
REF=data/genomes/human/Homo_sapiens.GRCh38.cdna.all.fa
ALN=data/genomes/human/cdna_asm_alignment/asmgene.grch38.cdna.paf
PAFTOOLS=tools/minimap2-2.23_x64-linux/paftools.js
PREFIX=output

rm -f ${PREFIX}*

export PATH="$PATH:$SCRATCH/tools/minimap2-2.23_x64-linux"
$EXE -cxsplice:hq -t 128 $ASM $REF > asmgene.asm.cdna.paf

$PAFTOOLS asmgene -a $ALN asmgene.asm.cdna.paf  > asmgene.asm.stats