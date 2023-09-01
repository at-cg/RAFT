# Running asmgene to get haplotype resolved multicopy genes retained and gene completeness

HAP1=asm/output.bp.hap1.p_ctg.fa
HAP2=asm/output.bp.hap2.p_ctg.fa
EXE=tools/minimap2-2.23_x64-linux/minimap2
REF=data/genomes/human/Homo_sapiens.GRCh38.cdna.all.fa
ALN=data/genomes/human/cdna_asm_alignment/asmgene.grch38.cdna.paf
PAFTOOLS=tools/minimap2-2.23_x64-linux/paftools.js
PREFIX=output

rm -f ${PREFIX}*

export PATH="$PATH:$SCRATCH/tools/minimap2-2.23_x64-linux"
$EXE -cxsplice:hq -t 128 $HAP1 $REF > asmgene.hap1.cdna.paf
$EXE -cxsplice:hq -t 128 $HAP2 $REF > asmgene.hap2.cdna.paf

$PAFTOOLS asmgene -a $ALN asmgene.hap1.cdna.paf  > asmgene.hap1.stats
$PAFTOOLS asmgene -a $ALN asmgene.hap2.cdna.paf  > asmgene.hap2.stats