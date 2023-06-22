EXE=$SCRATCH/tools/minigraph/minigraph
CONTIGS=../output.bp.p_ctg.fa
REF=$SCRATCH/data/genomes/human/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
PREFIX=output

rm -r ${PREFIX}

/usr/bin/time $EXE -xasm --show-unmap=yes -t 256 -o ${PREFIX}.paf $REF $CONTIGS
export PATH="$PATH:/pscratch/sd/m/mehak/tools/minimap2-2.23_x64-linux/"
/pscratch/sd/m/mehak/tools/minimap2-2.23_x64-linux/paftools.js asmstat ${REF}.fai  ${PREFIX}.paf