# Compute length of unresolved gaps

HAP1=data/genomes/human/HG002.hifiasm.trio.0.16.1.hap1.fa
HAP2=data/genomes/human/HG002.hifiasm.trio.0.16.1.hap2.fa
BEDTOOLS=tools/bedtools2/bedtools
READS=data/input/reads.fasta
EXE=tools/minimap2-2.23_x64-linux/minimap2
READIDS=asm/output.bp.r_utg.gfa.contained.headers
GENOMESIZE=data/genomes/human/HG002.hifiasm.trio.0.16.1.size
GAPS=../../map_mm_noncontained/mm2.exactmapped.nocov.merged.subtracted.noends.100.bed
SEQTK=$SCRATCH/tools/seqtk-1.3/seqtk

/usr/bin/time $SEQTK subseq $READS $READIDS > non-redundant.fasta
/usr/bin/time $EXE -t 32 -N 50 -cx map-ont $HAP1 non-redundant.fasta > mm2.paf
/usr/bin/time $EXE -t 32 -N 50 -cx map-ont $HAP2 non-redundant.fasta >> mm2.paf
cat mm2.paf | awk '{if ($3 == 0 && $2 == $4 && $2 == $10) print $6"\t"$8"\t"$9}' > tmp
cat mm2.paf | awk '{if ($8 == 0 && $7 == $9 && $7 == $10) print $6"\t"$8"\t"$9}' >> tmp
cat tmp | sort -k 1,1 -k2,2n -k3,3nr >  mm2.exactmapped.bed
$BEDTOOLS subtract -a $GAPS -b mm2.exactmapped.bed > unresolved_gaps.bed
cat unresolved_gaps.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' > unresolved_gaps.bed.sum
cat unresolved_gaps.bed | awk '{print ($3-$2)}' | sort -n > unresolved_gaps.bed.lengths.sorted
rm tmp