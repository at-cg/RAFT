READS=../reads.fasta
OVERLAPS=../minimapAllToAllCigar/overlaps.paf
SEQTK=$SCRATCH/tools/seqtk-1.3/seqtk
EXE=$SCRATCH/tools/minimap2-2.23_x64-linux/minimap2
BEDTOOLS=$SCRATCH/tools/bedtools2/bedtools
HAP1=$SCRATCH/data/genomes/human/HG002.hifiasm.trio.0.16.1.hap1.fa
HAP2=$SCRATCH/data/genomes/human/HG002.hifiasm.trio.0.16.1.hap2.fa
GENOMESIZE=$SCRATCH/data/genomes/human/HG002.hifiasm.trio.0.16.1.size

for MINIDENTITY in 100
do
  echo "MINIDENTITY=" ${MINIDENTITY}
  cat ${OVERLAPS} | awk -v minidnty="$MINIDENTITY" '{if ($3 == 0 && $2 == $4 && $2 < $7 && $10*100.0/$11 >= minidnty) print $0}' | cut -f1 > tmp1
  cat ${OVERLAPS} | awk -v minidnty="$MINIDENTITY" '{if ($8 == 0 && $7 == $9 && $7 < $2 && $10*100.0/$11 >= minidnty) print $0}' | cut -f6 > tmp2
  cat tmp1 tmp2 | sort | uniq > contained.${MINIDENTITY}.txt
  cat ../reads.fasta.headers | sort | sed 's/>//g' > reads.fasta.sorted.headers
  comm -23 reads.fasta.sorted.headers contained.${MINIDENTITY}.txt > non-contained.${MINIDENTITY}.txt
  wc -l contained.${MINIDENTITY}.txt
  wc -l non-contained.${MINIDENTITY}.txt
  $SEQTK subseq $READS non-contained.${MINIDENTITY}.txt > non-contained.fasta
  /usr/bin/time $EXE -t 32 -N 50 -cx map-ont $HAP1 non-contained.fasta > mm2.${MINIDENTITY}.paf
  /usr/bin/time $EXE -t 32 -N 50 -cx map-ont $HAP2 non-contained.fasta >> mm2.${MINIDENTITY}.paf
  cat mm2.${MINIDENTITY}.paf | awk '{if ($3 == 0 && $2 == $4 && $2 == $10) print $6"\t"$8"\t"$9}' | sort -k 1,1 -k2,2n -k3,3nr >  mm2.exactmapped.${MINIDENTITY}.bed
  $BEDTOOLS genomecov -i mm2.exactmapped.${MINIDENTITY}.bed -g $GENOMESIZE -bga | awk '{if ($4==0) print $0}' > mm2.exactmapped.nocov.${MINIDENTITY}.bed
  $BEDTOOLS merge -d 500 -i mm2.exactmapped.nocov.${MINIDENTITY}.bed >  mm2.exactmapped.nocov.merged.${MINIDENTITY}.bed
  $BEDTOOLS subtract -A -a mm2.exactmapped.nocov.merged.${MINIDENTITY}.bed -b ../reads.fasta.headers.nocov.bed > mm2.exactmapped.nocov.merged.subtracted.${MINIDENTITY}.bed
  $BEDTOOLS genomecov -i mm2.exactmapped.${MINIDENTITY}.bed -g $GENOMESIZE > mm2.exactmapped.cov.${MINIDENTITY}.hist
  python3 $SCRATCH/code/printEdgeBedIntervals.py $GENOMESIZE 25000 > genome_25kbp_ends.bed
  $BEDTOOLS subtract -A -a mm2.exactmapped.nocov.merged.subtracted.${MINIDENTITY}.bed -b genome_25kbp_ends.bed > mm2.exactmapped.nocov.merged.subtracted.noends.${MINIDENTITY}.bed
  cat mm2.exactmapped.nocov.merged.subtracted.noends.${MINIDENTITY}.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' > mm2.exactmapped.nocov.merged.subtracted.noends.${MINIDENTITY}.bed.sum
  cat mm2.exactmapped.nocov.merged.subtracted.noends.${MINIDENTITY}.bed | awk '{print ($3-$2)}' | sort -n >  mm2.exactmapped.nocov.merged.subtracted.noends.${MINIDENTITY}.bed.lengths.sorted
done

rm tmp1 tmp2 non-contained.fasta  