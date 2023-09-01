# Identify regions where long repeats are located

BEDTOOLS=tools/bedtools2/bedtools

cat chopper.long_repeats.bed | awk -v EX=20000 '{print $1,$2-EX,$3+EX}' OFS='\t' > chopper.long_repeats.extended.bed
cat chopper.long_repeats.extended.bed | awk '{if ($2<0) print $1,0,$3; else print $1,$2,$3;}' OFS='\t' > chopper.long_repeats.extended.pos.bed
$BEDTOOLS sort -i chopper.long_repeats.extended.pos.bed > chopper.long_repeats.extended.pos.sorted.bed
$BEDTOOLS merge -d 50 -i chopper.long_repeats.extended.pos.sorted.bed > chopper.long_repeats.extended.pos.sorted.merged.bed