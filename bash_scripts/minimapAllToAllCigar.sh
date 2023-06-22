EXE=$SCRATCH/tools/minimap2-2.23_x64-linux/minimap2
READS=../reads.fasta
PAFTOOLS=$SCRATCH/tools/minimap2-2.23_x64-linux/paftools.js

/usr/bin/time $EXE -t 256 -w 101 -k 27 -g 500 -B 8 -O 8,48 -E 4,2 -cx ava-ont $READS $READS > overlaps.paf