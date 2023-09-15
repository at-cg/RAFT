#!/bin/bash

# https://anaconda.org/bioconda/tidk
TIDK=/path/to/tidk
SAMTOOLS=/path/to/samtools
COUNTER=~/bash_scripts/T2T-count.pl

SEQ=/path/to/concatenated/assembly/contigs/output.bp.hap1hap2.p_ctg.fa

PREFIX=output
TELOMERE=TTAGGG

$SAMTOOLS --version
/usr/bin/time $SAMTOOLS faidx $SEQ

if [ ! -d t2tCounts/ ]; then
	mkdir t2tCounts
fi

#input 10000 for window parameter is default
$TIDK --version
$TIDK search --dir ./t2tCounts/ --extension tsv --output $PREFIX --string $TELOMERE --window 10000 $SEQ
sed 's/\t/,/g' ./t2tCounts/${PREFIX}_telomeric_repeat_windows.tsv > ./t2tCounts/${PREFIX}_telomeric_repeat_windows.csv

#100000 is the length of prefix and suffix of contig where $TELOMERE is searched; 100 is the lower bound cutoff;
perl $COUNTER ./t2tCounts/${PREFIX}_telomeric_repeat_windows.csv ${SEQ}.fai 100000 100 > ./t2tCounts/t2tSearch.out

#list out all contigs which satisfy telomere count criterion
awk '$3>100 && $4>100' ./t2tCounts/t2tSearch.out > ./t2tCounts/search.list

# move output of perlscript to t2tCounts folder
mv Plot_T2T.txt ./t2tCounts/Plot_T2T.txt
