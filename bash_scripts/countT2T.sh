#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --constraint=cpu
#SBATCH --output=BATCH_OUTPUT
#SBATCH --error=BATCH_OUTPUT

REF=/global/cfs/cdirs/m4320/references/human/T2T/chm13.draft_v2.0.fasta
HAP1=/path/to/asm.bp.hap1.p_ctg.fa
HAP2=/path/to/asm.bp.hap2.p_ctg.fa
HAP12=hap1hap2.p_ctg.fa
WORKFLOW=HPRC_assessAsemblyCompletness.wdl  #in RAFT/bash_scripts/
EXE=/path/to/cromwell-85.jar  #.jar file available at https://github.com/broadinstitute/cromwell/releases/tag/85
THREADS=32

#load dependencies
#Ensure that the executables "seqtk", "NCRF", "bioawk", "mashamp" are available in search path
#https://anaconda.org/bioconda/seqtk
#https://anaconda.org/bioconda/ncrf
#https://anaconda.org/bioconda/mashmap
#https://anaconda.org/bioconda/bioawk

#create input in JSON format
JSON_FMT='{"assessAssemblyCompletness.reference":"%s","assessAssemblyCompletness.assembly":"%s","assessAssemblyCompletness.threadCount":"%s"}\n'
printf "$JSON_FMT" "$REF" "$HAP12" "$THREADS" > input.json


cat $HAP1 $HAP2 > $HAP12
java -jar $EXE run $WORKFLOW --inputs input.json

#T2T contigs are saved in a output file
#To find the location of the output file, run:   grep -A2 "outputs" BATCH_OUTPUT
