## <a name="intro"></a>Introduction

Removal of contained reads has long been a weakness of overlap-layout-consensus (OLC) assemblers. RAFT (**R**epeat **A**ware **F**ragmentation **T**ool) is an algorithm designed to improve assembly quality by rescuing contained reads. RAFT breaks long reads into smaller sub-reads by following an algorithm described in our [preprint](#papers). The read fragmentation allows an OLC assembler to retain contained reads during string graph construction. When input reads have non-uniform lengths, retaining contained reads improves assembly contiguity and base-level accuracy. The inputs to RAFT include an error-corrected read file in FASTA/FASTQ format and an all-vs-all alignment file in PAF format. It performs read fragmentation and outputs the fragmented reads in FASTA format. 

We recommend users to use [hifiasm](https://github.com/chhylp123/hifiasm) for the initial steps (read error correction, all-vs-all overlap computation) and also for the final step (assembly of fragmented reads). The assembly output format of hifiasm is described [here](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output). 

_The RAFT-hifiasm workflow is recommended for long accurate reads with non-uniform length distribution (e.g., ONT Duplex, or a mixture of ONT Duplex and HiFi reads). ONT UL reads can optionally be [integrated](https://github.com/chhylp123/hifiasm#ul) during the final assembly step._

## <a name="started"></a>Try RAFT-hifiasm Workflow on Small Test Data
The entire test workflow below will take about 3-4 minutes. Users can either run the commands one by one or copy the commands into an executable script.

```sh
# Install RAFT 
git clone https://github.com/at-cg/RAFT.git
cd RAFT && make && cd ..

# Install hifiasm (requiring g++ and zlib)
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make -j4 && cd ..

mkdir -p assembly && cd assembly/

# Get small test data
wget https://github.com/chhylp123/hifiasm/releases/download/v0.7/chr11-2M.fa.gz

# First run of hifiasm with 4 threads to obtain error corrected reads and coverage estimate
../hifiasm/hifiasm -o errorcorrect -t4 --write-ec chr11-2M.fa.gz 2> errorcorrect.log
COVERAGE=$(grep "homozygous" errorcorrect.log | tail -1 | awk '{print $6}')

# Second run of hifiasm to obtain all-vs-all read overlaps as a paf file
../hifiasm/hifiasm -o getOverlaps -t4 --dbg-ovec errorcorrect.ec.fa 2> getOverlaps.log
# Merge cis and trans overlaps
cat getOverlaps.0.ovlp.paf getOverlaps.1.ovlp.paf > overlaps.paf

# RAFT fragments the error corrected reads
../RAFT/raft -e ${COVERAGE} -o fragmented errorcorrect.ec.fa overlaps.paf

# Final hifiasm run to obtain assembly of fragmented reads
# A single round of error correction (-r1) is enough here
../hifiasm/hifiasm -o finalasm -t4 -r1 fragmented.reads.fasta 2> finalasm.log
ls finalasm*p_ctg.gfa
```
For large inputs, users are recommended to increase the thread count depending on the number of the cores available for use. RAFT-hifiasm workflow takes about 9 hours and ~100 GB RAM using 128 threads on a multicore [Perlmutter CPU-based node](https://docs.nersc.gov/systems/perlmutter/architecture/) to process 32x ONT Duplex human data.

## <a name="use"></a>Usage Details

```sh
raft [options] <input-reads.fa> <in.paf>
```

The following options can be used to customize the behavior of the program. The default values are set if there is no custom requirement.

    -r INT: Set the resolution of local coverage [50]
    -e INT: Set the estimated coverage
    -m NUM: Set the coverage multiplier for high coverage [1.5]
    -l INT: Set the desired read length [20000]
    -p INT: Set the minimum repeat length to be preserved [5000]
    -f INT: Set the flanking length for repeats [1000]
    -v INT: Set the overlap length between fragmented reads [500]
    -o STR: Set the prefix of output files ["raft"]

## <a name="install"></a>Installation

1. Clone the repository:
```sh
git clone https://github.com/at-cg/RAFT.git
```

2. Compile the source code:
```sh
cd RAFT
make
```

## <a name="examples"></a>Examples

1. Run RAFT with estimated coverage 20:
```sh
raft -e 20 -m 1.3 -o output <input_reads> <input_overlaps>
```

2. Run RAFT with custom parameters:
```sh
raft -e 20 -m 1.3 -p 7000 -f 500 -v 500 -l 15000 -o output <input_reads> <input_overlaps>
```

## <a name="output"></a>Output Files
RAFT outputs the following files:
1. Coverage information for each read in \`prefix\`.coverage.txt
2. For input reads simulated using seqrequester, it outputs additional information for debugging 
    1. the positions of long repeats in reference contigs in \`prefix\`.long_repeats.bed
    2. the positions of long repeats in reads in \`prefix\`.long_repeats.txt

## <a name="papers"></a>Preprint
Sudhanva Shyam Kamath, Mehak Bindra, Debnath Pal, Chirag Jain. [Telomere-to-telomere assembly by preserving contained reads](https://doi.org/10.1101/2023.11.07.565066). Biorxiv (November 2023).
