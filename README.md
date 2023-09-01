## <a name="intro"></a>Introduction

Removal of contained reads has long been a weakness of overlap-layout-consensus (OLC) assemblers. RAFT is an algorithm designed to improve assembly quality by rescuing contained reads. RAFT breaks long reads into smaller fragments by following an algorithm described in our [preprint](#papers). The read fragmentation allows an OLC assembler to retain contained reads during string graph construction. The inputs to RAFT is an error-corrected read file in FASTA/FASTQ format and an all-vs-all alignment file in PAF format. It performs read fragmentation and outputs the fragmented reads in FASTA format. 

We recommend users to use [hifiasm](https://github.com/chhylp123/hifiasm) for the initial steps (read error correction, all-vs-all overlap computation) and also for the final step (assembly of fragmented reads). The assembly output format of hifiasm is described [here](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output). The RAFT-hifiasm workflow is designed to work with ONT Duplex, or a mixture of ONT Duplex and HiFI reads. ONT UL reads can optionally be [integrated](https://github.com/chhylp123/hifiasm#ul) during the final assembly step.

## <a name="started"></a>Try RAFT-hifiasm Workflow on Small Test Data

```sh
# Install hifiasm (requiring g++ and zlib)
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make -j4 && cd ..

# Install RAFT 
git clone https://github.com/at-cg/RAFT.git
cd RAFT && make && cd ..

mkdir -p assembly && cd assembly/

# Get small test data
wget https://github.com/chhylp123/hifiasm/releases/download/v0.7/chr11-2M.fa.gz

# First run of hifiasm to obtain error corrected reads and homozygous coverage estimate
../hifiasm/hifiasm -o errorcorrect -t4 -f0 --write-ec chr11-2M.fa.gz 2> errorcorrect.log
COVERAGE=$(grep "homozygous" errorcorrect.log | tail -1 | awk '{print $6}')

# Second run of hifiasm to obtain all-vs-all read overlaps as a paf file
../hifiasm/hifiasm -o getOverlaps -t4 -f0 --dbg-ovec errorcorrect.ec.fa 2> getOverlaps.log
# Merge cis and trans overlaps
cat getOverlaps.0.ovlp.paf getOverlaps.1.ovlp.paf > overlaps.paf

# RAFT fragments the error corrected reads
../RAFT/raft -e ${COVERAGE} -o fragmented errorcorrect.ec.fa overlaps.paf

# Final hifiasm run to obtain assembly of fragmented reads
# A single round of error correction (-r1) is enough here
../hifiasm/hifiasm -o finalasm -t4 -f0 -r1 fragmented.reads.fasta 2> finalasm.log
```

## <a name="use"></a>Usage Details

```sh
raft [options] <input-reads.fa> <in.paf>
```

The following options can be used to customize the behavior of the program:

    -r INT: Set the resolution of local coverage [50].
    -e INT: Set the estimated coverage [0].
    -m NUM: Set the coverage multiplier for high coverage [1.5].
    -l INT: Set the desired read length [20000].
    -p INT: Set the minimum repeat length to be preserved [5000].
    -f INT: Set the flanking length for repeats [1000].
    -v INT: Set the overlap length between fragmented reads [500].
    -o STR: Set the prefix of output files ["raft"].

## <a name="install"></a>Installation

1. Clone the repository:
```sh
git clone https://github.com/MehakBindra/RAFT.git
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
1. Coverage information for each read in `output.coverage.txt`
2. For input reads simulated using seqrequester, it outputs additional information for debugging 
    1. the positions of long repeats in reference contigs in `output.long_repeats.bed`
    2. the positions of long repeats in reads in `output.long_repeats.txt`

## <a name="papers"></a>Preprint
Coming soon
