## <a name="started"></a>Getting Started

```sh
Assumed directory structure
./
|____data/
|____hifiasm/
|____RAFT/
|____errorCorrect/
|____overlaps/
|____assembly/

# Install hifiasm (requiring g++ and zlib)
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make
cd ..

# Install raft 
git clone https://github.com/username/RAFT.git
cd raft && make
cd ..

# Run on test data (use -f0 for small datasets)
mkdir -p data && cd data
wget https://github.com/chhylp123/hifiasm/releases/download/v0.7/chr11-2M.fa.gz
cd ..

# First run of hifiasm to obtain error corrected reads and homozygous coverage estimate
cd errorCorrect/
../hifiasm/hifiasm -o output -t4 -f0 --write-ec ./data/chr11-2M.fa.gz 2> output.log
# Extract primary contigs in FASTA format
awk '/^S/{print ">"$2;print $3}' output.bp.p_ctg.gfa > output.p_ctg.fa
mv output.ec.fa ../data/error_free_reads.fa
cd ../

# Extract estimated coverage of dataset
COVERAGE=$(grep "homozygous" ./hifiasm/firstRun.log | tail -1 | awk '{print $6}')

# Second run of hifiasm to obtain cis and trans overlaps as a paf file
cd overlaps/
../hifiasm/hifiasm -o output -t4 -f0 --dbg-ovec ../data/error_free_reads.fa 2> output.log
# Merge cis and trans overlaps into a single PAF file
cat output.0.ovlp.paf output.1.ovlp.paf > ../data/overlaps.paf
cd ../

# RAFT fragments the error free reads using the overlap information
# Repeats longer than 5000 are preserved in the fragmented reads
./raft/raft -e ${COVERAGE} -p 5000 -o output ./data/errorFree.ec.fa ./data/overlaps.paf
mv output_reads.fa ./data/fragmented_reads.fa

# Final hifiasm run to obtain assembly from fragmented set of reads.
cd assembly/
../hifiasm/hifiasm -o asm -t4 -f0 -r1 ../data/fragmented_reads.fa 2> asm.log
# Extract primary contigs in FASTA format
awk '/^S/{print ">"$2;print $3}' ./hifiasm/asm.bp.p_ctg.gfa > asm.p_ctg.fa
```

## Table of Contents
- [Getting Started](#started)
- [Introduction](#intro)
- [Usage](#use)
    - [Options](#opt)
- [Installation](#install)
- [Examples](#examples)
- [Output Files](#output)

## <a name="intro"></a>Introduction

RAFT is a pipeline designed to improve assembly continuity. It contains a command-line tool, `raft`, designed to break long reads into smaller fragments based on specified parameters. RAFT takes input read files in FASTA/FASTQ format and alignment files in PAF format and performs read fragmentation according to the given parameters. The resulting fragments are outputted to `output_reads.fasta`.

## <a name="use"></a>Usage

```sh
raft [options] <input-reads.fa> <in.paf>
```

### <a name="opt"></a>Options

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
2. For simulated reads, it outputs 
    1. the positions of long repeats in reference contigs in `output.long_repeats.bed`
    2. the positions of long repeats in reads in `output.long_repeats.txt`