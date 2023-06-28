## <a name="intro"></a>Introduction

The chopper program is a command-line tool designed to break long reads into smaller fragments based on specified parameters. It takes input read files in FASTA/FASTQ format and alignment files in PAF format and performs read fragmentation according to the given parameters. The resulting fragments are outputted to output_reads.fasta. 

## <a name="use"></a>Usage

```sh
chopper [options] <input-reads.fa> <in.paf>
```

Options

The following options can be used to customize the behavior of the program:

    -r INT: Set the resolution of local coverage [50].
    -e INT: Set the estimated coverage [0].
    -m NUM: Set the coverage multiplier for high coverage [1.5].
    -l INT: Set the desirec read length [20000].
    -p INT: Set the minimum repeat length to be preserved [5000].
    -f INT: Set the flanking length for repets [1000].
    -v INT: Set the overlap length between fragmented reads [500].
    -o STR: Set the prefix of output files ["chopper"].

## <a name="install"></a>Installation

1. Clone the repository:
```sh
git clone https://github.com/username/chopper.git
```

2.Compile the source code:
```sh
cd chopper
make
```

## <a name="examples"></a>Examples

1. Run Chopper with estimated coverage 20:
```sh
chopper -e 20 -m 1.3 -o output <input_reads> <input_overlaps>
```

2. Run Chopper with custom parameters:
```sh
chopper -e 20 -m 1.3 -p 7000 -f 500 -v 500 -l 15000 -o output <input_reads> <input_overlaps>
```

## <a name="output"></a>Output Files
Chopper outputs coverage information for each read and positions of long repeats in case of simulated reads.