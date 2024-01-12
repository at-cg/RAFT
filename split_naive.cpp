#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "kseq.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

int splitSeq(const char *filename, std::string outputFilename, int subreadLength)
{
  gzFile fp;
  kseq_t *seq;
  int l;
  fp = gzopen(filename, "r");
  seq = kseq_init(fp);
  int num = 0;
  std::string header, sequence, line;
  std::vector<std::string> subreads;
  std::ofstream outputFile(outputFilename);

  while ((l = kseq_read(seq)) >= 0)
  { 
    header = seq->name.s;
    sequence = seq->seq.s;

    for (size_t i = 0; i < sequence.length(); i += subreadLength) {
      subreads.push_back(sequence.substr(i, subreadLength));
    }

    for (size_t i = 0; i < subreads.size(); ++i) {
      outputFile << ">" << header << "_" << i + 1 << '\n' << subreads[i] << '\n';
    }

    subreads.clear();
    num++;
  }

  kseq_destroy(seq);
  gzclose(fp);
  outputFile.close();

  return num;
}

void printHelp()
{
  std::cout << "Purpose: Split input reads naively into non-overlapping subreads. The output format is FASTA\n";
  std::cout << "Usage: split_naive <inputfilename> <outputfilename> SPLITLEN\n";
  std::cout << "Example: split_naive input.fastq output.fragmented.fasta 20000\n";
  exit(1);
}

int main(int argc, char** argv) 
{
  //Usage: Input file, output file, subread length
  if (argc < 4)
    printHelp();

  splitSeq(argv[1], argv[2], std::stoi(argv[3]));
  return 0;
}




