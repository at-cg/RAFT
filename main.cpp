#include <iostream>
#include <unistd.h>
#include <chrono>

#include "chop.hpp"

void printHelp(const algoParams &params)
{
    std::cout << "Usage: chopper [options] <input-reads.fa> <in.paf>\n";
    std::cout << "  -r NUM     resolution of coverage " << params.reso << "\n";
    std::cout << "  -k NUM     kmer length " << params.kmer_length << "\n";
    std::cout << "  -p NUM     kmer frac " << params.kmer_frac << "\n";
    std::cout << "  -l NUM     repeat_length " << params.repeat_length << "\n";
    std::cout << "  -i NUM     interval_length " << params.interval_length << "\n";
    std::cout << "  -f NUM     flanking_length " << params.flanking_length << "\n";
    std::cout << "  -o FILE    prefix of output files " << params.outputfilename << "\n";
    exit(1);
}

int main(int argc, char *argv[])
{
    algoParams params;

    params.initParams();
    int option;

    while ((option = getopt(argc, argv, "r:k:p:l:i:f:o:a:b:")) != -1)
    {
        switch (option)
        {
        case 'r':
            params.reso = atoi(optarg);
            break;
        case 'k':
            params.kmer_length = atoi(optarg);
            break;
        case 'p':
            params.kmer_frac = atof(optarg);
            break;
        case 'l':
            params.repeat_length = atoi(optarg);
            break;
        case 'i':
            params.interval_length = atoi(optarg);
            break;
        case 'f':
            params.flanking_length = atoi(optarg);
            break;
        case 'o':
            params.outputfilename = optarg;
            break;
        case 'a':
            params.additional_kmers = atoi(optarg);
            break;
        case 'b':
            params.bloom_filter_element_count = atoi(optarg);
            break;
        default:
            printHelp(params);
        }
    }

    //print usage
    if (argc <= optind + 1)
        printHelp(params);

    params.printParams();

    auto tStart = std::chrono::system_clock::now();
    std::cout << "INFO, main(), started timer\n";

    break_long_reads(argv[optind], argv[optind + 1], params);

    std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - tStart);
    std::cout << "INFO, main(), program completed after " << wctduration.count() << " seconds\n";

    // log complete command given by user
    fprintf(stdout, "INFO, %s(), CMD:", __func__);
    for (int i = 0; i < argc; ++i)
        fprintf(stdout, " %s", argv[i]);
    std::cout << "\n";

    return 0;
}