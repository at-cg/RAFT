#include <iostream>
#include <unistd.h>
#include <chrono>

#include "chop.hpp"

void printHelp(const algoParams &params)
{
    std::cout << "Usage: chopper [options] <input-reads.fa> <in.paf>\n";
    std::cout << "  -r NUM     resolution of coverage " << params.reso << "\n";
    std::cout << "  -e NUM     estimated coverage " << params.est_cov << "\n";
    std::cout << "  -m NUM     coverage multiplier " << params.cov_mul << "\n";
    std::cout << "  -l NUM     repeat_length " << params.repeat_length << "\n";
    std::cout << "  -v NUM     overlap_length " << params.overlap_length << "\n";
    std::cout << "  -u NUM     uniform_read_length " << params.uniform_read_length << "\n";
    std::cout << "  -t NUM     read_length_threshold " << params.read_length_threshold << "\n";
    std::cout << "  -o FILE    prefix of output files " << params.outputfilename << "\n";
    std::cout << "  -d NUM     real reads " << params.real_reads << "\n";
    std::cout << "  -p NUM     hifiasm overlaps " << params.hifiasm_overlaps << "\n";
    exit(1);
}

int main(int argc, char *argv[])
{
    algoParams params;

    params.initParams();
    int option;

    while ((option = getopt(argc, argv, "r:e:m:l:v:u:t:o:d:p:")) != -1)
    {
        switch (option)
        {
        case 'r':
            params.reso = atoi(optarg);
            break;
        case 'e':
            params.est_cov = atoi(optarg);
            break;
        case 'm':
            params.cov_mul = std::stod(optarg);
            break;
        case 'l':
            params.repeat_length = atoi(optarg);
            break;
        case 'v':
            params.overlap_length = atoi(optarg);
            break;
        case 'u':
            params.uniform_read_length = atoi(optarg);
            break;
        case 't':
            params.read_length_threshold = atoi(optarg);
            break;
        case 'o':
            params.outputfilename = optarg;
            break;
        case 'd':
            params.real_reads = atoi(optarg);
            break;
        case 'p':
            params.hifiasm_overlaps = atoi(optarg);
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