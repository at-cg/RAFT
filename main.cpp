#include <iostream>
#include <unistd.h>

#include "chop.hpp"

void printHelp(const algoParams &params)
{
    std::cout << "Usage: chopper [options] <input-reads.fa> <in.paf>\n";
    std::cout << "  -t NUM     thread count, default " << params.threads << "\n";
    std::cout << "  -r NUM     resolution of masks, repeat annotation, coverage " << params.reso << "\n";
    std::cout << "  -e NUM     estimated coverage " << params.est_cov << "\n";
    std::cout << "  -m NUM     minimum coverage " << params.min_cov << "\n";
    std::cout << "  -f NUM     coverage fraction" << params.cov_frac << "\n";
    std::cout << "  -g NUM     repeat_annotation_gap_thres " << params.repeat_annotation_gap_thres << "\n";
    std::cout << "  -o FILE    prefix of output files " << params.outputfilename << "\n";
    std::cout << "  -L FILE    dump algorithm log\n";
    exit(1);
}

int main(int argc, char *argv[])
{
    algoParams params;

    params.initParams();
    int option;

    while ((option = getopt(argc, argv, "t:r:e:m:f:g:l:o:L:")) != -1)
    {
        switch (option)
        {
        case 't':
            params.threads = atoi(optarg);
            break;
        case 'r':
            params.reso = atoi(optarg);
            break;
        case 'e':
            params.est_cov = atoi(optarg);
            break;
        case 'm':
            params.min_cov = atoi(optarg);
            break;
        case 'f':
            params.cov_frac = atoi(optarg);
            break;
        case 'g':
            params.repeat_annotation_gap_thres = atoi(optarg);
            break;
        case 'l':
            params.repeat_length = atoi(optarg);
            break;
        case 'o':
            params.outputfilename = optarg;
            break;
        case 'L':
            params.logFileName = optarg;
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