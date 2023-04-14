#ifndef PARAM_CHOPPER_H
#define PARAM_CHOPPER_H

    struct algoParams
{
    int reso;
    int est_cov;
    double cov_mul;
    int repeat_length;
    int interval_length;
    int read_length;
    int overlap_length;
    int flanking_length;
    std::string outputfilename;
    int real_reads;
    int symmetric_overlaps;

    void initParams()
    {
        reso = 50;
        est_cov = 0;
        cov_mul = 1.5;
        repeat_length = 1000;
        interval_length = 10000;
        read_length = 20000;
        overlap_length = 500;
        flanking_length = 2000;
        outputfilename = "chopper";
        real_reads = 1;
        symmetric_overlaps = 0;
    }

    void printParams()
    {
        std::cout << "INFO, printParams(), reso = " << reso << "\n";
        std::cout << "INFO, printParams(), est_cov = " << est_cov << "\n";
        std::cout << "INFO, printParams(), cov_mul = " << cov_mul << "\n";
        std::cout << "INFO, printParams(), repeat_length = " << repeat_length << "\n";
        std::cout << "INFO, printParams(), interval_length = " << interval_length << "\n";
        std::cout << "INFO, printParams(), read_length = " << read_length << "\n";
        std::cout << "INFO, printParams(), overlap_length = " << overlap_length << "\n";
        std::cout << "INFO, printParams(), flanking_length = " << flanking_length << "\n";
    }
};

#endif
