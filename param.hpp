#ifndef PARAM_CHOPPER_H
#define PARAM_CHOPPER_H

    struct algoParams
{
    int reso;
    int kmer_length;
    float kmer_frac;
    int repeat_length;
    int interval_length;
    int flanking_length;
    std::string outputfilename;
    int real_reads;
    int bloom_filter_element_count;
    int debug;

    void initParams()
    {
        reso = 500;
        kmer_length = 21;
        kmer_frac=0.8;
        repeat_length = 10000;
        interval_length = 10000;
        flanking_length = 2000;
        outputfilename = "chopper";
        real_reads = 1;
        bloom_filter_element_count = 1000;
        debug=0;
    }

    void printParams()
    {
        std::cout << "INFO, printParams(), reso = " << reso << "\n";
        std::cout << "INFO, printParams(), kmer_length = " << kmer_length << "\n";
        std::cout << "INFO, printParams(), kmer_frac = " << kmer_frac << "\n";
        std::cout << "INFO, printParams(), repeat_length = " << repeat_length << "\n";
        std::cout << "INFO, printParams(), interval_length = " << interval_length << "\n";
        std::cout << "INFO, printParams(), flanking_length = " << flanking_length << "\n";
        std::cout << "INFO, printParams(), bloom_filter_element_count = " << bloom_filter_element_count << "\n";
    }
};

#endif