#ifndef PARAM_CHOPPER_H
#define PARAM_CHOPPER_H

struct algoParams
{
    int algo;
    int reso;
    int est_cov;
    int min_cov;
    int cov_frac;
    int repeat_annotation_gap_thres;
    int repeat_length;
    int read_length_threshold;
    std::string logFileName;
    std::string outputfilename;

    void initParams()
    {
        algo = 1;
        reso = 50;
        est_cov = 0;
        min_cov = 0;
        cov_frac = 3;
        repeat_annotation_gap_thres = 200;
        repeat_length = 10000;
        read_length_threshold = 25000;
        outputfilename = "prefix";
    }

    void printParams()
    {
        std::cout << "INFO, printParams(), algorithm = " << algo << "\n";
        std::cout << "INFO, printParams(), reso = " << reso << "\n";
        std::cout << "INFO, printParams(), est_cov = " << est_cov << "\n";
        std::cout << "INFO, printParams(), min_cov = " << min_cov << "\n";
        std::cout << "INFO, printParams(), cov_frac = " << cov_frac << "\n";
        std::cout << "INFO, printParams(), repeat_annotation_gap_thres = " << repeat_annotation_gap_thres << "\n";
        std::cout << "INFO, printParams(), repeat_length = " << repeat_length << "\n";
        std::cout << "INFO, printParams(), read_length_threshold = " << read_length_threshold << "\n";
    }
};

#endif