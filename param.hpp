#ifndef PARAM_CHOPPER_H
#define PARAM_CHOPPER_H

struct algoParams
{
    int reso;
/*     int cov_frac;
    int repeat_annotation_gap_thres; */
    int est_cov;
    double cov_mul;
    int repeat_length;
    int overlap_length;
    int uniform_read_length;
    int read_length_threshold;
    std::string outputfilename;

    void initParams()
    {
        reso = 50;
/*         cov_frac = 3;
        repeat_annotation_gap_thres = 200; */
        est_cov = 0;
        cov_mul = 1.5;
        repeat_length = 10000;
        overlap_length = 10000;
        uniform_read_length = overlap_length * 2;
        read_length_threshold = 25000;
        outputfilename = "chopper";
    }

    void printParams()
    {
        std::cout << "INFO, printParams(), reso = " << reso << "\n";
/*         std::cout << "INFO, printParams(), cov_frac = " << cov_frac << "\n";
        std::cout << "INFO, printParams(), repeat_annotation_gap_thres = " << repeat_annotation_gap_thres << "\n"; */
        std::cout << "INFO, printParams(), est_cov = " << est_cov << "\n";
        std::cout << "INFO, printParams(), cov_mul = " << cov_mul << "\n";
        std::cout << "INFO, printParams(), repeat_length = " << repeat_length << "\n";
        std::cout << "INFO, printParams(), overlap_length = " << overlap_length << "\n";
        std::cout << "INFO, printParams(), uniform_read_length = " << uniform_read_length << "\n";
        std::cout << "INFO, printParams(), read_length_threshold = " << read_length_threshold << "\n";
    }
};

#endif