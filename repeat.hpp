#include <iostream>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <fstream>

#include "read.hpp"
#include "param.hpp"
#include "bloom_filter.hpp"

#ifndef ATCG
#define ATCG
unsigned char seq_nt4_table[256] = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

#endif

#ifndef COMPARE_START_REPEATS
#define COMPARE_START_REPEATS
bool compare_start_repeat(std::pair<int, int> repeat1, std::pair<int, int> repeat2)
{
    return std::get<0>(repeat1) < std::get<0>(repeat2);
}
#endif

void profileCoverage(bloom_filter *kmer_filter, std::vector<std::pair<int, float>> &coverage, Read *read,
                     const algoParams &param, int reso)
{

    int intervals = read->len / reso;

    if (read->len % reso)
    {
        intervals++;
    }

    for (int i = 0; i < intervals; i++)
    {
        coverage.push_back(std::pair<int, int>());
        coverage[i].first = i * reso;
    }

    int j=0, z=0;
    int k = param.kmer_length;
    int candidiate_kmers=0, high_freq_kmers = 0;
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1, kmer[2] = {0, 0};

    for (int i = 0; i < read->len; ++i)
    {
        int c = seq_nt4_table[(uint8_t)read->bases[i]];
        if (i >= ((j + 1) * reso)){
            coverage[j].second = (float)high_freq_kmers/candidiate_kmers;
            j++;
            high_freq_kmers=0;
            candidiate_kmers=0;
        }
        kmer[0] = (kmer[0] << 2 | c) & mask;         // forward k-mer
        kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1; // reverse k-mer
        z = kmer[0] < kmer[1] ? 0 : 1; // strand
        if (i >= k)
        {
            candidiate_kmers++;
             if (kmer_filter->contains(kmer[z]))
            {
                high_freq_kmers++;
            }
        }
    }

    coverage[j].second = (float)high_freq_kmers / candidiate_kmers;

    return;
}

void repeat_annotate1(std::vector<Read *> &reads, bloom_filter *kmer_filter, const algoParams &param)
{

    int n_read = reads.size();
    std::ofstream cov(param.outputfilename + ".coverage1.txt");
    std::ofstream long_repeats(param.outputfilename + ".long_repeats1.txt");
    std::ofstream long_repeats_bed(param.outputfilename + ".long_repeats1.bed");

    long long total_repeat_length = 0;
    long long total_read_length = 0;

    for (int i = 0; i < n_read; i++)
    {
        total_read_length = total_read_length + reads[i]->len;
        std::vector<std::pair<int, float>> coverage;
        profileCoverage(kmer_filter, coverage, reads[i], param, param.reso1);

        cov << "read " << i << " ";
        for (int j = 0; j < coverage.size(); j++)
            cov << coverage[j].first << "," << coverage[j].second << " ";
        cov << std::endl;

        // get the longest consecutive region that has high coverage, high coverage = estimated coverage * 1.5
        int start = 0;
        int end = start;
        int s, e = 0;
        for (int j = 0; j < coverage.size(); j++)
        {
            if (coverage[j].second >= param.kmer_frac)
            {
                end = coverage[j].first + param.reso1;
            }
            else
            {
                if ((end - start) >= param.repeat_length)
                {
                    int flanking_length = std::min(param.flanking_frac * (end - start), float(param.flanking_length));
                    total_repeat_length = total_repeat_length + end - start;

                    s = start - flanking_length;
                    e = end + flanking_length;

                    if (s <= 0)
                    {
                        s = 0;
                    }

                    if ( e >= reads[i]->len)
                    {
                        e = reads[i]->len;
                    }

                    reads[i]->long_repeats1.push_back(std::pair<int, int>(s, e));
                    
                }

                start = coverage[j + 1].first;
                end = start;
            }
        }

        if ((end - start) >= param.repeat_length)
        {
            int flanking_length = std::min(param.flanking_frac * (end - start), float(param.flanking_length));
            total_repeat_length = total_repeat_length + end - start;

            s = start - flanking_length;
            e = end + flanking_length;

            if (s <= 0)
            {
                s = 0;
            }

            if (e >= reads[i]->len)
            {
                e = reads[i]->len;
            }

            reads[i]->long_repeats1.push_back(std::pair<int, int>(s, e));
            
        }
    }

    double fraction_of_repeat_length = (double)total_repeat_length / total_read_length;

    fprintf(stdout, "fraction_of_repeat_length %f \n", fraction_of_repeat_length);

    for (int i = 0; i < n_read; i++)
    {
        long_repeats << "read " << i << ", ";
        for (int j = 0; j < reads[i]->long_repeats1.size(); j++)
        {
            long_repeats << reads[i]->long_repeats1[j].first << "," << reads[i]->long_repeats1[j].second << "    ";
            
            if(!param.real_reads){
                if (reads[i]->align.compare("forward") == 0)
                {
                    long_repeats_bed << reads[i]->chr << "\t" << reads[i]->start_pos + reads[i]->long_repeats1[j].first
                                    << "\t" << reads[i]->start_pos + reads[i]->long_repeats1[j].second << std::endl;
                }
                else if (reads[i]->align.compare("reverse") == 0)
                {
                    long_repeats_bed << reads[i]->chr << "\t" << reads[i]->end_pos - reads[i]->long_repeats1[j].second
                                     << "\t" << reads[i]->end_pos - reads[i]->long_repeats1[j].first << std::endl;
                }
            }
        }

        long_repeats << std::endl;
    }
}

void repeat_annotate2(std::vector<Read *> &reads, bloom_filter *kmer_filter, const algoParams &param)
{

    int n_read = reads.size();
    std::ofstream cov(param.outputfilename + ".coverage2.txt");
    std::ofstream long_repeats(param.outputfilename + ".long_repeats2.txt");
    std::ofstream long_repeats_bed(param.outputfilename + ".long_repeats2.bed");

    long long total_repeat_length = 0;
    long long total_read_length = 0;

    for (int i = 0; i < n_read; i++)
    {
        total_read_length = total_read_length + reads[i]->len;
        std::vector<std::pair<int, float>> coverage;
        profileCoverage(kmer_filter, coverage, reads[i], param, param.reso2);

        cov << "read " << i << " ";
        for (int j = 0; j < coverage.size(); j++)
            cov << coverage[j].first << "," << coverage[j].second << " ";
        cov << std::endl;

        // get the longest consecutive region that has high coverage, high coverage = estimated coverage * 1.5
        int start = 0;
        int end = start;
        int s, e = 0;
        for (int j = 0; j < coverage.size(); j++)
        {
            if (coverage[j].second >= param.kmer_frac)
            {
                end = coverage[j].first + param.reso2;
            }
            else
            {
                if ((end - start) >= param.repeat_length)
                {
                    int flanking_length = std::min(param.flanking_frac * (end - start), float(param.flanking_length));
                    total_repeat_length = total_repeat_length + end - start;

                    s = start - flanking_length;
                    e = end + flanking_length;

                    if (s <= 0)
                    {
                        s = 0;
                    }

                    if (e >= reads[i]->len)
                    {
                        e = reads[i]->len;
                    }

                    reads[i]->long_repeats2.push_back(std::pair<int, int>(s, e));
                    
                }

                start = coverage[j + 1].first;
                end = start;
            }
        }

        if ((end - start) >= param.repeat_length)
        {
            int flanking_length = std::min(param.flanking_frac * (end - start), float(param.flanking_length));
            total_repeat_length = total_repeat_length + end - start;

            s = start - flanking_length;
            e = end + flanking_length;

            if (s <= 0)
            {
                s = 0;
            }

            if (e >= reads[i]->len)
            {
                e = reads[i]->len;
            }

            reads[i]->long_repeats2.push_back(std::pair<int, int>(s, e));
            
        }
    }

    double fraction_of_repeat_length = (double)total_repeat_length / total_read_length;

    fprintf(stdout, "fraction_of_repeat_length %f \n", fraction_of_repeat_length);

    for (int i = 0; i < n_read; i++)
    {
        long_repeats << "read " << i << ", ";
        for (int j = 0; j < reads[i]->long_repeats2.size(); j++)
        {
            long_repeats << reads[i]->long_repeats2[j].first << "," << reads[i]->long_repeats2[j].second << "    ";

            if (!param.real_reads)
            {
                if (reads[i]->align.compare("forward") == 0)
                {
                    long_repeats_bed << reads[i]->chr << "\t" << reads[i]->start_pos + reads[i]->long_repeats2[j].first
                                     << "\t" << reads[i]->start_pos + reads[i]->long_repeats2[j].second << std::endl;
                }
                else if (reads[i]->align.compare("reverse") == 0)
                {
                    long_repeats_bed << reads[i]->chr << "\t" << reads[i]->end_pos - reads[i]->long_repeats2[j].second
                                     << "\t" << reads[i]->end_pos - reads[i]->long_repeats2[j].first << std::endl;
                }
            }
        }

        long_repeats << std::endl;
    }
}