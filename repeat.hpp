#include <iostream>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <fstream>

#include "read.hpp"
#include "param.hpp"
#include "overlap.hpp"
#include "bloom_filter.hpp"

#ifndef COMPARE_EVENT
#define COMPARE_EVENT
bool compare_event(std::pair<int, int> event1, std::pair<int, int> event2)
{
    return event1.first < event2.first;
}
#endif

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

void profileCoverage1(std::vector<Overlap *> &alignments, std::vector<std::pair<int, int>> &coverage, Read *read,
    const algoParams &param)
{
    int reso = param.reso;
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

    // Returns coverage, which is a pair of ints <i*reso, coverage at position i*reso of read a>
    std::vector<std::pair<int, int>> events;
    for (int i = 0; i < alignments.size(); i++)
    {
        if (alignments[i]->read_A_id_ == read->id)
        {
            events.push_back(std::pair<int, int>(alignments[i]->read_A_match_start_ , alignments[i]->read_A_match_end_ - 1));
        }
        else if (!param.symmetric_overlaps && alignments[i]->read_B_id_ == read->id)
        {
            events.push_back(std::pair<int, int>(alignments[i]->read_B_match_start_ ,  alignments[i]->read_B_match_end_ - 1));
        }
    }

    std::sort(events.begin(), events.end(), compare_event);

    int pos = 0;
    int i = 0;
    while (pos < events.size())
    {
        while ((events[pos].first < (i + 1) * reso) and (pos < events.size()))
        {
            int k = i;
            while (events[pos].second >= k * reso)
            {
                coverage[k].second++;
                k++;
            }
            pos++;
        }
        i++;
    }
    return;
}

void profileCoverage2(bloom_filter *kmer_filter, std::vector<std::pair<int, float>> &coverage, Read *read,
                     const algoParams &param)
{
    int reso = param.reso*10;
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

    int j = 0, z = 0;
    int k = param.kmer_length;
    int candidiate_kmers, high_freq_kmers = 0;
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1, kmer[2] = {0, 0};

    for (int i = 0; i < read->len; ++i)
    {
        int c = seq_nt4_table[(uint8_t)read->bases[i]];
        if (i >= ((j + 1) * reso))
        {
            coverage[j].second = (float)high_freq_kmers / candidiate_kmers;
            j++;
            high_freq_kmers = 0;
            candidiate_kmers = 0;
        }
        kmer[0] = (kmer[0] << 2 | c) & mask;             // forward k-mer
        kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1; // reverse k-mer
        z = kmer[0] < kmer[1] ? 0 : 1;                   // strand
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

void repeat_annotate1(std::vector<Read *> reads, std::vector<std::vector<Overlap *>> idx_pileup, const algoParams &param)
{

    int n_read = reads.size();
    std::ofstream cov(param.outputfilename + ".coverage1.txt");
    std::ofstream long_repeats(param.outputfilename + ".long_repeats1.txt");
    std::ofstream long_repeats_bed(param.outputfilename + ".long_repeats1.bed");
 
    int cov_est = param.est_cov;
    int high_cov = cov_est * param.cov_mul;

    long long total_coverage = 0;
    int total_windows = 0;
    long long total_repeat_length =0;
    long long total_read_length=0;

    for (int i = 0; i < n_read; i++)
    {
        total_read_length = total_read_length + reads[i]->len;
        std::vector<std::pair<int, int>> coverage;
        profileCoverage1(idx_pileup[i], coverage, reads[i], param);

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
            total_coverage = total_coverage + coverage[j].second;
            total_windows++;

            if (coverage[j].second >= high_cov)
            {
                end = coverage[j].first + param.reso;
            }
            else
            {
                if ((end - start) >= param.repeat_length)
                {
                    total_repeat_length = total_repeat_length + end - start;

                    s = start - param.flanking_length;
                    e = end + param.flanking_length;

                    if (s <= 0)
                    {
                        s = 0;
                    }

                    if ( e >= reads[i]->len)
                    {
                        e = reads[i]->len;
                    }

                    if (reads[i]->long_repeats1.size() && (s <= reads[i]->long_repeats1.back().second))
                    {

                        reads[i]->long_repeats1.back().second = e;
                    }
                    else
                    {
                        reads[i]->long_repeats1.push_back(std::pair<int, int>(s, e));
                    }
                }

                start = coverage[j + 1].first;
                end = start;
            }
        }

        if ((end - start) >= param.repeat_length)
        {
            total_repeat_length=total_repeat_length + end - start;

            s = start - param.flanking_length;
            e = end + param.flanking_length;

            if (s <= 0)
            {
                s = 0;
            }

            if (e >= reads[i]->len)
            {
                e = reads[i]->len;
            }

            if (reads[i]->long_repeats1.size() && (s <= reads[i]->long_repeats1.back().second))
            {

                reads[i]->long_repeats1.back().second = e;
            }
            else
            {
                reads[i]->long_repeats1.push_back(std::pair<int, int>(s, e));
            }
        }
    }

    double coverage_per_window = (double)total_coverage / total_windows;
    double fraction_of_repeat_length = (double)total_repeat_length / total_read_length;

    fprintf(stdout, "coverage per window is %f \n", coverage_per_window);
    fprintf(stdout, "coverage per window/average coverage is %f \n", coverage_per_window/cov_est);
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

void repeat_annotate2(std::vector<Read *> reads, bloom_filter *kmer_filter, const algoParams &param)
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
        profileCoverage2(kmer_filter, coverage, reads[i], param);

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
                end = coverage[j].first + param.reso;
            }
            else
            {
                if ((end - start) >= param.repeat_length)
                {
                    total_repeat_length = total_repeat_length + end - start;

                    s = start - param.flanking_length;
                    e = end + param.flanking_length;

                    if (s <= 0)
                    {
                        s = 0;
                    }

                    if (e >= reads[i]->len)
                    {
                        e = reads[i]->len;
                    }

                    if (reads[i]->long_repeats2.size() && (s <= reads[i]->long_repeats2.back().second))
                    {

                        reads[i]->long_repeats2.back().second = e;
                    }
                    else
                    {
                        reads[i]->long_repeats2.push_back(std::pair<int, int>(s, e));
                    }
                }

                start = coverage[j + 1].first;
                end = start;
            }
        }

        if ((end - start) >= param.repeat_length)
        {
            total_repeat_length = total_repeat_length + end - start;

            s = start - param.flanking_length;
            e = end + param.flanking_length;

            if (s <= 0)
            {
                s = 0;
            }

            if (e >= reads[i]->len)
            {
                e = reads[i]->len;
            }

            if (reads[i]->long_repeats2.size() && (s <= reads[i]->long_repeats2.back().second))
            {

                reads[i]->long_repeats2.back().second = e;
            }
            else
            {
                reads[i]->long_repeats2.push_back(std::pair<int, int>(s, e));
            }
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
