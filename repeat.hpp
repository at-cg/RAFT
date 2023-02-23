#include <iostream>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <fstream>

#include "paf.hpp"
#include "read.hpp"
#include "overlap.hpp"
#include "param.hpp"

#ifndef COMPARE_EVENT
#define COMPARE_EVENT
bool compare_event(std::pair<int, int> event1, std::pair<int, int> event2)
{
    return event1.first < event2.first;
}
#endif

#ifndef PAIR_ASCEND
#define PAIR_ASCEND
bool pairAscend(int &firstElem, int &secondElem)
{
    return firstElem < secondElem;
}
#endif

#ifndef PAIR_DESCEND
#define PAIR_DESCEND
bool pairDescend(int &firstElem, int &secondElem)
{
    return firstElem > secondElem;
}
#endif


void profileCoverage(std::vector<Overlap *> &alignments, std::vector<std::pair<int, int>> &coverage, int reso, Read *read)
{

    int intervals = read->len / reso;

    if (intervals % reso)
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
            events.push_back(std::pair<int, int>(alignments[i]->read_A_match_start_ - 1, alignments[i]->read_A_match_end_ - 1));
        }
        else if (alignments[i]->read_B_id_ == read->id)
        {
            events.push_back(std::pair<int, int>(alignments[i]->read_B_match_start_ - 1,  alignments[i]->read_B_match_end_ - 1));
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

void repeat_annotate(std::vector<Read *> reads, const algoParams &param, std::vector<std::vector<Overlap *>> idx_pileup)
{
    int n_read = reads.size();

    std::ofstream cov(param.outputfilename + ".coverage.txt");
    std::ofstream repeat_reg(param.outputfilename + ".repeats.txt");
    std::ofstream long_repeats(param.outputfilename + ".long_repeats.txt");
    std::ofstream long_repeats_bed(param.outputfilename + ".long_repeats.bed");

    int cov_est = param.est_cov;
    int high_cov = cov_est * param.cov_mul;
    int count_long_repeat_reads = 0;

    long long total_coverage = 0;
    int total_windows = 0;

    for (int i = 0; i < n_read; i++)
    {
        repeat_reg << "read " << i << ", ";

        std::vector<std::pair<int, int>> coverage;

        // profileCoverage: get the coverage based on pile-o-gram
        profileCoverage(idx_pileup[i], coverage, param.reso, reads[i]);

        cov << "read " << i << " ";
        for (int j = 0; j < coverage.size(); j++)
            cov << coverage[j].first << "," << coverage[j].second << " ";
        cov << std::endl;

        // get the longest consecutive region that has high coverage, high coverage = estimated coverage * 1.5
        int start = 0;
        int end = start;
        int maxlen = 0, maxstart = 0, maxend = 0;
        for (int j = 0; j < coverage.size(); j++)
        {
            total_coverage = total_coverage + coverage[j].second;
            total_windows++;
            if (coverage[j].second >= high_cov)
            {
                end = coverage[j].first + param.reso - 1;
            }
            else
            {
                if (end > start)
                {
                    if ((end - start) > maxlen)
                    {
                        maxlen = end - start;
                        maxstart = start;
                        maxend = end;
                    }
                }
                start = coverage[j + 1].first;
                end = start;
            }
        }

        if (end > start)
        {
            if ((end - start) > maxlen)
            {
                maxlen = end - start;
                maxstart = start;
                maxend = end;
            }
        }

        repeat_reg << "longest_repeat_reg " << maxstart << "," << maxend << "," << maxlen << std::endl;

        if (maxlen >= param.repeat_length)
        {
            reads[i]->preserve = 1;
            long_repeats << "read " << i << ", ";
            long_repeats << maxstart << "," << maxend << "," << maxlen << std::endl;
            count_long_repeat_reads++;
            if (!param.real_reads)
            {
                if (reads[i]->align.compare("forward") == 0)
                {
                    long_repeats_bed << reads[i]->chr << "\t" << reads[i]->start_pos + maxstart << "\t" << reads[i]->start_pos + maxend << std::endl;
                }
                else if (reads[i]->align.compare("reverse") == 0)
                {
                    long_repeats_bed << reads[i]->chr << "\t" << reads[i]->end_pos - maxend << "\t" << reads[i]->end_pos - maxstart << std::endl;
                }
            }
        }
    }

    double coverage_per_window = (double)total_coverage / total_windows;
    fprintf(stdout, "coverage per window is %f \n", coverage_per_window);
    fprintf(stdout, "coverage per window/average coverage is %f \n", coverage_per_window / cov_est);

    fprintf(stdout, "INFO, Number of reads with long repeats:  %d\n", count_long_repeat_reads);
}
