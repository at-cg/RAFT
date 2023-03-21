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

void profileCoverage(std::vector<Overlap *> &alignments, std::vector<std::pair<int, int>> &coverage, Read *read,
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

void repeat_annotate(std::vector<Read *> reads, std::vector<std::vector<Overlap *>> idx_pileup, const algoParams &param)
{

    int n_read = reads.size();
    std::ofstream cov(param.outputfilename + ".coverage.txt");
    std::ofstream long_repeats(param.outputfilename + ".long_repeats.txt");
    std::ofstream long_repeats_bed(param.outputfilename + ".long_repeats.bed");
 
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
        profileCoverage(idx_pileup[i], coverage, reads[i], param);

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

                    if (reads[i]->long_repeats.size() && (s <= reads[i]->long_repeats.back().second))
                    {

                        reads[i]->long_repeats.back().second = e;
                    }
                    else
                    {
                        reads[i]->long_repeats.push_back(std::pair<int, int>(s, e));
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

            if (reads[i]->long_repeats.size() && (s <= reads[i]->long_repeats.back().second))
            {

                reads[i]->long_repeats.back().second = e;
            }
            else
            {
                reads[i]->long_repeats.push_back(std::pair<int, int>(s, e));
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
        for (int j = 0; j < reads[i]->long_repeats.size(); j++)
        {
            long_repeats << reads[i]->long_repeats[j].first << "," << reads[i]->long_repeats[j].second << "    ";
            
            if(!param.real_reads){
                if (reads[i]->align.compare("forward") == 0)
                {
                    long_repeats_bed << reads[i]->chr << "\t" << reads[i]->start_pos + reads[i]->long_repeats[j].first
                                    << "\t" << reads[i]->start_pos + reads[i]->long_repeats[j].second << std::endl;
                }
                else if (reads[i]->align.compare("reverse") == 0)
                {
                    long_repeats_bed << reads[i]->chr << "\t" << reads[i]->end_pos - reads[i]->long_repeats[j].second
                                     << "\t" << reads[i]->end_pos - reads[i]->long_repeats[j].first << std::endl;
                }
            }
        }

        long_repeats << std::endl;
    }

}
