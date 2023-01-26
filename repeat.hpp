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


void profileCoverage1(std::vector<Overlap *> &alignments, std::vector<std::pair<int, int>> &coverage, int reso, Read *read)
{
    // Returns coverage, which is a pair of ints <i*reso, coverage at position i*reso of read a>
    std::vector<std::pair<int, int>> events;
    for (int i = 0; i < alignments.size(); i++)
    {
        if (alignments[i]->read_A_id_ == read->id)
        {
            events.push_back(std::pair<int, int>(alignments[i]->read_A_match_start_ - 1, 1));
            events.push_back(std::pair<int, int>(alignments[i]->read_A_match_end_ - 1, -1));
        }
        else if (alignments[i]->read_B_id_ == read->id)
        {
            events.push_back(std::pair<int, int>(alignments[i]->read_B_match_start_ - 1, 1));
            events.push_back(std::pair<int, int>(alignments[i]->read_B_match_end_ - 1, -1));
        }
    }

    std::sort(events.begin(), events.end(), compare_event);

    int pos = 0;
    int i = 0;
    int count = 0;
    while (pos < events.size())
    {
        while ((events[pos].first < (i + 1) * reso) and (pos < events.size()))
        {
            count += events[pos].second;
            pos++;
        }
        coverage.push_back(std::pair<int, int>(i * reso, count));
        i++;
    }
    return;
}

void repeat_annotate1(std::vector<Read *> reads, std::vector<Overlap *> aln, const algoParams &param, std::vector<std::vector<Overlap *>> idx_pileup)
{
    int n_read = reads.size();
    std::vector<std::vector<std::pair<int, int>>> coverages;
    std::vector<std::vector<std::pair<int, int>>> cgs; // coverage gradient;
    std::vector<std::pair<int, int>> maskvec;
    std::vector<std::vector<std::pair<int, int>>> repeat_annotation;
    std::unordered_map<int, std::vector<std::tuple<int, int, int>>> repeats;

    std::ofstream cov(param.outputfilename + ".coverage1.txt");
    std::ofstream cov_grad(param.outputfilename + ".cov_grad1.txt");
    std::ofstream mask(param.outputfilename + ".mask1.txt");
    std::ofstream repeat_anno(param.outputfilename + ".repeat_anno1.txt");
    std::ofstream long_repeats(param.outputfilename + ".long_repeats1.txt");
    std::ofstream long_repeats_bed(param.outputfilename + ".long_repeats1.bed");

    fprintf(stdout, "INFO, length of alignments  %lu()\n", aln.size());

    for (int i = 0; i < n_read; i++)
    {
        repeat_annotation.push_back(std::vector<std::pair<int, int>>());
        coverages.push_back(std::vector<std::pair<int, int>>());
        cgs.push_back(std::vector<std::pair<int, int>>());
        maskvec.push_back(std::pair<int, int>());
    }

    fprintf(stdout, "INFO, profile coverage\n");

    for (int i = 0; i < n_read; i++)
    {
        std::vector<std::pair<int, int>> coverage;
        std::vector<std::pair<int, int>> cg;

        // profileCoverage: get the coverage based on pile-o-gram
        profileCoverage1(idx_pileup[i], coverage, param.reso, reads[i]);

        cov << "read " << i << " ";
        for (int j = 0; j < coverage.size(); j++)
            cov << coverage[j].first << "," << coverage[j].second << " ";
        cov << std::endl;

        cov_grad << "read " << i << " ";
        // Computes coverage gradients.
        if (coverage.size() >= 2)
        {
            for (int j = 1; j < coverage.size(); j++)
            {
                cg.push_back(std::pair<int, int>(coverage[j].first, coverage[j].second - coverage[j - 1].second));
                cov_grad << cg[j-1].first << "," << cg[j-1].second << " ";
            }
        }

        cov_grad << std::endl;

        coverages[i] = coverage;
        cgs[i] = cg;
    }

    fprintf(stdout, "INFO, profile coverage done\n");

    int total_slot = 0;
    long int total_cov = 0;

    std::vector<int> read_coverage;
    long int read_cov = 0;
    int read_slot = 0;
    // Finding the average coverage

    for (int i = 0; i < n_read; i++)
    {
        read_cov = 0;
        read_slot = 0;
        for (int j = 0; j < coverages[i].size(); j++)
        {
            read_cov += coverages[i][j].second;
            read_slot++;
        }
        total_cov += read_cov;
        total_slot += read_slot;
        int mean_read_cov = read_cov / std::max(1, read_slot);
        read_coverage.push_back(mean_read_cov);
    }

    size_t median_id = read_coverage.size() / 2;
    if (median_id > 0)
        std::nth_element(read_coverage.begin(), read_coverage.begin() + median_id, read_coverage.end());

    int cov_est = read_coverage[median_id];

    int mean_cov_est = total_cov / total_slot;

    fprintf(stdout, "INFO, Estimated mean coverage:  %d\n", mean_cov_est);
    fprintf(stdout, "INFO, Estimated median coverage:  %d\n", cov_est);

    int min_cov = cov_est / param.cov_frac;

    fprintf(stdout, "INFO, Estimated min coverage:  %d\n", min_cov);

    fprintf(stdout, "INFO, mask vector\n");

    for (int i = 0; i < n_read; i++)
    {

        // get the longest consecutive region that has decent coverage, decent coverage = estimated coverage / 3
        int start = 0;
        int end = start;
        int maxlen = 0, maxstart = 0, maxend = 0;
        for (int j = 0; j < coverages[i].size(); j++)
        {
            if (coverages[i][j].second > min_cov)
            {
                end = coverages[i][j].first + param.reso;
            }
            else
            {
                if (end > start)
                {
                    if (end - start > maxlen)
                    {
                        maxlen = end - start;
                        maxstart = start;
                        maxend = end;
                    }
                }
                start = coverages[i][j + 1].first;
                end = start;
            }
        }

        if (end > start)
        {
            if (end - start > maxlen)
            {
                maxlen = end - start;
                maxstart = start;
                maxend = end;
            }
        }

        maskvec[i] = (std::pair<int, int>(maxstart, maxend));
        mask << "read " << i << " " << maxstart << " " << maxend << std::endl;
    }

    fprintf(stdout, "INFO, mask vector done\n");

    fprintf(stdout, "INFO, repeats\n");

    // detect repeats based on coverage gradient, mark it has rising (1) or falling (-1)
    for (int i = 0; i < n_read; i++)
    {
        std::vector<std::pair<int, int>> anno;
        for (int j = 0; j < cgs[i].size(); j++)
        {

            if ((cgs[i][j].first >= maskvec[i].first) and (cgs[i][j].first <= maskvec[i].second))
            {
                if (cgs[i][j].second > std::min(
                                           std::max((coverages[i][j].second + min_cov) / param.cov_frac, cov_est / 2),
                                           cov_est * 2))
                {
                    anno.push_back(std::pair<int, int>(cgs[i][j].first, 1));
                }
                else if (cgs[i][j].second < -std::min(
                                                std::max((coverages[i][j].second + min_cov) / param.cov_frac, cov_est / 2),
                                                cov_est * 2))
                {
                    anno.push_back(std::pair<int, int>(cgs[i][j].first, -1));
                }
            }
        }
        repeat_annotation[i] = (anno);
    }

    // // clean it a bit, merge consecutive 1 or consecutive -1 if their position is within gap_threshold
    for (int i = 0; i < n_read; i++)
    {
        for (std::vector<std::pair<int, int>>::iterator iter = repeat_annotation[i].begin(); iter < repeat_annotation[i].end();)
        {
            if (iter + 1 < repeat_annotation[i].end())
            {
                if (((iter->second == 1) and ((iter + 1)->second == 1)) and
                    ((iter + 1)->first - iter->first < param.repeat_annotation_gap_thres))
                {
                    repeat_annotation[i].erase((iter + 1));
                }
                else if (((iter->second == -1) and ((iter + 1)->second == -1)) and
                         ((iter + 1)->first - iter->first < param.repeat_annotation_gap_thres))
                {
                    iter = repeat_annotation[i].erase(iter);
                }
                else
                    iter++;
            }
            else
                iter++;
        }
    }

    fprintf(stdout, "INFO, repeats done\n");

    int count_repeat_reads = 0;
    int count_repeat_annos = 0;

    for (int i = 0; i < n_read; i++)
    {
        repeat_anno << "read " << i << ", ";
        repeat_anno << "repeat annos " << repeat_annotation[i].size() << ":";
        for (auto &r : repeat_annotation[i])
        {
            count_repeat_annos++;
            repeat_anno << r.first << ",";
            repeat_anno << r.second << " ";
        }
        repeat_anno << std::endl;
        if (repeat_annotation[i].size() > 0)
            count_repeat_reads++;
    }

    fprintf(stdout, "INFO, Number of reads with repeats:  %d\n", count_repeat_reads);
    fprintf(stdout, "INFO, Number of repeat annos:  %d\n", count_repeat_annos);

    for (int i = 0; i < n_read; i++)
    {
        repeats[i] = std::vector<std::tuple<int, int, int>>();

        for (int j = 0; j < repeat_annotation[i].size(); j++)
        {
            if (repeat_annotation[i][j].second == -1) // repeat ending, negative gradient
            {
                bool bridged = true;
                int support = 0;
                int num_reads_at_end = 1;

                std::vector<int> read_other_ends;

                for (int k = 0; k < idx_pileup[i].size(); k++)
                {
                    int temp_id;
                    temp_id = idx_pileup[i][k]->read_B_id_;



                    if ((idx_pileup[i][k]->read_A_id_==i) and (idx_pileup[i][k]->read_A_match_end_ >
                         repeat_annotation[i][j].first - param.repeat_annotation_gap_thres / 2) and
                        (idx_pileup[i][k]->read_A_match_end_ <
                         repeat_annotation[i][j].first + param.repeat_annotation_gap_thres / 2))
                        {
                            read_other_ends.push_back(idx_pileup[i][k]->read_A_match_start_ - 1);
                            support++;
                        }
                    else if ((idx_pileup[i][k]->read_B_id_ == i) and (idx_pileup[i][k]->read_B_match_end_ > 
                            repeat_annotation[i][j].first - param.repeat_annotation_gap_thres / 2) and
                            (idx_pileup[i][k]->read_B_match_end_ <
                            repeat_annotation[i][j].first + param.repeat_annotation_gap_thres / 2))
                            {
                                read_other_ends.push_back(idx_pileup[i][k]->read_B_match_start_ - 1);
                                support++;
                            }
                }

                if (support < (cov_est / 2))
                {
                    continue;
                }

                int bridge_threshold = 3 * (support / 4);

                std::sort(read_other_ends.begin(), read_other_ends.end(), pairAscend);

                int num_reads_considered = 0;
                int num_reads_extending_to_end = 0;
                int num_reads_with_internal_overlaps = 0;
                int repeat_length = 1;

                for (int id = 0; id < read_other_ends.size(); ++id)
                {
                    if (read_other_ends[id] - maskvec[i].first < 200)
                    {
                        num_reads_considered++;
                        num_reads_extending_to_end++;

                        if ((num_reads_extending_to_end > bridge_threshold) or
                            ((num_reads_considered > bridge_threshold) and
                             (read_other_ends[id] - read_other_ends[0] > param.repeat_annotation_gap_thres)))
                        {
                            bridged = false;
                            repeat_length = repeat_annotation[i][j].first - maskvec[i].first;
                            break;
                        }
                    }
                    else
                    {
                        num_reads_with_internal_overlaps++;
                        num_reads_considered++;
                        int id1 = id + 1;
                        int pileup_length = 1;

                        while (id1 < read_other_ends.size())
                        {
                            if (read_other_ends[id1] - read_other_ends[id] < param.repeat_annotation_gap_thres)
                            {
                                pileup_length++;
                                id1++;
                            }
                            else
                            {
                                break;
                            }
                        }

                        if (pileup_length > bridge_threshold)
                        {
                            bridged = true;
                            repeat_length = repeat_annotation[i][j].first - read_other_ends[id];
                            break;
                        }
                    }
                }
                if (repeat_length > param.repeat_length)
                {
                    repeats[i].push_back(std::tuple<int, int, int>(repeat_annotation[i][j].first, -1, repeat_length));
                    reads[i]->preserve = 1;
                    if (reads[i]->align.compare("forward") == 0)
                    {
                        long_repeats_bed << reads[i]->chr << "\t" << reads[i]->start_pos + repeat_annotation[i][j].first - repeat_length << "\t" << reads[i]->start_pos + repeat_annotation[i][j].first << std::endl;
                    }
                    else if (reads[i]->align.compare("reverse") == 0)
                    {
                        long_repeats_bed << reads[i]->chr << "\t" << reads[i]->end_pos - repeat_annotation[i][j].first << "\t" << reads[i]->end_pos - repeat_annotation[i][j].first + repeat_length << std::endl;
                    }
                }
            }
            else
            {

                bool bridged = true;
                int support = 0;
                int num_reads_at_end = 1;

                std::vector<int> read_other_ends;

                for (int k = 0; k < idx_pileup[i].size(); k++)
                {

                    if ((idx_pileup[i][k]->read_A_id_ == i) and (idx_pileup[i][k]->read_A_match_start_ > repeat_annotation[i][j].first - param.repeat_annotation_gap_thres / 2) and
                        (idx_pileup[i][k]->read_A_match_start_ <
                         repeat_annotation[i][j].first + param.repeat_annotation_gap_thres / 2))
                    {
                        read_other_ends.push_back(idx_pileup[i][k]->read_A_match_end_ - 1);
                        support++;
                    }
                    else if ((idx_pileup[i][k]->read_B_id_ == i) and (idx_pileup[i][k]->read_B_match_start_ >  
                            repeat_annotation[i][j].first - param.repeat_annotation_gap_thres / 2) and
                             (idx_pileup[i][k]->read_B_match_start_ <
                              repeat_annotation[i][j].first + param.repeat_annotation_gap_thres / 2))
                    {
                        read_other_ends.push_back(idx_pileup[i][k]->read_B_match_end_ - 1);
                        support++;
                    }
                }

                if (support < (cov_est / 2))
                {
                    continue;
                }

                int bridge_threshold = 3 * (support / 4);

                std::sort(read_other_ends.begin(), read_other_ends.end(), pairDescend); // Sort in descending order

                int num_reads_considered = 0;
                int num_reads_extending_to_end = 0;
                int num_reads_with_internal_overlaps = 0;
                int repeat_length = 0;

                for (int id = 0; id < read_other_ends.size(); ++id)
                {
                    if (maskvec[i].second - read_other_ends[id] < param.repeat_annotation_gap_thres)
                    {
                        num_reads_considered++;
                        num_reads_extending_to_end++;

                        if ((num_reads_extending_to_end > bridge_threshold) or
                            ((num_reads_considered > bridge_threshold) and
                             (read_other_ends[0] - read_other_ends[id] > param.repeat_annotation_gap_thres)))
                        {
                            bridged = false;
                            repeat_length = maskvec[i].second - repeat_annotation[i][j].first;
                            break;
                        }
                    }
                    else
                    {
                        num_reads_with_internal_overlaps++;
                        num_reads_considered++;
                        int id1 = id + 1;
                        int pileup_length = 1;

                        while (id1 < read_other_ends.size())
                        {
                            if (read_other_ends[id] - read_other_ends[id1] < param.repeat_annotation_gap_thres)
                            {
                                pileup_length++;
                                id1++;
                            }
                            else
                            {
                                break;
                            }
                        }

                        if (pileup_length > bridge_threshold)
                        {
                            bridged = true;
                            repeat_length = read_other_ends[id] - repeat_annotation[i][j].first;
                            break;
                        }
                    }
                }

                if (repeat_length > param.repeat_length)
                {
                    repeats[i].push_back(std::tuple<int, int, int>(repeat_annotation[i][j].first, -1, repeat_length));
                    reads[i]->preserve = 1;
                    if (reads[i]->align.compare("forward") == 0)
                    {
                        long_repeats_bed << reads[i]->chr << "\t" << reads[i]->start_pos + repeat_annotation[i][j].first << "\t" << reads[i]->start_pos + repeat_annotation[i][j].first + repeat_length << std::endl;
                    }
                    else if (reads[i]->align.compare("reverse") == 0)
                    {
                        long_repeats_bed << reads[i]->chr << "\t" << reads[i]->end_pos - repeat_annotation[i][j].first - repeat_length << "\t" << reads[i]->end_pos - repeat_annotation[i][j].first << std::endl;
                    }
                }
            }
        }
    }

    int count_long_repeat_reads = 0;
    int count_long_repeats = 0;

    for (int i = 0; i < n_read; i++)
    {
        long_repeats << "read " << i << ", ";
        long_repeats << "long repeats " << repeats[i].size() << ":";
        for (auto &r : repeats[i])
        {
            if (std::get<2>(r) > param.repeat_length)
            {
                count_long_repeats++;
                long_repeats << std::get<0>(r) << ",";
                long_repeats << std::get<1>(r) << ",";
                long_repeats << std::get<2>(r) << " ";
            }
        }
        long_repeats << std::endl;
        if (repeats[i].size() > 0)
            count_long_repeat_reads++;
    }

    fprintf(stdout, "INFO, Number of reads with long repeats:  %d\n", count_long_repeat_reads);
    fprintf(stdout, "INFO, Number of long repeats:  %d\n", count_long_repeats);
}


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
        if (alignments[i]->read_B_id_ == read->id)
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

void repeat_annotate(std::vector<Read *> reads, std::vector<Overlap *> aln, const algoParams &param, std::vector<std::vector<Overlap *>> idx_pileup)
{

    int n_read = reads.size();
    std::vector<std::vector<std::pair<int, int>>> coverages;
    std::vector<std::tuple<int, int, int>> repeats;

    std::ofstream cov(param.outputfilename + ".coverage2.txt");
    std::ofstream repeat_reg(param.outputfilename + ".repeats2.txt");
    std::ofstream long_repeats(param.outputfilename + ".long_repeats2.txt");
    std::ofstream long_repeats_bed(param.outputfilename + ".long_repeats2.bed");


    for (int i = 0; i < n_read; i++)
    {
        coverages.push_back(std::vector<std::pair<int, int>>());
        repeats.push_back(std::tuple<int, int, int>());
    }

    fprintf(stdout, "INFO, profile coverage\n");

    for (int i = 0; i < n_read; i++)
    {
        std::vector<std::pair<int, int>> coverage;

        // profileCoverage: get the coverage based on pile-o-gram
        profileCoverage(idx_pileup[i], coverage, param.reso, reads[i]);

        cov << "read " << i << " ";
        for (int j = 0; j < coverage.size(); j++)
            cov << coverage[j].first << "," << coverage[j].second << " ";
        cov << std::endl;

        coverages[i] = coverage;
    }

    fprintf(stdout, "INFO, profile coverage done\n");

    int cov_est = param.est_cov;
    int high_cov = cov_est * param.cov_mul;
    int count_long_repeat_reads = 0;

    for (int i = 0; i < n_read; i++)
    {
        repeat_reg << "read " << i << ", ";

        // get the longest consecutive region that has high coverage, high coverage = estimated coverage * 1.5
        int start = 0;
        int end = start;
        int maxlen = 0, maxstart = 0, maxend = 0;
        for (int j = 0; j < coverages[i].size(); j++)
        {
            if (coverages[i][j].second >= high_cov)
            {
                end = coverages[i][j].first + param.reso - 1;
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
                start = coverages[i][j + 1].first;
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

        if (maxlen > param.repeat_length)
        {
            repeats[i] = std::tuple<int, int, int>(maxstart, maxend, maxlen);
            reads[i]->preserve = 1;
            long_repeats << "read " << i << ", ";
            long_repeats << maxstart << "," << maxend << "," << maxlen << std::endl;
            count_long_repeat_reads++;
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

    fprintf(stdout, "INFO, Number of reads with long repeats:  %d\n", count_long_repeat_reads);
}
