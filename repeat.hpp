#include <iostream>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <fstream>

#include "paf.hpp"
#include "read.hpp"
#include "overlap.hpp"
#include "param.hpp"

#ifndef COMPARE_OVERLAP
#define COMPARE_OVERLAP
bool compare_overlap(Overlap *ovl1, Overlap *ovl2)
{
    // Returns True if the sum of the match lengths of the two reads in ovl1 > the sum of the  overlap lengths of the two reads in ovl2
    // Returns False otherwise.
    return ((ovl1->read_A_match_end_ - ovl1->read_A_match_start_ + ovl1->read_B_match_end_ - ovl1->read_B_match_start_) > (ovl2->read_A_match_end_ - ovl2->read_A_match_start_ + ovl2->read_B_match_end_ - ovl2->read_B_match_start_));
}
#endif

#ifndef COMPARE_EVENT
#define COMPARE_EVENT
bool compare_event(std::pair<int, int> event1, std::pair<int, int> event2)
{
    return event1.first < event2.first;
}
#endif

#ifndef PAIR_ASCEND
#define PAIR_ASCEND
bool pairAscend(const std::pair<int, int> &firstElem, const std::pair<int, int> &secondElem)
{
    return firstElem.first < secondElem.first;
}
#endif

#ifndef PAIR_DESCEND
#define PAIR_DESCEND
bool pairDescend(const std::pair<int, int> &firstElem, const std::pair<int, int> &secondElem)
{
    return firstElem.first > secondElem.first;
}
#endif

void profileCoverage(std::vector<Overlap *> &alignments, std::vector<std::pair<int, int>> &coverage, int reso)
{
    // Returns coverage, which is a pair of ints <i*reso, coverage at position i*reso of read a>
    std::vector<std::pair<int, int>> events;
    for (int i = 0; i < alignments.size(); i++)
    {
        events.push_back(std::pair<int, int>(alignments[i]->read_A_match_start_ , 1));
        events.push_back(std::pair<int, int>(alignments[i]->read_A_match_end_ , -1));
    }

    std::sort(events.begin(), events.end(), compare_event);

    int pos = 0;
    int i = 0;
    int count = 0;
    while (pos < events.size())
    {
        while ((events[pos].first <= i * reso) and (pos < events.size()))
        {
            count += events[pos].second;
            pos++;
        }
        coverage.push_back(std::pair<int, int>(i * reso, count));
        i++;
    }
    return;
}

void repeat_annotate(std::vector<Read *> reads, std::vector<Overlap *> aln, const algoParams &param)
{
    int n_read = reads.size();
    std::vector<std::vector<std::pair<int, int>>> coverages(n_read);
    std::vector<std::vector<std::pair<int, int>>> cgs(n_read); // coverage gradient;
    std::vector<std::pair<int, int>> maskvec;
    std::vector<std::vector<std::pair<int, int>>> repeat_annotation;
    std::unordered_map<int, std::vector<std::tuple<int, int, int>>> repeats;

    std::ofstream cov(param.outputfilename + ".coverage.txt");
    std::ofstream mask(param.outputfilename + ".mas");
    std::ofstream repeat_anno(param.outputfilename + ".repeat_anno.txt");
    std::ofstream bridged_repeats(param.outputfilename + ".repeats.txt");

    std::vector<std::vector<Overlap *>> idx_pileup; // this is the pileup
    std::unordered_map<int, std::vector<std::pair<int, int>>> self_aln_list;

    int r_begin = aln.front()->read_A_id_;
    int r_end = aln.back()->read_A_id_;


    fprintf(stdout, "INFO, begin %d(), end %d()\n",r_begin, r_end);
    fprintf(stdout, "INFO, length of alignments  %lu()\n", aln.size());

    for (int i = 0; i < n_read; i++)
    {
        idx_pileup.push_back(std::vector<Overlap *>());
        repeat_annotation.push_back(std::vector<std::pair<int, int>>());
        coverages.push_back(std::vector<std::pair<int, int>>());
        cgs.push_back(std::vector<std::pair<int, int>>());
        maskvec.push_back(std::pair<int, int>());
    }

    for (int i = 0; i < aln.size(); i++)
    {
        if (aln[i]->read_A_id_ == aln[i]->read_B_id_)
        {
            aln[i]->active = false;
            if (self_aln_list.find(aln[i]->read_A_id_) == self_aln_list.end())
                self_aln_list[aln[i]->read_A_id_] = std::vector<std::pair<int, int>>();

            self_aln_list[aln[i]->read_A_id_].push_back(std::pair<int, int>(aln[i]->read_A_match_start_, aln[i]->read_A_match_end_));
            self_aln_list[aln[i]->read_A_id_].push_back(std::pair<int, int>(aln[i]->read_B_match_start_, aln[i]->read_B_match_end_));
        }
        if (aln[i]->active)
        {
            idx_pileup[aln[i]->read_A_id_].push_back(aln[i]);
        }
    }

    for (int i = 0; i < n_read; i++)
    { // sort overlaps of a reads
        std::sort(idx_pileup[i].begin(), idx_pileup[i].end(), compare_overlap);
    }

    fprintf(stdout, "INFO, profile coverage\n");

    for (int i = r_begin; i <= r_end; i++)
    {
        std::vector<std::pair<int, int>> coverage;
        std::vector<std::pair<int, int>> cg;

        // profileCoverage: get the coverage based on pile-o-gram
        profileCoverage(idx_pileup[i], coverage, param.reso);

        cov << "read " << i << " ";
        for (int j = 0; j < coverage.size(); j++)
            cov << coverage[j].first << "," << coverage[j].second << " ";
        cov << std::endl;

        // Computes coverage gradients.
        if (coverage.size() >= 2)
            for (int j = 0; j < coverage.size() - 1; j++)
            {
                cg.push_back(std::pair<int, int>(coverage[j].first, coverage[j + 1].second - coverage[j].second));
            }
        else
            cg.push_back(std::pair<int, int>(0, 0));

        coverages[i]=coverage;
        cgs[i]=cg;
    }

    fprintf(stdout, "INFO, profile coverage done\n");

    int total_slot = 0;
    long int total_cov = 0;

    std::vector<int> read_coverage;
    long int read_cov = 0;
    int read_slot = 0;
    // Finding the average coverage

    for (int i = r_begin; i <= r_end; i++)
    {
        if (reads[i]->len < 5000)
            continue;
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

    if (param.est_cov != 0)
        cov_est = param.est_cov;

    fprintf(stdout, "INFO, Estimated mean coverage:  %d\n", mean_cov_est);
    fprintf(stdout, "INFO, Estimated median coverage:  %d\n", cov_est);

    int min_cov = cov_est / param.cov_frac;

    if (param.min_cov > min_cov)
        min_cov = param.min_cov;

    fprintf(stdout, "INFO, mask vector\n");



    for (int i = r_begin; i <= r_end; i++)
    {

        // get the longest consecutive region that has decent coverage, decent coverage = estimated coverage / 3
        int start = 0;
        int end = start;
        int maxlen = 0, maxstart = 0, maxend = 0;
        int start_pos = 0;
        int end_pos = start_pos;
        int max_start_pos = 0, max_end_pos = 0;
        for (int j = 0; j < coverages[i].size(); j++)
        {
            if (coverages[i][j].second > min_cov)
            {
                end = coverages[i][j].first;
                end_pos = j;
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
                        max_start_pos = start_pos;
                        max_end_pos = end_pos;
                    }
                }
                start = coverages[i][j].first + 1;
                start_pos = j + 1;
                end_pos = start_pos;
                end = start;
            }
        }

        maskvec[i] = (std::pair<int, int>(maxstart, maxend));
        mask << i << " " << maxstart << " " << maxend << std::endl;
    }

    fprintf(stdout, "INFO, mask vector done\n");

    fprintf(stdout, "INFO, repeats\n");

    // detect repeats based on coverage gradient, mark it has rising (1) or falling (-1)
    for (int i = r_begin; i <= r_end; i++)
    {
        std::vector<std::pair<int, int>> anno;
        for (int j = 0; j < cgs[i].size() - 1; j++)
        {

            if ((cgs[i][j].first >= maskvec[i].first) and (cgs[i][j].first <= maskvec[i].second))
            {
                if (cgs[i][j].second > std::min(
                                           std::max((coverages[i][j].second + min_cov) / param.cov_frac, cov_est),
                                           cov_est * 2))
                {
                    anno.push_back(std::pair<int, int>(cgs[i][j].first, 1));
                }
                else if (cgs[i][j].second < -std::min(
                                                std::max((coverages[i][j].second + min_cov) / param.cov_frac, cov_est),
                         cov_est * 2))
                {
                    anno.push_back(std::pair<int, int>(cgs[i][j].first, -1));
                }
            }
        }
        repeat_annotation[i] = (anno);
    }

    int count_repeat_reads = 0;

    for (int i = r_begin; i <= r_end; i++)
    {

        if (repeat_annotation[i].size() > 1)
            count_repeat_reads++;
    }

    fprintf(stdout, "INFO, Number of reads with 2 repeat annotations before merging:  %d\n", count_repeat_reads);

    // // clean it a bit, merge consecutive 1 or consecutive -1 if their position is within gap_threshold
    for (int i = r_begin; i <= r_end; i++)
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
    

    count_repeat_reads = 0;

    for (int i = r_begin; i <= r_end; i++){
        repeat_anno << "read " << i << " ";
        repeat_anno << repeat_annotation[i].size();
        repeat_anno << std::endl;
        if (repeat_annotation[i].size() > 1)
            count_repeat_reads++;
    }

    fprintf(stdout, "INFO, Number of reads with 2 repeat annotations:  %d\n", count_repeat_reads);

    for (int i = r_begin; i <= r_end; i++)
    {
        repeats[i] = std::vector<std::tuple<int, int, int>>();

        for (int j = 0; j < repeat_annotation[i].size(); j++)
        {
            if (repeat_annotation[i][j].second == -1) // repeat ending, negative gradient
            {
                bool bridged = true;
                int support = 0;
                int num_reads_at_end = 1;

                std::vector<std::pair<int, int>> read_other_ends;

                for (int k = 0; k < idx_pileup[i].size(); k++)
                {
                    int left_overhang, right_overhang;
                    int temp_id;
                    temp_id = idx_pileup[i][k]->read_B_id_;

                    if (idx_pileup[i][k]->reverse_complement_match_ == 0)
                    {
                        right_overhang = std::max(maskvec[temp_id].second - idx_pileup[i][k]->read_B_match_end_, 0);
                        left_overhang = std::max(idx_pileup[i][k]->read_B_match_start_ - maskvec[temp_id].first, 0);
                    }
                    else if (idx_pileup[i][k]->reverse_complement_match_ == 1)
                    {
                        right_overhang = std::max(idx_pileup[i][k]->read_B_match_start_ - maskvec[temp_id].first, 0);
                        left_overhang = std::max(maskvec[temp_id].second - idx_pileup[i][k]->read_B_match_end_, 0);
                    }

                    if (right_overhang > 0)
                    {
                        if ((idx_pileup[i][k]->read_A_match_end_ >
                             repeat_annotation[i][j].first - param.repeat_annotation_gap_thres / 2) and
                            (idx_pileup[i][k]->read_A_match_end_ <
                             repeat_annotation[i][j].first + param.repeat_annotation_gap_thres / 2))
                        {

                            std::pair<int, int> other_end;
                            other_end.first = idx_pileup[i][k]->read_A_match_start_;
                            other_end.second = left_overhang;
                            read_other_ends.push_back(other_end);
                            support++;
                        }
                    }
                }

                if (support < (min_cov * 2))
                {
                    continue;
                }

                int bridge_threshold = 3 * (support/4);

                std::sort(read_other_ends.begin(), read_other_ends.end(), pairAscend);

                int num_reads_considered = 0;
                int num_reads_extending_to_end = 0;
                int num_reads_with_internal_overlaps = 0;
                int repeat_length = 1;

                for (int id = 0; id < read_other_ends.size(); ++id)
                {
                    if (read_other_ends[id].first - maskvec[i].first < 200)
                    {
                        num_reads_considered++;
                        num_reads_extending_to_end++;

                        if ((num_reads_extending_to_end > bridge_threshold) or
                            ((num_reads_considered > bridge_threshold) and
                                (read_other_ends[id].first - read_other_ends[0].first > param.repeat_annotation_gap_thres)))
                        {
                            bridged = false;
                            break;
                        }
                    }
                    else if (read_other_ends[id].second < 0)
                    {
                        num_reads_considered++;

                        if ((num_reads_extending_to_end > bridge_threshold) or
                            ((num_reads_considered > bridge_threshold) and
                             (read_other_ends[id].first - read_other_ends[0].first > param.repeat_annotation_gap_thres)))
                        {
                            bridged = false;
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
                            if (read_other_ends[id1].first - read_other_ends[id].first < param.repeat_annotation_gap_thres)
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
                            int median_id = (id1-id)/2;
                            repeat_length = repeat_annotation[i][j].first - read_other_ends[median_id].first;
                            break;
                        }
                    }
                }
                if (bridged){
                    repeats[i].push_back(std::tuple<int, int, int>(repeat_annotation[i][j].first, -1, repeat_length));
                    if(repeat_length > param.repeat_length){
                        reads[i]->preserve = 1;
                    }
                }

            }
            else{

                bool bridged = true;
                int support = 0;
                int num_reads_at_end = 1;

                std::vector<std::pair<int, int>> read_other_ends;

                for (int k = 0; k < idx_pileup[i].size(); k++)
                {
                    int left_overhang, right_overhang;
                    int temp_id;
                    temp_id = idx_pileup[i][k]->read_B_id_;

                    if (idx_pileup[i][k]->reverse_complement_match_ == 0)
                    {
                        right_overhang = std::max(maskvec[temp_id].second - idx_pileup[i][k]->read_B_match_end_, 0);
                        left_overhang = std::max(idx_pileup[i][k]->read_B_match_start_ - maskvec[temp_id].first, 0);
                    }
                    else if (idx_pileup[i][k]->reverse_complement_match_ == 1)
                    {
                        right_overhang = std::max(idx_pileup[i][k]->read_B_match_start_ - maskvec[temp_id].first, 0);
                        left_overhang = std::max(maskvec[temp_id].second - idx_pileup[i][k]->read_B_match_end_, 0);
                    }

                    if (left_overhang > 0)
                    {
                        if ((idx_pileup[i][k]->read_A_match_start_ >
                             repeat_annotation[i][j].first - param.repeat_annotation_gap_thres / 2) and
                            (idx_pileup[i][k]->read_A_match_start_ <
                             repeat_annotation[i][j].first + param.repeat_annotation_gap_thres / 2))
                        {

                            std::pair<int, int> other_end;
                            other_end.first = idx_pileup[i][k]->read_A_match_end_;
                            other_end.second = right_overhang;
                            read_other_ends.push_back(other_end);
                            support++;
                        }
                    }
                }
                if (support < (min_cov * 2))
                {
                    continue;
                }

                int bridge_threshold = 3 * (support/4);

                    std::sort(read_other_ends.begin(), read_other_ends.end(), pairDescend); // Sort in descending order

                int num_reads_considered = 0;
                int num_reads_extending_to_end = 0;
                int num_reads_with_internal_overlaps = 0;
                int repeat_length = 0;
 
                for (int id = 0; id < read_other_ends.size(); ++id)
                {
                    if (maskvec[i].second - read_other_ends[id].first < param.repeat_annotation_gap_thres)
                    {
                        num_reads_considered++;
                        num_reads_extending_to_end++;

                        if ((num_reads_extending_to_end > bridge_threshold) or
                            ((num_reads_considered > bridge_threshold) and
                             (read_other_ends[0].first - read_other_ends[id].first > param.repeat_annotation_gap_thres)))
                        {
                            bridged = false;
                            break;
                        }
                    }
                    else if (read_other_ends[id].second < 0)
                    {
                        num_reads_considered++;

                        if ((num_reads_extending_to_end > bridge_threshold) or
                            ((num_reads_considered > bridge_threshold)) and
                                (read_other_ends[0].first - read_other_ends[id].first > param.repeat_annotation_gap_thres))
                        {
                            bridged = false;
                            break;
                        }
                    }
                    else if (read_other_ends[id].second > 0)
                    {
                        num_reads_with_internal_overlaps++;
                        num_reads_considered++;
                        int id1 = id + 1;
                        int pileup_length = 1;

                        while (id1 < read_other_ends.size())
                        {
                            if (read_other_ends[id].first - read_other_ends[id1].first < param.repeat_annotation_gap_thres)
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
                            int median_id = (id1 - id) / 2;
                            repeat_length = read_other_ends[median_id].first - repeat_annotation[i][j].first;
                            break;
                        }
                    }
                }

                if (bridged){
                    repeats[i].push_back(std::tuple<int, int, int>(repeat_annotation[i][j].first, 1, repeat_length));
                    if (repeat_length > param.repeat_length)
                    {
                        reads[i]->preserve = 1;
                    }
                }
            }
        }
    }

    int count_bridged_repeat_reads = 0;
    int count_long_bridged_repeats = 0;

    for (int i = r_begin; i <= r_end; i++)
    {
        bridged_repeats << "read " << i << " ";
        bridged_repeats << repeats[i].size() << " ";
        for (auto &r : repeats[i])
        {
            bridged_repeats << std::get<2>(r) << " ";
            if(std::get<2>(r) > param.repeat_length)
                count_long_bridged_repeats++;
        }
        bridged_repeats << std::endl;
        if (repeats[i].size() > 0)
            count_bridged_repeat_reads++;
    }

    fprintf(stdout, "INFO, Number of reads with bridged repeats:  %d\n", count_bridged_repeat_reads);
    fprintf(stdout, "INFO, Number of long bridged repeats:  %d\n", count_long_bridged_repeats);
}
