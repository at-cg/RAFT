#include "repeat.hpp"
#include <unistd.h>
#include <regex>
#include <unordered_map>

#ifndef COMPARE_READ
#define COMPARE_READ
bool compare_read(Read *read1, Read *read2)
{
    return read1->id < read2->id;
}
#endif

int get_id_from_string(const char *name_str)
{
    const char *sub0 = strchr(name_str, '=') + 1;
    const char *sub1 = strchr(sub0, ',');

    char substr[15];
    strncpy(substr, sub0, strlen(sub0) - strlen(sub1));
    substr[strlen(sub0) - strlen(sub1)] = 0;
    return atoi(substr);
}

int get_start_pos_from_string(const char *name_str)
{
    const char *sub0 = strchr(name_str, ',');
    const char *sub1 = strchr(sub0, '=') + 1;
    const char *sub2 = strchr(sub1, '-');

    char substr[15];
    strncpy(substr, sub1, strlen(sub1) - strlen(sub2));
    substr[strlen(sub1) - strlen(sub2)] = 0;
    return atoi(substr);
}

int get_end_pos_from_string(const char *name_str)
{

    const char *sub0 = strchr(name_str, '-') + 1;
    const char *sub1 = strchr(sub0, ',');

    char substr[15];
    strncpy(substr, sub0, strlen(sub0) - strlen(sub1));
    substr[strlen(sub0) - strlen(sub1)] = 0;
    return atoi(substr);
}

std::string get_alignment_from_string(const char *name_str)
{

    const char *sub0 = strchr(name_str, ',') + 1;
    const char *sub1 = strchr(sub0, ',');

    char substr[15];
    strncpy(substr, sub0, strlen(sub0) - strlen(sub1));
    substr[strlen(sub0) - strlen(sub1)] = 0;
    return std::string(substr);
}

std::string get_chr_from_string(const char *name_str)
{

    const char *sub0 = strrchr(name_str, ',') + 1;

    char substr[15];
    strncpy(substr, sub0, strlen(sub0));
    substr[strlen(sub0)] = 0;
    return std::string(substr);
}

// save fastq read identifier into hash table, and give it an integer id
int addStringToMap(const std::string &str, std::unordered_map<std::string, int> &umap)
{
    if (umap.find(str) != umap.end())
    {
        return umap[str];
    }
    else
    {
        int key = (int)umap.size();
        umap[str] = key;
        return key;
    }
}

// parse + save all reads
int loadFASTA(const char *fn, std::vector<Read *> &reads, std::unordered_map<std::string, int> &umap, struct algoParams &param)
{
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(fn, "r");
    seq = kseq_init(fp);
    int num = 0;

    while ((l = kseq_read(seq)) >= 0)
    {
        if (num == 0)
        {
            if (std::regex_match(seq->name.s, std::regex("^read=[0-9]+,[a-z]+,position=[0-9]+-[0-9]+,length=[0-9]+,(.*)")))
            {
                param.real_reads = 0;
            }
            fprintf(stdout, "Real Reads %d \n", param.real_reads);
        }

        int read_id = addStringToMap(std::string(seq->name.s), umap);

        if (param.real_reads)
        {
            Read *new_r = new Read(read_id, strlen(seq->seq.s), std::string(seq->name.s),
                                   std::string(seq->seq.s));
            reads.push_back(new_r);
        }
        else
        {
            Read *new_r = new Read(read_id, strlen(seq->seq.s), std::string(seq->name.s), std::string(seq->seq.s), 
                            get_start_pos_from_string(seq->name.s), get_end_pos_from_string(seq->name.s),
                            get_alignment_from_string(seq->name.s), get_chr_from_string(seq->name.s));
            reads.push_back(new_r);
        }

        num++;
    }

    kseq_destroy(seq);
    gzclose(fp);

    return num;
}

void create_pileup(const char *paffilename, std::vector<Read *> &reads, std::vector<std::vector<Overlap *>> &idx_pileup,
                   std::unordered_map<std::string, int> &umap, struct algoParams &param)
{
    paf_file_t *fp;
    paf_rec_t r;
    fp = paf_open(paffilename);
    int num = 0;
    int check_sym_ovlp = 0;
    if (!param.symmetric_overlaps)
        check_sym_ovlp = 1;

    Overlap *first_ovl = new Overlap();

    // int count_of_non_overlaps = 0;
    while (paf_read(fp, &r) >= 0)
    {
        // if (r.qe - r.qs == r.ql || r.te - r.ts == r.tl ||
        //     (r.rev == 0 && r.qs > 0 && r.qe == r.ql && r.ts == 0 && r.te < r.tl) ||
        //     (r.rev == 0 && r.qs == 0 && r.qe < r.ql && r.ts > 0 && r.te == r.tl) ||
        //     (r.rev == 1 && r.qs == 0 && r.qe < r.ql && r.ts == 0 && r.te < r.tl) ||
        //     (r.rev == 1 && r.qs > 0 && r.qe == r.ql && r.ts > 0 && r.te == r.tl))

        Overlap *new_ovl = new Overlap();

        new_ovl->read_A_match_start_ = r.qs;
        new_ovl->read_B_match_start_ = r.ts;
        new_ovl->read_A_match_end_ = r.qe;
        new_ovl->read_B_match_end_ = r.te;

        new_ovl->read_A_id_ = addStringToMap(std::string(r.qn), umap);
        new_ovl->read_B_id_ = addStringToMap(std::string(r.tn), umap);
        
        idx_pileup[new_ovl->read_A_id_].push_back(new_ovl);
        if (new_ovl->read_A_id_ != new_ovl->read_B_id_ && !param.symmetric_overlaps)
        {
            idx_pileup[new_ovl->read_B_id_].push_back(new_ovl);
        }

        if (num == 0)
        {
            first_ovl = new_ovl;
        }
        else if (check_sym_ovlp && first_ovl->read_A_id_ == new_ovl->read_B_id_ &&
                 first_ovl->read_B_id_ == new_ovl->read_A_id_ &&
                 first_ovl->read_A_match_start_ == new_ovl->read_B_match_start_ &&
                 first_ovl->read_A_match_end_ == new_ovl->read_B_match_end_ &&
                 first_ovl->read_B_match_start_ == new_ovl->read_A_match_start_ &&
                 first_ovl->read_B_match_end_ == new_ovl->read_A_match_end_)
        {
            param.symmetric_overlaps = 1;
            check_sym_ovlp = 0;
        }

        num++;
    }

    fprintf(stdout, "INFO, Symmetric overlaps %d \n", param.symmetric_overlaps);
    fprintf(stdout, "INFO, length of alignments  %d()\n", num);
}

void break_reads(const algoParams &param, int n_read, std::vector<Read *> &reads, std::ofstream &reads_final)
{
    int read_num = 1;
    int interval_length = param.interval_length;

    for (int i = 0; i < n_read; i++)
    {

        std::string read_name = reads[i]->name;
        std::string read_seq = reads[i]->bases;
        int read_length = reads[i]->len;
        int start_pos = reads[i]->start_pos;
        int end_pos = reads[i]->end_pos;
        std::string align = reads[i]->align;
        std::string chr = reads[i]->chr;

        int parts = read_length / interval_length;

        std::vector<int> initial_stars;
        std::vector<int> final_stars;

        initial_stars.push_back(0);

        for (int j = 1; j < parts; j++)
        {
            initial_stars.push_back(j * interval_length);
        }
        initial_stars.push_back(read_length);

        final_stars.push_back(initial_stars[0]);

        int pos = 1;

        for (int k = 0; k < reads[i]->long_repeats.size(); k++)
        {
            while (reads[i]->long_repeats[k].first > initial_stars[pos] and (pos < initial_stars.size() - 1))
            {
                final_stars.push_back(initial_stars[pos]);
                pos++;
            }
            while (reads[i]->long_repeats[k].second >= initial_stars[pos] and (pos < initial_stars.size() - 1))
            {
                pos++;
            }
        }

        while (pos < initial_stars.size())
        {
            final_stars.push_back(initial_stars[pos]);
            pos++;
        }

        int div = param.read_length / param.interval_length;

        if (final_stars.size() <= (div+1))
        {
            if (!param.real_reads)
            {
                reads_final << ">read=" << read_num << "," << align << ",position="
                            << start_pos << "-" << end_pos
                            << ",length=" << read_length
                            << read_name.substr(read_name.find_last_of(',')) << "\n";
            }
            else
            {
                reads_final << ">read=" << read_num << "," << read_name << "\n";
            }

            reads_final << read_seq << "\n";
            read_num++;
        }
        else
        {
            int fragments = 1 + (final_stars.size() - (div+1)) / div;

            int remaining_markers = (final_stars.size() - (div + 1)) % div;

            if ( remaining_markers ){
                fragments ++;
            }

            int pos = 0;

            for (int j = 1; j <= fragments-1; j++)
            {
                int overlap_length = param.overlap_length;
                if(j==1){
                    overlap_length=0;
                }

                int last_pos = final_stars[pos + div];
                if (j==fragments)
                {
                    last_pos = final_stars.back();
                }

                if (!param.real_reads)
                {
                    if (align.compare("forward") == 0)
                    {
                        reads_final << ">read=" << read_num << "," << align << ",position="
                                    << start_pos + final_stars[pos] - overlap_length << "-"
                                    << start_pos + last_pos
                                    << ",length=" << last_pos - final_stars[pos] + overlap_length
                                    << read_name.substr(read_name.find_last_of(',')) << "\n";
                    }
                    else if (align.compare("reverse") == 0)
                    {
                        reads_final << ">read=" << read_num << "," << align << ",position="
                                    << end_pos - last_pos << "-"
                                    << end_pos - final_stars[pos] + overlap_length
                                    << ",length=" << last_pos - final_stars[pos] + overlap_length
                                    << read_name.substr(read_name.find_last_of(',')) << "\n";
                    }
                }
                else
                {
                    reads_final << ">read=" << read_num << "," << read_name << "\n";
                }
                reads_final << read_seq.substr(final_stars[pos] - overlap_length, last_pos - final_stars[pos] + overlap_length) << "\n";
                read_num++;
                pos = pos + div;
            }

            int overlap_length = param.overlap_length;
            int starting_position = final_stars[pos];

            if (pos == 0)
            {
                overlap_length = 0;
            }

            if (pos != 0 && final_stars.back() - final_stars[pos] < param.read_length)
            {
                starting_position = final_stars.back()-param.read_length;
            }


            if (!param.real_reads)
            {
                if (align.compare("forward") == 0)
                {
                    reads_final << ">read=" << read_num << "," << align << ",position="
                                << start_pos + starting_position - overlap_length << "-"
                                << start_pos + final_stars.back()
                                << ",length=" << final_stars.back() - starting_position + overlap_length
                                << read_name.substr(read_name.find_last_of(',')) << "\n";
                }
                else if (align.compare("reverse") == 0)
                {
                    reads_final << ">read=" << read_num << "," << align << ",position="
                                << end_pos - final_stars.back() << "-"
                                << end_pos - starting_position + overlap_length
                                << ",length=" << final_stars.back() - starting_position + overlap_length
                                << read_name.substr(read_name.find_last_of(',')) << "\n";
                }
            }
            else
            {
                reads_final << ">read=" << read_num << "," << read_name << "\n";
            }
            reads_final << read_seq.substr(starting_position - overlap_length, final_stars.back() - starting_position + overlap_length) << "\n";
            read_num++;
        }
    }
}


void break_long_reads(const char *readfilename, const char *paffilename, const char *repeatreadsfilename, struct algoParams &param)
{

    std::ofstream reads_final("output_reads.fasta");

    int n_read;
    std::vector<Read *> reads;

    // hash: read id -> number
    std::unordered_map<std::string, int> umap; // size = count of reads

    n_read = loadFASTA(readfilename, reads, umap, param);

    std::vector<std::vector<Overlap *>> idx_pileup; // this is the pileup

    for (int i = 0; i < n_read; i++)
    {
        idx_pileup.push_back(std::vector<Overlap *>());
    }

    create_pileup(paffilename, reads, idx_pileup, umap, param);

    umap.clear();

    repeat_annotate(reads, idx_pileup, param);

    break_reads(param, n_read, reads, reads_final);
}
