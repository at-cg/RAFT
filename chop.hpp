#include "repeat.hpp"
#include <unistd.h>
#include <regex>
#include <unordered_map>

#ifndef COMPARE_READ
#define COMPARE_READ
bool compare_read(Read * read1, Read * read2)
{
    return read1->id < read2->id;
}
#endif

#ifndef COMPARE_OVERLAP
#define COMPARE_OVERLAP
bool compare_overlap(Overlap *ovl1, Overlap *ovl2)
{
    // Returns True if the sum of the match lengths of the two reads in ovl1 > the sum of the  overlap lengths of the two reads in ovl2
    // Returns False otherwise.
    return ((ovl1->read_A_match_end_ - ovl1->read_A_match_start_ + ovl1->read_B_match_end_ - ovl1->read_B_match_start_) > (ovl2->read_A_match_end_ - ovl2->read_A_match_start_ + ovl2->read_B_match_end_ - ovl2->read_B_match_start_));
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
            if (param.real_reads)
            {
                Read *new_r = new Read(addStringToMap(std::string(seq->name.s), umap),
                 strlen(seq->seq.s), std::string(seq->name.s),
                                       std::string(seq->seq.s));
                reads.push_back(new_r);
            }
            else
            {
                Read *new_r = new Read(get_id_from_string(seq->name.s) - 1, strlen(seq->seq.s), std::string(seq->name.s),
                                       std::string(seq->seq.s), get_start_pos_from_string(seq->name.s), get_end_pos_from_string(seq->name.s),
                                       get_alignment_from_string(seq->name.s), get_chr_from_string(seq->name.s));
                reads.push_back(new_r);
            }
            num++;
    }

    kseq_destroy(seq);
    gzclose(fp);

    if(!param.real_reads){
        std::sort(reads.begin(), reads.end(), compare_read);
    }

    return num;
}

void create_pileup(const char *paffilename, std::vector<std::vector<Overlap *>> &idx_pileup,
                   std::unordered_map<std::string, int> &umap, struct algoParams &param)
{
    paf_file_t *fp;
    paf_rec_t r;
    fp = paf_open(paffilename);
    int num = 0;
    int check_sym_ovlp = 1;

    Overlap *first_ovl = new Overlap();
    //int count_of_non_overlaps = 0;
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
            if (param.real_reads)
            {
                new_ovl->read_A_id_ = addStringToMap(std::string(r.qn), umap);
                new_ovl->read_B_id_ = addStringToMap(std::string(r.tn), umap);
            }
            else
            {
                new_ovl->read_A_id_ = get_id_from_string(r.qn) - 1;
                new_ovl->read_B_id_ = get_id_from_string(r.tn) - 1;
            }
            // } else{
            //     count_of_non_overlaps++;
            // }


            if (new_ovl->read_A_id_ == new_ovl->read_B_id_)
            {
                idx_pileup[new_ovl->read_A_id_].push_back(new_ovl);
            }
            else
            {
                idx_pileup[new_ovl->read_A_id_].push_back(new_ovl);
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

    fprintf(stdout, "Symmetric overlaps %d \n", param.symmetric_overlaps);
    fprintf(stdout, "INFO, length of alignments  %d()\n", num);
}

void break_real_reads(const algoParams &param, int n_read, std::vector<Read *> &reads, std::ofstream &reads_final)
{
    int read_num = 1;

    int overlap_length = param.overlap_length;
    int uniform_read_length = param.uniform_read_length;
    int distance = uniform_read_length - overlap_length;

    int count_of_eligible_preserved_reads = 0;
    int count_of_eligible_fragmented_reads = 0;

    for (int i = 0; i < n_read; i++)
    {

            std::string read_name = reads[i]->name;
            std::string read_seq = reads[i]->bases;
            int read_length = reads[i]->len;

            if (read_length <= param.read_length_threshold)
            {
                reads_final << ">read=" << read_num << ", " << read_name << "\n";
                reads_final << read_seq << "\n";
                read_num++;
            }
            else if (reads[i]->preserve)
            {
                count_of_eligible_preserved_reads++;
                reads_final << ">read=" << read_num << ", " << read_name << "\n";
                reads_final << read_seq << "\n";
                read_num++;
            }
            else
            {
                count_of_eligible_fragmented_reads++;
                int parts = read_length / distance;
                int j;

                for (j = 0; j < parts - 2; j++)
                {

                    reads_final << ">read=" << read_num << ", " << read_name << "\n";
                    reads_final << read_seq.substr(0 + j * distance, uniform_read_length) << "\n";
                    read_num++;
                }

                int last_length = read_length - ((parts - 2) * distance);

                reads_final << ">read=" << read_num << ", " << read_name << "\n";
                reads_final << read_seq.substr(0 + j * distance, last_length) << "\n";
                read_num++;
                
            }
    }
    fprintf(stdout, "fraction of eligible preserved reads %f \n", double(count_of_eligible_preserved_reads) / (count_of_eligible_fragmented_reads + count_of_eligible_preserved_reads));
}

void break_simulated_reads(const algoParams &param, int n_read, std::vector<Read *> &reads, std::ofstream &reads_final)
{
    std::ofstream bed_fragmented(param.outputfilename + ".fragmentation.bed");
    std::ofstream bed_preserved(param.outputfilename + ".preserved.bed");
    int read_num = 1;

    int overlap_length = param.overlap_length;
    int uniform_read_length = param.uniform_read_length;
    int distance = uniform_read_length - overlap_length;

    int count_of_eligible_preserved_reads = 0;
    int count_of_eligible_fragmented_reads = 0;

    for (int i = 0; i < n_read; i++)
    {

            std::string read_name = reads[i]->name;
            std::string read_seq = reads[i]->bases;
            int read_length = reads[i]->len;
            int start_pos = reads[i]->start_pos;
            int end_pos = reads[i]->end_pos;
            std::string align = reads[i]->align;
            std::string chr = reads[i]->chr;

            if (read_length <= param.read_length_threshold)
            {
                reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";
                reads_final << read_seq << "\n";
                read_num++;
            }
            else if (reads[i]->preserve)
            {
                count_of_eligible_preserved_reads++;
                reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";
                reads_final << read_seq << "\n";
                read_num++;
                bed_preserved << chr << "\t" << start_pos << "\t" << end_pos << std::endl;
            }
            else
            {
                count_of_eligible_fragmented_reads++;
                int parts = read_length / distance;
                int j;

                if (align.compare("forward") == 0)
                {

                    for (j = 0; j < parts - 2; j++)
                    {

                        reads_final << ">read=" << read_num << "," << align << ",position="
                                    << start_pos + j * distance << "-" << start_pos + j * distance + uniform_read_length
                                    << ",length=" << uniform_read_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                        reads_final << read_seq.substr(0 + j * distance, uniform_read_length) << "\n";
                        read_num++;
                    }

                    int last_length = read_length - ((parts - 2) * distance);

                    reads_final << ">read=" << read_num << "," << align << ",position="
                                << start_pos + j * distance << "-" << start_pos + j * distance + last_length
                                << ",length=" << last_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                    reads_final << read_seq.substr(0 + j * distance, last_length) << "\n";
                    read_num++;
                    
                }
                else if (align.compare("reverse") == 0)
                {
                    for (j = 0; j < parts - 2; j++)
                    {

                        reads_final << ">read=" << read_num << "," << align << ",position="
                                    << end_pos - j * distance - uniform_read_length << "-" << end_pos - j * distance
                                    << ",length=" << uniform_read_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                        reads_final << read_seq.substr(0 + j * distance, uniform_read_length) << "\n";
                        read_num++;
                    }

                    int last_length = read_length - ((parts - 2) * distance);

                    reads_final << ">read=" << read_num << "," << align << ",position="
                                << end_pos - j * distance - last_length << "-" << end_pos - j * distance
                                << ",length=" << last_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                    reads_final << read_seq.substr(0 + j * distance, last_length) << "\n";
                    read_num++;
                    
                }

                bed_fragmented << chr << "\t" << start_pos << "\t" << end_pos << std::endl;
            }
    }
    fprintf(stdout, "fraction of eligible preserved reads %f \n", double(count_of_eligible_preserved_reads) / (count_of_eligible_fragmented_reads + count_of_eligible_preserved_reads));
}

void break_long_reads(const char *readfilename, const char *paffilename, struct algoParams &param)
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

    create_pileup(paffilename, idx_pileup, umap, param);
    repeat_annotate(reads, idx_pileup, param);

    if (param.real_reads)
    {
            break_real_reads(param, n_read, reads, reads_final);
    }
    else
    {
            break_simulated_reads(param, n_read, reads, reads_final);
    }
}
