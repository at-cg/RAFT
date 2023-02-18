#include "repeat.hpp"
#include <unistd.h>

#ifndef COMPARE_READ
#define COMPARE_READ
bool compare_read(Read * read1, Read * read2)
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

void loadPAF(const char *fn, std::vector<Overlap *> &alns)
{
    paf_file_t *fp;
    paf_rec_t r;
    fp = paf_open(fn);
    int num = 0;
    int count_of_non_overlaps = 0;
    while (paf_read(fp, &r) >= 0)
    {
        // if (r.qe - r.qs == r.ql || r.te - r.ts == r.tl ||
        //     (r.rev == 0 && r.qs > 0 && r.qe == r.ql && r.ts == 0 && r.te < r.tl) ||
        //     (r.rev == 0 && r.qs == 0 && r.qe < r.ql && r.ts > 0 && r.te == r.tl) ||
        //     (r.rev == 1 && r.qs == 0 && r.qe < r.ql && r.ts == 0 && r.te < r.tl) ||
        //     (r.rev == 1 && r.qs > 0 && r.qe == r.ql && r.ts > 0 && r.te == r.tl))
        // {
            num++;
            Overlap *new_ovl = new Overlap();

            new_ovl->read_A_match_start_ = r.qs;
            new_ovl->read_B_match_start_ = r.ts;
            new_ovl->read_A_match_end_ = r.qe;
            new_ovl->read_B_match_end_ = r.te;
            new_ovl->read_A_id_ = get_id_from_string(r.qn) - 1;
            new_ovl->read_B_id_ = get_id_from_string(r.tn) - 1;
            alns.push_back(new_ovl);
        // } else{
        //     count_of_non_overlaps++;
        // }
    }
}

// parse + save all reads
int loadFASTA(const char *fn, std::vector<Read *> &reads, const algoParams &param)
{
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(fn, "r");
    seq = kseq_init(fp);
    int num = 0;
 
    while ((l = kseq_read(seq)) >= 0)
    {       
            if(param.real_reads){
                Read *new_r = new Read(get_id_from_string(seq->name.s) - 1, strlen(seq->seq.s), std::string(seq->name.s),
                                       std::string(seq->seq.s));
                reads.push_back(new_r);
            }else{
                Read *new_r = new Read(get_id_from_string(seq->name.s) - 1, strlen(seq->seq.s), std::string(seq->name.s),
                                    std::string(seq->seq.s), get_start_pos_from_string(seq->name.s), get_end_pos_from_string(seq->name.s),
                                    get_alignment_from_string(seq->name.s), get_chr_from_string(seq->name.s));
                reads.push_back(new_r);
            }
            
            num++;
    }

    kseq_destroy(seq);
    gzclose(fp);

   std::sort(reads.begin(), reads.end(), compare_read);

    return num;
}

void create_pileup(const char *paffilename, std::vector<std::vector<Overlap *>> &idx_pileup)
{
    std::vector<Overlap *> aln;
    loadPAF(paffilename, aln);

    fprintf(stdout, "INFO, length of alignments  %lu()\n", aln.size());

    for (int i = 0; i < aln.size(); i++)
    {
            if (aln[i]->read_A_id_ == aln[i]->read_B_id_)
            {
                idx_pileup[aln[i]->read_A_id_].push_back(aln[i]);
            }
            else
            {
                idx_pileup[aln[i]->read_A_id_].push_back(aln[i]);
                idx_pileup[aln[i]->read_B_id_].push_back(aln[i]);
            }
    }
}

void break_real_reads(const algoParams &param, int n_read, std::vector<Read *> &reads, std::ofstream &reads_final)
{
    int read_num = 1;

    int overlap_length = param.overlap_length;
    int uniform_read_length = param.uniform_read_length;
    int distance = uniform_read_length - overlap_length;

    for (int i = 0; i < n_read; i++)
    {

            std::string read_name = reads[i]->name;
            std::string read_seq = reads[i]->bases;
            int read_length = reads[i]->len;

            if (reads[i]->preserve)
            {
                int non_repeat_start = 0;
                int repeat_start = 0;
                int repeat_end = 0;
                for (int j = 0; j <= reads[i]->long_repeats.size(); j++)
                {

                    if (j != reads[i]->long_repeats.size())
                    {

                        repeat_start = std::get<0>(reads[i]->long_repeats[j]);
                        repeat_end = std::get<1>(reads[i]->long_repeats[j]);
                    }
                    else
                    {
                        repeat_start = read_length;
                        repeat_end = read_length;
                    }

                    if (repeat_start == 0 && repeat_end == read_length)
                    {
                        reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";
                        reads_final << read_seq << "\n";
                        read_num++;
                    }
                    else if (repeat_start == 0)
                    {

                        reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";
                        reads_final << read_seq.substr(0, repeat_end + overlap_length) << "\n";
                        read_num++;
                    }
                    else
                    {

                        int parts = (repeat_start - non_repeat_start) / distance;
                        int k = 0;

                        int overlap_length2 = overlap_length;
                        if (repeat_end == read_length)
                        {
                            overlap_length2 = 0;
                        }

                        for (k = 0; k < parts - 1; k++)
                        {

                            reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";

                            reads_final << read_seq.substr(non_repeat_start + k * distance, uniform_read_length) << "\n";
                            read_num++;
                        }

                        int last_length = (repeat_start - non_repeat_start) - ((parts - 1) * distance);

                        if (last_length > overlap_length)
                        {
                            reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";

                            reads_final << read_seq.substr(non_repeat_start + k * distance, last_length) << "\n";
                            read_num++;
                        }

                        if (repeat_start != read_length)
                        {

                            reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";

                            reads_final << read_seq.substr(repeat_start - overlap_length, repeat_end - repeat_start + overlap_length + overlap_length2) << "\n";
                            read_num++;
                        }
                    }

                    if (repeat_end == read_length)
                    {
                        break;
                    }
                    else
                    {
                        non_repeat_start = repeat_end + 1;
                    }
                }
            }
            else if (read_length <= param.read_length_threshold)
            {
                reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";
                reads_final << read_seq << "\n";
                read_num++;
            }
            else
            {
                int parts = read_length / distance;
                int j;

                for (j = 0; j < parts - 1; j++)
                {

                    reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";

                    reads_final << read_seq.substr(0 + j * distance, uniform_read_length) << "\n";
                    read_num++;
                }

                int last_length = read_length - ((parts - 1) * distance);

                if (last_length > overlap_length)
                {
                    reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";

                    reads_final << read_seq.substr(0 + j * distance, last_length) << "\n";
                    read_num++;
                }
            }
    }
}

void break_simulated_reads(const algoParams &param, int n_read, std::vector<Read *> &reads, std::ofstream &reads_final)
{
    std::ofstream bed_fragmented(param.outputfilename + ".fragmentation.bed");
    std::ofstream bed_preserved(param.outputfilename + ".preserved.bed");
    int read_num = 1;

    int overlap_length = param.overlap_length;
    int uniform_read_length = param.uniform_read_length;
    int distance = uniform_read_length - overlap_length;

    for (int i = 0; i < n_read; i++)
    {

            std::string read_name = reads[i]->name;
            std::string read_seq = reads[i]->bases;
            int read_length = reads[i]->len;
            int start_pos = reads[i]->start_pos;
            int end_pos = reads[i]->end_pos;
            std::string align = reads[i]->align;
            std::string chr = reads[i]->chr;

            if (reads[i]->preserve)
            {
                int non_repeat_start = 0;
                int repeat_start = 0;
                int repeat_end = 0;
                for (int j = 0; j <= reads[i]->long_repeats.size(); j++)
                {

                    if (j != reads[i]->long_repeats.size())
                    {

                        repeat_start = std::get<0>(reads[i]->long_repeats[j]);
                        repeat_end = std::get<1>(reads[i]->long_repeats[j]);
                    }
                    else
                    {
                        repeat_start = read_length;
                        repeat_end = read_length;
                    }

                    if (repeat_start == 0 && repeat_end == read_length)
                    {
                        reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";
                        reads_final << read_seq << "\n";
                        read_num++;
                        bed_preserved << chr << "\t" << start_pos << "\t" << end_pos << std::endl;
                    }
                    else if (repeat_start == 0)
                    {
                        if (align.compare("forward") == 0)
                        {
                            reads_final << ">read=" << read_num << "," << align << ",position="
                                        << start_pos << "-" << start_pos + repeat_end + overlap_length
                                        << ",length=" << repeat_end + overlap_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                            reads_final << read_seq.substr(0, repeat_end + overlap_length) << "\n";
                            read_num++;
                            bed_preserved << chr << "\t" << start_pos << "\t" << start_pos + repeat_end + overlap_length << std::endl;
                        }
                        else if (align.compare("reverse") == 0)
                        {
                            reads_final << ">read=" << read_num << "," << align << ",position="
                                        << end_pos - repeat_end - overlap_length << "-" << end_pos
                                        << ",length=" << repeat_end + overlap_length
                                        << read_name.substr(read_name.find_last_of(',')) << "\n";

                            reads_final << read_seq.substr(0, repeat_end + overlap_length) << "\n";
                            read_num++;
                            bed_preserved << chr << "\t" << end_pos - repeat_end - overlap_length << "\t" << end_pos << std::endl;
                        }
                    }
                    else
                    {

                        int parts = (repeat_start - non_repeat_start) / distance;
                        int k = 0;

                        int overlap_length2 = overlap_length;
                        if (repeat_end == read_length)
                        {
                            overlap_length2 = 0;
                        }

                        if (align.compare("forward") == 0)
                        {

                            for (k = 0; k < parts - 1; k++)
                            {

                                reads_final << ">read=" << read_num << "," << align << ",position="
                                            << start_pos + non_repeat_start + k * distance << "-"
                                            << start_pos + non_repeat_start + k * distance + uniform_read_length
                                            << ",length=" << uniform_read_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                                reads_final << read_seq.substr(non_repeat_start + k * distance, uniform_read_length) << "\n";
                                read_num++;
                            }

                            int last_length = (repeat_start - non_repeat_start) - ((parts - 1) * distance);

                            if (last_length > overlap_length)
                            {
                                reads_final << ">read=" << read_num << "," << align << ",position="
                                            << start_pos + non_repeat_start + k * distance << "-"
                                            << start_pos + non_repeat_start + k * distance + last_length
                                            << ",length=" << last_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                                reads_final << read_seq.substr(non_repeat_start + k * distance, last_length) << "\n";
                                read_num++;
                            }

                            bed_fragmented << chr << "\t" << start_pos + non_repeat_start
                                           << "\t" << start_pos + non_repeat_start + k * distance + last_length << std::endl;

                            if (repeat_start != read_length)
                            {

                                reads_final << ">read=" << read_num << "," << align << ",position="
                                            << start_pos + repeat_start - overlap_length << "-"
                                            << start_pos + repeat_end + overlap_length2
                                            << ",length=" << repeat_end - repeat_start + overlap_length + overlap_length2
                                            << read_name.substr(read_name.find_last_of(',')) << "\n";

                                reads_final << read_seq.substr(repeat_start - overlap_length, repeat_end - repeat_start + overlap_length + overlap_length2) << "\n";
                                read_num++;
                                bed_preserved << chr << "\t" << start_pos + repeat_start - overlap_length
                                              << "\t" << start_pos + repeat_end + overlap_length2 << std::endl;
                            }
                        }
                        else if (align.compare("reverse") == 0)
                        {
                            for (k = 0; k < parts - 1; k++)
                            {

                                reads_final << ">read=" << read_num << "," << align << ",position="
                                            << end_pos - non_repeat_start - k * distance - uniform_read_length << "-"
                                            << end_pos - non_repeat_start - k * distance
                                            << ",length=" << uniform_read_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                                reads_final << read_seq.substr(non_repeat_start + k * distance, uniform_read_length) << "\n";
                                read_num++;
                            }

                            int last_length = (repeat_start - non_repeat_start) - ((parts - 1) * distance);

                            if (last_length > overlap_length)
                            {

                                reads_final << ">read=" << read_num << "," << align << ",position="
                                            << end_pos - non_repeat_start - k * distance - last_length << "-"
                                            << end_pos - non_repeat_start - k * distance
                                            << ",length=" << last_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                                reads_final << read_seq.substr(non_repeat_start + k * distance, last_length) << "\n";
                                read_num++;
                            }

                            bed_fragmented << chr << "\t" << end_pos - non_repeat_start - k * distance - last_length
                                           << "\t" << end_pos - non_repeat_start << std::endl;

                            if (repeat_start != read_length)
                            {
                                reads_final << ">read=" << read_num << "," << align << ",position="
                                            << end_pos - repeat_end - overlap_length2 << "-"
                                            << end_pos - repeat_start + overlap_length
                                            << ",length=" << repeat_end - repeat_start + overlap_length + overlap_length2
                                            << read_name.substr(read_name.find_last_of(',')) << "\n";

                                reads_final << read_seq.substr(repeat_start - overlap_length, repeat_end - repeat_start + overlap_length + overlap_length2) << "\n";
                                read_num++;

                                bed_preserved << chr << "\t" << end_pos - repeat_end - overlap_length2
                                              << "\t" << end_pos - repeat_start + overlap_length << std::endl;
                            }
                        }
                    }

                    if (repeat_end == read_length)
                    {
                        break;
                    }
                    else
                    {
                        non_repeat_start = repeat_end + 1;
                    }
                }
            }
            else if (read_length <= param.read_length_threshold)
            {
                reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";
                reads_final << read_seq << "\n";
                read_num++;
            }
            else
            {
                int parts = read_length / distance;
                int j;

                if (align.compare("forward") == 0)
                {

                    for (j = 0; j < parts - 1; j++)
                    {

                        reads_final << ">read=" << read_num << "," << align << ",position="
                                    << start_pos + j * distance << "-" << start_pos + j * distance + uniform_read_length
                                    << ",length=" << uniform_read_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                        reads_final << read_seq.substr(0 + j * distance, uniform_read_length) << "\n";
                        read_num++;
                    }

                    int last_length = read_length - ((parts - 1) * distance);

                    if (last_length > overlap_length)
                    {
                        reads_final << ">read=" << read_num << "," << align << ",position="
                                    << start_pos + j * distance << "-" << start_pos + j * distance + last_length
                                    << ",length=" << last_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                        reads_final << read_seq.substr(0 + j * distance, last_length) << "\n";
                        read_num++;
                    }
                }
                else if (align.compare("reverse") == 0)
                {
                    for (j = 0; j < parts - 1; j++)
                    {

                        reads_final << ">read=" << read_num << "," << align << ",position="
                                    << end_pos - j * distance - uniform_read_length << "-" << end_pos - j * distance
                                    << ",length=" << uniform_read_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                        reads_final << read_seq.substr(0 + j * distance, uniform_read_length) << "\n";
                        read_num++;
                    }

                    int last_length = read_length - ((parts - 1) * distance);

                    if (last_length > overlap_length)
                    {

                        reads_final << ">read=" << read_num << "," << align << ",position="
                                    << end_pos - j * distance - last_length << "-" << end_pos - j * distance
                                    << ",length=" << last_length << read_name.substr(read_name.find_last_of(',')) << "\n";

                        reads_final << read_seq.substr(0 + j * distance, last_length) << "\n";
                        read_num++;
                    }
                }

                bed_fragmented << chr << "\t" << start_pos << "\t" << end_pos << std::endl;
            }
    }
}

void break_long_reads(const char *readfilename, const char *paffilename, const algoParams &param)
{

    std::ofstream reads_final("output_reads.fasta");

    int n_read;
    std::vector<Read *> reads;

    n_read = loadFASTA(readfilename, reads, param);
    std::vector<std::vector<Overlap *>> idx_pileup; // this is the pileup

    for (int i = 0; i < n_read; i++)
    {
            idx_pileup.push_back(std::vector<Overlap *>());
    }

    create_pileup(paffilename, idx_pileup);

    repeat_annotate(reads, param, idx_pileup);

    if(param.real_reads){
        break_real_reads(param, n_read, reads, reads_final);
    }else{
        break_simulated_reads(param, n_read, reads, reads_final);
    }
 
}