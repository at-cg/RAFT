#include "repeat.hpp"

int get_id_from_string(const char *name_str)
{

    const char *sub0 = strchr(name_str, '=');
    const char *sub1 = sub0 + 1;
    const char *sub2 = strchr(sub1, ',');

    char substr[15];
    strncpy(substr, sub1, strlen(sub1) - strlen(sub2));
    substr[strlen(sub1) - strlen(sub2)] = 0;
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

std::string get_alignment_from_string(const char *name_str)
{

    const char *sub0 = strchr(name_str, ',');
    const char *sub1 = sub0 + 1;
    const char *sub2 = strchr(sub1, ',');

    char substr[15];
    strncpy(substr, sub1, strlen(sub1) - strlen(sub2));
    substr[strlen(sub1) - strlen(sub2)] = 0;
    return std::string(substr);
}

int loadPAF(const char *fn, std::vector<Overlap *> &alns)
{
    paf_file_t *fp;
    paf_rec_t r;
    fp = paf_open(fn);
    int num = 0;
    while (paf_read(fp, &r) >= 0)
    {
        num++;
        Overlap *new_ovl = new Overlap();

        new_ovl->read_A_match_start_ = r.qs;
        new_ovl->read_B_match_start_ = r.ts;
        new_ovl->read_A_match_end_ = r.qe;
        new_ovl->read_B_match_end_ = r.te;
        new_ovl->alen = r.ql;
        new_ovl->blen = r.tl;
        new_ovl->reverse_complement_match_ = r.rev;
        new_ovl->diffs = 0;
        new_ovl->read_A_id_ = get_id_from_string(r.qn) - 1;
        new_ovl->read_B_id_ = get_id_from_string(r.tn) - 1;
        alns.push_back(new_ovl);
    }
    return num;
}

// parse + save all reads
int loadFASTA(const char *fn, std::vector<Read *> &reads)
{
    gzFile fp;
    kseq_t *seq;
    int l;
    int num = 0;

    fp = gzopen(fn, "r");
    seq = kseq_init(fp);

    while ((l = kseq_read(seq)) >= 0)
    {
        Read *new_r = new Read(num, strlen(seq->seq.s), std::string(seq->name.s), std::string(seq->seq.s), get_start_pos_from_string(seq->name.s), get_alignment_from_string(seq->name.s));
        reads.push_back(new_r);
        num++;
    }

    kseq_destroy(seq);
    gzclose(fp);

    return num;
}

void break_long_reads(const char *readfilename, const char *paffilename, const algoParams &param)
{
    std::ofstream reads_final("output_reads.fasta");

    int n_read;
    int64_t n_aln = 0;
    std::vector<Read *> reads;
    std::vector<Overlap *> aln;

    n_read = loadFASTA(readfilename, reads);
    n_aln = loadPAF(paffilename, aln);
    
    repeat_annotate(reads, aln, param);
    int read_num = 1;

    int overlap_length = param.repeat_length;
    int uniform_read_length = overlap_length * 2;

    for (int i = 0; i < n_read; i++){

        std::string read_name = reads[i]->name;
        std::string read_seq = reads[i]->bases;
        int read_length = reads[i]->len;
        int start_pos = reads[i]->start_pos;
        std::string align = reads[i]->align;

        if (reads[i]->preserve)
        {
            reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";
            reads_final << read_seq << "\n";
            read_num++;
        }
        else if (read_length <= 20000)
        {
            reads_final << ">read=" << read_num << read_name.substr(read_name.find(',')) << "\n";
            reads_final << read_seq << "\n";
            read_num++;
        }
        else
        {   
            int parts = read_length / overlap_length;
            int j;

            for (j=0; j < parts-1; j++){

                reads_final << ">read=" << read_num << "," << align << ",position=" << start_pos + j * overlap_length << "-" << start_pos + uniform_read_length + j * overlap_length
                            << ",length=" << 10000 << read_name.substr(read_name.find_last_of(',')) << "\n";

                reads_final << read_seq.substr(0 + j * overlap_length, uniform_read_length) << "\n";
                read_num++;
            }

            int last_length = read_length - (parts - 1) * overlap_length;

            reads_final << "read=" << read_num << "," << align << ",position=" << start_pos + j * overlap_length << "-" << start_pos + j * overlap_length + last_length
                        << ",length=" << last_length << read_name.substr(read_name.find_last_of(',')) << "\n";

            reads_final << read_seq.substr(0 + j * overlap_length, uniform_read_length) << "\n";
            read_num++;

        }
    }
}