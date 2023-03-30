#include "repeat.hpp"
#include "paf.hpp"
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
                    param.real_reads=0;                
                }
                fprintf(stdout, "Real Reads %d \n", param.real_reads);
            }
            if(param.real_reads){
                int read_id = addStringToMap(std::string(seq->name.s), umap);
                Read *new_r = new Read(read_id, strlen(seq->seq.s), std::string(seq->name.s),
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

    if(!parama.real_reads)
        std::sort(reads.begin(), reads.end(), compare_read);

    return num;
}

void save_repeat_reads(const char *fn, std::vector<Read *> &reads, 
    std::unordered_map<std::string, int> &umap, struct algoParams &param)
{
    int saved_reads=0;
    
    std::ifstream idt(fn);
    std::string read_name;

    while (idt >> read_name){
        int read_id = addStringToMap(std::string(read_name), umap);
        reads[read_id]->save=1;
        saved_reads++;
    }

    fprintf(stdout, "INFO, Number of saved reads from file %d \n", saved_reads);
}

void create_pileup(const char *paffilename, std::vector<Read *> &reads, std::vector<std::vector<Overlap *>> &idx_pileup,
                   std::unordered_map<std::string, int> &umap, struct algoParams &param)
{
    paf_file_t *fp;
    paf_rec_t r;
    fp = paf_open(paffilename);
    int num = 0;
    int check_sym_ovlp = 1;
    int saved_reads = 0;

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
            if(param.real_reads){
                new_ovl->read_A_id_ = addStringToMap(std::string(r.qn), umap);
                new_ovl->read_B_id_ = addStringToMap(std::string(r.tn), umap);
            }else{
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
                first_ovl =  new_ovl;
            }
            else if (check_sym_ovlp && first_ovl->read_A_id_ == new_ovl->read_B_id_ &&
                     first_ovl->read_B_id_ == new_ovl->read_A_id_ &&
                     first_ovl->read_A_match_start_ == new_ovl->read_B_match_start_ &&
                     first_ovl->read_A_match_end_ == new_ovl->read_B_match_end_ &&
                     first_ovl->read_B_match_start_ == new_ovl->read_A_match_start_ &&
                     first_ovl->read_B_match_end_ == new_ovl->read_A_match_end_)
            {
                param.symmetric_overlaps=1;
                check_sym_ovlp=0;
            }

            if (reads[new_ovl->read_A_id_]->save && !reads[new_ovl->read_B_id_]->save && !reads[new_ovl->read_B_id_]->save_overlap)
            {
                reads[new_ovl->read_B_id_]->save_overlap=1;
                saved_reads++;
            }
            else if (!reads[new_ovl->read_A_id_]->save && !reads[new_ovl->read_A_id_]->save_overlap && reads[new_ovl->read_B_id_]->save)
            {
                reads[new_ovl->read_A_id_]->save_overlap = 1;
                saved_reads++;
            }

            num++;
    }

    fprintf(stdout, "INFO, Symmetric overlaps %d \n", param.symmetric_overlaps);
    fprintf(stdout, "INFO, length of alignments  %d()\n", num);
    fprintf(stdout, "INFO, Number of saved reads from overlaps %d \n", saved_reads);
}

uint64_t encodeKmer(const std::string &str)
{
    uint64_t kmer[2] = {0, 0};
    int k = str.length();
    uint64_t shift1 = 2 * (k - 1);

    for (int i = 0; i < k; ++i)
    {
            int c = seq_nt4_table[(uint8_t)str[i]];
            kmer[0] = kmer[0] << 2 | c;
            kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1;
    }

    return kmer[0] < kmer[1] ? kmer[0] : kmer[1];
}

bloom_filter *loadHighFreqKMers(const char *kmer_freq_filename, struct algoParams &param)
{

    std::ifstream idt(kmer_freq_filename);

    std::string kmer;
    int num = 0;
    int freq;
    while (idt >> kmer >> freq)
            num++;

    // kmer length used for kmer counting and mapping must be consistent
    if (num > 0)
    {
            if (kmer.length() != param.kmer_length)
            {
                fprintf(stderr, "ERROR: input list of k-mers and parameter k are inconsistent\n");
                abort();
            }
    }

    // set up bloom filter
    bloom_parameters parameters;
    parameters.false_positive_probability = 0.001;
    parameters.projected_element_count = std::max(num, param.bloom_filter_element_count);
    parameters.compute_optimal_parameters();

    bloom_filter *kmer_filter = new bloom_filter(parameters);

    // read the file again
    idt.clear();
    idt.seekg(0);
    while (idt >> kmer >> freq)
    {
            kmer_filter->insert(encodeKmer(kmer));
    }

    fprintf(stdout, "Number of high freq k-mers %d \n", num);

    return kmer_filter;
}

void create_kmer_from_repetitive_reads(bloom_filter *kmer_filter, const char *repetitive_reads_filename, const algoParams &param)
{
    int k = param.kmer_length;
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1;
    int z = 0;
    std::ifstream idt(repetitive_reads_filename);
    std::string line1, line2;
    int added_kmers = 0;

    while (getline(idt, line1) && getline(idt, line2))
    {
            uint64_t kmer[2] = {0, 0};
            line2.erase(std::remove(line2.begin(), line2.end(), '\n'), line2.end());

            for (int i = 0; i < line2.length(); ++i)
            {
                int c = seq_nt4_table[(uint8_t)line2[i]];
                kmer[0] = (kmer[0] << 2 | c) & mask;             // forward k-mer
                kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1; // reverse k-mer
                z = kmer[0] < kmer[1] ? 0 : 1;                   // strand
                if (i >= k)
                {
                    if (!kmer_filter->contains(kmer[z]))
                    {
                        added_kmers++;
                        kmer_filter->insert(kmer[z]);
                    }
                }
            }
    }

    fprintf(stdout, "Number of additional high freq k-mers in repetitive reads %d \n", added_kmers);
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
        if (reads[i]->save || reads[i]->save_overlap)
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

        } else{

            int parts = read_length / interval_length;

            std::vector<int> initial_stars;
            std::vector<int> final_stars1;
            std::vector<int> final_stars2;

            initial_stars.push_back(0);

            for (int j = 1; j < parts; j++)
            {
                initial_stars.push_back(j*interval_length);
            }
            initial_stars.push_back(read_length);

            final_stars1.push_back(initial_stars[0]);

            int pos = 1;

            for (int k = 0; k < reads[i]->long_repeats1.size(); k++)
            {
                while (reads[i]->long_repeats1[k].first > initial_stars[pos] and (pos < initial_stars.size()-1))
                {
                    final_stars1.push_back(initial_stars[pos]);
                    pos++;
                }
                while (reads[i]->long_repeats1[k].second >= initial_stars[pos] and (pos < initial_stars.size() - 1))
                {
                    pos++;
                }
            }

            while(pos<initial_stars.size()){
                final_stars1.push_back(initial_stars[pos]);
                pos++;
            }

            final_stars2.push_back(final_stars1[0]);

            pos = 1;

            for (int k = 0; k < reads[i]->long_repeats2.size(); k++)
            {
                while (reads[i]->long_repeats2[k].first > final_stars1[pos] and (pos < final_stars1.size() - 1))
                {
                    final_stars2.push_back(initial_stars[pos]);
                    pos++;
                }
                while (reads[i]->long_repeats2[k].second >= final_stars1[pos] and (pos < final_stars1.size() - 1))
                {
                    pos++;
                }
            }

            while (pos < final_stars1.size())
            {
                final_stars2.push_back(final_stars1[pos]);
                pos++;
            }

            if (final_stars2.size() == 2)
            {
                if (!param.real_reads)
                {
                    reads_final << ">read=" << read_num << "," << align << ",position="
                                << start_pos << "-" << end_pos
                                << ",length=" << read_length
                                << read_name.substr(read_name.find_last_of(',')) << "\n";
                }else{
                    reads_final << ">read=" << read_num << "," << read_name << "\n";
                }

                    reads_final << read_seq << "\n";
                    read_num++;
            }

            else {
                    for (int j = 0; j < final_stars2.size() - 2; j++)
                    {

                    if (!param.real_reads)
                    {
                        if (align.compare("forward") == 0)
                        {
                            reads_final << ">read=" << read_num << "," << align << ",position="
                                        << start_pos + final_stars2[j] << "-"
                                        << start_pos + final_stars2[j + 2]
                                        << ",length=" << final_stars2[j + 2] - final_stars2[j]
                                        << read_name.substr(read_name.find_last_of(',')) << "\n";
                        }
                        else if (align.compare("reverse") == 0)
                        {
                            reads_final << ">read=" << read_num << "," << align << ",position="
                                        << end_pos - final_stars2[j + 2] << "-"
                                        << end_pos - final_stars2[j]
                                        << ",length=" << final_stars2[j + 2] - final_stars2[j]
                                        << read_name.substr(read_name.find_last_of(',')) << "\n";
                        }
                    }else{
                        reads_final << ">read=" << read_num << "," << read_name << "\n";
                    }
                    reads_final << read_seq.substr(final_stars2[j], final_stars2[j + 2] - final_stars2[j]) << "\n";
                    read_num++;
                    }
            }
        }
    }
}

void break_long_reads(const char *readfilename, const char *paffilename, const char *kmer_freq_filename,
                      const char *repeatreadsfilename, struct algoParams &param)
{

    std::ofstream reads_final("output_reads.fasta");

    int n_read;
    std::vector<Read *> reads;

    // hash: read id -> number
    std::unordered_map<std::string, int> umap; // size = count of reads

    n_read = loadFASTA(readfilename, reads, umap, param);
    if (repeatreadsfilename)
        save_repeat_reads(repeatreadsfilename, reads, umap, param);

    std::vector<std::vector<Overlap *>> idx_pileup; // this is the pileup

    for (int i = 0; i < n_read; i++)
    {
        idx_pileup.push_back(std::vector<Overlap *>());
    }

    create_pileup(paffilename, reads, idx_pileup, umap, param);

    bloom_filter *kmer_filter = loadHighFreqKMers(kmer_freq_filename, param);

    if (repeatreadsfilename)
        create_kmer_from_repetitive_reads(kmer_filter, "repetitive.fasta", param);

    repeat_annotate1(reads, idx_pileup, param);
    repeat_annotate2(reads, kmer_filter, param);

    break_reads(param, n_read, reads, reads_final);
}
