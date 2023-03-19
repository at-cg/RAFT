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

   std::sort(reads.begin(), reads.end(), compare_read);

    return num;
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

bloom_filter* loadHighFreqKMers(const char *kmer_freq_filename, struct algoParams &param)
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
    parameters.projected_element_count = std::max(num, 50000000);
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
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1 ;
    int z = 0;
    std::ifstream idt(repetitive_reads_filename);
    std::string line1, line2;
    int added_kmers = 0;

    while (getline(idt, line1) && getline(idt, line2)){
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

        int parts = read_length / interval_length;

        std::vector<int> initial_stars;
        std::vector<int> final_stars;

        initial_stars.push_back(0);

        for (int j = 1; j < parts; j++)
        {
            initial_stars.push_back(j*interval_length);
        }
        initial_stars.push_back(read_length);

        final_stars.push_back(initial_stars[0]);

        int pos = 1;

        for (int k = 0; k < reads[i]->long_repeats.size(); k++)
        {
            while (reads[i]->long_repeats[k].first > initial_stars[pos] and (pos < initial_stars.size()-1))
            {
                final_stars.push_back(initial_stars[pos]);
                pos++;
            }
            while (reads[i]->long_repeats[k].second >= initial_stars[pos] and (pos < initial_stars.size() - 1))
            {
                pos++;
            }
        }

        while(pos<initial_stars.size()){
            final_stars.push_back(initial_stars[pos]);
            pos++;
        }

        if(final_stars.size()==2){
            if (!param.real_reads)
            {
                reads_final << ">read=" << read_num << "," << align << ",position="
                            << start_pos << "-" << end_pos
                            << ",length=" << read_length
                            << read_name.substr(read_name.find_last_of(',')) << "\n";
            }else{
                reads_final << ">read=" << read_num << ", " << read_name << "\n";
            }

                reads_final << read_seq << "\n";
                read_num++;
        }
        else {
            for (int j=0; j < final_stars.size()-2; j++){

                if (!param.real_reads)
                {
                    if (align.compare("forward") == 0)
                    {
                        reads_final << ">read=" << read_num << "," << align << ",position="
                                    << start_pos + final_stars[j] << "-"
                                    << start_pos + final_stars[j + 2]
                                    << ",length=" << final_stars[j + 2] - final_stars[j]
                                    << read_name.substr(read_name.find_last_of(',')) << "\n";
                    }
                    else if (align.compare("reverse") == 0)
                    {
                        reads_final << ">read=" << read_num << "," << align << ",position="
                                    << end_pos - final_stars[j + 2] << "-"
                                    << end_pos - final_stars[j]
                                    << ",length=" << final_stars[j + 2] - final_stars[j]
                                    << read_name.substr(read_name.find_last_of(',')) << "\n";
                    }
                }else{
                    reads_final << ">read=" << read_num << ", " << read_name << "\n";
                }
                    reads_final << read_seq.substr(final_stars[j], final_stars[j + 2] - final_stars[j]) << "\n";
                    read_num++;
            }
        }
    }
}

void break_long_reads(const char *readfilename, const char *kmer_freq_filename, struct algoParams &param)
{

    std::ofstream reads_final("output_reads.fasta");

    int n_read;
    std::vector<Read *> reads;

    // hash: read id -> number
    std::unordered_map<std::string, int> umap; // size = count of reads

    n_read = loadFASTA(readfilename, reads, umap, param);

    bloom_filter* kmer_filter = loadHighFreqKMers(kmer_freq_filename, param);

    if(param.additional_kmers) create_kmer_from_repetitive_reads(kmer_filter, "repetitive.fasta", param);

    repeat_annotate(reads, kmer_filter, param);

    break_reads(param, n_read, reads, reads_final);
}
