#include <map>
#include "chop.hpp"

std::map<uint64_t, int> loadAllKMers(const char *kmer_freq_filename)
{

    std::map<uint64_t, int> map;

    std::ifstream idt(kmer_freq_filename);

    std::string kmer;
    int freq;

    while (idt >> kmer >> freq)
    {
        map[encodeKmer(kmer)] = freq;
    }

    return map;
}

void create_kmer_hist_from_reads(const char *reads_filename, const char *kmer_freq_filename, const algoParams &param)
{

    std::map<uint64_t, int> map = loadAllKMers(kmer_freq_filename);

    int k = param.kmer_length;
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1;
    int z = 0;
    std::ifstream idt(reads_filename);
    std::string line1, line2;
    int added_kmers = 0;

    while (getline(idt, line1) && getline(idt, line2))
    {
        fprintf(stdout, "%s\n", line1.c_str());
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
                int freq = map.at(kmer[z]);
                fprintf(stdout, "%d: ", freq);
                for (int j=0; j < std::min(freq, 100) ; j++){
                    fprintf(stdout, "*");
                }
                fprintf(stdout, "\n");
            }
        }
    }
}