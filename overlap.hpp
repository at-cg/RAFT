#ifndef OVERLAP_H
#define OVERLAP_H

#include <vector>
#include <iostream>


class Overlap
{
public:
    Overlap(){};

    ~Overlap(){};

    void show() { printf("%d %d %d [%d...%d]/%d x [%d...%d]/%d %d diffs\n", read_A_id_, read_B_id_,
                         reverse_complement_match_,
                         read_A_match_start_, read_A_match_end_, alen, read_B_match_start_, read_B_match_end_, blen, diffs); };
    int read_A_id_, read_B_id_;
    int alen; // length of read a
    int blen; // length of read b
     int diffs; // differences
    int read_A_match_start_, read_B_match_start_; // starting position alignment in read a and read b
    int read_A_match_end_, read_B_match_end_; // ending position of alignment in read a and read b
    int reverse_complement_match_; // reverse_complement_match_, reverse complement = 1, same direction = 0
    bool active = true;
    int length;
};

#endif