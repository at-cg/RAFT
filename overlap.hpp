#ifndef OVERLAP_H
#define OVERLAP_H

#include <vector>
#include <iostream>


class Overlap
{
public:
    Overlap(){};

    ~Overlap(){};
    
    int read_A_id_, read_B_id_;
    int read_A_match_start_, read_B_match_start_; // starting position alignment in read a and read b
    int read_A_match_end_, read_B_match_end_; // ending position of alignment in read a and read b
    double identity; 
};

#endif