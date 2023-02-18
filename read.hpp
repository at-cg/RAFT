#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include <stdio.h>

typedef std::pair<int, int> Interval;

class Read
{ // read class
public:

    int id;            // id, start from 0
    std::string name;  // read name
    std::string bases; // read bases
    int len;
    int start_pos;
    int end_pos;
    std::string align;
    std::string chr;
    int preserve = 0;

    Read(int id, int length, std::string name, std::string bases, int start_pos, int end_pos, std::string align, std::string chr) : id(id), bases(bases), name(name), len(length), start_pos(start_pos), end_pos(end_pos), align(align), chr(chr){};
    Read(int id, int length, std::string name, std::string bases) : id(id), bases(bases), name(name), len(length){};


    void showRead()
    {
        std::cout << "read " << id << "\n";
        std::cout << ">" << name << "\n";
        std::cout << bases << "\n";
    }
};

#endif