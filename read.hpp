#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include <stdio.h>

typedef std::pair<int, int> Interval;

enum AlignType
{
    REVERSE,
    FORWARD
};

class Read
{ // read class
public:

    int id;            // id, start from 0
    std::string name;  // read name
    std::string bases; // read bases
    int len;
    int start_pos;
    std::string align;
    int preserve = 0;

    Read(int id, int length, std::string name, std::string bases, int start_pos, std::string align) : id(id), bases(bases), name(name), len(length), start_pos(start_pos), align(align){};
    Read(int id, int length, std::string name, std::string bases) : id(id), bases(bases), name(name), len(length){};
    Read(int id, std::string name, std::string bases) : id(id), bases(bases), name(name){};

    bool active = true;

    void showRead()
    {
        std::cout << "read " << id << "\n";
        std::cout << ">" << name << "\n";
        std::cout << bases << "\n";
    }
};

#endif