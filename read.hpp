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
    int preserve = 0;
    std::vector<std::tuple<int, int, int>> long_repeats;

    Read(int id, std::string name, std::string bases) : id(id), bases(bases), name(name){};

    void showRead()
    {
        std::cout << "read " << id << "\n";
        std::cout << ">" << name << "\n";
        std::cout << bases << "\n";
    }
};

#endif