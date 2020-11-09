//
// Created by Malysheva, Nadezhda on 10.07.20.
//
#include "Utility.h"
#include <fstream>


double sampleRandUni(std::mt19937_64 &generator)
{
    std::uniform_real_distribution<> randuni;
    double r = randuni(generator);
    while (r == 0)
    {
        r = randuni(generator);
    }
    return r;
}

void printBenFile(std::string fileName, const std::vector<BenStructure> &benToFile)
{
    std::ofstream benFile;
    benFile.open(fileName);
    for (auto &it: benToFile)
    {
        std::string stateStr = "False";
        if (it.state)
        {
            stateStr = "True";
        }
        benFile << it.t << " " << it.u << " " << it.v << " " << stateStr << std::endl;
    }
    benFile.close();

}
