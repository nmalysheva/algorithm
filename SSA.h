//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#ifndef ALGO_SSA_H
#define ALGO_SSA_H



#include "ContactNetwork.h"
#include <random>
#include <vector>

class SSA
{
public:
    SSA();
    void exe();
    //time assumed to be in years
    void execute(double tStart, double tEnd, ContactNetwork &contNetwork,
                 std::vector<double> &tSteps, std::vector<uint32_t> &nInfected,
                 std::vector<std::vector<size_t>> &degreeDistr/*, std::vector<BenStructure> &benToFile*/);
    ~SSA() {};

private:

    double sampleRandUni();
    void   executeReaction(ContactNetwork & contNetwork, std::string reactId, double rStart,
                           double rBound, double time, uint32_t &nInf/*, std::vector<BenStructure> &benToFile*/);
    double recycleRandUni();

    lemon::ListGraph::Edge binarySearch(std::vector<std::pair<double, lemon::ListGraph::Edge>> propCumSum, size_t indL, size_t indR, double rStart, double rBound);

private:
    std::random_device rDev;
    std::mt19937_64 generator;
    std::uniform_real_distribution<> randuni;
};


#endif //ALGO_SSA_H
