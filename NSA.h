//
// Created by Malysheva, Nadezhda on 2019-07-30.
//

#ifndef ALGO_NSA_H
#define ALGO_NSA_H

#include "ContactNetwork.h"
#include <random>


class NSA
{
public:
    NSA();
    void execute(double tStart, double tEnd, ContactNetwork &contNetwork, std::vector<double> &tSteps,
            std::vector<uint32_t> &nInfected, std::vector<std::vector<size_t>> &degreeDistr,
            double epsilon, size_t &nRejections, size_t &nAcceptance, size_t &nThin/*,
                 std::vector<BenStructure> &benToFile*/);
    ~NSA() {};

private:

    double  getPropUpperLimit (double lookAheadTime, ContactNetwork & contNetwork) const;
    double  recycleRandUni(double r); //recycle random number so we dont need to sample it again
    double sampleRandUni();
private:
    void executeReaction(ContactNetwork & contNetwork, const std::string & reactId,
                              double rStart, double rBound, double time, uint32_t &nInf/*, std::vector<BenStructure> &benToFile*/);

    std::random_device rDev;
    std::mt19937_64 generator;
    std::uniform_real_distribution<> randuni;
};


#endif //ALGO_NSA_H
