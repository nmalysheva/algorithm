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
                 std::vector<double> &tInfect, std::unordered_map<Specie::State, std::vector<uint32_t>> &populationState,
                 std::vector<uint32_t> &numberOfTransmitEdges,
                 std::vector<std::vector<size_t>> &degreeDistr, const std::string &saveDegreeDistMode,
            double epsilon, size_t &nRejections, size_t &nAcceptance, size_t &nThin/*,
                 std::vector<BenStructure> &benToFile*/);
    ~NSA() {};

private:

    double  getPropUpperLimit (double lookAheadTime, ContactNetwork & contNetwork,
                               double dignosisUpperLimit, double deathUpperLimit) const;

    double sampleRandUni();
private:
    void executeReaction(ContactNetwork & contNetwork, const std::string &reactId, double rStart,
                         double rBound, double time, uint32_t &nInf,
                         std::vector<std::pair<double, lemon::ListGraph::Edge>> &propTransmit,
                         std::vector<std::pair<double, lemon::ListGraph::Node>> &propDiagnos,
                         std::vector<std::pair<double, lemon::ListGraph::Node>> &propDeath,
                         std::vector<double> &tSteps, std::vector<double> &tInfect,
                         std::unordered_map<Specie::State, std::vector<uint32_t>> &populationState,
                         std::vector<uint32_t> &numberOfTransmitEdges,
                         std::vector<std::vector<size_t>> &degreeDistr, const std::string &saveDegreeDistMode);

    std::random_device rDev;
    std::mt19937_64 generator;
    std::uniform_real_distribution<> randuni;
};


#endif //ALGO_NSA_H
