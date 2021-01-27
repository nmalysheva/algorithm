/* Created by Malysheva, Nadezhda on 2019-07-28.
 *
 *
 * Gillespie algorithm (SSA)
 * TODO:: make it STATIC
 */


#ifndef ALGO_SSA_H
#define ALGO_SSA_H



#include "ContactNetwork.h"
#include "types.h"
#include <random>
#include <vector>

class SSA
{
public:
    SSA();
    void exe(); //
    void execute(double tStart, double tEnd, ContactNetwork &contNetwork,
                 NetworkStorage &nwStorage,
                 /*std::vector<double> &tSteps,*/ std::vector<double> &tInfect,
                 /*std::unordered_map<Specie::State, std::vector<uint32_t>> &populationState,*/
                 std::vector<uint32_t> &numberOfTransmitEdges,
                 //std::vector<std::vector<size_t>> &degreeDistr,
                 const std::string &saveDegreeDist/*, std::vector<BenStructure> &benToFile*/);
    ~SSA() {};

private:

    void   executeReaction(ContactNetwork & contNetwork, const std::string &reactId, double rStart,
                           double rBound, double time, uint32_t &nInf,
                            std::vector<std::pair<double, lemon::ListGraph::Edge>> &propDel,
                            std::vector<std::pair<double, lemon::ListGraph::Edge>> &propAdd,
                            std::vector<std::pair<double, lemon::ListGraph::Edge>> &propTransmit,
                            std::vector<std::pair<double, lemon::ListGraph::Node>> &propDiagnos,
                            std::vector<std::pair<double, lemon::ListGraph::Node>> &propDeath,
                            /*std::vector<double> &tSteps, */ std::vector<double> &tInfect,
                           //std::unordered_map<Specie::State, std::vector<uint32_t>> &populationState,
                           NetworkStorage &nwStorage,
                           std::vector<uint32_t> &numberOfTransmitEdges,
                            /*std::vector<std::vector<size_t>> &degreeDistr,*/ const std::string &saveDegreeDist
                           /*, std::vector<BenStructure> &benToFile*/);
    double recycleRandUni();

private:
    std::random_device rDev;
    std::mt19937_64 generator;
    std::uniform_real_distribution<> randuni;
};


#endif //ALGO_SSA_H
