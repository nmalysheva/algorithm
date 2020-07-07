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
    void PoissonTauleap(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork, double epsilon,
                        std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                        bool updateDegreeDistr/*, std::vector<BenStructure> &benToFile*/);
    void RKF45Approximation(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork,
                            double dtMax, double dtMin, double errorMax, double errorMin, std::vector<double> &timeSteps,
                            std::vector<std::vector<size_t>> &degreeDistr, bool updateDegreeDistr/*, std::vector<BenStructure> &*/);
    void MidpointApproximation(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork,
                            double dtMax, double dtMin, double errorMax, double errorMin, std::vector<double> &timeSteps,
                            std::vector<std::vector<size_t>> &degreeDistr, bool updateDegreeDistr = true);
    ~NSA() {};

private:

    double  getPropUpperLimit (double lookAheadTime, ContactNetwork & contNetwork) const;
    double  recycleRandUni(double r); //recycle random number so we dont need to sample it again
    void selectTimeStepAndK(double &tau, const std::unordered_map<std::string, double> &propensities, const std::vector<int> &k,
            size_t &kDel, size_t &kAdd, double tau1, double tau2, double aCrit);

    double proposeTau1(size_t lDel, size_t lAdd, size_t nAdd, double epsilon, std::vector<double> mu, std::vector<double> sigmaSq);
    double proposeTau2(double aCrit);
    double sampleRandUni();
    lemon::ListGraph::Edge binarySearch(std::vector<std::pair<double, lemon::ListGraph::Edge>> &propCumSum,
            size_t indL, size_t indR, double rStart, double rBound);
private:
    void executeReaction(ContactNetwork & contNetwork, std::string reactId,
                              double rStart, double rBound, double time, uint32_t &nInf/*, std::vector<BenStructure> &benToFile*/);

    double getNadd(double X, double Y, double lam, double theta, std::vector<std::vector<double>> a, std::vector<double> d, double dt);
    double getNdel(double X, double Y, double lam, double theta, std::vector<std::vector<double>> a, std::vector<double> d, double dt);
    std::random_device rDev;
    std::mt19937_64 generator;
    std::uniform_real_distribution<> randuni;
};


#endif //ALGO_NSA_H
