//
// Created by Malysheva, Nadezhda on 13.07.20.
//

#ifndef ALGO_POISSONTAULEAP_H
#define ALGO_POISSONTAULEAP_H

#include "ContactNetwork.h"

void PoissonTauleap(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork, double epsilon,
                    std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                    bool updateDegreeDistr, std::mt19937_64 &generator/*, std::vector<BenStructure> &benToFile*/);

void selectTimeStepAndK(double &tau, const std::unordered_map<std::string, double> &propensities, const std::vector<int> &k,
                        size_t &kDel, size_t &kAdd, double tau1, double tau2, double aCrit, std::mt19937_64 &generator);

double proposeTau1(size_t lDel, size_t lAdd, size_t nAdd, double epsilon, std::vector<double> mu, std::vector<double> sigmaSq);
double proposeTau2(double aCrit, std::mt19937_64 &generator);

void executeSSA(size_t n, double tEnd, ContactNetwork & contNetwork, double &t,
                double &tLastNetworkUpdate, std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                bool updateDegreeDistr, std::mt19937_64 &generator);
void updateDegreeDistribution(bool updateDegreeDistr, double t, std::vector<double> &timeSteps,
                              std::vector<std::vector<size_t>> &degreeDistr, /*const*/ ContactNetwork  &contNetwork);

int splitRandomNumber(int kDel, int kAdd, int &N, std::mt19937_64 &generator);

void updateNetwork(std::vector<BenStructure> &benToFile, std::vector<int> k, int nDel, std::mt19937_64 &generator,
                   std::vector<std::pair<double, lemon::ListGraph::Edge>> &propAdd,
                   std::vector<std::pair<double, lemon::ListGraph::Edge>> &propDel, double t,
                   ContactNetwork & contNetwork,
                   std::unordered_map<std::string, double> &propensities);

void updateNetwork2(std::vector<BenStructure> &benToFile, std::vector<int> k, int nDel, std::mt19937_64 &generator,
                   std::vector<std::pair<double, lemon::ListGraph::Edge>> &propAdd,
                   std::vector<std::pair<double, lemon::ListGraph::Edge>> &propDel, double t,
                   ContactNetwork & contNetwork,
                   std::unordered_map<std::string, double> &propensities,
                    size_t nnn);

#endif //ALGO_POISSONTAULEAP_H
