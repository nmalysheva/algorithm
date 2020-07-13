//
// Created by Malysheva, Nadezhda on 13.07.20.
//

#ifndef ALGO_POISSONTAULEAP_H
#define ALGO_POISSONTAULEAP_H

#include "ContactNetwork.h"

void PoissonTauleap(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork, double epsilon,
                    std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                    bool updateDegreeDistr, std::mt19937_64 generator/*, std::vector<BenStructure> &benToFile*/);

void selectTimeStepAndK(double &tau, const std::unordered_map<std::string, double> &propensities, const std::vector<int> &k,
                        size_t &kDel, size_t &kAdd, double tau1, double tau2, double aCrit, std::mt19937_64 generator);

double proposeTau1(size_t lDel, size_t lAdd, size_t nAdd, double epsilon, std::vector<double> mu, std::vector<double> sigmaSq);
double proposeTau2(double aCrit, std::mt19937_64 generator);

void executeSSA(size_t n, double tEnd, ContactNetwork & contNetwork, double &t,
                double &tLastNetworkUpdate, std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                bool updateDegreeDistr, std::mt19937_64 generator);

#endif //ALGO_POISSONTAULEAP_H
