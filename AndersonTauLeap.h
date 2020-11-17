//
// Created by Malysheva, Nadezhda on 21.07.20.
//

#ifndef ALGO_ANDERSONTAULEAP_H
#define ALGO_ANDERSONTAULEAP_H

#endif //ALGO_ANDERSONTAULEAP_H

#include "ContactNetwork.h"

void AndersonTauLeap(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork, double epsilon,
                    std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                    bool updateDegreeDistr, std::mt19937_64 &generator/*, std::vector<BenStructure> &benToFile*/);

void updateNetwork(std::vector<BenStructure> &benToFile, std::vector<int> k, int nDel, std::mt19937_64 &generator,
                    std::vector<std::pair<double, lemon::ListGraph::Edge>> &propAdd,
                    std::vector<std::pair<double, lemon::ListGraph::Edge>> &propDel, double t,
                    ContactNetwork & contNetwork,
                    std::vector<double> &props);

void updateDegreeDistributio(bool updateDegreeDistr, double t, std::vector<double> &timeSteps,
                              std::vector<std::vector<size_t>> &degreeDistr, const ContactNetwork  &contNetwork);

/*void exactAlgorithm(size_t n, double tEnd, ContactNetwork & contNetwork, double &t,
                    double &tLastNetworkUpdate, std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                    bool updateDegreeDistr,  std::mt19937_64 &generator,
                    std::vector<double> &T, std::vector<double> &C,
                    std::vector<std::vector<std::pair<double, int>>> &S);*/

void executeSSA(size_t n, double tEnd, ContactNetwork & contNetwork, double &t,
                double &tLastNetworkUpdate, std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                bool updateDegreeDistr, std::mt19937_64 &generator,
                std::vector<double> &T, std::vector<double> &C,
                std::vector<std::vector<std::pair<double, int>>> &S);

double getTau(size_t nParts, std::vector<std::vector<int>> nu, std::vector<double> props,  double epsilon,
              std::vector<size_t> X);