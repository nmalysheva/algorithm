//
// Created by Malysheva, Nadezhda on 13.07.20.
//

#ifndef ALGO_RKF45_H
#define ALGO_RKF45_H

#include "ContactNetwork.h"

void RKF45Approximation(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork,
                        double dtMax, double dtMin, double errorMax, double errorMin, std::vector<double> &timeSteps,
                        std::vector<std::vector<size_t>> &degreeDistr, bool updateDegreeDistr/*, std::vector<BenStructure> &*/, std::mt19937_64 generator);

double getNadd(double X, double Y, double lam, double theta, std::vector<std::vector<double>> a, std::vector<double> d, double dt);
double getNdel(double X, double Y, double lam, double theta, std::vector<std::vector<double>> a, std::vector<double> d, double dt);

#endif //ALGO_RKF45_H
