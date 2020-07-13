//
// Created by Malysheva, Nadezhda on 13.07.20.
//

#ifndef ALGO_MIDPOINTAPPROXIMATION_H
#define ALGO_MIDPOINTAPPROXIMATION_H

#include "ContactNetwork.h"

void MidpointApproximation(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork,
                           double dtMax, double dtMin, double errorMax, double errorMin, std::vector<double> &timeSteps,
                           std::vector<std::vector<size_t>> &degreeDistr, bool updateDegreeDistr = true);

#endif //ALGO_MIDPOINTAPPROXIMATION_H
