//
// Created by Malysheva, Nadezhda on 10.07.20.
//

#ifndef ALGO_UTILITY_H
#define ALGO_UTILITY_H


#include <lemon/list_graph.h>
#include <random>
#include "ContactNetwork.h"

size_t binarySearch(std::vector<std::pair<double, lemon::ListGraph::Edge>> propCumSum,
        size_t indL, size_t indR, double rStart, double rBound);
void printBenFile(std::string fileName, const std::vector<BenStructure> &benToFile);

double sampleRandUni(std::mt19937_64 &generator);
#endif //ALGO_UTILITY_H
