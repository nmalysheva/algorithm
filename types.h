//
// Created by Malysheva, Nadezhda on 26.01.21.
//

#ifndef ALGO_TYPES_H
#define ALGO_TYPES_H

#include "Specie.h"
#include <vector>

struct specieState
{
    int id;
    Specie::State state;
    std::vector<int> contacts;
};
using NetworkStorage = std::vector<std::pair<double, std::vector<specieState>>>;
#endif //ALGO_TYPES_H
