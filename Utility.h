//
// Created by Malysheva, Nadezhda on 10.07.20.
//

#ifndef ALGO_UTILITY_H
#define ALGO_UTILITY_H


#include <lemon/list_graph.h>
#include <random>
#include "ContactNetwork.h"

template <typename T> size_t binarySearch(std::vector<std::pair<double, T>> propCumSum,
        size_t indL, size_t indR, double rStart, double rBound)
{
    size_t result = size_t(-1);

    if (indR == indL)
    {
        if (propCumSum.at(indR).first + rStart  >= rBound)
        {
            result = indR;
        }

        return result;

    }
    else if (indR > indL)
    {
        int mid = indL + (indR - indL) / 2;

        // If the element is present at the middle
        // itself
        if (propCumSum.at(mid).first + rStart < rBound)
        {
            return binarySearch(propCumSum, mid + 1, indR, rStart, rBound);
        }

            // If element is smaller than mid, then
            // it can only be present in left subarray
        else
        {
            return binarySearch(propCumSum, indL, mid, rStart, rBound);
        }

        // Else the element can only be present
        // in right subarray
        //return binarySearch(propCumSum, mid + 1, indR, rBound);
    }

    // We reach here when element is not
    // present in array

}



void printBenFile(std::string fileName, const std::vector<BenStructure> &benToFile);

double sampleRandUni(std::mt19937_64 &generator);
#endif //ALGO_UTILITY_H
