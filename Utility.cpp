//
// Created by Malysheva, Nadezhda on 10.07.20.
//
#include "Utility.h"

lemon::ListGraph::Edge binarySearch(std::vector<std::pair<double, lemon::ListGraph::Edge>> propCumSum,
                                         size_t indL, size_t indR, double rStart, double rBound)
{
    lemon::ListGraph::Edge result(lemon::INVALID);

    if (indR == indL)
    {
        return propCumSum.at(indR).second;
    }
    else if (indR >= indL)
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
    //return result;
}
