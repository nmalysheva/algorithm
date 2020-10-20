//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include "SSA.h"
#include <string>
#include <iostream>
#include <unistd.h>
#include "Utility.h"

SSA::SSA()
{
    generator = std::mt19937_64(rDev());
    generator.seed(::time(nullptr) * getpid()); //to change the seed for every run
    randuni  = std::uniform_real_distribution<>(0.0, 1.0);

}

double SSA::sampleRandUni()
{
    double r = randuni(generator);
    while (r == 0)
    {
        r = randuni(generator);
    }
    return r;
}

void SSA::exe()
{
    std::cout << "pybind works";
}

void SSA::execute(double tStart, double tEnd, ContactNetwork &contNetwork,
                  std::vector<double> &tSteps, std::vector<uint32_t> &nInfected,
                  std::vector<std::vector<size_t>> &degreeDistr/*, std::vector<BenStructure> &benToFile*/)
{
    std::vector<BenStructure> benToFile = contNetwork.getBenStructure(0);

    uint32_t  nInf = contNetwork.countByState(Specie::State::I);
    double time = tStart;
    tSteps.push_back(time);
    nInfected.push_back(nInf);
    degreeDistr.push_back(contNetwork.getDegreeDistribution());

    std::vector<std::pair<double, lemon::ListGraph::Edge>> propDel;
    propDel.reserve(1e6 + 1);
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propAdd;
    propAdd.reserve(1e6 + 1);

    std::unordered_map<std::string, double >propensities {
            {"edge_del", 0},
            {"edge_add", 0},
            {"transmission", 0},
            {"death", 0},
            {"birth", 0},
    };
    while (time < tEnd)
    {
        propDel = contNetwork.getEdgeDeletionRateSum();
        propAdd = contNetwork.getEdgeAdditionRateSum();

        propensities.at("edge_del") = propDel.at(propDel.size() - 1).first;
        propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;

        propensities.at("transmission") = contNetwork.getTransmissionRateSum();
        propensities.at("death") = contNetwork.getDeathRateSum();
        propensities.at("birth") = contNetwork.getBirthRateSum();

        double propensitieSum = 0;
        for (auto &it: propensities)
        {
            propensitieSum += it.second;
        }

        if (propensitieSum == 0)
        {
            time = tEnd;

            tSteps.push_back(time);
            nInfected.push_back(nInf);
            degreeDistr.push_back(contNetwork.getDegreeDistribution());

            break;
        }
        double r = sampleRandUni();
        double proposedTime = 1 / propensitieSum * std::log(1/r);
        if (time + proposedTime > tEnd )
        {
            time = tEnd;

            tSteps.push_back(time);
            nInfected.push_back(nInf);
            degreeDistr.push_back(contNetwork.getDegreeDistribution());

            break;
        }
        else
        {
            time += proposedTime;
            r = sampleRandUni();

            double pSum = 0;
            for (auto &it: propensities)
            {

                if (pSum + it.second >= propensitieSum * r)
                {
                   // executeReaction(contNetwork, it.first, pSum, propensitieSum * r, time, nInf/*, benToFile*/);

                    /// TEMPORAL SOLUTION

                    if (it.first == "edge_del")
                    {
                        size_t index = binarySearch(propDel, 0, propDel.size() - 1, pSum, propensitieSum * r);
                        std::pair<int, int> b = contNetwork.removeEdge(propDel.at(index).second);
                        //propDel.erase(propDel.begin() + index);
                        benToFile.emplace_back(time, b.first, b.second, false);
                    }
                    else if (it.first == "edge_add")
                    {
                        size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, pSum, propensitieSum * r);
                        std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);
                        //propAdd.erase(propAdd.begin() + index);
                        benToFile.emplace_back(time, b.first, b.second, true);
                    }

                    /// TEMPORAL SOLUTION
                    tSteps.push_back(time);
                    nInfected.push_back(nInf);
                    degreeDistr.push_back(contNetwork.getDegreeDistribution());

                    break;
                }
                pSum += it.second;

            }

        }
    }
}

void SSA::executeReaction(ContactNetwork & contNetwork, const std::string &reactId,
                          double rStart, double rBound, double time, uint32_t &nInf/*,
                          std::vector<BenStructure> &benToFile*/)
{
    if (reactId == "edge_del")
    {
        /*BenStructure b (time, -1, -1, false);
        contNetwork.executeEdgeDeletion(rStart, rBound, b);
        benToFile.push_back(b);*/

        //contNetwork.executeEdgeDeletion(rStart, rBound);

    }

    else if (reactId == "edge_add")
    {
        /*BenStructure b (time, -1, -1, true);
        contNetwork.executeEdgeAddition(rStart, rBound, b);
        benToFile.push_back(b);*/

        //contNetwork.executeEdgeAddition(rStart, rBound);
        //std::cout  << "edge_add "<<std::endl;
    }
    else if (reactId == "transmission")
    {
        contNetwork.executeTransmission(rStart, rBound, time);
        nInf++;
        //std::cout << "transmission " << time << " " << contNetwork.countByState(Specie::I) << " " << contNetwork.size()<<std::endl;
    }

    else if (reactId == "death")
    {
        contNetwork.executeDeath(rStart, rBound);
        nInf = contNetwork.countByState(Specie::State::I);
        //std::cout  << "death " << time << " " << contNetwork.countByState(Specie::I) << " " << contNetwork.size()  <<std::endl;
    }

    else if (reactId == "birth")
    {
        contNetwork.executeBirth(rStart, rBound);
        //std::cout << "birth " << time << " " << contNetwork.countByState(Specie::I)  << " " << contNetwork.size() <<std::endl;
    }
    //std::cout << "------------------" <<std::endl;
}