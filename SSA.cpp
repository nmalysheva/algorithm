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


void SSA::exe()
{
    std::cout << "pybind works";
}

void SSA::execute(double tStart, double tEnd, ContactNetwork &contNetwork,
                  std::vector<double> &tSteps, std::vector<uint32_t> &nInfected,
                  std::vector<std::vector<size_t>> &degreeDistr/*, std::vector<BenStructure> &benToFile*/)
{
    size_t cntDel = 0;
    size_t cntAdd = 0;
    std::vector<BenStructure> benToFile = contNetwork.getBenStructure(0);

    uint32_t  nInf = contNetwork.countByState(Specie::State::I);
    double time = tStart;
    tSteps.push_back(time);
    nInfected.push_back(nInf);
    degreeDistr.push_back(contNetwork.getDegreeDistribution());

    std::vector<std::pair<double, lemon::ListGraph::Edge>> propDel;
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propAdd;
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propTransmit;
    std::vector<std::pair<double, lemon::ListGraph::Node>> propDiagnos;
    std::vector<std::pair<double, lemon::ListGraph::Node>> propDeath;

    std::unordered_map<std::string, double >propensities {
            {"edge_del", 0},
            {"edge_add", 0},
            {"transmission", 0},
            {"diagnosis", 0},
            {"death", 0},
            {"birth", 0},
    };
    while (time < tEnd)
    {
        propDel = contNetwork.getEdgeDeletionRateSum();
        propAdd = contNetwork.getEdgeAdditionRateSum();
        propTransmit = contNetwork.getTransmissionRateSum();
        propDiagnos = contNetwork.getDiagnosisRateSum();
        propDeath = contNetwork.getDeathRateSum();

        propensities.at("edge_del") = propDel.at(propDel.size() - 1).first;
        propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;

        propensities.at("transmission") = propTransmit.at(propTransmit.size() - 1).first;
        propensities.at("diagnosis") = propDiagnos.at(propDiagnos.size() - 1).first;
        propensities.at("death") = propDeath.at(propDeath.size() - 1).first;

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
        double r = sampleRandUni(generator);
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
            r = sampleRandUni(generator);

            double pSum = 0;
            for (auto &it: propensities)
            {

                if (pSum + it.second >= propensitieSum * r)
                {
                    //std::cout << it.first << "; ";

                    if (it.first == "edge_del")
                    {
                        cntDel++;
                    }
                    if (it.first == "edge_add")
                    {
                        cntAdd++;
                    }
                    executeReaction(contNetwork, it.first, pSum, propensitieSum * r, time, nInf,
                                    propDel, propAdd, propTransmit, propDiagnos, propDeath);
                    //std::cout << std::endl;
                    tSteps.push_back(time);
                    nInfected.push_back(nInf);
                    degreeDistr.push_back(contNetwork.getDegreeDistribution());
                    break;
                }
                pSum += it.second;

            }

        }
    }
    std::cout << "del = " << cntDel << ", add = "<< cntAdd <<", edges = " << contNetwork.countEdges()<<std::endl;
}

void SSA::executeReaction(ContactNetwork & contNetwork, const std::string &reactId, double rStart,
                          double rBound, double time, uint32_t &nInf,
                           std::vector<std::pair<double, lemon::ListGraph::Edge>> &propDel,
                           std::vector<std::pair<double, lemon::ListGraph::Edge>> &propAdd,
                           std::vector<std::pair<double, lemon::ListGraph::Edge>> &propTransmit,
                           std::vector<std::pair<double, lemon::ListGraph::Node>> &propDiagnos,
                           std::vector<std::pair<double, lemon::ListGraph::Node>> &propDeath
                          /*,std::vector<BenStructure> &benToFile*/)
{
    if (reactId == "edge_del")
    {
        size_t index = binarySearch(propDel, 0, propDel.size() - 1, rStart, rBound);
        std::pair<int, int> b = contNetwork.removeEdge(propDel.at(index).second);

        //std::cout << b.first << ", " << b.second;
        //benToFile.emplace_back(time, b.first, b.second, false);

    }

    else if (reactId == "edge_add")
    {
        size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, rStart, rBound);
        std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);
        //benToFile.emplace_back(time, b.first, b.second, true);
        //std::cout << b.first << ", " << b.second;

    }
    else if (reactId == "transmission")
    {
        size_t index = binarySearch(propTransmit, 0, propTransmit.size() - 1, rStart, rBound);
        contNetwork.executeTransmission(propTransmit.at(index).second, time);
        nInf++;
    }

    else if (reactId == "diagnosis")
    {
        size_t index = binarySearch(propDiagnos, 0, propDiagnos.size() - 1, rStart, rBound);
        //std::pair<int, int> b = contNetwork.
        contNetwork.executeDiagnosis(propDiagnos.at(index).second, time);
    }

    else if (reactId == "death")
    {
        size_t index = binarySearch(propDeath, 0, propDeath.size() - 1, rStart, rBound);
        contNetwork.executeDeath(propDeath.at(index).second);
        nInf = contNetwork.countByState(Specie::State::I) + contNetwork.countByState(Specie::State::D);
    }

    else if (reactId == "birth")
    {
        //contNetwork.executeBirth(rStart, rBound);
    }
}