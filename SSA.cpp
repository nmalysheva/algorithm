//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include "SSA.h"
#include <string>
#include <iostream>
#include <unistd.h>
#include "Utility.h"

//for ben structure
#include <fstream>

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
                  NetworkStorage &nwStorage, std::vector<double> &tInfect,
                  std::vector<uint32_t> &numberOfTransmitEdges,
                  const std::string & saveDegreeDistMode)
{
    std::vector<BenStructure> benToFile = contNetwork.getBenStructure(0);

    double time = tStart;


    //used to store not ever contact update in "c" mode, but each contactSaveFactor-th
    size_t contactSaveFactor = 20;
    size_t contactSaveCounter = 0;

    /*if (saveDegreeDistMode == "c" || saveDegreeDistMode == "v")
    {
        nwStorage.emplace_back(time, contNetwork.getNetworkState());
    }*/

    uint32_t  nInf = contNetwork.countByState(Specie::State::I);
    /*if (saveDegreeDistMode == "v")
    {
        populationState.at(Specie::I).push_back(nInf);
        populationState.at(Specie::D).push_back(contNetwork.countByState(Specie::D));
        populationState.at(Specie::S).push_back(contNetwork.countByState(Specie::S));
        numberOfTransmitEdges.push_back(contNetwork.getTransmissionRateSum().size() - 1);
    }*/

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

            /*if (saveDegreeDistMode == "c" || saveDegreeDistMode == "v")
            {
                nwStorage.emplace_back(time, contNetwork.getNetworkState());
            }*/
            break;
        }
        double r = sampleRandUni(generator);
        double proposedTime = 1 / propensitieSum * std::log(1/r);
        if (time + proposedTime > tEnd )
        {
            time = tEnd;

            /*if (saveDegreeDistMode == "c" || saveDegreeDistMode == "v")
            {
                nwStorage.emplace_back(time, contNetwork.getNetworkState());
            }*/
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
                    executeReaction(contNetwork, it.first, pSum, propensitieSum * r, time, nInf,
                                    propDel, propAdd, propTransmit, propDiagnos, propDeath,
                                    tInfect,nwStorage,numberOfTransmitEdges, saveDegreeDistMode, contactSaveFactor, contactSaveCounter);
                    break;
                }
                pSum += it.second;

            }

        }
    }
}

void SSA::executeReaction(ContactNetwork & contNetwork, const std::string &reactId, double rStart,
                          double rBound, double time, uint32_t &nInf,
                           std::vector<std::pair<double, lemon::ListGraph::Edge>> &propDel,
                           std::vector<std::pair<double, lemon::ListGraph::Edge>> &propAdd,
                           std::vector<std::pair<double, lemon::ListGraph::Edge>> &propTransmit,
                           std::vector<std::pair<double, lemon::ListGraph::Node>> &propDiagnos,
                           std::vector<std::pair<double, lemon::ListGraph::Node>> &propDeath,
                            std::vector<double> &tInfect,
                          NetworkStorage &nwStorage,
                          std::vector<uint32_t> &numberOfTransmitEdges,
                          const std::string &saveDegreeDistMode,
                          const size_t contactSaveFactor,
                          size_t &contactSaveCounter
                           /*,std::vector<BenStructure> &benToFile*/)
{
    if (reactId == "edge_del")
    {
        size_t index = binarySearch(propDel, 0, propDel.size() - 1, rStart, rBound);
        std::pair<int, int> b = contNetwork.removeEdge(propDel.at(index).second);

        /*if (saveDegreeDistMode == "c")
        {
            contactSaveCounter ++;
            if (contactSaveCounter % contactSaveFactor == 0)
            {
                nwStorage.emplace_back(time, contNetwork.getNetworkState());
            }
        }*/


        //
        /*std::ofstream out;
        out.open("Ben_Cont_Dyn_SSA.txt", std::ios::app);
        out << time << " " << b.first << " " << b.second << " False" << std::endl ;
        out.close();*/

    }

    else if (reactId == "edge_add")
    {
        size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, rStart, rBound);
        std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);

        /*if (saveDegreeDistMode == "c")
        {
            contactSaveCounter ++;
            if (contactSaveCounter % contactSaveFactor == 0)
            {
                nwStorage.emplace_back(time, contNetwork.getNetworkState());
            }
        }*/

        /*std::ofstream out;
        out.open("Ben_Cont_Dyn_SSA.txt", std::ios::app);
        out << time << " " << b.first << " " << b.second << " True" << std::endl ;
        out.close();*/

    }
    else if (reactId == "transmission")
    {
        size_t index = binarySearch(propTransmit, 0, propTransmit.size() - 1, rStart, rBound);
        contNetwork.executeTransmission(propTransmit.at(index).second, time);
        nInf++;

        /*if (saveDegreeDistMode == "v")
        {
            nwStorage.emplace_back(time, contNetwork.getNetworkState());
            tInfect.push_back(time);
        }*/

    }

    else if (reactId == "diagnosis")
    {
        size_t index = binarySearch(propDiagnos, 0, propDiagnos.size() - 1, rStart, rBound);
        contNetwork.executeDiagnosis(propDiagnos.at(index).second, time);
        /*if (saveDegreeDistMode == "v")
        {
            nwStorage.emplace_back(time, contNetwork.getNetworkState());
        }*/
    }

    else if (reactId == "death")
    {
        size_t index = binarySearch(propDeath, 0, propDeath.size() - 1, rStart, rBound);
        contNetwork.executeDeath(propDeath.at(index).second);
        nInf = contNetwork.countByState(Specie::State::I) + contNetwork.countByState(Specie::State::D);

        /*if (saveDegreeDistMode == "v")
        {
            nwStorage.emplace_back(time, contNetwork.getNetworkState());
        }*/
    }

    else if (reactId == "birth")
    {
        //contNetwork.executeBirth(rStart, rBound);
    }
}