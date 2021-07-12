//
// Created by Malysheva, Nadezhda on 2019-07-30.
//

#include "NSA.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <unistd.h>
#include "Utility.h"

#include "AndersonTauLeap.h"
NSA::NSA()
{
    generator = std::mt19937_64(rDev());
    generator.seed(::time(nullptr) * getpid()); //to change the seed for every run
    randuni  = std::uniform_real_distribution<>(0.0, 1.0);

}

void NSA::execute(double tStart, double tEnd, ContactNetwork &contNetwork, NetworkStorage &nwStorage,
                  std::vector<double> &tInfect, const std::string &saveDegreeDistMode,
                  double epsilon, size_t &nRejections, size_t &nAcceptance, size_t &nThin/*,
                  std::vector<BenStructure> &benToFile*/)
{
    uint32_t  nInf = contNetwork.countByState(Specie::State::I) + + contNetwork.countByState(Specie::State::D);
    double time = tStart;


   /* if (saveDegreeDistMode == "c" || saveDegreeDistMode == "v")
    {
        nwStorage.emplace_back(time, contNetwork.getNetworkState());
    }*/

    double lookAheadTime  =  0; //init look-ahead time
    double propUpperLimit = -1; //init upper limit for propensitie sum
    double networkLastUpdate = tStart;

    double proposedTime = -1;
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propTransmit;
    propTransmit = contNetwork.getTransmissionRateSum();

    std::vector<std::pair<double, lemon::ListGraph::Node>> propDeath;
    propDeath = contNetwork.getDeathRateSum();

    std::vector<std::pair<double, lemon::ListGraph::Node>> propDiagnos;
    propDiagnos = contNetwork.getDiagnosisRateSum();

    std::unordered_map<std::string, double >propensities{
            {"transmission", propTransmit.at(propTransmit.size() - 1).first},
            {"diagnosis",propDiagnos.at(propDiagnos.size() - 1).first},
            {"death", propDeath.at(propDeath.size() - 1).first},
            {"birth", contNetwork.getBirthRateSum()},
    };

    while (time < tEnd)
    {
        //choose look-ahead time
        lookAheadTime = tEnd - time;
        propUpperLimit = getPropUpperLimit(lookAheadTime, contNetwork,
                                           propDiagnos.at(propDiagnos.size() - 1).first,
                                           propDeath.at(propDeath.size() - 1).first);

        // if there is no transmition possible anymore
        if (propUpperLimit == 0)
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
            double r = sampleRandUni();

            proposedTime = 1 / propUpperLimit * std::log(1/r);
            //std::cout << "proposedTime = " << proposedTime << std::endl;
            if (proposedTime > lookAheadTime)
            {
                nRejections ++;
                time += lookAheadTime;
                /*if (saveDegreeDistMode == "c" || saveDegreeDistMode == "v")
                {
                    nwStorage.emplace_back(time, contNetwork.getNetworkState());
                }*/
            }
            else
            {
                time += proposedTime;

                double tmpUpd = networkLastUpdate;
                AndersonTauLeap(networkLastUpdate, time, contNetwork, epsilon, nwStorage, "v"/*, benToFile*/, generator);

                if (tmpUpd < networkLastUpdate)
                {
                    networkLastUpdate = time;
                }
                //networkLastUpdate = time;

                //if (tmpUpd < networkLastUpdate)
                //{
                    //std::cout << "time = " << time <<", lastUpdate = " << networkLastUpdate <<std::endl;
                    propTransmit = contNetwork.getTransmissionRateSum();
                    propensities.at("transmission") = propTransmit.at(propTransmit.size() - 1).first;
                    propDiagnos = contNetwork.getDiagnosisRateSum();
                    propensities.at("diagnosis") = propDiagnos.at(propDiagnos.size() - 1).first;
                    propDeath = contNetwork.getDeathRateSum();
                    propensities.at("death") = propDeath.at(propDeath.size() - 1).first;
                    propensities.at("birth") = contNetwork.getBirthRateSum();

                    /*if (saveDegreeDistMode == "c")
                    {
                        nwStorage.emplace_back(time, contNetwork.getNetworkState());
                    }*/
                //}

                double propensitieSum = 0;
                for (auto &it: propensities)
                {
                    propensitieSum += it.second;
                }

                r = sampleRandUni(); //TODO recycle random number, not sample it

                if (propensitieSum >= propUpperLimit * r)
                {
                    nAcceptance ++;
                    double pSum = 0;

                    for (auto &it: propensities)
                    {

                        if (pSum + it.second >= propUpperLimit * r)
                        {
                            if (it.first == "diagnosis")
                            {
                                //networkLastUpdate = time;
                                /*std::vector<double> T(2, 0);
                                std::vector<int> C(2, 0);

                                std::vector<std::vector<std::pair<double, int>>> S;
                                S.reserve(1e3 + 1);
                                std::vector<std::pair<double, int>> temp = {{0.0, 0}};
                                S.push_back(temp);
                                S.push_back(temp);
                                double t = time - proposedTime;
                                executeSSA(1e6, tEnd, contNetwork, t, networkLastUpdate, nwStorage,
                                           saveDegreeDistMode, generator, T, C, S);*/
                                networkLastUpdate = time;
                            }
                            executeReaction(contNetwork, it.first, pSum,propUpperLimit * r, time, nInf,
                                    propTransmit,propDiagnos,propDeath, tInfect, nwStorage, saveDegreeDistMode);

                            propTransmit = contNetwork.getTransmissionRateSum();
                            propensities.at("transmission") = propTransmit.at(propTransmit.size() - 1).first;

                            propDeath = contNetwork.getDeathRateSum();
                            propensities.at("death") = propDeath.at(propDeath.size() - 1).first;

                            propDiagnos = contNetwork.getDiagnosisRateSum();
                            propensities.at("diagnosis") = propDiagnos.at(propDiagnos.size() - 1).first;

                            propensities.at("birth") = contNetwork.getBirthRateSum();

                            break;
                        }
                        pSum += it.second;

                    }
                }
                else
                {
                   nThin++;
                   //std::cout << "thin" << std::endl;
                }

            }
        }
    }

   // std::cout << std::endl;
}

double  NSA::getPropUpperLimit (double lookAheadTime, ContactNetwork & contNetwork, double dignosisUpperLimit, double deathUpperLimit) const
{
    double maxContInf = contNetwork.getMaxContactsLimitOfInfected(lookAheadTime);
    //std::cout << "max.cont.inf: " << maxContInf <<  std::endl;
    double maxContSusc = contNetwork.getMaxContactsLimitOfSusceptible(lookAheadTime);
    //std::cout << "max.cont.susc: " << maxContSusc <<  std::endl;
   // std::cout << "n.inf : " << contNetwork.countByState(Specie::I) + contNetwork.countByState(Specie::D)  <<  std::endl;

    double rmc = std::min(maxContInf, maxContSusc);

    size_t tmp = contNetwork.countByState(Specie::I) + contNetwork.countByState(Specie::D);
    size_t tmp2 = tmp * contNetwork.countByState(Specie::S);
    rmc = std::min(rmc, static_cast<double>(tmp2));

    //std::cout << "rmc = " << rmc <<  std::endl;
    double result = rmc * contNetwork.getTransmissionRateLimit() + dignosisUpperLimit + deathUpperLimit;//+

    //temporary!! For test
     //result = static_cast<double>(tmp2) * contNetwork.getTransmissionRateLimit() + dignosisUpperLimit + deathUpperLimit;

    return result;
}

double NSA::sampleRandUni()
{
    double r = randuni(generator);
    while (r == 0)
    {
        r = randuni(generator);
    }
    return r;
}


void NSA::executeReaction(ContactNetwork & contNetwork, const std::string &reactId, double rStart,
                          double rBound, double time, uint32_t &nInf,
                          std::vector<std::pair<double, lemon::ListGraph::Edge>> &propTransmit,
                          std::vector<std::pair<double, lemon::ListGraph::Node>> &propDiagnos,
                          std::vector<std::pair<double, lemon::ListGraph::Node>> &propDeath,
                          std::vector<double> &tInfect,NetworkStorage &nwStorage, const std::string &saveDegreeDistMode)
{

    if (reactId == "transmission")
    {
        size_t index = binarySearch(propTransmit, 0, propTransmit.size() - 1, rStart, rBound);
        contNetwork.executeTransmission(propTransmit.at(index).second, time);
        nInf++;

       /* if (saveDegreeDistMode == "v")
        {
            nwStorage.emplace_back(time, contNetwork.getNetworkState());
            tInfect.push_back(time);
        }*/
    }
    else if (reactId == "diagnosis")
    {
        size_t index = binarySearch(propDiagnos, 0, propDiagnos.size() - 1, rStart, rBound);
        contNetwork.executeDiagnosis(propDiagnos.at(index).second, time);

       /* if (saveDegreeDistMode == "v")
        {
            nwStorage.emplace_back(time, contNetwork.getNetworkState());
        }*/
    }

    else if (reactId  == "death" )
    {
        size_t index = binarySearch(propDeath, 0, propDeath.size() - 1, rStart, rBound);
        contNetwork.executeDeath(propDeath.at(index).second);
        nInf = contNetwork.countByState(Specie::State::I) + contNetwork.countByState(Specie::State::D);

        /*if (saveDegreeDistMode == "v")
        {
            nwStorage.emplace_back(time, contNetwork.getNetworkState());
        }*/
    }
}


