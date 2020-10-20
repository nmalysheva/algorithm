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
#include "PoissonTauLeap.h"

NSA::NSA()
{
    generator = std::mt19937_64(rDev());
    generator.seed(::time(nullptr) * getpid()); //to change the seed for every run
    randuni  = std::uniform_real_distribution<>(0.0, 1.0);

}

void NSA::execute(double tStart, double tEnd, ContactNetwork &contNetwork,
                  std::vector<double> &tSteps, std::vector<uint32_t> &nInfected, std::vector<std::vector<size_t>> &degreeDistr,
                  double epsilon, size_t &nRejections, size_t &nAcceptance, size_t &nThin/*,
                  std::vector<BenStructure> &benToFile*/)
{
    uint32_t  nInf = contNetwork.countByState(Specie::State::I);
    double time = tStart;

    std::vector<double> tUpds;
    tSteps.push_back(time);
    nInfected.push_back(nInf);
    degreeDistr.push_back(contNetwork.getDegreeDistribution());

    double lookAheadTime  =  0; //init look-ahead time
    double propUpperLimit = -1; //init upper limit for propensitie sum
    double networkLastUpdate = tStart;

    tUpds.push_back(networkLastUpdate);

    double proposedTime = -1;

    std::unordered_map<std::string, double >propensities{
            {"transmission", contNetwork.getTransmissionRateSum()},
            {"death", contNetwork.getDeathRateSum()},
            {"birth", contNetwork.getBirthRateSum()},
    };

    while (time < tEnd)
    {
        //choose look-ahead time
        lookAheadTime = tEnd - time;
        propUpperLimit = getPropUpperLimit(lookAheadTime, contNetwork);

        // if there is no transmition possible anymore
        if (propUpperLimit == 0)
        {
            std::cout << "B=0" << std::endl;
            time = tEnd;
            tSteps.push_back(time);
            nInfected.push_back(nInf);
            degreeDistr.push_back(contNetwork.getDegreeDistribution());
            break;
        }
        else
        {
            double r = sampleRandUni();

            proposedTime = 1 / propUpperLimit * std::log(1/r);
            if (proposedTime > lookAheadTime)
            {
                nRejections ++;
                //std::cout<<"reject" <<std::endl;
                time += lookAheadTime;
                tSteps.push_back(time);
                nInfected.push_back(nInf);
                degreeDistr.push_back(contNetwork.getDegreeDistribution());
            }
            else
            {
                time += proposedTime;

                //tUpds.push_back(networkLastUpdate);
                PoissonTauleap(networkLastUpdate, time, contNetwork, epsilon, tSteps, degreeDistr, false/*, benToFile*/, generator);
                //tUpds.push_back(networkLastUpdate);
                //RKF45Approximation(networkLastUpdate, time, contNetwork, (time - networkLastUpdate)/2, (time - networkLastUpdate) * 1e-4, 0.1, 1e-4, tSteps, degreeDistr, false);
                std::cout << "_______________-________________" <<std::endl;

                propensities.at("transmission") = contNetwork.getTransmissionRateSum();
                propensities.at("death") = contNetwork.getDeathRateSum();
                propensities.at("birth") = contNetwork.getBirthRateSum();

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
                            std::cout << "accepted" << std::endl;
                            executeReaction(contNetwork, it.first, pSum, propUpperLimit * r, time, nInf/*, benToFile*/);

                            if (it.first  == "death" )
                            {
                                nInf = contNetwork.countByState(Specie::State::I);
                            }
                            tSteps.push_back(time);
                            nInfected.push_back(nInf);
                            degreeDistr.push_back(contNetwork.getDegreeDistribution());

                            propensities.at("transmission") = contNetwork.getTransmissionRateSum();
                            propensities.at("death") = contNetwork.getDeathRateSum();
                            propensities.at("birth") = contNetwork.getBirthRateSum();

                            break;
                        }
                        pSum += it.second;

                    }
                }
                else
                {
                   nThin++;
                   std::cout << "thin" << std::endl;
                }

            }
        }
    }

    std::cout << "network updates: ";
    for (auto &it: tUpds)
    {
        std::cout << it << ", ";
    }
    std::cout << std::endl;
}

double  NSA::getPropUpperLimit (double lookAheadTime, ContactNetwork & contNetwork) const
{
    /*double result = contNetwork.getMaxContactsLimitOfInfected()  * contNetwork.getTransmissionRateLimit() +
              contNetwork.getDeathRateSum() + contNetwork.getBirthRateSum();*/

    /*size_t edgesLimit = contNetwork.countByState(Specie::I) * 3.5 * lookAheadTime + contNetwork.getNumberContactsOfInfected();
    edgesLimit = std::min(edgesLimit, contNetwork.getMaxContactsLimitOfInfected());
    double result = contNetwork.getDeathRateSum() + contNetwork.getBirthRateSum() + edgesLimit * contNetwork.getTransmissionRateLimit();
     */
    //std::cout << "-----------------------------" <<  std::endl;

    double maxContInf = contNetwork.getMaxContactsLimitOfInfected(lookAheadTime);
    //std::cout << "max.cont.inf: " << maxContInf <<  std::endl;
    //std::cout << "****************************************" <<  std::endl;
    double maxContSusc = contNetwork.getMaxContactsLimitOfSusceptible(lookAheadTime);
    //std::cout << "****************************************" <<  std::endl;
    //std::cout << "max.cont.susc: " << maxContSusc <<  std::endl;
    //std::cout << "death rate: " << contNetwork.getDeathRateSum() <<  std::endl;
    //std::cout << "birth rate: " << contNetwork.getBirthRateSum() <<  std::endl;
    std::cout << "n.inf : " << contNetwork.countByState(Specie::I) <<"; "; //<<  std::endl;

    double rmc = std::min(maxContInf, maxContSusc);
    //std::cout << "rmc : " << rmc <<  std::endl;
    double result = rmc * contNetwork.getTransmissionRateLimit() +
                    contNetwork.getDeathRateSum() + contNetwork.getBirthRateSum();
    //std::cout << "-----------------------------" <<  std::endl;
    return result;
}

double  NSA::recycleRandUni(double r)
{
    //TODO proper recycling
    return r;
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


void NSA::executeReaction(ContactNetwork & contNetwork, const std::string &reactId,
                          double rStart, double rBound, double time, uint32_t &nInf/*, std::vector<BenStructure> &benToFile*/)
{
    //std::cout << reactId <<std::endl;
    if (reactId == "edge_del")
    {
        /*BenStructure b (time, -1, -1, false);
        contNetwork.executeEdgeDeletion(rStart, rBound, b);
        benToFile.push_back(b);*/

        contNetwork.executeEdgeDeletion(rStart, rBound);
    }

    else if (reactId == "edge_add")
    {
        /*BenStructure b (time, -1, -1, true);
        contNetwork.executeEdgeAddition(rStart, rBound, b);
        benToFile.push_back(b);*/

        contNetwork.executeEdgeAddition(rStart, rBound);
    }

    else if (reactId == "transmission")
    {
        contNetwork.executeTransmission(rStart, rBound, time);
        nInf++;
    }

    else if (reactId == "death")
    {
        contNetwork.executeDeath(rStart, rBound);
    }

    else if (reactId == "birth")
    {
        contNetwork.executeBirth(rStart, rBound);
    }
}


