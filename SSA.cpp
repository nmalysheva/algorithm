//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include "SSA.h"
#include <string>
#include <iostream>

SSA::SSA()
{
    generator = std::mt19937_64(rDev());
    generator.seed(::time(NULL)); //to change the seed for every run
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
                  std::vector<double> &tSteps, std::vector<uint32_t> &nInfected, std::vector<std::vector<size_t>> &degreeDistr)
{
    //time_t programStart, programEnd;
    //std::time(&programStart);
    uint32_t  nInf = contNetwork.countByState(Specie::State::I);
    double time = tStart;
    tSteps.push_back(time);
    nInfected.push_back(nInf);
    degreeDistr.push_back(contNetwork.getDegreeDistribution());

    size_t nDel = 0;
    size_t nAdd = 0;

    std::unordered_map<std::string, double >propensities {
            {"edge_del", contNetwork.getEdgeDeletionRateSum(nDel)},
            {"edge_add", contNetwork.getEdgeAdditionRateSum(nAdd)},
            {"transmission", contNetwork.getTransmissionRateSum()},
            {"death", contNetwork.getDeathRateSum()},
            {"birth", contNetwork.getBirthRateSum()},
    };
    while (time < tEnd)
    {
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
                    executeReaction(contNetwork, it.first, pSum, propensitieSum * r, time, nInf);

                    /*if (it.first == "transmission")
                    {*/
                        tSteps.push_back(time);
                        nInfected.push_back(nInf);
                        degreeDistr.push_back(contNetwork.getDegreeDistribution());


                    //}

                    /*if (it.first  == "death" )
                    {
                        size_t newNInf = contNetwork.countByState(Specie::State::I);
                        if (newNInf != nInf)
                        {
                            nInf = newNInf;
                            tSteps.push_back(time);
                            nInfected.push_back(nInf);
                            degreeDistr.push_back(contNetwork.getDegreeDistribution());
                        }
                    }*/

                    break;
                }
                pSum += it.second;

            }

        }
        //update propensities
        propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum(nDel);
        propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum(nAdd);
        propensities.at("transmission") = contNetwork.getTransmissionRateSum();
        propensities.at("death") = contNetwork.getDeathRateSum();
        propensities.at("birth") = contNetwork.getBirthRateSum();

    }

    //std::time(&programEnd);
    //std::cout << "time of execution: " << programEnd - programStart << std::endl;

}

void SSA::executeReaction(ContactNetwork & contNetwork, std::string reactId,
                          double rStart, double rBound, double time, uint32_t &nInf)
{
    //std::cout << reactId <<std::endl;
    if (reactId == "edge_del")
    {
        contNetwork.executeEdgeDeletion(rStart, rBound);
       // std::cout  << "edge_del "<<std::endl;
    }

    else if (reactId == "edge_add")
    {
        contNetwork.executeEdgeAddition(rStart, rBound);
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
        //std::cout  << "death " << time << " " << contNetwork.countByState(Specie::I) << " " << contNetwork.size()  <<std::endl;
    }

    else if (reactId == "birth")
    {
        contNetwork.executeBirth(rStart, rBound);
        //std::cout << "birth " << time << " " << contNetwork.countByState(Specie::I)  << " " << contNetwork.size() <<std::endl;
    }
    //std::cout << "------------------" <<std::endl;
}