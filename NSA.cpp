//
// Created by Malysheva, Nadezhda on 2019-07-30.
//

#include "NSA.h"
#include "math.h"
#include <algorithm>

NSA::NSA()
{
    generator = std::mt19937_64(rDev());
    generator.seed(::time(NULL)); //to change the seed for every run
    randuni  = std::uniform_real_distribution<>(0.0, 1.0);

}

void NSA::execute(double tStart, double tEnd, ContactNetwork &contNetwork,
                  std::vector<double> &tSteps, std::vector<uint32_t> &nInfected, std::vector<std::vector<size_t>> &degreeDistr,
                  double epsilon, size_t &nRejections, size_t &nAcceptance, size_t &nThin)
{
    //time_t programStart, programEnd;
    //std::time(&programStart);
    uint32_t  nInf = contNetwork.countByState(Specie::State::I);
    double time = tStart;
    tSteps.push_back(time);
    nInfected.push_back(nInf);
    degreeDistr.push_back(contNetwork.getDegreeDistribution());

    double lookAheadTime  =  0; //init look-ahead time
    double propUpperLimit = -1; //init upper limit for propensitie sum
    double networkLastUpdate = tStart;

    double proposedTime = -1;

    std::unordered_map<std::string, double >propensities {
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
                //std::cout<<"rehject" <<std::endl;
                time += lookAheadTime;
                tSteps.push_back(time);
                nInfected.push_back(nInf);
                degreeDistr.push_back(contNetwork.getDegreeDistribution());
            }
            else
            {
                time += proposedTime;

                updateContacts(networkLastUpdate, time, contNetwork, epsilon);
                //networkLastUpdate = time;

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
                            executeReaction(contNetwork, it.first, pSum, propUpperLimit * r, time, nInf);

                            if (it.first == "transmission")
                            {
                                tSteps.push_back(time);
                                nInfected.push_back(nInf);
                                degreeDistr.push_back(contNetwork.getDegreeDistribution());
                            }

                            if (it.first  == "death" )
                            {
                                size_t newNInf = contNetwork.countByState(Specie::State::I);
                                if (newNInf != nInf)
                                {
                                    nInf = newNInf;
                                    tSteps.push_back(time);
                                    nInfected.push_back(nInf);
                                    degreeDistr.push_back(contNetwork.getDegreeDistribution());
                                }
                            }

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
                }

            }

        }

    }
    //std::time(&programEnd);
    //std::cout << "time of execution: " << programEnd - programStart << std::endl;

}

double  NSA::getPropUpperLimit (double lookAheadTime, ContactNetwork & contNetwork) const
{
    double result = contNetwork.getMaxContactsLimitOfInfected() * lookAheadTime * contNetwork.getTransmissionRateLimit() +
              contNetwork.getDeathRateSum() + contNetwork.getBirthRateSum();
    return result;
}

double  NSA::recycleRandUni(double r)
{
    //TODO proper recycling
    return r;
}


void NSA::updateContacts(double &tStart, double tEnd, ContactNetwork & contNetwork, double epsilon)
{// here we consider only two reactions: assembling and disassembling of contact (edge)
    // only deletion can be critical since it can lead to the negative result
    /*std::vector<int> nu = {-1, 1}; // -1 for edge disappearance, 1 for adding
    double time = tStart;

    std::vector<double> propensities(nu.size());*/

    double time = tStart;
    std::unordered_map<std::string, double >propensities {
            {"edge_del", contNetwork.getEdgeDeletionRateSum()},
            {"edge_add", contNetwork.getEdgeAdditionRateSum()}
    };

    //std:: cout << "add: "<<propensities.at("edge_add") <<std::endl;

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


            break;
        }

        double tau = proposeTimestep(epsilon, contNetwork);

        if (time + tau > tEnd)
        {
            //time = tEnd;
            tau = tEnd - time;
            //std:: cout <<"tau="<< tau <<" time="<< time <<" tend="<< tEnd<<std::endl;
            //break;
        }

        time += tau;

        tStart = time;
        std::unordered_map<std::string, size_t >k {
                {"edge_del", 0},
                {"edge_add", 0}
        };

        std::unordered_map<std::string, size_t >kMax {
                {"edge_del", contNetwork.countEdges()},
                {"edge_add", contNetwork.getAmountOfEdgesPossibleToAdd()}
        };

        for (auto &it: propensities )
        {
            double p = 1;

            if (it.second * tau <= kMax.at(it.first))
            {
                p = it.second * tau / kMax.at(it.first);
            }

            k.at(it.first) = sampleRandBinomial(kMax.at(it.first), p);
        }

        std::vector<std::string> order;

        //size_t i = 0;
        for (auto &it: k )
        {
           // i++;
            order.insert(order.end(), it.second, it.first);
        }

        std::shuffle(order.begin(), order.end(), generator);
        for (size_t i = 0; i < order.size(); i++)
        {
            double r = sampleRandUni();

            uint32_t  ni = 0;
            executeReaction(contNetwork, order.at(i), 0, r * propensities.at(order.at(i)), time, ni);

            propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum();
            propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum();
        }
        //std:: cout <<"---------------"<<std::endl;

    }
}

double  NSA::proposeTimestep(double epsilon, ContactNetwork & contNetwork) const
{

    //std:: cout <<"edges="<< contNetwork.countEdges() <<" delrsum="<< contNetwork.getEdgeDeletionRateSum()<<std::endl;

    size_t nEdgesExist = contNetwork.countEdges();
    size_t nEdgesAdd = contNetwork.getAmountOfEdgesPossibleToAdd();
    double tau = std::min(epsilon * nEdgesExist /contNetwork.getEdgeDeletionRateSum(),
                          epsilon * nEdgesAdd / contNetwork.getEdgeAdditionRateSum());
    std::cout << "edge del.r.s.: " << contNetwork.getEdgeDeletionRateSum() << std::endl;
    std::cout << "edge add.r.s.: " << contNetwork.getEdgeAdditionRateSum() << std::endl;
    std::cout << "ed. exist: " << nEdgesExist << std::endl;
    std::cout << "ed. to_add:" << nEdgesAdd << std::endl;
    std::cout << "rel1: " << nEdgesExist /contNetwork.getEdgeDeletionRateSum() << std::endl;
    std::cout << "rel2: " << nEdgesAdd / contNetwork.getEdgeAdditionRateSum() << std::endl;

    std::cout << "tau: " << tau << std::endl;
    return tau;
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

int NSA::sampleRandBinomial(int nTrials, double p)
{
    std::binomial_distribution<> randBinomial(nTrials, p);
    int result  = randBinomial(generator);
    return result;
}


void NSA::executeReaction(ContactNetwork & contNetwork, std::string reactId,
                          double rStart, double rBound, double time, uint32_t &nInf)
{
    //std::cout << reactId <<std::endl;
    if (reactId == "edge_del")
    {
        contNetwork.executeEdgeDelition(rStart, rBound);
    }

    else if (reactId == "edge_add")
    {
        contNetwork.executeEdgeAddition(rStart, rBound);
    }
    else if (reactId == "transmission")
    {
        contNetwork.executeTransmission(rStart, rBound, time);
        nInf++;
        //std::cout << "transmission " << time << " " << contNetwork.countByState(Specie::I) << " " << contNetwork.size() <<std::endl;
    }

    else if (reactId == "death")
    {
        contNetwork.executeDeath(rStart, rBound);
         //std::cout  << "death " << time << " " << contNetwork.countByState(Specie::I)  << " " << contNetwork.size() <<std::endl;
    }

    else if (reactId == "birth")
    {
        contNetwork.executeBirth(rStart, rBound);
        //std::cout << "birth " << time << " " << contNetwork.countByState(Specie::I)  << " " << contNetwork.size() <<std::endl;
    }
    //std::cout << "------------------" <<std::endl;
}


void NSA::BDtauleap(double tStart, double tEnd, ContactNetwork & contNetwork, double epsilon,
               std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr)
{
    double time = tStart;
    timeSteps.push_back(time);
    degreeDistr.push_back(contNetwork.getDegreeDistribution());
    std::unordered_map<std::string, double >propensities {
            {"edge_del", contNetwork.getEdgeDeletionRateSum()},
            {"edge_add", contNetwork.getEdgeAdditionRateSum()}
    };

    //std:: cout << "add: "<<propensities.at("edge_add") <<std::endl;

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
            timeSteps.push_back(time);
            degreeDistr.push_back(contNetwork.getDegreeDistribution());

            break;
        }

        double tau = proposeTimestep(epsilon, contNetwork);
        //std::cout << "tau=" <<tau << std::endl;

        if (time + tau > tEnd)
        {
            tau = tEnd - time;
        }

        time += tau;

        std::unordered_map<std::string, size_t >k {
                {"edge_del", 0},
                {"edge_add", 0}
        };

        std::unordered_map<std::string, size_t >kMax {
                {"edge_del", contNetwork.countEdges()},
                {"edge_add", contNetwork.getAmountOfEdgesPossibleToAdd()}
        };

        for (auto &it: propensities )
        {
            double p = 1;

            if (it.second * tau <= kMax.at(it.first))
            {
                p = it.second * tau / kMax.at(it.first);
            }

            k.at(it.first) = sampleRandBinomial(kMax.at(it.first), p);
        }

        std::vector<std::string> order;

        //size_t i = 0;
        for (auto &it: k )
        {
            // i++;
            order.insert(order.end(), it.second, it.first);
        }

        std::shuffle(order.begin(), order.end(), generator);
        std:: cout << "reactions fired:" << order.size() << std::endl;
        std:: cout << "deletion fired:" << k.at("edge_del") << std::endl;
        std:: cout << "addition fired:" << k.at("edge_add") << std::endl;
        for (size_t i = 0; i < order.size(); i++)
        {
            std::cout <<order.at(i) << " ";

            double r = sampleRandUni();

            uint32_t  ni = 0;
            executeReaction(contNetwork, order.at(i), 0, r * propensities.at(order.at(i)), time, ni);

            propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum();
            propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum();
        }
        std::cout <<std::endl;
        std:: cout << "------------" << std::endl;
        timeSteps.push_back(time);

        degreeDistr.push_back(contNetwork.getDegreeDistribution());
        //std:: cout <<"---------------"<<std::endl;

    }

}