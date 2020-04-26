//
// Created by Malysheva, Nadezhda on 2019-07-30.
//

#include "NSA.h"
#include "math.h"
#include <algorithm>
#include <vector>
#include <unistd.h>

NSA::NSA()
{
    generator = std::mt19937_64(rDev());
    generator.seed(::time(NULL) * getpid()); //to change the seed for every run
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

                PoissonTauleap(networkLastUpdate, time, contNetwork, epsilon, tSteps, degreeDistr, false);
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
                            executeReaction(contNetwork, it.first, pSum, propUpperLimit * r, time, nInf);

                            if (it.first == "transmission")
                            {
                                tSteps.push_back(time);
                                nInfected.push_back(nInf);
                                degreeDistr.push_back(contNetwork.getDegreeDistribution());
                            }

                            if (it.first  == "death" )
                            {
                                nInf = contNetwork.countByState(Specie::State::I);
                                tSteps.push_back(time);
                                nInfected.push_back(nInf);
                                degreeDistr.push_back(contNetwork.getDegreeDistribution());
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
                   //time += proposedTime;//!!!!!!!
                   std::cout << "thin" << std::endl;
                }

            }
        }
    }
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

double NSA::proposeTau1(size_t lDel, size_t lAdd, size_t nAdd, double epsilon, std::vector<double> mu, std::vector<double> sigmaSq)
{
    double timeStep = std::numeric_limits<double>::infinity();


    if (lDel >= 10)
    {
        timeStep = std::min(timeStep, std::max(epsilon * lDel, 1.0) / std::abs(mu.at(0)));
        timeStep = std::min(timeStep, std::max(epsilon * lDel, 1.0) * std::max(epsilon * lDel, 1.0) / sigmaSq.at(0));
    }

    if (lAdd >= 10)
    {
        timeStep = std::min(timeStep, std::max(epsilon * /*lAdd*/ nAdd, 1.0) / std::abs(mu.at(1)));
        timeStep = std::min(timeStep, std::max(epsilon * /*lAdd*/ nAdd, 1.0) * std::max(epsilon * /*lAdd*/ nAdd, 1.0) / sigmaSq.at(1));
    }

    return timeStep;

}

double NSA::proposeTau2(double aCrit)
{
    double tau2 = std::numeric_limits<double>::infinity();
    if (aCrit > 0)
    {
        double r = sampleRandUni();
        tau2 = 1 / aCrit * std::log(1 / r);
    }
    return tau2;
}

void NSA::selectTimeStepAndK(double &tau, const std::unordered_map<std::string, double> &propensities, const std::vector<int> &k,
        size_t &kDel, size_t &kAdd,double tau1, double tau2, double aCrit)
{
    tau = 0;

    if (tau1 < tau2)
    {
        std::cout << "tau1: " << tau1 << std::endl;
        tau = tau1;

    }
    else
    {
        std::cout << "tau2: " << tau2 << std::endl;
        tau = tau2;
       // size_t randK1, randK2;
        double r = sampleRandUni();

        if (k.at(0) == 0)
        {
            if (propensities.at("edge_del") > r * aCrit)
            {
                kDel = 1;
                //std::cout << "kDel: " << kDel << std::endl;
            }
            else if (k.at(1) == 0)
            {
                if (propensities.at("edge_del") + propensities.at("edge_add") > r * aCrit)
                {
                    //k.at(1) = 1;
                    kAdd = 1;
                    //std::cout << "kAdd: " << kAdd << std::endl;
                }

            }
        }
        else if (k.at(1) == 0)
        {
            if (propensities.at("edge_del") + propensities.at("edge_add") > r * aCrit)
            {
                kAdd = 1;
                //k.at(1) = 1;
                //std::cout << "kAdd: " << kAdd << std::endl;
            }
        }
    }

    if (k.at(0) == -1)
    {
        std::poisson_distribution<size_t> poiss(propensities.at("edge_del") * tau);
        kDel = poiss(generator);
        //std::cout<<"kDel_poiss: " <<kDel <<std::endl;
    }

    if (k.at(1) == -1)
    {
        std::poisson_distribution<size_t> poiss(propensities.at("edge_add") * tau);
        kAdd = poiss(generator);
        //std::cout<<"kAdd_poiss: " <<kAdd <<std::endl;
    }

    if (propensities.at("edge_del") == 0)
    {
        kDel = 0;
    }

    if (propensities.at("edge_add") == 0)
    {
        kAdd = 0;
    }
}

void NSA::PoissonTauleap(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork, double epsilon,
                    std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr, bool updateDegreeDistr)
{

    int N = 2;
    int M = 2;

    std::unordered_map<std::string, double> propensities {
            {"edge_del", 0},
            {"edge_add", 0} };

    double t = tLastNetworkUpdate;

    if (updateDegreeDistr)
    {
        timeSteps.push_back(t);
        degreeDistr.push_back(contNetwork.getDegreeDistribution());
    }

    while (t < tEnd)
    {

        //std::cout << "new iteration; t = " << t << std::endl;

        size_t nAdd = 0;
        size_t nDel = 0;

        propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum(nDel);
        propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum(nAdd);
        std::cout << "edge_del_prop " << propensities.at("edge_del") << std::endl;
        std::cout << "edge_add_prop " << propensities.at("edge_add") << std::endl;

        //bool equil = abs(propensities.at("edge_del") - propensities.at("edge_add")) < 0.05 * std::min(propensities.at("edge_del"), propensities.at("edge_add"));
        //std::cout << "equilibrium: " << equil << std::endl;

        if (propensities.at("edge_del") + propensities.at("edge_add") == 0)
        {
            t = tEnd;
            if (updateDegreeDistr)
            {
                timeSteps.push_back(t);
                degreeDistr.push_back(contNetwork.getDegreeDistribution());
            }
            break;
        }

        //std::cout<<"edge del prop: " << propensities.at("edge_del") <<std::endl;
        //std::cout<<"edge add prop: " << propensities.at("edge_add") <<std::endl;
        //std::cout<<"nAdd: " << nAdd <<std::endl;


        std::vector<double> mu(N, 0);
        std::vector<double> sigmaSq(N, 0);
        std::vector<int> k(M, 0);


        size_t lDel = nDel;
        std::cout << "lDel = " << lDel << std::endl;

        double aCrit = 0;

        if (propensities.at("edge_del") > 0)
        {
            if (lDel < 10)
            {
                aCrit += propensities.at("edge_del");
            }
            else
            {
                mu.at(0) -=  propensities.at("edge_del");
                sigmaSq.at(0) +=  propensities.at("edge_del");

                mu.at(1) += propensities.at("edge_del");
                sigmaSq.at(1) +=  propensities.at("edge_del");

                k.at(0) = -1;//static_cast<size_t> (-1);

            }
        }

        size_t lAdd = 0;

        if (propensities.at("edge_add") > 0)
        {
            lAdd = contNetwork.subgraph();
            std::cout << "lAdd = " << lAdd << std::endl;

            if (lAdd < 10)
            {
                aCrit += propensities.at("edge_add");
            }
            else
            {
                mu.at(0) +=  propensities.at("edge_add");
                sigmaSq.at(0) +=   propensities.at("edge_add");

                mu.at(1) -=  propensities.at("edge_add");
                sigmaSq.at(1) +=  propensities.at("edge_add");
                k.at(1) = -1;//static_cast<size_t> (-1);
            }
        }

        double tau1 = proposeTau1(lDel, lAdd, nAdd, epsilon, mu, sigmaSq);

        bool flag = true;
        while (flag)
        {
            //std::cout <<"tau1: " << tau1 <<std::endl;
            if (tau1 < 10.0 / (propensities.at("edge_del") + propensities.at("edge_add") ) )
            {
                //SSA
                std::cout << "prop.sum: " << (propensities.at("edge_del") + propensities.at("edge_add") ) << std::endl;
                std::cout << "lAdd: " << lAdd << "; lDel: " << lDel << std::endl;
                std::cout << "ssa. start time =  " << t <<"; ";
                for (size_t i = 0; i < 100; i++)
                {
                    propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum(nDel);
                    propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum(nAdd);

                    double propensitiesSum = propensities.at("edge_del") + propensities.at("edge_add");

                    if (propensitiesSum == 0)
                    {
                        tLastNetworkUpdate = tEnd; //used to update netw.Upd.Time
                        t = tEnd;
                        if (updateDegreeDistr)
                        {
                            timeSteps.push_back(t);
                            degreeDistr.push_back(contNetwork.getDegreeDistribution());
                        }
                        break;
                    }

                    double r = sampleRandUni();

                    double proposedTime = 1 / propensitiesSum * std::log(1 / r);

                    if (t + proposedTime > tEnd)
                    {
                        tLastNetworkUpdate = t;
                        t = tEnd;
                        if (updateDegreeDistr)
                        {
                            timeSteps.push_back(t);
                            degreeDistr.push_back(contNetwork.getDegreeDistribution());
                        }
                        break;
                    }

                    t += proposedTime;
                    tLastNetworkUpdate = t; //used to update netw.Upd.Time

                    r = sampleRandUni();
                    //deletion
                    if (propensities.at("edge_del") > r * propensitiesSum)
                    {
                        std::cout << i << "=del, ";
                        contNetwork.executeEdgeDeletion(0, r * propensitiesSum);
                    }
                    else
                    {
                        std::cout << i << "=add, ";
                        contNetwork.executeEdgeAddition(propensities.at("edge_del"), r * propensitiesSum);
                    }
                }
                std::cout << std::endl;
                flag = false;
            }
            else
            {

                double tau2 = proposeTau2(aCrit);

                double tau = 0;

                size_t kAdd = 0;
                size_t kDel = 0;

                selectTimeStepAndK(tau, propensities, k, kDel, kAdd, tau1, tau2, aCrit);

                //std::cout << "tau = "
                if (kDel > lDel || kAdd > lAdd)
                {
                    tau1 = tau1 / 2;
                }

                else if (t + tau > tEnd)
                {
                    tLastNetworkUpdate = t;
                    t = tEnd;
                    flag = false;
                }

                else
                {
                    t = t + tau;
                    std::cout << "t = " << t << std::endl;
                    tLastNetworkUpdate = t;
                    k.at(0) = kDel;
                    k.at(1) = kAdd;

                    //TODO get EdgeidList of network and Complement. May be EdgeUIDs?


                    std::vector<int> order;

                    std::cout << "k_Del = " << k.at(0) << "; k_Add = " << k.at(1) << std::endl;


                    for (size_t ind = 0; ind < k.size(); ind++)
                    {
                        order.insert(order.end(), k.at(ind), ind);
                    }
                    std::shuffle(order.begin(), order.end(), generator);

                    int maxEdgesDelete = nDel;
                    for (size_t ind = 0; ind < order.size(); ind++)
                    {
                        //std::cout << order.at(ind) << " ";
                        if (order.at(ind) == 1)
                        {
                            int maxEdgesAdd    = contNetwork.getAmountOfEdgesToAdd();
                            std::uniform_int_distribution<size_t> dist(0, maxEdgesAdd - 1);
                            size_t edgeNum = dist(generator);
                            contNetwork.executeEdgeAddition(edgeNum, maxEdgesAdd);
                            maxEdgesDelete++;
                        }
                        else if (order.at(ind) == 0)

                        {
                            std::uniform_int_distribution<size_t> dist(0, maxEdgesDelete - 1);

                            size_t ind = dist(generator);
                            contNetwork.executeEdgeDeletion(ind, maxEdgesDelete);
                            maxEdgesDelete--;
                        }
                    }
                    flag = false;
                }
            }
        }

        if (updateDegreeDistr)
        {
            timeSteps.push_back(t);
            degreeDistr.push_back(contNetwork.getDegreeDistribution());
        }
    }
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


void NSA::executeReaction(ContactNetwork & contNetwork, std::string reactId,
                          double rStart, double rBound, double time, uint32_t &nInf)
{
    //std::cout << reactId <<std::endl;
    if (reactId == "edge_del")
    {
        contNetwork.executeEdgeDeletion(rStart, rBound);
    }

    else if (reactId == "edge_add")
    {
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


void NSA::RKF45Approximation(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork,
        double dtMax, double dtMin, double errorMax, double errorMin, std::vector<double> &timeSteps,
        std::vector<std::vector<size_t>> &degreeDistr, bool updateDegreeDistr)
{
    double t = tLastNetworkUpdate;

    if (updateDegreeDistr)
    {
        timeSteps.push_back(t);
        degreeDistr.push_back(contNetwork.getDegreeDistribution());
    }

    double dt = dtMax;
    while (t < tEnd)
    {
        dt = std::min(std::max(dt, dtMin), dtMax);

        double nEdgesToAddTotal = static_cast<double>(contNetwork.getAmountOfEdgesToAdd());

        /*
         *
         */
    }
}