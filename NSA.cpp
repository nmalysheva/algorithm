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

                //tUpds.push_back(networkLastUpdate);
                PoissonTauleap(networkLastUpdate, time, contNetwork, epsilon, tSteps, degreeDistr, false/*, benToFile*/);
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
        double r = sampleRandUni();

        if (k.at(0) == 0)
        {
            if (propensities.at("edge_del") >= r * aCrit)
            {
                kDel = 1;
            }
            else if (k.at(1) == 0)
            {
                if (propensities.at("edge_del") + propensities.at("edge_add") >= r * aCrit)
                {
                    kAdd = 1;
                }
            }
        }
        else if (k.at(1) == 0)
        {
            if (propensities.at("edge_del") + propensities.at("edge_add") >= r * aCrit)
            {
                kAdd = 1;
            }
        }
    }

    if (k.at(0) == -1)
    {
        std::poisson_distribution<size_t> poiss(propensities.at("edge_del") * tau);
        kDel = poiss(generator);
    }

    if (k.at(1) == -1)
    {
        std::poisson_distribution<size_t> poiss(propensities.at("edge_add") * tau);
        kAdd = poiss(generator);
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
                    std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                    bool updateDegreeDistr/*, std::vector<BenStructure> &benToFile*/)
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

    std::vector<std::pair<double, lemon::ListGraph::Edge>> propDel;
    propDel.reserve(1e6 + 1);
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propAdd;
    propAdd.reserve(1e6 + 1);

    while (t < tEnd)
    {
        size_t nAdd = 0;
        size_t nDel = 0;


        //propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum(nDel);
        //propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum(nAdd);

        propDel = contNetwork.getEdgeDeletionRateSum();
        propAdd = contNetwork.getEdgeAdditionRateSum();
        propensities.at("edge_del") = propDel.at(propDel.size() - 1).first;
        propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;

        nDel = propDel.size() - 1;
        nAdd = propAdd.size() - 1;

        //std::cout << "edge_del_prop " << propensities.at("edge_del") << std::endl;
        //std::cout << "edge_add_prop " << propensities.at("edge_add") << std::endl;

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

        std::vector<double> mu(N, 0);
        std::vector<double> sigmaSq(N, 0);
        std::vector<int> k(M, 0);


        size_t lDel = nDel;
        //std::cout << "lDel = " << lDel << std::endl;

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
            lAdd = contNetwork.getAmountOfEdgesToAddSafe();
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
                    propDel = contNetwork.getEdgeDeletionRateSum();
                    propAdd = contNetwork.getEdgeAdditionRateSum();
                    propensities.at("edge_del") = propDel.at(propDel.size() - 1).first;
                    propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;
                    //propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum(nDel);
                    //propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum(nAdd);

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
                    if (propensities.at("edge_del") >= r * propensitiesSum)
                    {
                        std::cout << i << "=del, ";
                        //BenStructure b(t, -1, -1, false);
                        //contNetwork.executeEdgeDeletion(0, r * propensitiesSum/*, b*/);
                        lemon::ListGraph::Edge e = binarySearch(propDel, 0, propDel.size() - 1, 0, r * propensitiesSum);
                        std::pair<int, int> b = contNetwork.removeEdge(e);
                        std::cout <<  b.first << ", " << b.second << std::endl;

                        //benToFile.push_back(b);
                    }
                    else
                    {
                        std::cout << i << "=add, ";
                        //BenStructure b(t, -1, -1, true);
                        //contNetwork.executeEdgeAddition(propensities.at("edge_del"), r * propensitiesSum/*, b*/);
                        lemon::ListGraph::Edge e = binarySearch(propAdd, 0, propAdd.size() - 1, propensities.at("edge_del"), r * propensitiesSum);
                        contNetwork.addEdge(e);
                        //benToFile.push_back(b);
                    }
                }
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
                        if (order.at(ind) == 1)
                        {
                            double r = randuni(generator);
                            while (r == 0)
                            {
                                r = randuni(generator);
                            }
                            //BenStructure b(t, -1, -1, true);
                            //contNetwork.executeEdgeAddition(0, r * propensities.at("edge_add")/*, b*/);
                            lemon::ListGraph::Edge e = binarySearch(propAdd, 0, propAdd.size() - 1, 0, r * propensities.at("edge_add"));
                            contNetwork.addEdge(e);
                            //benToFile.push_back(b);
                            maxEdgesDelete++;
                        }
                        else if (order.at(ind) == 0)
                        {
                            //std::uniform_int_distribution<size_t> dist(0, maxEdgesDelete - 1);
                            /*NOTE: we can delete uniformely when all edge del. rates are the same. When we have adaptivity case
                             * - we can not do this anymore!!!*/
                            //size_t ind = dist(generator);
                            //contNetwork.executeEdgeDeletion(ind);

                            double r = randuni(generator);
                            while (r == 0)
                            {
                                r = randuni(generator);
                            }
                            //BenStructure b(t, -1, -1, false);
                            //contNetwork.executeEdgeDeletion(0, r * propensities.at("edge_del")/*, b*/);
                            lemon::ListGraph::Edge e = binarySearch(propDel, 0, propDel.size() - 1, 0, r * propensities.at("edge_del"));
                            contNetwork.removeEdge(e);
                            //std::cout <<  b.first << ", " << b.second << std::endl;
                            //benToFile.push_back(b);
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


void NSA::RKF45Approximation(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork,
        double dtMax, double dtMin, double errorMax, double errorMin, std::vector<double> &timeSteps,
        std::vector<std::vector<size_t>> &degreeDistr, bool updateDegreeDistr/*, std::vector<BenStructure> &benToFile*/)
{
    std::vector<std::vector<double>> a({{1.0/4.0,        0.0,           0.0,           0.0,          0.0},
                                           {3.0/32.0,       9.0/32.0,      0.0,           0.0,          0.0},
                                           {1932.0/2197.0, 7200.0/2197.0, 7296.0/2197.0, 0.0,          0.0},
                                           {439.0/216.0,   8.0,           3680.0/513.0, 845.0/4104.0, 0.0},
                                           {8.0/27.0,      2.0,          3544.0/2565.0, 1859.0/4104.0, 11.0/40.0}});


    std::vector<double> b({25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0});
    std::vector<double> d({16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, 9.0/50.0, 2.0/55.0});


    //double theta = contNetwork.getExpectedEdgeDeletionRate();
    //std::cout << "theta: " << theta << std::endl;

    double t = tLastNetworkUpdate;

    if (updateDegreeDistr)
    {
        timeSteps.push_back(t);
        degreeDistr.push_back(contNetwork.getDegreeDistribution());
    }

    dtMax = (tEnd - tLastNetworkUpdate);
    double dt = dtMax;
    while (t < tEnd)
    {
        dt = std::min(std::max(dt, dtMin), dtMax);

        size_t nX = 0;
        double lam = contNetwork.getExpectedEdgeAdditionRate(nX);
        double theta = contNetwork.getExpectedEdgeDeletionRate();
        double X_ = static_cast<double> (nX);
        double X = static_cast<double> (contNetwork.countEdges());
        double Xmax = static_cast<double>(contNetwork.getAmountOfEdgesToAddSafe());

        std::cout <<"X = " << X << "X_ = " << X_<< std::endl;
        double k1 = dt * (X_ * lam - theta * X);

        double d_X_k2 = a.at(0).at(0) * k1;
        std::cout <<"d_X_k2  = " << d_X_k2 << std::endl;
        double X_k2 = X + d_X_k2;
        //double X_appr = (Xmax + X_ - 2 * d_X_k2) / 2;
        double X_appr = X_ - d_X_k2;
        std::cout <<"X_appr2  = " << X_appr << ", ";
        double k2 = dt * (X_appr * lam - X_k2 * theta);

        double d_X_k3 = a.at(1).at(0) * k1 + a.at(1).at(1) * k2;

        double X_k3 = X + d_X_k3;
        //X_appr = (Xmax + X_ - 2 * d_X_k3) / 2;
        X_appr = X_ - d_X_k3;
        std::cout <<"X_appr3  = " << X_appr << ", ";
        double k3 = dt * (X_appr * lam - X_k3 * theta);

        double d_X_k4 = a.at(2).at(0) * k1 - a.at(2).at(1) * k2 + a.at(2).at(2) * k3;

        double X_k4 = X + d_X_k4;
        //X_appr = (Xmax + X_ - 2 * d_X_k4) / 2;
        X_appr = X_ - d_X_k4;
        std::cout <<"X_appr4  = " << X_appr << ", ";

        double k4 = dt * (X_appr * lam - X_k4 * theta);

        double d_X_k5 = a.at(3).at(0) * k1 - a.at(3).at(1) * k2 + a.at(3).at(2) * k3 - a.at(3).at(3) * k4;

        double X_k5 = X + d_X_k5;
        //X_appr = (Xmax + X_ - 2 * d_X_k5) / 2;
        X_appr = X_ - d_X_k5;
        std::cout <<"X_appr5  = " << X_appr << ", ";
        double k5 = dt * (X_appr * lam - X_k5 * theta);

        double d_X_k6 = - a.at(4).at(0) * k1 + a.at(4).at(1) * k2 - a.at(4).at(2) * k3 + a.at(4).at(3) * k4 - a.at(4).at(4) * k5;

        double X_k6 = X + d_X_k6;
        X_appr = X_ - d_X_k6;
        std::cout <<"X_appr6  = " << X_appr << std::endl;
        //X_appr = (Xmax + X_ - 2 * d_X_k6) / 2;

        double k6 = dt * (X_appr * lam - X_k6 * theta);

        double n1 = b.at(0) * k1 + b.at(1) * k2 + b.at(2) * k3 + b.at(3) * k4 + b.at(4) * k5;
        double n2 = d.at(0) * k1 + d.at(1) * k2 + d.at(2) * k3 + d.at(3) * k4 - d.at(4) * k5 + d.at(5) * k6;

        double nDel = getNdel(X, X_, lam, theta, a, d, dt);
        double nAdd = getNadd(X, X_, lam, theta, a, d, dt);


        double err = abs(n1 - n2);

        if (err > errorMax && dt > dtMin)
        {
            dt = dt / 2;
        }
        else
        {
            if (t + dt > tEnd)
            {
                tLastNetworkUpdate = t;
                if (updateDegreeDistr)
                {
                    timeSteps.push_back(tEnd);
                    degreeDistr.push_back(contNetwork.getDegreeDistribution());
                }
                break;
            }
            t = t + dt;

            //std::poisson_distribution<int> poiss(abs(n2));
            std::poisson_distribution<int> poiss(nAdd);
            int mAdd  = poiss(generator);
            std::poisson_distribution<int> poiss2(nDel);
            int mDel  = poiss2(generator);
            int m  = mAdd - mDel;
            std::cout <<"n2 = " << n2 << "; nDel = " << nDel << "; nAdd = " << nAdd  << "; nDiff = " << nAdd - nDel << std::endl;
            int nEdges =  contNetwork.countEdges();

            if (X + m < 0 || m > Xmax)
            {
                dt = dt/2;
                throw std::domain_error("LOOOOL");
            }
            else
            {
                size_t nSurv = contNetwork.updateSurvivalProbability(mDel, mAdd/*, benToFile, t*/);

                std::cout << "nSurv = " << nSurv << std::endl;
                double nnn = contNetwork.getEdgeAdditionRateSum(nX);
                std::cout << "nEdges = " << contNetwork.countEdges() << std::endl;
                std::cout << "to add = " << abs(X + m - nSurv) << std::endl;
                int change = X + m - nSurv;

                if (change > 0)
                {
                    for (int i = 0; i < change; i ++)
                    {

                        double r = randuni(generator);
                        while (r == 0)
                        {
                            r = randuni(generator);
                        }
                        /*BenStructure b (t, -1, -1, true);
                        contNetwork.executeEdgeAddition(0, r * nnn, b);
                        benToFile.push_back(b);*/

                        contNetwork.executeEdgeAddition(0, r * nnn);
                    }

                }
                else
                {
                    int nDel =  contNetwork.countEdges();
                    for (int i = 0; i < abs(change); i ++)
                    {
                        std::uniform_int_distribution<size_t> dist(0, nDel - 1);

                        size_t ind = dist(generator);
                        contNetwork.executeEdgeDeletion(ind);
                        nDel--;
                    }
                }

                if (updateDegreeDistr)
                {
                    timeSteps.push_back(t);
                    degreeDistr.push_back(contNetwork.getDegreeDistribution());
                }

            }

            if (err < errorMin)
            {
                dt = dt * 2;
            }
        }
    }
}

/*void NSA::MidpointApproximation(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork,
                             double dtMax, double dtMin, double errorMax, double errorMin, std::vector<double> &timeSteps,
                             std::vector<std::vector<size_t>> &degreeDistr, bool updateDegreeDistr)
{
    double t = tLastNetworkUpdate;
    double dt = 1e-3;
   // double theta = contNetwork.getExpectedEdgeDeletionRate();

    timeSteps.push_back(t);
    degreeDistr.push_back(contNetwork.getDegreeDistribution());

    while (t < tEnd)
    {
        size_t nX = 0; //num of edges to add
        double theta = contNetwork.getExpectedEdgeDeletionRate();
        double lam = contNetwork.getExpectedEdgeAdditionRate(nX);
        std::cout << "lam = " << lam << std::endl;
        double X_ = static_cast<double> (nX);
        double X = static_cast<double> (contNetwork.countEdges());
        double Xmax = static_cast<double>(contNetwork.getAmountOfEdgesToAddSafe());
        double dX2 = dt/2 * (X_ * lam - theta * X);
        //double X_appr = (Xmax + X_ - 2 * dX2) / 2;
        double X_appr = X_ - dX2;
        std::cout << "X_appr = " << X_appr << std::endl;
        std::cout << "X_ = " << X_ << std::endl;
        std::cout << "X_max = " << Xmax << std::endl;
        std::cout << "X = " << X << std::endl;

        double addition = X_appr * lam;
        double deletion = (X + dX2) * theta;

        std::poisson_distribution<int> poiss1(dt * addition);
        int add  = poiss1(generator);
        std::poisson_distribution<int> poiss2(dt * deletion);
        int del  = poiss2(generator);

        int m = add - del;
        std::cout << "add = " << add << "; del =" << del << "; m=" << m  << std::endl;

        if (X + m < 0 || m > Xmax)
        //if (false)
        {

            std::cout << "add = " << add << "; del =" << del << "; m=" << m <<"  ERROR!!!" << std::endl;
            dt = dt/2;
            throw std::domain_error("LOOOOL");
        }
        else
        {
            t += dt;
            timeSteps.push_back(t);
            int nDel =  contNetwork.countEdges();
            size_t nSurv = contNetwork.updateSurvivalProbability(del, add);

            std::cout << "nSurv = " << nSurv << std::endl;
            double nnn = contNetwork.getEdgeAdditionRateSum(nX);
            std::cout << "nEdges = " << contNetwork.countEdges() << std::endl;
            std::cout << "to add = " << abs(X + m - nSurv) << std::endl;
            int change = X + m - nSurv;

            if (change > 0)
            {
                for (int i = 0; i < change; i ++)
                {

                    double r = randuni(generator);
                    while (r == 0)
                    {
                        r = randuni(generator);
                    }
                    contNetwork.executeEdgeAddition(0, r * nnn);
                }

            }
            else
            {
                int nDel =  contNetwork.countEdges();
                for (int i = 0; i < abs(change); i ++)
                {
                    std::uniform_int_distribution<size_t> dist(0, nDel - 1);

                    size_t ind = dist(generator);
                    contNetwork.executeEdgeDeletion(ind);
                    nDel--;
                }
            }

            degreeDistr.push_back(contNetwork.getDegreeDistribution());

        }

    }
}
*/
double NSA::getNadd(double X, double Y, double lam, double theta, std::vector<std::vector<double>> a, std::vector<double> d, double dt)
{
    double nAdd = X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 5)*std::pow(dt, 6)*theta
            + 4*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 6)*std::pow(theta, 2)
            + 6*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 6)*std::pow(theta, 3)
            + 4*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 6)*std::pow(theta, 4)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 6)*std::pow(theta, 5)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            - 3*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            - 3*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            - 3*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            - 3*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            + 3*X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            + 3*X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*X*a.at(0).at(0)*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(0).at(0)*a.at(1).at(1)*d.at(2)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - X*a.at(0).at(0)*a.at(1).at(1)*d.at(2)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            + 3*X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            + 3*X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            + X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            - X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(0).at(0)*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*X*a.at(0).at(0)*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - X*a.at(0).at(0)*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(0).at(0)*a.at(2).at(1)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + X*a.at(0).at(0)*a.at(2).at(1)*d.at(3)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(0).at(0)*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*X*a.at(0).at(0)*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + X*a.at(0).at(0)*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(0).at(0)*a.at(3).at(1)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - X*a.at(0).at(0)*a.at(3).at(1)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(0).at(0)*a.at(4).at(1)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - X*a.at(0).at(0)*a.at(4).at(1)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            - 3*X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            - 3*X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            + X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(1).at(0)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*X*a.at(1).at(0)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + X*a.at(1).at(0)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(1).at(0)*a.at(2).at(2)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - X*a.at(1).at(0)*a.at(2).at(2)*d.at(3)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(1).at(0)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*X*a.at(1).at(0)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - X*a.at(1).at(0)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(1).at(0)*a.at(3).at(2)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + X*a.at(1).at(0)*a.at(3).at(2)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(1).at(0)*a.at(4).at(2)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + X*a.at(1).at(0)*a.at(4).at(2)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(1).at(0)*d.at(2)*lam*std::pow(dt, 2)*theta
            - X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            - 3*X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            - 3*X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            + X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*X*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + X*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - X*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*X*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - X*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + X*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + X*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(1).at(1)*d.at(2)*lam*std::pow(dt, 2)*theta
            + X*a.at(2).at(0)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*X*a.at(2).at(0)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + X*a.at(2).at(0)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(2).at(0)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - X*a.at(2).at(0)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(2).at(0)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - X*a.at(2).at(0)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(2).at(0)*d.at(3)*lam*std::pow(dt, 2)*theta
            - X*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*X*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - X*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + X*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + X*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(2).at(1)*d.at(3)*lam*std::pow(dt, 2)*theta
            + X*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*X*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + X*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - X*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - X*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(2).at(2)*d.at(3)*lam*std::pow(dt, 2)*theta
            + X*a.at(3).at(0)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + X*a.at(3).at(0)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(3).at(0)*d.at(4)*lam*std::pow(dt, 2)*theta
            - X*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - X*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(3).at(1)*d.at(4)*lam*std::pow(dt, 2)*theta
            + X*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + X*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(3).at(2)*d.at(4)*lam*std::pow(dt, 2)*theta
            - X*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - X*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 2)*theta
            - X*a.at(4).at(0)*d.at(5)*lam*std::pow(dt, 2)*theta
            + X*a.at(4).at(1)*d.at(5)*lam*std::pow(dt, 2)*theta
            - X*a.at(4).at(2)*d.at(5)*lam*std::pow(dt, 2)*theta
            + X*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 2)*theta
            - X*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 2)*theta
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 6)*std::pow(dt, 6)
            - 4*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 5)*std::pow(dt, 6)*theta
            - 6*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 6)*std::pow(theta, 2)
            - 4*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 6)*std::pow(theta, 3)
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 6)*std::pow(theta, 4)
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 5)*std::pow(dt, 5)
            + 3*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            + 3*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 5)*std::pow(dt, 5)
            + 3*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            + 3*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(lam, 4)*std::pow(dt, 4)
            - 2*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 5)*std::pow(dt, 5)
            - 3*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            - 3*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(lam, 4)*std::pow(dt, 4)
            + 2*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 4)
            + 2*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(0).at(0)*a.at(1).at(1)*d.at(2)*std::pow(lam, 3)*std::pow(dt, 3)
            + Y*a.at(0).at(0)*a.at(1).at(1)*d.at(2)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 5)*std::pow(dt, 5)
            - 3*Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            - 3*Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            - Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            + Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(lam, 4)*std::pow(dt, 4)
            + 2*Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(0).at(0)*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 4)
            + 2*Y*a.at(0).at(0)*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + Y*a.at(0).at(0)*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(0).at(0)*a.at(2).at(1)*d.at(3)*std::pow(lam, 3)*std::pow(dt, 3)
            - Y*a.at(0).at(0)*a.at(2).at(1)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(0).at(0)*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 4)
            - 2*Y*a.at(0).at(0)*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - Y*a.at(0).at(0)*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(0).at(0)*a.at(3).at(1)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 3)
            + Y*a.at(0).at(0)*a.at(3).at(1)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(0).at(0)*a.at(4).at(1)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 3)
            + Y*a.at(0).at(0)*a.at(4).at(1)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 5)*std::pow(dt, 5)
            + 3*Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            + 3*Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            + Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 4)*std::pow(dt, 4)
            - 2*Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(1).at(0)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 4)
            - 2*Y*a.at(1).at(0)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - Y*a.at(1).at(0)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(1).at(0)*a.at(2).at(2)*d.at(3)*std::pow(lam, 3)*std::pow(dt, 3)
            + Y*a.at(1).at(0)*a.at(2).at(2)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(1).at(0)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 4)
            + 2*Y*a.at(1).at(0)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + Y*a.at(1).at(0)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(1).at(0)*a.at(3).at(2)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 3)
            - Y*a.at(1).at(0)*a.at(3).at(2)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(1).at(0)*a.at(4).at(2)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 3)
            - Y*a.at(1).at(0)*a.at(4).at(2)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(1).at(0)*d.at(2)*std::pow(lam, 2)*std::pow(dt, 2)
            + Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 5)*std::pow(dt, 5)
            + 3*Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            + 3*Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            + Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 4)*std::pow(dt, 4)
            - 2*Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 4)
            - 2*Y*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - Y*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(lam, 3)*std::pow(dt, 3)
            + Y*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 4)
            + 2*Y*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + Y*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 3)
            - Y*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 3)
            - Y*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(1).at(1)*d.at(2)*std::pow(lam, 2)*std::pow(dt, 2)
            - Y*a.at(2).at(0)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 4)
            - 2*Y*a.at(2).at(0)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - Y*a.at(2).at(0)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(2).at(0)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 3)
            + Y*a.at(2).at(0)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(2).at(0)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 3)
            + Y*a.at(2).at(0)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(2).at(0)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 2)
            + Y*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 4)
            + 2*Y*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + Y*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 3)
            - Y*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 3)
            - Y*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(2).at(1)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 2)
            - Y*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 4)
            - 2*Y*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - Y*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 3)
            + Y*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 3)
            + Y*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(2).at(2)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 2)
            - Y*a.at(3).at(0)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 3)
            - Y*a.at(3).at(0)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(3).at(0)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 2)
            + Y*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 3)
            + Y*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(3).at(1)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 2)
            - Y*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 3)
            - Y*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(3).at(2)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 2)
            + Y*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 3)
            + Y*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 2)
            + Y*a.at(4).at(0)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 2)
            - Y*a.at(4).at(1)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 2)
            + Y*a.at(4).at(2)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 2)
            - Y*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 2)
            + Y*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 2)
            + Y*d.at(0)*lam*dt
            + Y*d.at(2)*lam*dt
            + Y*d.at(3)*lam*dt
            - Y*d.at(4)*lam*dt
            + Y*d.at(5)*lam*dt;
    return nAdd;

}

double NSA::getNdel(double X, double Y, double lam, double theta, std::vector<std::vector<double>> a, std::vector<double> d, double dt)
{
    double ndel =  -X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 6)*std::pow(theta, 2)
            - 4*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 6)*std::pow(theta, 3)
            - 6*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 6)*std::pow(theta, 4)
            - 4*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 6)*std::pow(theta, 5)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(dt, 6)*std::pow(theta, 6)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            + 3*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            + 3*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(dt, 5)*std::pow(theta, 5)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            + 3*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            + 3*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(dt, 5)*std::pow(theta, 5)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - 2*X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(dt, 4)*std::pow(theta, 4)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            - 3*X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - 3*X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            - X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(dt, 5)*std::pow(theta, 5)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + 2*X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(dt, 4)*std::pow(theta, 4)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + 2*X*a.at(0).at(0)*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(0).at(0)*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(dt, 4)*std::pow(theta, 4)
            + X*a.at(0).at(0)*a.at(1).at(1)*d.at(2)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(0).at(0)*a.at(1).at(1)*d.at(2)*std::pow(dt, 3)*std::pow(theta, 3)
            - X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            - 3*X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - 3*X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            - X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(dt, 5)*std::pow(theta, 5)
            + X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + 2*X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(dt, 4)*std::pow(theta, 4)
            + X*a.at(0).at(0)*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + 2*X*a.at(0).at(0)*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(0).at(0)*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(dt, 4)*std::pow(theta, 4)
            - X*a.at(0).at(0)*a.at(2).at(1)*d.at(3)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(0).at(0)*a.at(2).at(1)*d.at(3)*std::pow(dt, 3)*std::pow(theta, 3)
            - X*a.at(0).at(0)*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - 2*X*a.at(0).at(0)*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(0).at(0)*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(dt, 4)*std::pow(theta, 4)
            + X*a.at(0).at(0)*a.at(3).at(1)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(0).at(0)*a.at(3).at(1)*d.at(4)*std::pow(dt, 3)*std::pow(theta, 3)
            + X*a.at(0).at(0)*a.at(4).at(1)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(0).at(0)*a.at(4).at(1)*d.at(5)*std::pow(dt, 3)*std::pow(theta, 3)
            + X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            + 3*X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            + 3*X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            + X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(dt, 5)*std::pow(theta, 5)
            - X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - 2*X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(dt, 4)*std::pow(theta, 4)
            - X*a.at(1).at(0)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - 2*X*a.at(1).at(0)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(1).at(0)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(dt, 4)*std::pow(theta, 4)
            + X*a.at(1).at(0)*a.at(2).at(2)*d.at(3)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(1).at(0)*a.at(2).at(2)*d.at(3)*std::pow(dt, 3)*std::pow(theta, 3)
            + X*a.at(1).at(0)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + 2*X*a.at(1).at(0)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(1).at(0)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(dt, 4)*std::pow(theta, 4)
            - X*a.at(1).at(0)*a.at(3).at(2)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(1).at(0)*a.at(3).at(2)*d.at(4)*std::pow(dt, 3)*std::pow(theta, 3)
            - X*a.at(1).at(0)*a.at(4).at(2)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(1).at(0)*a.at(4).at(2)*d.at(5)*std::pow(dt, 3)*std::pow(theta, 3)
            - X*a.at(1).at(0)*d.at(2)*std::pow(dt, 2)*std::pow(theta, 2)
            + X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            + 3*X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            + 3*X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            + X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(dt, 5)*std::pow(theta, 5)
            - X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - 2*X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(dt, 4)*std::pow(theta, 4)
            - X*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - 2*X*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(dt, 4)*std::pow(theta, 4)
            + X*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(dt, 3)*std::pow(theta, 3)
            + X*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + 2*X*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(dt, 4)*std::pow(theta, 4)
            - X*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(dt, 3)*std::pow(theta, 3)
            - X*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(dt, 3)*std::pow(theta, 3)
            - X*a.at(1).at(1)*d.at(2)*std::pow(dt, 2)*std::pow(theta, 2)
            - X*a.at(2).at(0)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - 2*X*a.at(2).at(0)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(2).at(0)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(dt, 4)*std::pow(theta, 4)
            + X*a.at(2).at(0)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(2).at(0)*a.at(3).at(3)*d.at(4)*std::pow(dt, 3)*std::pow(theta, 3)
            + X*a.at(2).at(0)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(2).at(0)*a.at(4).at(3)*d.at(5)*std::pow(dt, 3)*std::pow(theta, 3)
            - X*a.at(2).at(0)*d.at(3)*std::pow(dt, 2)*std::pow(theta, 2)
            + X*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + 2*X*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + X*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(dt, 4)*std::pow(theta, 4)
            - X*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(dt, 3)*std::pow(theta, 3)
            - X*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(dt, 3)*std::pow(theta, 3)
            + X*a.at(2).at(1)*d.at(3)*std::pow(dt, 2)*std::pow(theta, 2)
            - X*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - 2*X*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - X*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(dt, 4)*std::pow(theta, 4)
            + X*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(dt, 3)*std::pow(theta, 3)
            + X*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(dt, 3)*std::pow(theta, 3)
            - X*a.at(2).at(2)*d.at(3)*std::pow(dt, 2)*std::pow(theta, 2)
            - X*a.at(3).at(0)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(3).at(0)*a.at(4).at(4)*d.at(5)*std::pow(dt, 3)*std::pow(theta, 3)
            + X*a.at(3).at(0)*d.at(4)*std::pow(dt, 2)*std::pow(theta, 2)
            + X*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(dt, 3)*std::pow(theta, 3)
            - X*a.at(3).at(1)*d.at(4)*std::pow(dt, 2)*std::pow(theta, 2)
            - X*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - X*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(dt, 3)*std::pow(theta, 3)
            + X*a.at(3).at(2)*d.at(4)*std::pow(dt, 2)*std::pow(theta, 2)
            + X*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + X*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(dt, 3)*std::pow(theta, 3)
            - X*a.at(3).at(3)*d.at(4)*std::pow(dt, 2)*std::pow(theta, 2)
            + X*a.at(4).at(0)*d.at(5)*std::pow(dt, 2)*std::pow(theta, 2)
            - X*a.at(4).at(1)*d.at(5)*std::pow(dt, 2)*std::pow(theta, 2)
            + X*a.at(4).at(2)*d.at(5)*std::pow(dt, 2)*std::pow(theta, 2)
            - X*a.at(4).at(3)*d.at(5)*std::pow(dt, 2)*std::pow(theta, 2)
            + X*a.at(4).at(4)*d.at(5)*std::pow(dt, 2)*std::pow(theta, 2)
            + X*d.at(0)*dt*theta
            + X*d.at(2)*dt*theta
            + X*d.at(3)*dt*theta
            - X*d.at(4)*dt*theta
            + X*d.at(5)*dt*theta
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 5)*std::pow(dt, 6)*theta
            + 4*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 6)*std::pow(theta, 2)
            + 6*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 6)*std::pow(theta, 3)
            + 4*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 6)*std::pow(theta, 4)
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 6)*std::pow(theta, 5)
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            - 3*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            - 3*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            - 3*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            - 3*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            + 3*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            + 3*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            + Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*Y*a.at(0).at(0)*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(0).at(0)*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - Y*a.at(0).at(0)*a.at(1).at(1)*d.at(2)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(0).at(0)*a.at(1).at(1)*d.at(2)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            + 3*Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            + 3*Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            + Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            - Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(0).at(0)*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - Y*a.at(0).at(0)*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*Y*a.at(0).at(0)*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(0).at(0)*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + Y*a.at(0).at(0)*a.at(2).at(1)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(0).at(0)*a.at(2).at(1)*d.at(3)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + Y*a.at(0).at(0)*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*Y*a.at(0).at(0)*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(0).at(0)*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - Y*a.at(0).at(0)*a.at(3).at(1)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(0).at(0)*a.at(3).at(1)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - Y*a.at(0).at(0)*a.at(4).at(1)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(0).at(0)*a.at(4).at(1)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            - 3*Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            - 3*Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            + Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(1).at(0)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + Y*a.at(1).at(0)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*Y*a.at(1).at(0)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(1).at(0)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - Y*a.at(1).at(0)*a.at(2).at(2)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(1).at(0)*a.at(2).at(2)*d.at(3)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - Y*a.at(1).at(0)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*Y*a.at(1).at(0)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(1).at(0)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + Y*a.at(1).at(0)*a.at(3).at(2)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(1).at(0)*a.at(3).at(2)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + Y*a.at(1).at(0)*a.at(4).at(2)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(1).at(0)*a.at(4).at(2)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + Y*a.at(1).at(0)*d.at(2)*lam*std::pow(dt, 2)*theta
            - Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 4)*std::pow(dt, 5)*theta
            - 3*Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 5)*std::pow(theta, 2)
            - 3*Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 5)*std::pow(theta, 3)
            - Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 5)*std::pow(theta, 4)
            + Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(1).at(1)*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + Y*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*Y*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(1).at(1)*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - Y*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(1).at(1)*a.at(2).at(2)*d.at(3)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - Y*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*Y*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(1).at(1)*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + Y*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(1).at(1)*a.at(3).at(2)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + Y*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(1).at(1)*a.at(4).at(2)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + Y*a.at(1).at(1)*d.at(2)*lam*std::pow(dt, 2)*theta
            + Y*a.at(2).at(0)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*Y*a.at(2).at(0)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(2).at(0)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - Y*a.at(2).at(0)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(2).at(0)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - Y*a.at(2).at(0)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(2).at(0)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + Y*a.at(2).at(0)*d.at(3)*lam*std::pow(dt, 2)*theta
            - Y*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            - 2*Y*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            - Y*a.at(2).at(1)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            + Y*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(2).at(1)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + Y*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(2).at(1)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - Y*a.at(2).at(1)*d.at(3)*lam*std::pow(dt, 2)*theta
            + Y*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 3)*std::pow(dt, 4)*theta
            + 2*Y*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 4)*std::pow(theta, 2)
            + Y*a.at(2).at(2)*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 4)*std::pow(theta, 3)
            - Y*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(2).at(2)*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - Y*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(2).at(2)*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + Y*a.at(2).at(2)*d.at(3)*lam*std::pow(dt, 2)*theta
            + Y*a.at(3).at(0)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(3).at(0)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - Y*a.at(3).at(0)*d.at(4)*lam*std::pow(dt, 2)*theta
            - Y*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(3).at(1)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + Y*a.at(3).at(1)*d.at(4)*lam*std::pow(dt, 2)*theta
            + Y*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            + Y*a.at(3).at(2)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            - Y*a.at(3).at(2)*d.at(4)*lam*std::pow(dt, 2)*theta
            - Y*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*std::pow(lam, 2)*std::pow(dt, 3)*theta
            - Y*a.at(3).at(3)*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 3)*std::pow(theta, 2)
            + Y*a.at(3).at(3)*d.at(4)*lam*std::pow(dt, 2)*theta
            - Y*a.at(4).at(0)*d.at(5)*lam*std::pow(dt, 2)*theta
            + Y*a.at(4).at(1)*d.at(5)*lam*std::pow(dt, 2)*theta
            - Y*a.at(4).at(2)*d.at(5)*lam*std::pow(dt, 2)*theta
            + Y*a.at(4).at(3)*d.at(5)*lam*std::pow(dt, 2)*theta
            - Y*a.at(4).at(4)*d.at(5)*lam*std::pow(dt, 2)*theta ;
    return ndel;

}

lemon::ListGraph::Edge NSA::binarySearch(std::vector<std::pair<double, lemon::ListGraph::Edge>> &propCumSum,
                                         size_t indL, size_t indR, double rStart, double rBound)
{
    lemon::ListGraph::Edge result(lemon::INVALID);

    if (indR == indL)
    {
        lemon::ListGraph::Edge result = propCumSum.at(indR).second;
        propCumSum.erase(propCumSum.begin() + indR);
        return result;

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