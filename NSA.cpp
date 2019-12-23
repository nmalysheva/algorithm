//
// Created by Malysheva, Nadezhda on 2019-07-30.
//

#include "NSA.h"
#include "math.h"
#include <algorithm>
#include <vector>

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

                //updateContacts(networkLastUpdate, time, contNetwork, epsilon);
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

double NSA::proposeTau1(size_t lDel, size_t lAdd, size_t nAdd, double epsilon, std::vector<double> mu, std::vector<double> sigmaSq)
{
    double timeStep = std::numeric_limits<double>::infinity();

    if (lDel >= 10)
    {
        timeStep = std::min(timeStep, std::max(epsilon * lDel, 1.0) / std::abs(mu.at(0)));
        timeStep = std::min(timeStep, std::max(epsilon * lDel, 1.0) * std::max(epsilon * lDel, 1.0) / sigmaSq.at(0));
        //std::cout <<" del tau: " << tau1 <<std::endl;
    }

    if (lAdd >= 10)
    {
        std::cout <<"nAdd: " << nAdd <<std::endl;
        std::cout <<"lAdd: " << lAdd <<std::endl;
        std::cout <<"lDel: " << lDel <<std::endl;

        timeStep = std::min(timeStep, std::max(epsilon * nAdd, 1.0) / std::abs(mu.at(1)));
        timeStep = std::min(timeStep, std::max(epsilon * nAdd, 1.0) * std::max(epsilon * nAdd, 1.0) / sigmaSq.at(1));
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

    //size_t randK1, randK2;

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
            if (propensities.at("edge_del") > r * aCrit) {
                //k.at(0) = 1;
                kDel = 1;
                std::cout << "kDel: " << kDel << std::endl;
            }
            else if (k.at(1) == 0)
            {
                if (propensities.at("edge_del") + propensities.at("edge_add") > r * aCrit)
                {
                    //k.at(1) = 1;
                    kAdd = 1;
                    std::cout << "kAdd: " << kAdd << std::endl;
                }

            }
        }
        else if (k.at(1) == 0)
        {
            if (propensities.at("edge_del") + propensities.at("edge_add") > r * aCrit)
            {
                kAdd = 1;
                //k.at(1) = 1;
                std::cout << "kAdd: " << kAdd << std::endl;
            }
        }
    }

    if (k.at(0) == -1)
    {
        std::poisson_distribution<size_t> poiss(propensities.at("edge_del") * tau);
        kDel = poiss(generator);
        std::cout<<"kDel: " <<kDel <<std::endl;
    }

    if (k.at(1) == -1)
    {
        std::poisson_distribution<size_t> poiss(propensities.at("edge_add") * tau);
        kAdd = poiss(generator);
        std::cout<<"kAdd: " <<kAdd <<std::endl;
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

void NSA::PoissonTauleap(double tStart, double tEnd, ContactNetwork & contNetwork, double epsilon,
                    std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr)
{
    //double eps = 0.03;

    int N = 2;
    int M = 2;

    std::unordered_map<std::string, double> propensities {
            {"edge_del", 0},
            {"edge_add", 0} };

    double t = tStart;

    timeSteps.push_back(t);
    degreeDistr.push_back(contNetwork.getDegreeDistribution());

    while (t < tEnd)
    {

        std::cout << "new iteration; t = " << t << std::endl;

        size_t nAdd = 0;
        size_t nDel = 0;

        propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum(nDel);
        propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum(nAdd);

        if (propensities.at("edge_del") + propensities.at("edge_add") == 0)
        {
            t = tEnd;
            timeSteps.push_back(t);
            degreeDistr.push_back(contNetwork.getDegreeDistribution());
            break;
        }

        //std::cout<<"edge del prop: " << propensities.at("edge_del") <<std::endl;
        //std::cout<<"edge add prop: " << propensities.at("edge_add") <<std::endl;
        //std::cout<<"nAdd: " << nAdd <<std::endl;


        std::vector<double> mu(N, 0);
        std::vector<double> sigmaSq(N, 0);
        std::vector<int> k(M, 0);


        size_t lDel = nDel;

        //std::cout<<"lDel: " << lDel <<std::endl;

        double aCrit = 0;

        if (propensities.at("edge_del") > 0)
        {
            std::cout<<"pdel>0 " <<std::endl;
            if (lDel < 10)
            {
                aCrit += propensities.at("edge_del");
            }
            else
            {
                mu.at(0) += - propensities.at("edge_del");
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
            if (lAdd < 10)
            {
                aCrit += propensities.at("edge_add");
            }
            else
            {
                mu.at(0) += + propensities.at("edge_add");
                sigmaSq.at(0) +=  + propensities.at("edge_add");

                mu.at(1) += - propensities.at("edge_add");
                sigmaSq.at(1) +=  propensities.at("edge_add");
                k.at(1) = -1;//static_cast<size_t> (-1);
            }
        }
        std::cout<<"lAdd: " << lAdd <<std::endl;

        double tau1 = proposeTau1(lDel, lAdd, nAdd, epsilon, mu, sigmaSq);

        //std::cout <<"tau: " << tau1 <<std::endl;


        bool flag = true;

        while (flag)
        {
            std::cout <<"tau1: " << tau1 <<std::endl;
            if (tau1 < 10.0 / (propensities.at("edge_del") + propensities.at("edge_add") ) )
            {
                //SSA
                std::cout << "ssa: " << std::endl;
                for (size_t i = 0; i < 100; i++)
                {
                    propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum(nDel);
                    propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum(nAdd);

                    double propensitiesSum = propensities.at("edge_del") + propensities.at("edge_add");

                    if (propensitiesSum == 0)
                    {
                        t = tEnd;
                        timeSteps.push_back(t);
                        degreeDistr.push_back(contNetwork.getDegreeDistribution());
                        break;
                    }

                    double r = sampleRandUni();

                    double proposedTime = 1 / propensitiesSum * std::log(1 / r);

                    if (t + proposedTime > tEnd)
                    {
                        t = tEnd;
                        timeSteps.push_back(t);
                        degreeDistr.push_back(contNetwork.getDegreeDistribution());
                        break;
                    }

                    t += proposedTime;


                    r = sampleRandUni();
                    //deletion
                    if (propensities.at("edge_del") > r * propensitiesSum)
                    {
                        contNetwork.executeEdgeDeletion(0, r * propensitiesSum);
                    }
                    else
                    {
                        contNetwork.executeEdgeAddition(propensities.at("edge_del"), r * propensitiesSum);
                    }

                    //std::cout << "t=" << t << std::endl;

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

                if (kDel > lDel || kAdd > lAdd)
                {
                    tau1 = tau1 / 2;
                }

                else if (t + tau > tEnd)
                {
                   t = tEnd;
                   flag = false;
                }

                else
                {
                    t = t + tau;
                    k.at(0) = kDel;
                    k.at(1) = kAdd;

                    //TODO get EdgeidList of network and Complement. May be EdgeUIDs?


                    std::vector<int> order;

                    std::cout << "k1 = " << k.at(0) << "; k2 = " << k.at(1) << std::endl;

                    //size_t i = 0;
                    for (size_t ind = 0; ind < k.size(); ind++)
                    {
                        // i++;
                        order.insert(order.end(), k.at(ind), ind);

                    }
                    std::shuffle(order.begin(), order.end(), generator);

                    //std::cout <<"check1" <<std::endl;
                    int maxEdgesDelete = nDel;
                    for (size_t ind = 0; ind < order.size(); ind++)
                    {
                        //std::cout << order.at(ind) << " ";
                        if (order.at(ind) == 1)
                        {
                            int maxEdgesAdd    = contNetwork.getAmountOfEdgesToAdd();
                            //std::cout <<"check2" <<std::endl;
                            //  std::cout <<"Unif. add " <<std::endl;
                            //contNetwork.executeEdgeAdditionUniform();
                            //int maxEdges = lemon::countEdges(complement);
                            std::uniform_int_distribution<size_t> dist(0, maxEdgesAdd - 1);
                            size_t edgeNum = dist(generator);
                            contNetwork.executeEdgeAddition(edgeNum, maxEdgesAdd);
                            maxEdgesDelete++;
                            //std::cout <<"check3" <<std::endl;

                        }
                        else if (order.at(ind) == 0)

                        {

                            //  std::cout <<"Unif. del " <<std::endl;
                            //int maxEdges = lemon::countEdges(network);

                            //std::cout <<"maxEdgeId " << maxEdgeId << " " << lemon::countEdges(network) <<std::endl;
                            std::uniform_int_distribution<size_t> dist(0, maxEdgesDelete - 1);


                            size_t ind = dist(generator);
                            contNetwork.executeEdgeDeletion(ind, maxEdgesDelete);
                            maxEdgesDelete--;
                            //std::cout <<"check4" <<std::endl;


                               //std::cout<<lemon::countEdges(network)<<std::endl;
                        }
                    }
                    std::cout <<"executed" <<std::endl;
                    flag = false;
                }
            }
        }
        timeSteps.push_back(t);
        degreeDistr.push_back(contNetwork.getDegreeDistribution());

    }

}

/*void NSA::updateContacts(double &tStart, double tEnd, ContactNetwork & contNetwork, double epsilon)
{// here we consider only two reactions: assembling and disassembling of contact (edge)
    // only deletion can be critical since it can lead to the negative result
    /*std::vector<int> nu = {-1, 1}; // -1 for edge disappearance, 1 for adding
    double time = tStart;

    std::vector<double> propensities(nu.size());*/
    /*std::cout << "check1  " << std::endl;
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
                {"edge_add", contNetwork.getAmountOfEdgesToAddSafe()}
        };
        std::cout << "loooollll  " << std::endl;
        for (auto &it: propensities )
        {
            double p = 1;

            std::cout << "prop  " << it.second<<std::endl;
            if (it.second * tau <= kMax.at(it.first))
            {
                //std:: cout <<kMax.at(it.first)<<std::endl;
                p = it.second * tau / kMax.at(it.first);
                //std::cout << "loololololololol p " << p<<std::endl;
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
}*/

/*double  NSA::proposeTimestep(double epsilon, ContactNetwork & contNetwork) const
{

    //std:: cout <<"edges="<< contNetwork.countEdges() <<" delrsum="<< contNetwork.getEdgeDeletionRateSum()<<std::endl;

    size_t nEdgesExist = contNetwork.countEdges();
    size_t nEdgesAdd = contNetwork.getAmountOfEdgesToAddSafe();
    size_t nNodes = contNetwork.size();
    //double tau = 1.0/ (contNetwork.getEdgeDeletionRateSum() + contNetwork.getEdgeAdditionRateSum());
    double tau = std::min(epsilon * nEdgesExist /contNetwork.getEdgeDeletionRateSum(),
                         // epsilon * (nNodes * (nNodes - 1) / 2 - nEdgesExist)  / contNetwork.getEdgeAdditionRateSum());
                          epsilon * contNetwork.getAmountOfEdgesToAdd()  / contNetwork.getEdgeAdditionRateSum());
    std::cout << "edge del.r.s.: " << contNetwork.getEdgeDeletionRateSum() << std::endl;
    std::cout << "edge add.r.s.: " << contNetwork.getEdgeAdditionRateSum() << std::endl;
    std::cout << "ed. exist: " << nEdgesExist << std::endl;
    std::cout << "ed. to_add safe:" << nEdgesAdd << std::endl;
    std::cout << "ed. to_add:" << contNetwork.getAmountOfEdgesToAdd() << std::endl;
   /* std::cout << "rel1: " << nEdgesExist /contNetwork.getEdgeDeletionRateSum() << std::endl;
    std::cout << "rel2: " << nEdgesAdd / contNetwork.getEdgeAdditionRateSum() << std::endl;*/

    /*std::cout << "tau: " << tau << std::endl;
    return tau;
}*/

double NSA::sampleRandUni()
{
    double r = randuni(generator);
    while (r == 0)
    {
        r = randuni(generator);
    }
    return r;
}

/*int NSA::sampleRandBinomial(int nTrials, double p)
{
    std::binomial_distribution<> randBinomial(nTrials, p);
    int result  = randBinomial(generator);
    return result;
}*/


void NSA::executeReaction(ContactNetwork & contNetwork, std::string reactId,
                          double rStart, double rBound, double time, uint32_t &nInf)
{
    //std::cout << reactId <<std::endl;
    if (reactId == "edge_del")
    {
        contNetwork.executeEdgeDeletion(rStart, rBound);
        //contNetwork.executeEdgeDelitionUniform();
    }

    else if (reactId == "edge_add")
    {
        contNetwork.executeEdgeAddition(rStart, rBound);
        //contNetwork.executeEdgeAdditionUniform();

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


/*void NSA::BDtauleap(double tStart, double tEnd, ContactNetwork & contNetwork, double epsilon,
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
                {"edge_add", contNetwork.getAmountOfEdgesToAddSafe()}
        };

        for (auto &it: propensities )
        {
            double p = 1;

            std::cout << "prop  " <<it.first<<" "<< it.second<<std::endl;
            if (it.second * tau <= kMax.at(it.first))
            {
                p = it.second * tau / kMax.at(it.first);
                //std::cout << "p  " << p<<std::endl;
            }
            std::cout << "p  " << p<<std::endl;

            k.at(it.first) = sampleRandBinomial(kMax.at(it.first), p);
        }

        std::vector<std::string> order;

        //size_t i = 0;
        for (auto &it: k )
        {
            // i++;
            order.insert(order.end(), it.second, it.first);
        }

        //std::cout << order.size() <<std::endl;

        std::shuffle(order.begin(), order.end(), generator);
        std:: cout << "reactions fired:" << order.size() << std::endl;
        std:: cout << "deletion fired:" << k.at("edge_del") << std::endl;
        std:: cout << "addition fired:" << k.at("edge_add") << std::endl;
        std::cout << "flag" <<std::endl;
        for (size_t i = 0; i < order.size(); i++)
        {
            std::cout <<order.at(i) << " ";

            double r = sampleRandUni();

            uint32_t  ni = 0;

            //executeReaction(contNetwork, order.at(i), 0, r * propensities.at(order.at(i)), time, ni);
            if (order.at(i) == "edge_add")
            {
               // contNetwork.executeEdgeAdditionUniform();
            }
            if (order.at(i) == "edge_del")
            {

               // contNetwork.executeEdgeDeletionUniform();
            }




            //propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum();
            //propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum();
        }

        /*propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum();
        propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum();

        std::cout <<std::endl;
        std:: cout << "------------" << std::endl;
        timeSteps.push_back(time);

        degreeDistr.push_back(contNetwork.getDegreeDistribution());
        //std:: cout <<"---------------"<<std::endl;

    }

}*/