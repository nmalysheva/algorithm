//
// Created by Malysheva, Nadezhda on 13.07.20.
//

#include "PoissonTauLeap.h"
#include "Utility.h"


void PoissonTauleap(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork, double epsilon,
                         std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                         bool updateDegreeDistr, std::mt19937_64 &generator/*, std::vector<BenStructure> &benToFile*/)
{
    std::vector<BenStructure> benToFile = contNetwork.getBenStructure(0);

    int N = 2; // number of reactants -
    int M = 2; //number of reactions, size of propensity vector

    std::unordered_map<std::string, double> propensities {
            {"edge_del", 0},
            {"edge_add", 0} };

    double t = tLastNetworkUpdate;

    updateDegreeDistribution(updateDegreeDistr, t, timeSteps, degreeDistr, contNetwork);

    std::vector<std::pair<double, lemon::ListGraph::Edge>> propDel;
    propDel.reserve(1e6 + 1);
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propAdd;
    propAdd.reserve(1e6 + 1);

    while (t < tEnd)
    {
        //size_t nAdd = 0;
        //size_t nDel = 0;


        //propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum(nDel);
        //propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum(nAdd);

        propDel = contNetwork.getEdgeDeletionRateSum();
        propAdd = contNetwork.getEdgeAdditionRateSum();
        propensities.at("edge_del") = propDel.at(propDel.size() - 1).first;
        propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;

        size_t nDel = propDel.size() - 1;
        size_t nAdd = propAdd.size() - 1;

        if (propensities.at("edge_del") + propensities.at("edge_add") == 0)
        {
            t = tEnd;
            updateDegreeDistribution(updateDegreeDistr, t, timeSteps, degreeDistr, contNetwork);
            break;
        }

        std::vector<double> mu(N, 0);
        std::vector<double> sigmaSq(N, 0);
        std::vector<int> k(M, 0);


        size_t lDel = nDel;

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
            std::cout << "lAdd = " << lAdd << "; lDel = "<< lDel<< std::endl;

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
                //std::cout << "prop.sum: " << (propensities.at("edge_del") + propensities.at("edge_add") ) << std::endl;
                //std::cout << "lAdd: " << lAdd << "; lDel: " << lDel << std::endl;
                std::cout << "SSA. start time =  " << t <<"; ";

                executeSSA(100,  tEnd, contNetwork, t, tLastNetworkUpdate, timeSteps, degreeDistr, updateDegreeDistr, generator);
                flag = false;
            }
            else
            {

                double tau2 = proposeTau2(aCrit, generator);

                double tau = 0;

                size_t kAdd = 0;
                size_t kDel = 0;

                selectTimeStepAndK(tau, propensities, k, kDel, kAdd, tau1, tau2, aCrit, generator);

                size_t nnn = contNetwork.size();
                //If events are captured correctly
                /*if (kDel > lDel || kAdd > lAdd)
                {
                    tau1 = tau1 / 2;
                    // TODO! SPLIT
                }*/

                //size_t nnn = contNetwork.size();
                if (nDel - kDel + kAdd < 0 || nDel - kDel + kAdd > nnn*(nnn-1)/2 )
                {
                    tau1 = tau1 / 2;
                    std::cout <<"TAU/2" << std::endl;
                    // TODO! SPLIT
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

                    std::cout << "k_Del = " << k.at(0) << "; k_Add = " << k.at(1) << std::endl;

                    //updateNetwork(benToFile, k, nDel, generator, propAdd, propDel, t, contNetwork, propensities);

                    updateNetwork2(benToFile, k, nDel, generator, propAdd, propDel, t, contNetwork, propensities, nnn);
                    flag = false;
                }
            }
        }
        updateDegreeDistribution(updateDegreeDistr, t, timeSteps, degreeDistr, contNetwork);
    }
}

double proposeTau1(size_t lDel, size_t lAdd, size_t nAdd, double epsilon, std::vector<double> mu, std::vector<double> sigmaSq)
{
    double timeStep = std::numeric_limits<double>::infinity();


    if (lDel >= 10 || lAdd >= 10)
    {
        timeStep = std::min(timeStep, std::max(epsilon * lDel, 1.0) / std::abs(mu.at(0)));
        timeStep = std::min(timeStep, std::max(epsilon * lDel, 1.0) * std::max(epsilon * lDel, 1.0) / sigmaSq.at(0));
        timeStep = std::min(timeStep, std::max(epsilon * /*lAdd*/ nAdd, 1.0) / std::abs(mu.at(1)));
        timeStep = std::min(timeStep, std::max(epsilon * /*lAdd*/ nAdd, 1.0) * std::max(epsilon * /*lAdd*/ nAdd, 1.0) / sigmaSq.at(1));
    }

    //if (lAdd >= 10)
    // {
    //    timeStep = std::min(timeStep, std::max(epsilon * /*lAdd*/ nAdd, 1.0) / std::abs(mu.at(1)));
    //    timeStep = std::min(timeStep, std::max(epsilon * /*lAdd*/ nAdd, 1.0) * std::max(epsilon * /*lAdd*/ nAdd, 1.0) / sigmaSq.at(1));
    //}

    return timeStep;

}

double proposeTau2(double aCrit, std::mt19937_64 &generator)
{
    double tau2 = std::numeric_limits<double>::infinity();
    if (aCrit > 0)
    {
        double r = sampleRandUni(generator);
        tau2 = 1 / aCrit * std::log(1 / r);
    }
    return tau2;
}

void selectTimeStepAndK(double &tau, const std::unordered_map<std::string, double> &propensities, const std::vector<int> &k,
                             size_t &kDel, size_t &kAdd,double tau1, double tau2, double aCrit, std::mt19937_64 &generator)
{
    tau = 0;

    if (tau1 < tau2)
    {
        //std::cout << "tau1: " << tau1 << std::endl;
        tau = tau1;

    }
    else
    {
        //std::cout << "tau2: " << tau2 << std::endl;
        tau = tau2;
        double r = sampleRandUni(generator);

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
        std::cout << "samplek1" << (propensities.at("edge_del") * tau) << std::endl;

    }

    if (k.at(1) == -1)
    {
        std::poisson_distribution<size_t> poiss(propensities.at("edge_add") * tau);
        kAdd = poiss(generator);
        std::cout << "samplek2" << (propensities.at("edge_add") * tau) << std::endl;
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

void executeSSA(size_t n, double tEnd, ContactNetwork & contNetwork, double &t,
                     double &tLastNetworkUpdate, std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                     bool updateDegreeDistr,  std::mt19937_64 &generator)
{
    std::unordered_map<std::string, double> propensities {
            {"edge_del", 0},
            {"edge_add", 0} };

    std::vector<std::pair<double, lemon::ListGraph::Edge>> propDel;
    propDel.reserve(1e6 + 1);
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propAdd;
    propAdd.reserve(1e6 + 1);

    for (size_t i = 0; i < n; i++)
    {
        propDel = contNetwork.getEdgeDeletionRateSum();
        propAdd = contNetwork.getEdgeAdditionRateSum();

        propensities.at("edge_del") = propDel.at(propDel.size() - 1).first;
        propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;

        if (propensities.at("edge_add") > 1e+6)
        {
            int a = 100;
        }
        std::cout << "propensities: " << propensities.at("edge_del") << ", " << propensities.at("edge_add") <<"; " << propAdd.size() <<std::endl;
        //propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum(nDel);
        //propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum(nAdd);

        double propensitiesSum = propensities.at("edge_del") + propensities.at("edge_add");

        if (propensitiesSum == 0)
        {
            tLastNetworkUpdate = tEnd; //used to update netw.Upd.Time
            t = tEnd;
            updateDegreeDistribution(updateDegreeDistr, t, timeSteps, degreeDistr, contNetwork);
            break;
        }

        double r = sampleRandUni(generator);

        double proposedTime = 1 / propensitiesSum * std::log(1 / r);

        if (t + proposedTime > tEnd)
        {
            tLastNetworkUpdate = t;
            t = tEnd;
            updateDegreeDistribution(updateDegreeDistr, t, timeSteps, degreeDistr, contNetwork);
            break;
        }

        t += proposedTime;
        tLastNetworkUpdate = t; //used to update netw. Upd.Time

        r = sampleRandUni(generator);
        //deletion
        if (propensities.at("edge_del") >= r * propensitiesSum)
        {
            std::cout << i << "=del, ";
            //BenStructure b(t, -1, -1, false);
            //contNetwork.executeEdgeDeletion(0, r * propensitiesSum/*, b*/);
            size_t index = binarySearch(propDel, 0, propDel.size() - 1, 0, r * propensitiesSum);

            if (propDel.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID SSA del!";
                throw std::domain_error(msg);
            }
            std::pair<int, int> b = contNetwork.removeEdge(propDel.at(index).second);
            propDel.erase(propDel.begin() + index);
            //benToFile.push_back(BenStructure(t, b.first, b.second, false));
        }
        else
        {
            std::cout << i << "=add, ";

            size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, propensities.at("edge_del"), propensitiesSum * r);

            if (propAdd.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID SSA add!";
                throw std::domain_error(msg);
            }
            std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);
            propAdd.erase(propAdd.begin() + index);
            //benToFile.push_back(BenStructure(t, b.first, b.second, true));
        }
        updateDegreeDistribution(updateDegreeDistr, t, timeSteps, degreeDistr, contNetwork);

    }
}

//TODO change contNetwork refernce to const
void updateDegreeDistribution(bool updateDegreeDistr, double t, std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr, /*const*/ ContactNetwork  &contNetwork)
{
    if (updateDegreeDistr)
    {
        timeSteps.push_back(t);
        degreeDistr.push_back(contNetwork.getDegreeDistribution());
    }
}

bool isCapturedCorrectly( size_t var1, size_t var2, size_t limit1, size_t limit2)
{
    bool result = ((var1 <=limit1) && (var2 <=limit2));
    return result;
}

int splitRandomNumber(int kDel, int kAdd, int &N, std::mt19937_64 &generator)
{
    if (kDel > N + kAdd)
    {
        std::uniform_int_distribution<> distrib(0, kAdd);
        int nAddNew = distrib(generator);
        int nAddNew2 = kAdd - nAddNew;

        std::uniform_int_distribution<> distrib2(0, kDel);
        int nDelNew = distrib2(generator);
        int nDelNew2 = kDel - nDelNew;
        int b1 = splitRandomNumber(nDelNew, nAddNew, N, generator);
        N = N +  b1;
        int b2 = splitRandomNumber(nDelNew2, nAddNew2, N, generator);
        N = N +  b2;
        return b1 + b2;

    }

    else
    {
        //N = N + kAdd - kDel;
        return kAdd - kDel;

    }
}

void updateNetwork(std::vector<BenStructure> &benToFile, std::vector<int> k, int nDel, std::mt19937_64 &generator,
                   std::vector<std::pair<double, lemon::ListGraph::Edge>> &propAdd,
                   std::vector<std::pair<double, lemon::ListGraph::Edge>> &propDel,
                   double t,
                   ContactNetwork & contNetwork,
                   std::unordered_map<std::string, double> &propensities)

{
    std::vector<int> order;
    for (size_t ind = 0; ind < k.size(); ind++)
    {
        order.insert(order.end(), k.at(ind), ind);
    }
    std::shuffle(order.begin(), order.end(), generator);

    //int maxEdgesDelete = nDel;
    for (auto i : order)
    {
        if (i == 1)
        {
            double r = sampleRandUni(generator);
            //BenStructure b(t, -1, -1, true);
            //contNetwork.executeEdgeAddition(0, r * propensities.at("edge_add")/*, b*/);
            propAdd = contNetwork.getEdgeAdditionRateSum();
            propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;
            size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, 0, r * propensities.at("edge_add"));

            if (propAdd.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID SSA del!";
                throw std::domain_error(msg);
            }

            std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);
            //propAdd.erase(propAdd.begin() + index);
            //benToFile.push_back(BenStructure(t, b.first, b.second, true));
            benToFile.emplace_back(t, b.first, b.second, true);
            //maxEdgesDelete++;
        }
        else if (i == 0)
        {
            //std::uniform_int_distribution<size_t> dist(0, maxEdgesDelete - 1);
            /*NOTE: we can delete uniformely when all edge del. rates are the same. When we have adaptivity case
             * - we can not do this anymore!!!*/
            //size_t ind = dist(generator);
            //contNetwork.executeEdgeDeletion(ind);
            propDel = contNetwork.getEdgeDeletionRateSum();
            propensities.at("edge_del") = propDel.at(propDel.size() - 1).first;
            double r = sampleRandUni(generator);
            size_t index = binarySearch(propDel, 0, propDel.size() - 1, 0, r * propensities.at("edge_del"));

            if (propDel.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID SSA del!";
                throw std::domain_error(msg);
            }
            propensities.at("edge_del") =
                    propensities.at("edge_del") - contNetwork.getEdgeDeletionRate(propDel.at(index).second);

            std::pair<int, int> b = contNetwork.removeEdge(propDel.at(index).second);
            propDel.erase(propDel.begin() + index);
            //benToFile.push_back(BenStructure(t, b.first, b.second, false));
            benToFile.emplace_back(t, b.first, b.second, false);
            //std::cout <<  b.first << ", " << b.second << std::endl;
            //benToFile.push_back(b);
            //maxEdgesDelete--;
        }
    }
}

void updateNetwork2(std::vector<BenStructure> &benToFile, std::vector<int> k, int nDel, std::mt19937_64 &generator,
                    std::vector<std::pair<double, lemon::ListGraph::Edge>> &propAdd,
                    std::vector<std::pair<double, lemon::ListGraph::Edge>> &propDel, double t,
                    ContactNetwork & contNetwork,
                    std::unordered_map<std::string, double> &propensities,
                    size_t nnn)
{
    int maxEdgesDelete = nDel;
    if (maxEdgesDelete < k.at(0))
    {
        std::cout << "EXceed NDEL!!!" << std::endl;
        for (auto i : propDel)
        {
            contNetwork.removeEdge(i.second);
        }

        /*propAdd = contNetwork.getEdgeAdditionRateSum();
        propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;*/
        for (size_t i = 0; i < k.at(1) - k.at(0) + maxEdgesDelete; i ++)
        {
            double r = sampleRandUni(generator);
            propAdd = contNetwork.getEdgeAdditionRateSum();
            propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;
            size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, 0, r * propensities.at("edge_add"));

            if (propAdd.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID SSA del!";
                throw std::domain_error(msg);
            }

            propensities.at("edge_add") =
                    propensities.at("edge_add") - contNetwork.getEdgeAdditionRate(propAdd.at(index).second);
            std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);
            propAdd.erase(propAdd.begin() + index);
        }
    }

    else
    {
        std::cout << "NORMAL!!!" << std::endl;
        //propDel = contNetwork.getEdgeDeletionRateSum();
        //propensities.at("edge_del") = propDel.at(propDel.size() - 1).first;
        for (size_t i = 0; i < k.at(0) ; i ++)
        {
            //std::cout << "DEL!!!" << std::endl;
            propDel = contNetwork.getEdgeDeletionRateSum();
            propensities.at("edge_del") = propDel.at(propDel.size() - 1).first;
            double r = sampleRandUni(generator);
            size_t index = binarySearch(propDel, 0, propDel.size() - 1, 0, r * propensities.at("edge_del"));

            //std::cout << "test1!!!" << std::endl;
            if (propDel.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID SSA del!";
                throw std::domain_error(msg);
            }

            //std::cout << "index = " << index << ",  size = " << propDel.size() - 1 << std::endl;
            //std::cout << propensities.at("edge_del") << " , ";
            propensities.at("edge_del") =
                    propensities.at("edge_del") - contNetwork.getEdgeDeletionRate(propDel.at(index).second);
            //std::cout << propensities.at("edge_del") << std::endl;
            //std::cout << "test2!!!  " << std::endl;
            std::pair<int, int> b = contNetwork.removeEdge(propDel.at(index).second);
            //std::cout << "b1 = " << b.first << ", b2 = " << b.second<< std::endl;
            //std::cout << "test3!!!" << std::endl;
            propDel.erase(propDel.begin() + index);
            //std::cout << "test4!!!" << std::endl;
            benToFile.emplace_back(t, b.first, b.second, false);
            maxEdgesDelete--;
        }


        //propAdd = contNetwork.getEdgeAdditionRateSum();
        //propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;
        for (size_t i = 0; i < k.at(1); i ++)
        {
            //std::cout << "ADD!!!" << std::endl;
            double r = sampleRandUni(generator);
            propAdd = contNetwork.getEdgeAdditionRateSum();
            propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;
            size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, 0, r * propensities.at("edge_add"));

            if (propAdd.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID SSA add!";
                throw std::domain_error(msg);
            }
            propensities.at("edge_add") =
                    propensities.at("edge_add") - contNetwork.getEdgeAdditionRate(propAdd.at(index).second);

            std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);
            //propDel.erase(propAdd.begin() + index);
            //benToFile.push_back(BenStructure(t, b.first, b.second, true));
            benToFile.emplace_back(t, b.first, b.second, true);
            maxEdgesDelete++;
        }

    }
}