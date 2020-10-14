//
// Created by Malysheva, Nadezhda on 21.07.20.
//

#include "AndersonTauLeap.h"
#include "Utility.h"
#include <queue>

void AndersonTauLeap(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork, double epsilon,
                     std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                     bool updateDegreeDistr, std::mt19937_64 &generator/*, std::vector<BenStructure> &benToFile*/)
{
    std::vector<BenStructure> benToFile = contNetwork.getBenStructure(0);

    int N = 2; // number of reactants -
    int M = 2; //number of reactions, size of propensity vector

    double eps = epsilon;
    std::vector<double> T(M, 0);
    std::vector<double> C(M, 0);

    std::vector<std::vector<std::pair<double, int>>> S;
    S.reserve(1e3 + 1);
    std::vector<std::pair<double, int>> temp = {{0.0, 0}};
    S.push_back(temp);
    S.push_back(temp);

    std::vector<double> propensities(M, 0);
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propDel;
    propDel.reserve(1e6 + 1);
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propAdd;
    propAdd.reserve(1e6 + 1);
    propDel = contNetwork.getEdgeDeletionRateSum();
    propAdd = contNetwork.getEdgeAdditionRateSum();

    propensities.at(0) = propDel.at(propDel.size() - 1).first;
    propensities.at(1) = propAdd.at(propAdd.size() - 1).first;

    std::vector<int> row(M, 0);

    std::vector<int> NN(M, 0);

    double t = tLastNetworkUpdate;

    double p = 0.75;
    double p1 = 0.9;
    double q = 0.98;
    updateDegreeDistributio(updateDegreeDistr, t, timeSteps, degreeDistr, contNetwork);

    double tau = 0.0023; //select tau acc. to proc.
    while (t < tEnd)
    {
        std::cout << "-------------" << std::endl;
        std::cout << "tau: " << tau << std::endl;
        std::cout << "props: " << propensities.at(0)  << ", " << propensities.at(1) << std::endl;

        if (fmod(t, 0.2) == 0)
        {
            executeSSA(50, tEnd, contNetwork, t, tLastNetworkUpdate, timeSteps, degreeDistr,
                       updateDegreeDistr, generator, T, C, S);
            continue;
        }
        for (int i = 0; i < M; i ++)
        {
            std::vector<std::pair<double, int>> Sk = S.at(i);
            int B = Sk.size() - 1;
            if (propensities.at(i) * tau + T.at(i) >= Sk.at(B).first)
            {
                std::cout << "lol1: " << propensities.at(i) * tau + T.at(i) - Sk.at(B).first << std::endl;
                std::poisson_distribution<int> poiss(propensities.at(i) * tau + T.at(i) - Sk.at(B).first);


                NN.at(i) = poiss(generator);
                std::cout << "poiss: " << NN.at(i) << std::endl;
                NN.at(i) = NN.at(i) + Sk.at(B).second - C.at(i);
                row.at(i) = B;
            }
            else
            {
                std::cout << "lol2: ";
                int index = 0;
                for (int k = 1; k < Sk.size(); k++)
                {
                    if (Sk.at(k - 1).first <= propensities.at(i) * tau + T.at(i) &&
                        propensities.at(i) * tau + T.at(i) < Sk.at(k).first)
                    {
                        index = k;
                        break;
                    }
                }
                //double r1 = static_cast<double> (T.at(i) + propensities.at(i) * tau - Sk.at(index - 1).first);
                std::cout << "index = " << index << std::endl;


                double r = static_cast<double>((T.at(i) + propensities.at(i) * tau - Sk.at(index - 1).first) /
                                               (Sk.at(index).first - Sk.at(index - 1).first));

                std::cout << "r = " << r << std::endl;
                std::binomial_distribution<> binom(Sk.at(index).second - Sk.at(index - 1).second, r);
                NN.at(i) = binom(generator) + Sk.at(index - 1).second -C.at(i);
                row.at(i) = index - 1;
            }
        }

        std::cout << "NN: " << NN.at(0) << "; " << NN.at(1) << std::endl;
        int change = NN.at(1) - NN.at(0);
        bool pass = std::abs(change) <= std::max(eps * contNetwork.countEdges(), 1.0) &&
                std::abs(change) <= std::max(eps * propAdd.size(), 1.0);
        //std::cout << pass << "; " << std::abs(change) << "; " << propAdd.size() << std::endl;
        if (pass)
        {
            std::cout << "pass" << "; " << std::abs(change) << "; " << propAdd.size() << std::endl;
            //std::cout << "NN: " << NN.at(0) << "; " << NN.at(1) << std::endl;
            t = t + tau;
            tLastNetworkUpdate = t;
            for (int i = 0; i < M; i ++)
            {
                S.at(i).erase(S.at(i).begin() + row.at(i), S.at(i).end());
                S.at(i).emplace(S.at(i).begin(), std::make_pair(propensities.at(i) * tau + T.at(i), C.at(i) + NN.at(i)));

                T.at(i) = propensities.at(i) * tau + T.at(i);
                C.at(i) = C.at(i) + NN.at(i);

            }

            bool pass2 = std::abs(change) <= std::max(0.75 * eps * contNetwork.countEdges(), 1.0) &&
                         std::abs(change) <= std::max(0.75 * eps * propAdd.size(), 1.0);
            if (pass2)
            {
                tau = std::pow(tau, q);
            }
            else
            {
                tau = tau * p1;
            }

            size_t nnn = contNetwork.size();
            updateNetwork2(benToFile, NN, contNetwork.countEdges(), generator, propAdd, propDel, t, contNetwork, propensities, nnn);

            propDel = contNetwork.getEdgeDeletionRateSum();
            propAdd = contNetwork.getEdgeAdditionRateSum();

            propensities.at(0) = propDel.at(propDel.size() - 1).first;
            propensities.at(1) = propAdd.at(propAdd.size() - 1).first;
            updateDegreeDistributio(updateDegreeDistr, t, timeSteps, degreeDistr, contNetwork);

        }
        else
        {
            for (int i = 0; i < M; i ++)
            {
                if (row.at(i) == S.at(i).size() - 1)
                {
                    S.at(i).emplace_back(std::make_pair(propensities.at(i) * tau + T.at(i), C.at(i) + NN.at(i)));
                }
                else
                {
                    S.at(i).emplace(S.at(i).begin() + row.at(i) + 1, std::make_pair(propensities.at(i) * tau + T.at(i), C.at(i) + NN.at(i)));
                }

            }

            tau = tau * p;
        }

    }
}


void updateNetwork2(std::vector<BenStructure> &benToFile, std::vector<int> k, int nDel, std::mt19937_64 &generator,
                    std::vector<std::pair<double, lemon::ListGraph::Edge>> &propAdd,
                    std::vector<std::pair<double, lemon::ListGraph::Edge>> &propDel, double t,
                    ContactNetwork & contNetwork,
                    std::vector<double> &props,
                    size_t nnn)
{
    std::unordered_map<std::string, double> propensities {
            {"edge_del", 0},
            {"edge_add", 0} };


    propensities.at("edge_del") = props.at(0);
    propensities.at("edge_add") = props.at(1);

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


//TODO change contNetwork refernce to const
void updateDegreeDistributio(bool updateDegreeDistr, double t, std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr, /*const*/ ContactNetwork  &contNetwork)
{
    if (updateDegreeDistr)
    {
        timeSteps.push_back(t);
        degreeDistr.push_back(contNetwork.getDegreeDistribution());
    }
}

void exactAlgorithm(size_t n, double tEnd, ContactNetwork & contNetwork, double &t,
                    double &tLastNetworkUpdate, std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                    bool updateDegreeDistr,  std::mt19937_64 &generator,
                    std::vector<double> &T, std::vector<double> &C,
                    std::vector<std::vector<std::pair<double, int>>> &S)
{
    std::vector<double> propensities(2, 0);
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propDel;
    propDel.reserve(1e6 + 1);
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propAdd;
    propAdd.reserve(1e6 + 1);
    propDel = contNetwork.getEdgeDeletionRateSum();
    propAdd = contNetwork.getEdgeAdditionRateSum();


    std::vector<double> P(2, 0);
    for (size_t i = 0; i < P.size(); i++)
    {
        if (S.at(i).empty())
        {
            double r = sampleRandUni(generator);
            P.at(i) = T.at(i) +  std::log(1 / r);
        }

        else
        {
            std::vector<std::pair<double, int>> Sk = S.at(i);
            int B = Sk.size() - 1;
            int index = -1;
            for (int l = 0; l < Sk.size() - 1; l++)
            {
                if (C.at(i) == Sk.at(l).second)
                {
                    if (C.at(i) < Sk.at(l + 1).second)
                    {
                        index = l;
                    }
                }
            }

            if (index == -1)
            {
                double r = sampleRandUni(generator);
                P.at(i) = Sk.at(B).first +  std::log(1 / r);
            }
            else
            {
                double r = sampleRandUni(generator);
                double tk = Sk.at(index + 1).first - Sk.at(index).first;
                int Nk = Sk.at(index + 1).second - Sk.at(index).second;
                P.at(i) = Sk.at(index).first +  tk * (1 - std::pow(r, 1/static_cast<double>(Nk)));
            }

        }


    }


    std::priority_queue<std::pair<double, int>>delta;
    for (size_t i = 0; i < 100; i++)
    {
        for (int i = 0; i < 2; i++)
        {
            auto pr = std::make_pair(((P.at(i) - T.at(i)) / propensities.at(i)), i);
            delta.push(pr);
        }
        auto pr = delta.top();


    }



}

void executeSSA(size_t n, double tEnd, ContactNetwork & contNetwork, double &t,
                double &tLastNetworkUpdate, std::vector<double> &timeSteps, std::vector<std::vector<size_t>> &degreeDistr,
                bool updateDegreeDistr, std::mt19937_64 &generator,
                std::vector<double> &T, std::vector<double> &C,
                std::vector<std::vector<std::pair<double, int>>> &S)

{
    std::vector<double> propensities(2, 0);
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propDel;
    propDel.reserve(1e6 + 1);
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propAdd;
    propAdd.reserve(1e6 + 1);

    for (size_t i = 0; i < n; i++)
    {
        propDel = contNetwork.getEdgeDeletionRateSum();
        propAdd = contNetwork.getEdgeAdditionRateSum();

        propensities.at(0) = propDel.at(propDel.size() - 1).first;
        propensities.at(1) = propAdd.at(propAdd.size() - 1).first;

        if (propensities.at(1) > 1e+6)
        {
            int a = 100;
        }
        std::cout << "propensities: " << propensities.at(0) << ", " << propensities.at(1) <<"; " << propAdd.size() <<std::endl;
        //propensities.at("edge_del") = contNetwork.getEdgeDeletionRateSum(nDel);
        //propensities.at("edge_add") = contNetwork.getEdgeAdditionRateSum(nAdd);

        double propensitiesSum = propensities.at(0) + propensities.at(1);

        if (propensitiesSum == 0)
        {
            tLastNetworkUpdate = tEnd; //used to update netw.Upd.Time
            t = tEnd;
            updateDegreeDistributio(updateDegreeDistr, t, timeSteps, degreeDistr, contNetwork);
            break;
        }

        double r = sampleRandUni(generator);

        double proposedTime = 1 / propensitiesSum * std::log(1 / r);

        if (t + proposedTime > tEnd)
        {
            tLastNetworkUpdate = t;
            t = tEnd;
            updateDegreeDistributio(updateDegreeDistr, t, timeSteps, degreeDistr, contNetwork);
            break;
        }

        t += proposedTime;
        tLastNetworkUpdate = t; //used to update netw. Upd.Time

        r = sampleRandUni(generator);
        //deletion
        if (propensities.at(0) >= r * propensitiesSum)
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
            C.at(0)++;
        }
        else
        {
            std::cout << i << "=add, ";

            size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, propensities.at(0), propensitiesSum * r);

            if (propAdd.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID SSA add!";
                throw std::domain_error(msg);
            }
            std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);
            propAdd.erase(propAdd.begin() + index);
            //benToFile.push_back(BenStructure(t, b.first, b.second, true));
            C.at(1)++;
        }

        T.at(0) = T.at(0) + propensities.at(0) * proposedTime;
        T.at(1) = T.at(1) + propensities.at(1) * proposedTime;

        for (int i = 0; i < 2; i ++)
        {
            for (int l = 0; l < S.at(i).size(); l ++)
            {
                if (S.at(i).at(l).first <= T.at(i))
                {

                    S.at(i).erase(S.at(i).begin() + l);
                    l--;
                }
            }
            S.at(i).emplace(S.at(i).begin(), std::make_pair(T.at(i), C.at(i)));
        }

        updateDegreeDistributio(updateDegreeDistr, t, timeSteps, degreeDistr, contNetwork);

    }
}