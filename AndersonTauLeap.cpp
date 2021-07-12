//
// Created by Malysheva, Nadezhda on 21.07.20.
//

#include "AndersonTauLeap.h"
#include "Utility.h"
#include "types.h"
#include <fstream>

void AndersonTauLeap(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork, double epsilon,
                     NetworkStorage &nwStorage,
                     const std::string &saveDegreeDistMode, std::mt19937_64 &generator/*, std::vector<BenStructure> &benToFile*/)
{
    std::vector<BenStructure> benToFile = contNetwork.getBenStructure(0);

    int N = 2; // number of reactants -
    int M = 2; //number of reactions, size of propensity vector

    std::vector<double> T(M, 0);
    std::vector<int> C(M, 0);

    std::vector<std::vector<std::pair<double, int>>> S;
    S.reserve(1e3 + 1);
    std::vector<std::pair<double, int>> temp = {{0.0, 0}};
    S.push_back(temp);
    S.push_back(temp);

    std::vector<double> propensities(M, 0);
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propDel = contNetwork.getEdgeDeletionRateSum();
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propAdd = contNetwork.getEdgeAdditionRateSum();

    propensities.at(0) = propDel.at(propDel.size() - 1).first;
    propensities.at(1) = propAdd.at(propAdd.size() - 1).first;

    std::vector<int> row(M, 0);

    std::vector<int> NN(M, 0);

    double t = tLastNetworkUpdate;

    double p = 0.75;
    double p1 = 0.9;
    double q = 0.98;

    /*if (saveDegreeDistMode == "c")
    {
        nwStorage.emplace_back(t, contNetwork.getNetworkState());
    }*/
    std::vector<std::vector<int>> nu = {{-1, 1}, {1, -1}};
    std::vector<size_t> X = {propDel.size() - 1, propAdd.size() - 1};
    double tau = getTau(N, nu, propensities, epsilon, X);
    while (t < tEnd)
    {
        if (propensities.at(0) + propensities.at(1) == 0)
        {
            t = tEnd;
            //tLastNetworkUpdate = t;
            /*if (saveDegreeDistMode == "c")
            {
                nwStorage.emplace_back(t, contNetwork.getNetworkState());
            }*/
            break;
        }

        /*if (t + tau > tEnd)
        {
            tau = tEnd - t;
        }*/

        //TODO !! CHECK THE CONDITION
        if (tau < 10.0 / (propensities.at(0) + propensities.at(1)))
        {
            executeSSA(100, tEnd, contNetwork, t, tLastNetworkUpdate, nwStorage,
                       saveDegreeDistMode, generator, T, C, S);
            propDel = contNetwork.getEdgeDeletionRateSum();
            propAdd = contNetwork.getEdgeAdditionRateSum();

            propensities.at(0) = propDel.at(propDel.size() - 1).first;
            propensities.at(1) = propAdd.at(propAdd.size() - 1).first;
            X = {propDel.size() - 1, propAdd.size() - 1};
            tau = getTau(2, nu, propensities, epsilon, X);
        }

        else
        {
            for (int i = 0; i < M; i ++)
            {
                std::vector<std::pair<double, int>> Sk = S.at(i);
                size_t B = Sk.size() - 1;
                if (propensities.at(i) * tau + T.at(i) >= Sk.at(B).first)
                {
                    //std::cout << "poisson" << std::endl;
                    std::poisson_distribution<int> poiss(propensities.at(i) * tau + T.at(i) - Sk.at(B).first);

                    int aaa = poiss(generator);
                    NN.at(i) = aaa + (Sk.at(B).second - C.at(i));
                    row.at(i) = B;
                }
                else
                {
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
                    if (index ==0)
                    {
                        std::cout << "index 0" << std::endl;
                    }

                    double r = (T.at(i) + propensities.at(i) * tau - Sk.at(index - 1).first) / (Sk.at(index).first - Sk.at(index - 1).first);


                    std::binomial_distribution<> binom(Sk.at(index).second - Sk.at(index - 1).second, r);

                    int aaa = binom(generator);
                    NN.at(i) = aaa + Sk.at(index - 1).second - C.at(i);
                    row.at(i) = index - 1;

                    if (NN.at(i) < 0)
                    {
                        //std::cout << "less 0; aaa =" << aaa << " SK2 minus " << Sk.at(index).second - Sk.at(index - 1).second << ", NN =" << NN.at(i) << "; Ci = " << C.at(i)<<std::endl;
                        throw std::domain_error("NN less than zero2");
                    }

                }
            }

            /*int change = NN.at(1) - NN.at(0);
            bool pass = std::abs(change) <= std::max(epsilon * (propDel.size() - 1), 1.0) &&
                        std::abs(change) <= std::max(epsilon * (propAdd.size() - 1), 1.0);*/

            bool pass = std::abs(NN.at(0)) <= std::max(epsilon * (propDel.size() - 1), 1.0) &&
                        std::abs(NN.at(1)) <= std::max(epsilon * (propAdd.size() - 1), 1.0);

            if (pass)
            {
                if (t + tau > tEnd)
                {
                    //std::cout << "exceed" << std::endl;
                    t = tEnd;
                    //tLastNetworkUpdate = t;
                    /*if (saveDegreeDistMode == "c")
                    {
                        nwStorage.emplace_back(t, contNetwork.getNetworkState());
                    }*/
                    break;
                }
                //std::cout << "accept" << std::endl;
                t = t + tau;
                for (int i = 0; i < M; i ++)
                {
                    S.at(i).erase(S.at(i).begin(),  S.at(i).begin() + row.at(i) + 1);
                    S.at(i).emplace(S.at(i).begin(), std::make_pair(propensities.at(i) * tau + T.at(i), C.at(i) + NN.at(i)));

                    T.at(i) = propensities.at(i) * tau + T.at(i);
                    C.at(i) = C.at(i) + NN.at(i);

                }

                /*bool pass2 = std::abs(change) <= std::max(0.75 * epsilon * contNetwork.countEdges(), 1.0) &&
                             std::abs(change) <= std::max(0.75 * epsilon * (propAdd.size() - 1), 1.0);*/
                bool pass2 = std::abs(NN.at(0)) <= std::max(0.75 * epsilon * contNetwork.countEdges(), 1.0) &&
                             std::abs(NN.at(1)) <= std::max(0.75 * epsilon * (propAdd.size() - 1), 1.0);

                if (pass2)
                {
                    tau = std::pow(tau, q);
                }
                else
                {
                    tau = tau * p1;
                }

                if (contNetwork.countEdges() + NN.at(1) < NN.at(0))
                {
                    std::string msg = "ERROR: INVALID del amnt!";
                     throw std::domain_error(msg);
                }

                //std::cout << "add: "<< NN.at(1) << ", del = " << NN.at(0) << std::endl;
                updateNetwork2(benToFile, NN, contNetwork.countEdges(), generator, propAdd, propDel, t, contNetwork, propensities);
                tLastNetworkUpdate = t;
                propDel = contNetwork.getEdgeDeletionRateSum();
                propAdd = contNetwork.getEdgeAdditionRateSum();

                propensities.at(0) = propDel.at(propDel.size() - 1).first;
                propensities.at(1) = propAdd.at(propAdd.size() - 1).first;

                /*if (saveDegreeDistMode == "c")
                {
                    nwStorage.emplace_back(t, contNetwork.getNetworkState());
                }*/

            }
            else
            {
                //std::cout << "reject" << std::endl;
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
}


/*void updateNetwork(std::vector<BenStructure> &benToFile, std::vector<int> k, int nDel, std::mt19937_64 &generator,
                    std::vector<std::pair<double, lemon::ListGraph::Edge>> &propAdd,
                    std::vector<std::pair<double, lemon::ListGraph::Edge>> &propDel, double t,
                    ContactNetwork & contNetwork,
                    std::vector<double> &props)
{
    //std::cout << "start = " << nDel << ", add: "<< k.at(1) << ", del = " << k.at(0) << ", ";
    std::unordered_map<std::string, double> propensities {
            {"edge_del", 0},
            {"edge_add", 0} };


    propensities.at("edge_del") = props.at(0);
    propensities.at("edge_add") = props.at(1);

    int maxEdgesDelete = nDel;
    if (maxEdgesDelete < k.at(0))
    {
        std::cout << "EXceed NDEL!!!" << std::endl;
        for (size_t i = 1; i < propDel.size(); i ++)
        {
            contNetwork.removeEdge(propDel.at(i).second);
        }

        propAdd = contNetwork.getEdgeAdditionRateSum();

        for (size_t i = 0; i < k.at(1) - k.at(0) + maxEdgesDelete; i ++)
        {
            //propAdd = contNetwork.getEdgeAdditionRateSum();
            propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;
            double r = sampleRandUni(generator);
            size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, 0, r * propensities.at("edge_add"));

            if (propAdd.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID update add 1!";
                throw std::domain_error(msg);
            }

            std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);
            propAdd.erase(propAdd.begin() + index);
        }
    }

    else
    {
        propDel = contNetwork.getEdgeDeletionRateSum();

        for (size_t i = 0; i < k.at(0) ; i ++)
        {
            propensities.at("edge_del") = propDel.at(propDel.size() - 1).first;
            double r = sampleRandUni(generator);
            size_t index = binarySearch(propDel, 0, propDel.size() - 1, 0, r * propensities.at("edge_del"));

            if (propDel.at(index).second == lemon::INVALID)
            {
                std::cout << "inv index = " << index << std::endl;
                std::cout << "prop inv = " << propensities.at("edge_del") << std::endl;
                std::string msg = "ERROR: INVALID update del!";
                throw std::domain_error(msg);
            }

            std::pair<int, int> b = contNetwork.removeEdge(propDel.at(index).second);
            propDel.erase(propDel.begin() + index);
            benToFile.emplace_back(t, b.first, b.second, false);
        }

        propAdd = contNetwork.getEdgeAdditionRateSum();
        for (size_t i = 0; i < k.at(1); i ++)
        {
            propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;
            double r = sampleRandUni(generator);

            size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, 0, r * propensities.at("edge_add"));

            if (propAdd.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID update add 2!";
                throw std::domain_error(msg);
            }

            std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);
            propAdd.erase(propAdd.begin() + index);
            //benToFile.push_back(BenStructure(t, b.first, b.second, true));
            benToFile.emplace_back(t, b.first, b.second, true);
        }
    }

}
*/



void executeSSA(size_t n, double tEnd, ContactNetwork & contNetwork, double &t,
                double &tLastNetworkUpdate, NetworkStorage &nwStorage,
                const std::string & saveDegreeDistMode, std::mt19937_64 &generator,
                std::vector<double> &T, std::vector<int> &C,
                std::vector<std::vector<std::pair<double, int>>> &S)

{
    //std::cout << "SSA " << n << std::endl;
    std::vector<double> propensities(2, 0);
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propDel;
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propAdd;

    for (size_t ind = 0; ind < n; ind++)
    //while (t < tEnd)
    {
        propDel = contNetwork.getEdgeDeletionRateSum();
        propAdd = contNetwork.getEdgeAdditionRateSum();

        propensities.at(0) = propDel.at(propDel.size() - 1).first;
        propensities.at(1) = propAdd.at(propAdd.size() - 1).first;

        double propensitiesSum = propensities.at(0) + propensities.at(1);

        if (propensitiesSum == 0)
        {
            t = tEnd;
            tLastNetworkUpdate = tEnd; //used to update netw.Upd.Time
            /*if (saveDegreeDistMode == "c")
            {
                nwStorage.emplace_back(t, contNetwork.getNetworkState());
            }*/
            break;
        }

        double r = sampleRandUni(generator);

        double proposedTime = std::log(1 / r) * 1 / propensitiesSum ;

        if (t + proposedTime > tEnd)
        {
            t = tEnd;
            tLastNetworkUpdate = t;
            /*if (saveDegreeDistMode == "c")
            {
                nwStorage.emplace_back(t, contNetwork.getNetworkState());
            }*/
            break;
        }

        t += proposedTime;
        tLastNetworkUpdate = t; //used to update netw. Upd.Time

        r = sampleRandUni(generator);
        //deletion
        if (propensities.at(0) >= r * propensitiesSum)
        {
            size_t index = binarySearch(propDel, 0, propDel.size() - 1, 0, r * propensitiesSum);

            if (propDel.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID SSA del!";
                throw std::domain_error(msg);
            }
            std::pair<int, int> b = contNetwork.removeEdge(propDel.at(index).second);
            //propDel.erase(propDel.begin() + index);
            //benToFile.push_back(BenStructure(t, b.first, b.second, false));
            C.at(0)++;
        }
        else
        {
            size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, propensities.at(0), propensitiesSum * r);

            if (propAdd.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID SSA add!";
                throw std::domain_error(msg);
            }
            std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);
            //propAdd.erase(propAdd.begin() + index);
            //benToFile.push_back(BenStructure(t, b.first, b.second, true));
            C.at(1)++;
        }

        T.at(0) = T.at(0) + propensities.at(0) * proposedTime;
        T.at(1) = T.at(1) + propensities.at(1) * proposedTime;


        for (int i = 0; i < 2; i ++)
        {
            for (size_t l = 0; l < S.at(i).size(); l ++)
            {
                if (S.at(i).at(l).first <= T.at(i) ||  S.at(i).at(l).second <= C.at(i))
                {

                    S.at(i).erase(S.at(i).begin() + l);
                    l--;
                }


            }
            S.at(i).emplace(S.at(i).begin(), std::make_pair(T.at(i), C.at(i)));
        }
        /*if (saveDegreeDistMode == "c")
        {
            nwStorage.emplace_back(t, contNetwork.getNetworkState());
        }*/

    }
}
double getTau(size_t nParts, std::vector<std::vector<int>> nu, std::vector<double> props,  double epsilon,
              std::vector<size_t> X)
{
    double tau = std::numeric_limits<double>::infinity();
    double gi = 1.0;
    for (size_t i = 0; i < nParts; i ++)
    {
        double mu = 0;
        double sigmaSq = 0;
        for (size_t j = 0; j < 2; j ++)
        {
            mu      += nu.at(j).at(i) * props.at(j);
            sigmaSq += nu.at(j).at(i) * nu.at(j).at(i) * props.at(j);
        }
        double a1 = std::max(epsilon * X.at(i) / gi, 1.0);
        double a2 = std::max(epsilon * X.at(i) / gi, 1.0);
        a2 = a2 * a2;

        double loc_min = std::min(a1 / abs(mu), a2 / sigmaSq);
        tau = std::min(tau, loc_min);

    }
    return tau;
}


void updateNetwork2(std::vector<BenStructure> &benToFile, std::vector<int> k, int nDel, std::mt19937_64 &generator,
                   std::vector<std::pair<double, lemon::ListGraph::Edge>> &propAdd,
                   std::vector<std::pair<double, lemon::ListGraph::Edge>> &propDel, double t,
                   ContactNetwork & contNetwork,
                   std::vector<double> &props)
{
    //std::cout << "start = " << nDel << ", add: "<< k.at(1) << ", del = " << k.at(0) << ", ";
    //std::cout << "NW update!!" << std::endl;
    //std::cout << "add: "<< k.at(1) << ", del = " << k.at(0) << std::endl;
    std::unordered_map<std::string, double> propensities {
            {"edge_del", props.at(0)},
            {"edge_add", props.at(1)} };

    int maxEdgesDelete = nDel;
    if (maxEdgesDelete < k.at(0))
    {
       // std::cout << "EXceed NDEL!!!" << std::endl;
        //std::cout << "del: " << k.at(0) << ", add: " << k.at(1) << ", total: " << maxEdgesDelete << std::endl;
        for (size_t i = 1; i < propDel.size(); i ++)
        {
            std::pair<int, int> b = contNetwork.removeEdge(propDel.at(i).second);

            //for ben structure
            /*std::ofstream out;
            out.open("Ben_Cont_Dyn_NSA.txt", std::ios::app);
            out << t << " " << b.first << " " << b.second << " False" << std::endl ;
            out.close();*/
        }

        propAdd = contNetwork.getEdgeAdditionRateSum();

        for (size_t i = 0; i < k.at(1) - k.at(0) + maxEdgesDelete; i ++)
        {
            //propAdd = contNetwork.getEdgeAdditionRateSum();
            propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;
            double r = sampleRandUni(generator);
            size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, 0, r * propensities.at("edge_add"));

            if (propAdd.at(index).second == lemon::INVALID)
            {
                std::string msg = "ERROR: INVALID update add 1!";
                throw std::domain_error(msg);
            }

            std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);
            propAdd.erase(propAdd.begin() + index);

            //for ben structure
            /*std::ofstream out;
            out.open("Ben_Cont_Dyn_NSA.txt", std::ios::app);
            out << t << " " << b.first << " " << b.second << " True" << std::endl ;
            out.close();*/
        }
    }

    else
    {
        //std::cout << "order" << std::endl;
        std::vector<int> order;
        for (size_t ind = 0; ind < k.size(); ind++)
        {
            order.insert(order.end(), k.at(ind), ind);
        }
        std::shuffle(order.begin(), order.end(), generator);

        /*propDel = contNetwork.getEdgeDeletionRateSum();
        propAdd = contNetwork.getEdgeAdditionRateSum();
        propensities.at("edge_del") = propDel.at(propDel.size() - 1).first;
        propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;
*/

        for (auto i : order)
        {
            //propDel = contNetwork.getEdgeDeletionRateSum();
            //propAdd = contNetwork.getEdgeAdditionRateSum();
            if (i == 0)
            {
               // std::cout << "order: " << i << std::endl;
                propensities.at("edge_del") = propDel.at(propDel.size() - 1).first;
                double r = sampleRandUni(generator);
                size_t index = binarySearch(propDel, 0, propDel.size() - 1, 0, r * propensities.at("edge_del"));

                if (propDel.at(index).second == lemon::INVALID)
                {
                    std::string msg = "ERROR: INVALID update del!";
                    throw std::domain_error(msg);
                }

                std::pair<int, int> b = contNetwork.removeEdge(propDel.at(index).second);

                //for ben structure
                /*std::ofstream out;
                out.open("Ben_Cont_Dyn_NSA.txt", std::ios::app);
                out << t << " " << b.first << " " << b.second << " False" << std::endl ;
                out.close();*/


                propDel.erase(propDel.begin() + index);

                lemon::ListGraph::Edge e = contNetwork.getComplementEdge(b.first, b.second);
                propAdd.emplace_back(propAdd.at(propAdd.size() - 1).first + contNetwork.getEdgeAdditionRate(e), e);
                benToFile.emplace_back(t, b.first, b.second, false);
            }
            else if (i == 1)
            {
                //std::cout << "order: " << i << std::endl;
                propensities.at("edge_add") = propAdd.at(propAdd.size() - 1).first;
                //std::cout << "prop.edge add:" << propensities.at("edge_add") << std::endl;
                double r = sampleRandUni(generator);

                size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, 0, r * propensities.at("edge_add"));

                if (propAdd.at(index).second == lemon::INVALID)
                {
                    std::string msg = "ERROR: INVALID update add 2!";
                    throw std::domain_error(msg);
                }

                std::pair<int, int> b = contNetwork.addEdge(propAdd.at(index).second);
                propAdd.erase(propAdd.begin() + index);

                //for ben structure
                /*std::ofstream out;
                out.open("Ben_Cont_Dyn_NSA.txt", std::ios::app);
                out << t << " " << b.first << " " << b.second << " True" << std::endl ;
                out.close();*/

                lemon::ListGraph::Edge e = contNetwork.getEdge(b.first, b.second);
                //std::cout << "del.rate emplace:" << contNetwork.getEdgeDeletionRate(e) << std::endl;
                propDel.emplace_back(propDel.at(propDel.size() - 1).first + contNetwork.getEdgeDeletionRate(e), e);
                //std::cout << "del.rate emplace:" << propDel.at(propDel.size() - 2).first <<"; "<< propDel.at(propDel.size() - 1).first << std::endl;
                //benToFile.emplace_back(t, b.first, b.second, true);

            }
        }

    }
}