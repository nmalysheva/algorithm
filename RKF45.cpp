//
// Created by Malysheva, Nadezhda on 13.07.20.
//

#include "RKF45.h"
#include "Utility.h"


void RKF45Approximation(double &tLastNetworkUpdate, double tEnd, ContactNetwork & contNetwork,
                             double dtMax, double dtMin, double errorMax, double errorMin, std::vector<double> &timeSteps,
                             std::vector<std::vector<size_t>> &degreeDistr, bool updateDegreeDistr/*, std::vector<BenStructure> &benToFile*/, std::mt19937_64 generator)
{
    std::vector<std::vector<double>> a({{1.0/4.0,        0.0,           0.0,           0.0,          0.0},
                                        {3.0/32.0,       9.0/32.0,      0.0,           0.0,          0.0},
                                        {1932.0/2197.0, 7200.0/2197.0, 7296.0/2197.0, 0.0,          0.0},
                                        {439.0/216.0,   8.0,           3680.0/513.0, 845.0/4104.0, 0.0},
                                        {8.0/27.0,      2.0,          3544.0/2565.0, 1859.0/4104.0, 11.0/40.0}});


    std::vector<double> b({25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0});
    std::vector<double> d({16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, 9.0/50.0, 2.0/55.0});


    std::vector<BenStructure> benToFile = contNetwork.getBenStructure(0);

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
            //int nEdges =  contNetwork.countEdges();

            if (X + m < 0 || m > Xmax)
            {
                dt = dt/2;
                throw std::domain_error("LOOOOL");
            }
            else
            {
                size_t nSurv = contNetwork.updateSurvivalProbability(mDel, mAdd, benToFile, t);
                //std::cout << "nSurv = " << nSurv << std::endl;

                //std::cout << "nEdges = " << contNetwork.countEdges() << std::endl;
                //std::cout << "to add = " << abs(X + m - nSurv) << std::endl;
                int change = X + m - nSurv;

                if (change > 0)
                {
                    propAdd = contNetwork.getEdgeAdditionRateSum();
                    nAdd = propAdd.size() - 1;
                    double nnn = propAdd.at(propAdd.size() - 1).first;
                    for (int i = 0; i < change; i ++)
                    {
                        //propAdd = contNetwork.getEdgeAdditionRateSum();

                        //double nnn = propAdd.at(propAdd.size() - 1).first;

                        double r = sampleRandUni(generator);
                        /*BenStructure b (t, -1, -1, true);
                        contNetwork.executeEdgeAddition(0, r * nnn, b);
                        benToFile.push_back(b);*/

                        //contNetwork.executeEdgeAddition(0, r * nnn);

                        //lemon::ListGraph::Edge e = binarySearch(propAdd, 0, propAdd.size() - 1, 0, r * nnn);
                        size_t index = binarySearch(propAdd, 0, propAdd.size() - 1, 0, r * nnn);

                        if (propAdd.at(index).second == lemon::INVALID)
                        {
                            std::string msg = "ERROR: INVALID SSA del!";
                            throw std::domain_error(msg);
                        }

                        std::pair<int, int> b = contNetwork.removeEdge(propAdd.at(index).second);
                        propDel.erase(propAdd.begin() + index);
                        benToFile.push_back(BenStructure(t, b.first, b.second, true));
                    }

                }
                else
                {
                    //int nDel =  contNetwork.countEdges();
                    propDel = contNetwork.getEdgeDeletionRateSum();
                    nDel = propDel.size() - 1;
                    double nnn = propDel.at(propDel.size() - 1).first;
                    for (int i = 0; i < abs(change); i ++)
                    {
                        double r = sampleRandUni(generator);
                        //lemon::ListGraph::Edge e = binarySearch(propDel, 0, propDel.size() - 1, 0, r * nnn);
                        size_t index = binarySearch(propDel, 0, propDel.size() - 1, 0, r * nnn);

                        if (propDel.at(index).second == lemon::INVALID)
                        {
                            std::string msg = "ERROR: INVALID SSA del!";
                            throw std::domain_error(msg);
                        }
                        std::pair<int, int> b = contNetwork.removeEdge(propDel.at(index).second);
                        propDel.erase(propDel.begin() + index);
                        benToFile.push_back(BenStructure(t, b.first, b.second, false));
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

    //-----------ben format
    /* std::string fileNameBen = "Ben_ContDyn_RKF45.txt";
     std::ofstream benFile;
     benFile.open(fileNameBen);
     for (auto &it: benToFile)
     {
         std::string stateStr = "";
         if (it.state)
         {
             stateStr = "True";
         }
         else
         {
             stateStr = "False";
         }
         benFile << it.t << " " << it.u << " " << it.v << " " << stateStr << std::endl;
     }
     benFile.close();*/
}

double getNadd(double X, double Y, double lam, double theta, std::vector<std::vector<double>> a, std::vector<double> d, double dt)
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

double getNdel(double X, double Y, double lam, double theta, std::vector<std::vector<double>> a, std::vector<double> d, double dt)
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