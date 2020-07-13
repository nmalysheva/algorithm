//
// Created by Malysheva, Nadezhda on 13.07.20.
//

#include "MidpointApproximation.h"

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