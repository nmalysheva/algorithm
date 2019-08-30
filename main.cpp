#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/adaptors.h>
#include "UniqueID.h"
#include <lemon/maps.h>
#include "ContactNetwork.h"
#include "SSA.h"
#include "NSA.h"
#include <chrono>
#include <unistd.h>
#include <fstream>
#include <string>

using namespace lemon;
using namespace std;

void executeSSA(double tStart, double tEnd, size_t nInfected, size_t nSusceptible, size_t nEdges,int maxContactsA, int MaxContactsB,
        double transmRate, double newContRate, double looseContRate, double dRate, double bRate, size_t simulationNumber)
{
    UniqueID().reset();
    ContactNetwork contNetwork(nInfected,
                               nSusceptible,
                               nEdges,
                               maxContactsA,
                               MaxContactsB,
                               transmRate,
                               newContRate,
                               looseContRate,
                               dRate,
                               bRate);

    size_t nPopulation = contNetwork.size();
    size_t startEdges = contNetwork.countEdges();
    std::cout <<"Nodes: "  << nPopulation <<std::endl;
    std::cout <<"Edges: "  << startEdges <<std::endl;

    std::vector<double> timeSteps;
    std::vector<uint32_t> infectedSteps;
    std::vector<std::vector<size_t>> degreeDistr; // vector stores degree dist. per time step

    SSA ssa;
    auto start_time = std::chrono::high_resolution_clock::now();
    ssa.execute(tStart, tEnd, contNetwork, timeSteps, infectedSteps, degreeDistr);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    ofstream newFile;
    newFile.open("SSA_" + std::to_string(nPopulation) + "_" + std::to_string(simulationNumber) + ".txt");

    //newFile << "duration in CPU time: " << time;
    newFile << "duration in milliseconds: " << chrono::duration <double, milli> (time).count() << std::endl;
    newFile << "nodes: " << nPopulation << std::endl;
    newFile << "edges: " << startEdges << std::endl;
    newFile << "infected: " << nInfected << std::endl;

    newFile << "max contact range: " << maxContactsA << " "  << MaxContactsB << std::endl;
    newFile << "transmission rate: " << transmRate << std::endl;
    newFile << "rate of make a new contact: " << newContRate << std::endl;
    newFile << "rate of loose a contact: " << looseContRate << std::endl;
    newFile << "death rate: " << dRate << std::endl;
    newFile << "birth rate: " << bRate << std::endl;

    for (size_t i = 0; i < timeSteps.size(); i++)
    {
        newFile << timeSteps.at(i) << ' ';
    }
    newFile << std::endl;

    for (size_t i = 0; i < infectedSteps.size(); i++)
    {
        newFile << infectedSteps.at(i) << ' ';
    }
    newFile << std::endl;
    newFile << "degree distribution:" << std::endl;
    for (size_t i = 0; i < degreeDistr.size(); i++)
    {
        std::vector<size_t> degrreDatTime = degreeDistr.at(i);
        for (size_t j = 0;  j< degrreDatTime.size(); j++)
        {
            newFile << degrreDatTime.at(j) << ' ';
        }
        newFile << std::endl;
    }

    newFile.close();

}

void executeNSA(double tStart, double tEnd, size_t nInfected, size_t nSusceptible, size_t nEdges,int maxContactsA, int MaxContactsB,
                double transmRate, double newContRate, double looseContRate, double dRate, double bRate, double epsilon, size_t simulationNumber)
{
    UniqueID().reset();
    ContactNetwork contNetwork(nInfected,
                               nSusceptible,
                               nEdges,
                               maxContactsA,
                               MaxContactsB,
                               transmRate,
                               newContRate,
                               looseContRate,
                               dRate,
                               bRate);

    size_t nPopulation = contNetwork.size();
    size_t startEdges = contNetwork.countEdges();
    std::cout <<"Nodes: "  << nPopulation <<std::endl;
    std::cout <<"Edges: "  << startEdges <<std::endl;

    std::vector<double> timeSteps;
    std::vector<uint32_t> infectedSteps;
    std::vector<std::vector<size_t>> degreeDistr; // vector stores degree dist. per time step

    size_t nRejections = 0;
    size_t nAcceptance = 0;
    size_t nThin = 0;


    NSA nsa;
    //ContactNetwork cNw2(20, amount.at(i) - 20, nEdges, 1, 20, 0.03, 1.2, 1.25, 0.0004, 2);
    auto start_time = std::chrono::high_resolution_clock::now();
    nsa.execute(tStart, tEnd, contNetwork, timeSteps, infectedSteps, degreeDistr, epsilon,nRejections, nAcceptance, nThin);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    ofstream newFile;
    newFile.open("NSA_" + std::to_string(nPopulation) + "_" + std::to_string(simulationNumber) + ".txt");

    newFile << "duration in milliseconds: " << chrono::duration <double, milli> (time).count() << std::endl;
    newFile << "nodes: " << nPopulation << std::endl;
    newFile << "edges: " << startEdges << std::endl;
    newFile << "infected: " << nInfected << std::endl;

    newFile << "max contact range: " << maxContactsA << " "  << MaxContactsB << std::endl;
    newFile << "transmission rate: " << transmRate << std::endl;
    newFile << "rate of make a new contact: " << newContRate << std::endl;
    newFile << "rate of loose a contact: " << looseContRate << std::endl;
    newFile << "death rate: " << dRate << std::endl;
    newFile << "birth rate: " << bRate << std::endl;
    newFile << "epsilon: " << epsilon << std::endl;

    for (size_t i = 0; i < timeSteps.size(); i++)
    {
        newFile << timeSteps.at(i) << ' ';
    }
    newFile << std::endl;

    for (size_t i = 0; i < infectedSteps.size(); i++)
    {
        newFile << infectedSteps.at(i) << ' ';
    }
    newFile <<  std::endl;
    newFile << "rejected: " << nRejections << std::endl;
    newFile << "accepted: " << nAcceptance << std::endl;
    newFile << "thined: "   <<nThin << std::endl;
    newFile << "degree distribution:" << std::endl;
    for (size_t i = 0; i < degreeDistr.size(); i++)
    {
        std::vector<size_t> degrreDatTime = degreeDistr.at(i);
        for (size_t j = 0;  j< degrreDatTime.size(); j++)
        {
            newFile << degrreDatTime.at(j) << ' ';
        }
        newFile << std::endl;
    }
    newFile.close();

}

void executeSSAOnlyContactUpdate(double tStart, double tEnd, size_t nInfected, size_t nSusceptible, size_t nEdges,int maxContactsA, int MaxContactsB,
                 double newContRate, double looseContRate,  size_t simulationNumber)
{
    UniqueID().reset();
    ContactNetwork contNetwork(nInfected,
                               nSusceptible,
                               nEdges,
                               maxContactsA,
                               MaxContactsB,
                               0,
                               newContRate,
                               looseContRate,
                               0,
                               0);

    size_t nPopulation = contNetwork.size();
    size_t startEdges = contNetwork.countEdges();
    std::cout <<"Nodes: "  << nPopulation <<std::endl;
    std::cout <<"Edges: "  << startEdges <<std::endl;

    std::vector<double> timeSteps;
    std::vector<uint32_t> infectedSteps;
    std::vector<std::vector<size_t>> degreeDistr; // vector stores degree dist. per time step

    SSA ssa;
    auto start_time = std::chrono::high_resolution_clock::now();
    ssa.execute(tStart, tEnd, contNetwork, timeSteps, infectedSteps, degreeDistr);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    ofstream newFile;
    newFile.open("ConntDyn_SSA_" + std::to_string(nPopulation) + "_" + std::to_string(simulationNumber) + ".txt");

    //newFile << "duration in CPU time: " << time;
    newFile << "duration in milliseconds: " << chrono::duration <double, milli> (time).count() << std::endl;
    newFile << "nodes: " << nPopulation << std::endl;
    newFile << "edges: " << startEdges << std::endl;
    newFile << "infected: " << nInfected << std::endl;

    newFile << "max contact range: " << maxContactsA << " "  << MaxContactsB << std::endl;
    newFile << "rate of make a new contact: " << newContRate << std::endl;
    newFile << "rate of loose a contact: " << looseContRate << std::endl;


    for (size_t i = 0; i < timeSteps.size(); i++)
    {
        newFile << timeSteps.at(i) << ' ';
    }
    newFile << std::endl;

    newFile << std::endl;
    newFile << "degree distribution:" << std::endl;
    for (size_t i = 0; i < degreeDistr.size(); i++)
    {
        std::vector<size_t> degrreDatTime = degreeDistr.at(i);
        for (size_t j = 0;  j< degrreDatTime.size(); j++)
        {
            newFile << degrreDatTime.at(j) << ' ';
        }
        newFile << std::endl;
    }

    newFile.close();

}


void executeNSAOnlyContactUpdate(double tStart, double tEnd, size_t nInfected, size_t nSusceptible, size_t nEdges,int maxContactsA, int MaxContactsB,
                double newContRate, double looseContRate, double epsilon, size_t simulationNumber)
{
    UniqueID().reset();
    ContactNetwork contNetwork(nInfected,
                               nSusceptible,
                               nEdges,
                               maxContactsA,
                               MaxContactsB,
                               0,
                               newContRate,
                               0,
                               0,
                               0);

    size_t nPopulation = contNetwork.size();
    size_t startEdges = contNetwork.countEdges();
    std::cout <<"Nodes: "  << nPopulation <<std::endl;
    std::cout <<"Edges: "  << startEdges <<std::endl;

    std::vector<double> timeSteps;
    std::vector<std::vector<size_t>> degreeDistr; // vector stores degree dist. per time step

    NSA nsa;
    //ContactNetwork cNw2(20, amount.at(i) - 20, nEdges, 1, 20, 0.03, 1.2, 1.25, 0.0004, 2);
    auto start_time = std::chrono::high_resolution_clock::now();
    nsa.BDtauleap(tStart, tEnd, contNetwork, epsilon, timeSteps, degreeDistr);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    ofstream newFile;
    newFile.open("ConntDyn_NSA_" + std::to_string(nPopulation) + "_" + std::to_string(simulationNumber) + ".txt");

    newFile << "duration in milliseconds: " << chrono::duration <double, milli> (time).count() << std::endl;
    newFile << "nodes: " << nPopulation << std::endl;
    newFile << "edges: " << startEdges << std::endl;
    newFile << "infected: " << nInfected << std::endl;

    newFile << "max contact range: " << maxContactsA << " "  << MaxContactsB << std::endl;
    newFile << "rate of make a new contact: " << newContRate << std::endl;
    newFile << "rate of loose a contact: " << looseContRate << std::endl;
    newFile << "epsilon: " << epsilon << std::endl;

    for (size_t i = 0; i < timeSteps.size(); i++)
    {
        newFile << timeSteps.at(i) << ' ';
    }
    newFile << std::endl;

    newFile << "degree distribution:" << std::endl;
    for (size_t i = 0; i < degreeDistr.size(); i++)
    {
        std::vector<size_t> degrreDatTime = degreeDistr.at(i);
        for (size_t j = 0;  j< degrreDatTime.size(); j++)
        {
            newFile << degrreDatTime.at(j) << ' ';
        }
        newFile << std::endl;
    }
    newFile.close();

}



int main(int argc, char* argv[])
{
    if (argc != 15)
    {
        std::cout << "wrong parameters" <<std::endl;
    }
    else
    {

        double tStart = std::strtod(argv[1], 0);
        std::cout<<tStart << std::endl;
        double tEnd   = std::strtod(argv[2], 0);
        std::cout<<tEnd << std::endl;
        size_t nInfected = std::stoi(argv[3]);
        std::cout<<nInfected << std::endl;
        size_t nSusceptible = std::stoi(argv[4]);
        std::cout<<nSusceptible << std::endl;
        size_t nEdges = std::stoi(argv[5]);
        std::cout<<nEdges << std::endl;
        int maxContactsA = std::stoi(argv[6]);
        std::cout<<maxContactsA << std::endl;
        int MaxContactsB = std::stoi(argv[7]);
        std::cout<<MaxContactsB << std::endl;
        double newContRate = std::strtod(argv[8], 0);
        std::cout<<newContRate << std::endl;
        double transmRate = std::strtod(argv[9], 0);
        std::cout<<transmRate << std::endl;
        double looseContRate = std::strtod(argv[10], 0);
        std::cout<<looseContRate << std::endl;
        double dRate = std::strtod(argv[11], 0);
        std::cout<<dRate << std::endl;
        double bRate = std::strtod(argv[12], 0);
        std::cout<<bRate << std::endl;
        double epsilon = std::strtod(argv[13], 0);
        std::cout<<epsilon << std::endl;
        size_t simulationNumber = std::stoi(argv[14]);
        std::cout<<simulationNumber << std::endl;

        executeNSA(tStart, tEnd, nInfected, nSusceptible, nEdges, maxContactsA, MaxContactsB,
                transmRate, newContRate, looseContRate, dRate, bRate, epsilon, simulationNumber);

    }
    /*std::vector<size_t> amount = {100, 500, 1000, 2000, 5000};

    std::vector <double> timesSSA;
    std::vector <double> timesNSA;
    for (size_t i = 0; i < amount.size() - 3; i++)
    {
        size_t nEdges = static_cast<size_t> (amount.at(i) * (amount.at(i) - 1) / 2 * 0.1);
        for (size_t j = 0; j < 100; j++)
        {
            //executeSSA(0, 2, 20, amount.at(i) - 20, nEdges, 1,  20, 0.03, 1.2, 1.25, 0.0004, 2, j);
            //executeNSA(0, 2, 20, amount.at(i) - 20, nEdges, 1,  20, 0.03, 1.2, 1.25, 0.0004, 2, j);
            executeSSAOnlyContactUpdate(0, 2, 20, amount.at(i) - 20, nEdges, 1,  20,  1.2, 1.25, j);
            executeNSAOnlyContactUpdate(0, 2, 20, amount.at(i) - 20, nEdges, 1,  20,  1.2, 1.25,0.03, j);

        }
    }

        std::cout << "------------------" <<std::endl;

*/
    return 0;
}