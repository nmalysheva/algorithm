#include <iostream>
#include "UniqueID.h"
#include "ContactNetwork.h"
#include "SSA.h"
#include "NSA.h"
#include <chrono>
#include <unistd.h>
#include <fstream>
#include <string>
#include <vector>
#include "PoissonTauLeap.h"
#include "AndersonTauLeap.h"
#include "RKF45.h"

using namespace lemon;
using namespace std;

void executeSSA(double tStart, double tEnd, size_t nInfected, size_t nSusceptible, size_t nEdges,int maxContactsA, int maxContactsB,
        double transmRate, double newContRate, double looseContRate, double dRate, double bRate, size_t simulationNumber)
{
    UniqueID().reset();
    ContactNetwork contNetwork(nInfected,
                               nSusceptible,
                               nEdges,
                               maxContactsA,
                               maxContactsB,
                               transmRate,
                               newContRate,
                               looseContRate,
                               dRate,
                               bRate);

    //std::vector<BenStructure> benToFile = contNetwork.getBenStructure(tStart);
    size_t nPopulation = contNetwork.size();
    size_t startEdges = contNetwork.countEdges();
    std::cout <<"Nodes: "  << nPopulation <<std::endl;
    std::cout <<"Edges: "  << startEdges <<std::endl;

    std::vector<double> timeSteps;
    std::vector<uint32_t> infectedSteps;
    std::vector<std::vector<size_t>> degreeDistr; // vector stores degree dist. per time step

    SSA ssa;
    auto start_time = std::chrono::high_resolution_clock::now();
    ssa.execute(tStart, tEnd, contNetwork, timeSteps, infectedSteps, degreeDistr/*, benToFile*/);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    char buffer[50];

    std::ostringstream strs;
    sprintf( buffer, "%.1e", newContRate);
    std::string newContRateStr(buffer);

    sprintf( buffer, "%.1e", looseContRate);
    std::string looseContRateStr(buffer);

    sprintf( buffer, "%.1e", bRate);
    std::string bRateStr(buffer);

    sprintf( buffer, "%.1e", dRate);
    std::string dRateStr(buffer);

    sprintf( buffer, "%.1e", transmRate);
    std::string transmRateStr(buffer);

    string fileName = "SSA_"  + std::to_string(nPopulation) + "_nInf_" + std::to_string(nInfected) +
                      "_MaxCont_" + std::to_string(maxContactsA) + "-" + std::to_string(maxContactsB) +
                      "_addR_" + newContRateStr + "_delR_" + looseContRateStr +
                      "_birthR_" + bRateStr + "_deathR_" + dRateStr +
                      "_trR_" + transmRateStr + "_" +
                      std::to_string(simulationNumber) + ".txt";

    ofstream newFile;
    newFile.open(fileName);

    //newFile << "duration in CPU time: " << time;
    newFile << "duration in milliseconds: " << chrono::duration <double, milli> (time).count() << std::endl;
    newFile << "nodes: " << nPopulation << std::endl;
    newFile << "edges: " << startEdges << std::endl;
    newFile << "infected: " << nInfected << std::endl;

    newFile << "max contact range: " << maxContactsA << " "  << maxContactsB << std::endl;
    newFile << "transmission rate: " << transmRate << std::endl;
    newFile << "rate of make a new contact: " << newContRate << std::endl;
    newFile << "rate of loose a contact: " << looseContRate << std::endl;
    newFile << "death rate: " << dRate << std::endl;
    newFile << "birth rate: " << bRate << std::endl;
    newFile << "transmission rate: " << transmRate << std::endl;

    for (const auto &tStep : timeSteps)
    {
        newFile << tStep << ' ';
    }
    newFile << std::endl;

    for (const auto &infectedStep : infectedSteps)
    {
        newFile << infectedStep << ' ';
    }
    newFile << std::endl;
    newFile << "degree distribution:" << std::endl;
    for (size_t i = 0; i < degreeDistr.size(); i++)
    {
        std::vector<size_t> degreeDatTime = degreeDistr.at(i);
        for (size_t j = 0;  j< degreeDatTime.size(); j++)
        {
            newFile << degreeDatTime.at(j) << ' ';
        }
        newFile << std::endl;
    }

    newFile.close();

}

void executeNSA(double tStart, double tEnd, size_t nInfected, size_t nSusceptible, size_t nEdges,int maxContactsA, int maxContactsB,
                double transmRate, double newContRate, double looseContRate, double dRate, double bRate, double epsilon, size_t simulationNumber)
{
    UniqueID().reset();
    ContactNetwork contNetwork(nInfected,
                               nSusceptible,
                               nEdges,
                               maxContactsA,
                               maxContactsB,
                               transmRate,
                               newContRate,
                               looseContRate,
                               dRate,
                               bRate);

    //std::vector<BenStructure> benToFile = contNetwork.getBenStructure(tStart);
    size_t nPopulation = contNetwork.size();
    size_t startEdges = contNetwork.countEdges();
    //std::cout <<"Nodes: "  << nPopulation <<std::endl;
    //std::cout <<"Edges: "  << startEdges <<std::endl;

    std::vector<double> timeSteps;
    std::vector<uint32_t> infectedSteps;
    std::vector<std::vector<size_t>> degreeDistr; // vector stores degree dist. per time step

    size_t nRejections = 0;
    size_t nAcceptance = 0;
    size_t nThin = 0;


    NSA nsa;
    //ContactNetwork cNw2(20, amount.at(i) - 20, nEdges, 1, 20, 0.03, 1.2, 1.25, 0.0004, 2);
    auto start_time = std::chrono::high_resolution_clock::now();
    nsa.execute(tStart, tEnd, contNetwork, timeSteps, infectedSteps, degreeDistr, epsilon,nRejections, nAcceptance, nThin/*, benToFile*/);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    char buffer[50];

    std::ostringstream strs;
    sprintf( buffer, "%.1e", newContRate);
    std::string newContRateStr(buffer);

    sprintf( buffer, "%.1e", looseContRate);
    std::string looseContRateStr(buffer);

    sprintf( buffer, "%.1e", bRate);
    std::string bRateStr(buffer);

    sprintf( buffer, "%.1e", dRate);
    std::string dRateStr(buffer);

    sprintf( buffer, "%.1e", transmRate);
    std::string transmRateStr(buffer);

    string fileName = "NSA_"  + std::to_string(nPopulation) + "_nInf_" + std::to_string(nInfected) +
                      "_MaxCont_" + std::to_string(maxContactsA) + "-" + std::to_string(maxContactsB) +
                      "_addR_" + newContRateStr + "_delR_" + looseContRateStr +
                      "_birthR_" + bRateStr + "_deathR_" + dRateStr +
                      "_trR_" + transmRateStr + "_" +
                      std::to_string(simulationNumber) + ".txt";

    ofstream newFile;
    newFile.open(fileName);

    newFile << "duration in milliseconds: " << chrono::duration <double, milli> (time).count() << std::endl;
    newFile << "nodes: " << nPopulation << std::endl;
    newFile << "edges: " << startEdges << std::endl;
    newFile << "infected: " << nInfected << std::endl;

    newFile << "max contact range: " << maxContactsA << " "  << maxContactsB << std::endl;
    newFile << "transmission rate: " << transmRate << std::endl;
    newFile << "rate of make a new contact: " << newContRate << std::endl;
    newFile << "rate of loose a contact: " << looseContRate << std::endl;
    newFile << "death rate: " << dRate << std::endl;
    newFile << "birth rate: " << bRate << std::endl;
    newFile << "transmission rate: " << transmRate << std::endl;
    newFile << "epsilon: " << epsilon << std::endl;

    for (const auto &tStep: timeSteps)
    {
        newFile << tStep << ' ';
    }
    newFile << std::endl;

    for (size_t i = 0; i < infectedSteps.size(); i++)
    {
        newFile << infectedSteps.at(i) << ' ';
    }
    newFile <<  std::endl;
    newFile << "rejected: " << nRejections << std::endl;
    newFile << "accepted: " << nAcceptance << std::endl;
    newFile << "thined: "   << nThin << std::endl;
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

/*void executeRKF(double tStart, double tEnd, size_t nInfected, size_t nSusceptible, size_t nEdges,int maxContactsA, int maxContactsB,
                double transmRate, double newContRate, double looseContRate, double dRate, double bRate, double epsilon, size_t simulationNumber)
{
    UniqueID().reset();
    ContactNetwork contNetwork(nInfected,
                               nSusceptible,
                               nEdges,
                               maxContactsA,
                               maxContactsB,
                               transmRate,
                               newContRate,
                               looseContRate,
                               dRate,
                               bRate);

    size_t nPopulation = contNetwork.size();
    size_t startEdges = contNetwork.countEdges();
    //std::cout <<"Nodes: "  << nPopulation <<std::endl;
    //std::cout <<"Edges: "  << startEdges <<std::endl;

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

    char buffer[50];

    std::ostringstream strs;
    sprintf( buffer, "%.1e", newContRate);
    std::string newContRateStr(buffer);

    sprintf( buffer, "%.1e", looseContRate);
    std::string looseContRateStr(buffer);

    sprintf( buffer, "%.1e", bRate);
    std::string bRateStr(buffer);

    sprintf( buffer, "%.1e", dRate);
    std::string dRateStr(buffer);

    sprintf( buffer, "%.1e", transmRate);
    std::string transmRateStr(buffer);

    string fileName = "NSA_"  + std::to_string(nPopulation) + "_nInf_" + std::to_string(nInfected) +
                      "_MaxCont_" + std::to_string(maxContactsA) + "-" + std::to_string(maxContactsB) +
                      "_addR_" + newContRateStr + "_delR_" + looseContRateStr +
                      "_birthR_" + bRateStr + "_deathR_" + dRateStr +
                      "_trR_" + transmRateStr + "_" +
                      std::to_string(simulationNumber) + ".txt";

    ofstream newFile;
    newFile.open(fileName);

    newFile << "duration in milliseconds: " << chrono::duration <double, milli> (time).count() << std::endl;
    newFile << "nodes: " << nPopulation << std::endl;
    newFile << "edges: " << startEdges << std::endl;
    newFile << "infected: " << nInfected << std::endl;

    newFile << "max contact range: " << maxContactsA << " "  << maxContactsB << std::endl;
    newFile << "transmission rate: " << transmRate << std::endl;
    newFile << "rate of make a new contact: " << newContRate << std::endl;
    newFile << "rate of loose a contact: " << looseContRate << std::endl;
    newFile << "death rate: " << dRate << std::endl;
    newFile << "birth rate: " << bRate << std::endl;
    newFile << "transmission rate: " << transmRate << std::endl;
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
    newFile << "thined: "   << nThin << std::endl;
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
}*/

void executeSSAOnlyContactUpdate(double tStart, double tEnd, size_t nPopulation, size_t nEdges,int maxContactsA, int maxContactsB,
                 double newContRate, double looseContRate,  size_t simulationNumber)
{
    UniqueID().reset();
    ContactNetwork contNetwork(0,
                               nPopulation,
                               nEdges,
                               maxContactsA,
                               maxContactsB,
                               0,
                               newContRate,
                               looseContRate,
                               0,
                               0);

    //std::vector<BenStructure> benToFile = contNetwork.getBenStructure(tStart);


    size_t startEdges = contNetwork.countEdges();
    std::cout <<"Nodes: "  << nPopulation <<std::endl;
    std::cout <<"Edges: "  << startEdges <<std::endl;

    std::vector<double> timeSteps;
    timeSteps.reserve(1e4 + 1);
    std::vector<uint32_t> infectedSteps;
    infectedSteps.reserve(1e4 + 1);
    std::vector<std::vector<size_t>> degreeDistr; // vector stores degree dist. per time step
    degreeDistr.reserve(1e4 + 1);
    SSA ssa;
    auto start_time = std::chrono::high_resolution_clock::now();
    ssa.execute(tStart, tEnd, contNetwork, timeSteps, infectedSteps, degreeDistr/*, benToFile*/);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    char buffer[50];

    std::ostringstream strs;
    sprintf( buffer, "%.1e", newContRate);
    std::string newContRateStr(buffer);

    sprintf( buffer, "%.1e", looseContRate);
    std::string looseContRateStr(buffer);

    string fileName = "ContDyn_SSA_" + std::to_string(nPopulation) + "_MaxCont_" +
                      std::to_string(maxContactsA) + "-" + std::to_string(maxContactsB) + "_addR_" +
                      newContRateStr +"_delR_" + looseContRateStr + "_" +
                      std::to_string(simulationNumber) + ".txt";


    ofstream newFile;
    newFile.open(fileName);

    //newFile << "duration in CPU time: " << time;
    newFile << "duration in milliseconds: " << chrono::duration <double, milli> (time).count() << std::endl;
    newFile << "nodes: " << nPopulation << std::endl;
    newFile << "edges: " << startEdges << std::endl;
    newFile << "population: " << nPopulation << std::endl;

    newFile << "max contact range: " << maxContactsA << " "  << maxContactsB << std::endl;
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
    //-----------ben format
    /*string fileNameBen = "Ben_ContDyn_SSA_" + std::to_string(nPopulation) + "_MaxCont_" +
                         std::to_string(maxContactsA) + "-" + std::to_string(maxContactsB) + "_addR_" +
                         newContRateStr +"_delR_" + looseContRateStr + "_" +
                         std::to_string(simulationNumber) + ".txt";
    ofstream benFile;
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
    benFile.close();
*/
}


void executeNSAOnlyContactUpdate(double tStart, double tEnd, size_t nPopulation, size_t nEdges,int maxContactsA, int maxContactsB,
                double newContRate, double looseContRate, double epsilon, size_t simulationNumber)
{
    UniqueID().reset();
    ContactNetwork contNetwork(0,
                               nPopulation,
                               nEdges,
                               maxContactsA,
                               maxContactsB,
                               0,
                               newContRate,
                               looseContRate,
                               0,
                               0);
    //size_t nPopulation = contNetwork.size();

    //std::vector<BenStructure> benToFile = contNetwork.getBenStructure(tStart);
    size_t startEdges = contNetwork.countEdges();
    std::cout <<"Nodes: "  << nPopulation <<std::endl;
    std::cout <<"Edges: "  << startEdges <<std::endl;

    std::vector<double> timeSteps;
    timeSteps.reserve(1e4 + 1);
    std::vector<std::vector<size_t>> degreeDistr; // vector stores degree dist. per time step
    degreeDistr.reserve(1e4 + 1);
    NSA nsa;
    auto start_time = std::chrono::high_resolution_clock::now();

    std::random_device rDev;
    std::mt19937_64 generator = std::mt19937_64(rDev());
    generator.seed(::time(nullptr) * getpid()); //to change the seed for every run

    //PoissonTauleap(tStart, tEnd, contNetwork, epsilon, timeSteps, degreeDistr, true/*, benToFile*/, generator);
    AndersonTauLeap(tStart, tEnd, contNetwork, epsilon, timeSteps, degreeDistr, true/*, benToFile*/, generator);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    char buffer[50];

    std::ostringstream strs;
    sprintf( buffer, "%.1e", newContRate);
    std::string newContRateStr(buffer);

    sprintf( buffer, "%.1e", looseContRate);
    std::string looseContRateStr(buffer);

    string fileName = "ContDyn_NSA_" + std::to_string(nPopulation) + "_MaxCont_" +
                      std::to_string(maxContactsA) + "-" + std::to_string(maxContactsB) + "_addR_" +
                      newContRateStr +"_delR_" + looseContRateStr + "_" +
                      std::to_string(simulationNumber) + ".txt";


    ofstream newFile;
    newFile.open(fileName);

    newFile << "duration in milliseconds: " << chrono::duration <double, milli> (time).count() << std::endl;
    newFile << "nodes: " << nPopulation << std::endl;
    newFile << "edges: " << startEdges << std::endl;
    newFile << "population: " << nPopulation << std::endl;

    newFile << "max contact range: " << "[" << maxContactsA << ", "  << maxContactsB << "]" << std::endl;
    newFile << "rate of make a new contact: " << newContRate << std::endl;
    newFile << "rate of loose a contact: " << looseContRate << std::endl;
    newFile << "epsilon: " << epsilon << std::endl;

    std::cout << " timeSteps.size: " <<  timeSteps.size() <<std::endl;
    for (const auto &tStep: timeSteps)
    {
        newFile << tStep << ' ';
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
    std::cout<<"save"<< std::endl;

    //-----------ben format
   /* string fileNameBen = "Ben_ContDyn_NSA_" + std::to_string(nPopulation) + "_MaxCont_" +
                         std::to_string(maxContactsA) + "-" + std::to_string(maxContactsB) + "_addR_" +
                         newContRateStr +"_delR_" + looseContRateStr + "_" +
                         std::to_string(simulationNumber) + ".txt";
    ofstream benFile;
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
    benFile.close();
*/
}

void executeRKF45OnlyContactUpdate(double tStart, double tEnd, size_t nPopulation, size_t nEdges,int maxContactsA, int maxContactsB,
                                 double newContRate, double looseContRate, double epsilon, size_t simulationNumber)
{
    UniqueID().reset();
    ContactNetwork contNetwork(0,
                               nPopulation,
                               nEdges,
                               maxContactsA,
                               maxContactsB,
                               0,
                               newContRate,
                               looseContRate,
                               0,
                               0);

    //std::vector<BenStructure> benToFile = contNetwork.getBenStructure(tStart);

    size_t startEdges = contNetwork.countEdges();
    std::cout <<"Nodes: "  << nPopulation <<std::endl;
    std::cout <<"Edges: "  << startEdges <<std::endl;

    std::vector<double> timeSteps;
    std::vector<std::vector<size_t>> degreeDistr; // vector stores degree dist. per time step
    NSA nsa;
    //ContactNetwork cNw2(20, amount.at(i) - 20, nEdges, 1, 20, 0.03, 1.2, 1.25, 0.0004, 2);

    std::random_device rDev;
    std::mt19937_64 generator = std::mt19937_64(rDev());
    generator.seed(::time(NULL) * getpid());

    double dtMax = (tEnd - tStart) / 2;
    double dtMin = (tEnd - tStart) * 1e-4;
    double errorMax = 1e-1;
    double errorMin = 1e-3;
    auto start_time = std::chrono::high_resolution_clock::now();
    RKF45Approximation(tStart, tEnd, contNetwork, dtMax, dtMin, errorMax, errorMin, timeSteps,
            degreeDistr, true/*, benToFile*/, generator);
    /*nsa.MidpointApproximation(tStart, tEnd, contNetwork, dtMax, dtMin, errorMax, errorMin, timeSteps,
                           degreeDistr);*/
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    char buffer[50];

    std::ostringstream strs;
    sprintf( buffer, "%.1e", newContRate);
    std::string newContRateStr(buffer);

    sprintf( buffer, "%.1e", looseContRate);
    std::string looseContRateStr(buffer);

    string fileName = "ContDyn_RKF45_" + std::to_string(nPopulation) + "_MaxCont_" +
                      std::to_string(maxContactsA) + "-" + std::to_string(maxContactsB) + "_addR_" +
                      newContRateStr +"_delR_" + looseContRateStr + "_" +
                      std::to_string(simulationNumber) + ".txt";
    /*string fileName = "ContDyn_Midpoint_" + std::to_string(nPopulation) + "_MaxCont_" +
                      std::to_string(maxContactsA) + "-" + std::to_string(maxContactsB) + "_addR_" +
                      newContRateStr +"_delR_" + looseContRateStr + "_" +
                      std::to_string(simulationNumber) + ".txt";*/


    ofstream newFile;
    newFile.open(fileName);

    newFile << "duration in milliseconds: " << chrono::duration <double, milli> (time).count() << std::endl;
    newFile << "nodes: " << nPopulation << std::endl;
    newFile << "edges: " << startEdges << std::endl;
    newFile << "population: " << nPopulation << std::endl;

    newFile << "max contact range: " << maxContactsA << " "  << maxContactsB << std::endl;
    newFile << "rate of make a new contact: " << newContRate << std::endl;
    newFile << "rate of loose a contact: " << looseContRate << std::endl;
    newFile << "epsilon: " << epsilon << std::endl;

    std::cout << " timeSteps.size: " <<  timeSteps.size() <<std::endl;
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
    std::cout<<"save"<< std::endl;


    //-----------ben format
    /*string fileNameBen = "Ben_ContDyn_RKF45_" + std::to_string(nPopulation) + "_MaxCont_" +
                      std::to_string(maxContactsA) + "-" + std::to_string(maxContactsB) + "_addR_" +
                      newContRateStr +"_delR_" + looseContRateStr + "_" +
                      std::to_string(simulationNumber) + ".txt";
    ofstream benFile;
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
    benFile.close();
*/
}

std::pair<double, double> getAndPrintSimulationParameters(char* timeBoundaries[], size_t &simulationNumber)
{
    simulationNumber = std::stoi(timeBoundaries[2]);
    std::cout << "Simulation number: " << simulationNumber << std::endl;

    double tStart = std::strtod(timeBoundaries[0], 0);
    double tEnd   = std::strtod(timeBoundaries[1], 0);
    std::cout << "Simulation time t=[" << tStart << ", " << tEnd << "]" << std::endl;
    return std::make_pair(tStart, tEnd);


}

void getAndPrintNetworkParameters(size_t &nPopulation, size_t &nEdges, int &maxContactsA, int &maxContactsB,
                                  double &newContRate, double &looseContRate,
                                  char* arr[])
{
    nPopulation = std::stoi(arr[0]);
    std::cout << "Population size: " << nPopulation << std::endl;

    nEdges = std::stoi(arr[1]);
    std::cout << "Initial edges: " << nEdges << std::endl;

    maxContactsA = std::stoi(arr[2]);
    maxContactsB = std::stoi(arr[3]);
    std::cout << "Max. contact boundaries =[" << maxContactsA << ", " << maxContactsB << "]" << std::endl;

    newContRate = std::strtod(arr[4], 0);
    std::cout << "New Contact Rate: " << newContRate << std::endl;

    looseContRate = std::strtod(arr[5], 0);
    std::cout << "Loose Contact Rate: " << looseContRate << std::endl;


}

void contactDynamics(int argc, char* argv[])
{
    if (argc != 10)
    {
        std::cout << "wrong parameters" <<std::endl;
    }
    else
    {
        size_t simulationNumber = 0;
        char* arr[] = {argv[1], argv[2], argv[9]};
        std::pair<double, double> simulationTime = getAndPrintSimulationParameters(arr, simulationNumber);

        size_t nPopulation = 0;
        size_t nEdges = 0;
        int maxContactsA = 0;
        int maxContactsB = 0;
        double newContRate = 0;
        double looseContRate = 0;

        char* arr1[]  = {argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]};
        getAndPrintNetworkParameters(nPopulation, nEdges, maxContactsA, maxContactsB,
                                      newContRate, looseContRate, arr1);
        double epsilon = 0.03;

        executeNSAOnlyContactUpdate(simulationTime.first, simulationTime.second, nPopulation, nEdges, maxContactsA, maxContactsB,
                         newContRate, looseContRate, epsilon, simulationNumber);
        executeSSAOnlyContactUpdate(simulationTime.first, simulationTime.second, nPopulation, nEdges, maxContactsA, maxContactsB,
                                    newContRate, looseContRate, simulationNumber);
        //executeRKF45OnlyContactUpdate(tStart, tEnd, nPopulation, nEdges, maxContactsA, MaxContactsB,
        //                              newContRate, looseContRate, epsilon, simulationNumber);


    }

}


void viralDynamics(int argc, char* argv[])
{
    if (argc != 14)
    {
        std::cout << "wrong parameters" <<std::endl;
    }
    else
    {

        char* arr[] = {argv[1], argv[2], argv[13]};
        size_t simulationNumber = 0;
        std::pair<double, double> simulationTime = getAndPrintSimulationParameters(arr, simulationNumber);

        size_t nPopulation = 0;
        size_t nEdges = 0;
        int maxContactsA = 0;
        int maxContactsB = 0;
        double newContRate = 0;
        double looseContRate = 0;

        char* arr1[]  = {argv[3], argv[5], argv[6], argv[7], argv[8], argv[9], argv[9]};
        getAndPrintNetworkParameters(nPopulation, nEdges, maxContactsA, maxContactsB,
                                     newContRate, looseContRate, arr1);


        size_t nInfected = std::stoi(argv[4]);
        std::cout << "Number of infected: " << nInfected << std::endl;


        double birthRate = std::strtod(argv[10], 0);
        std::cout << "Birth rate: " << birthRate << std::endl;

        double deathRate = std::strtod(argv[11], 0);
        std::cout << "Death rate: " << deathRate << std::endl;

        double transmitRate = std::strtod(argv[12], 0);
        std::cout << "Transmission rate: " << transmitRate << std::endl;


        double epsilon = 0.03;

        executeNSA(simulationTime.first, simulationTime.second, nInfected, nPopulation - nInfected, nEdges, maxContactsA, maxContactsB,
                   transmitRate, newContRate, looseContRate, deathRate, birthRate, epsilon, simulationNumber);
        executeSSA(simulationTime.first, simulationTime.second, nInfected, nPopulation - nInfected, nEdges, maxContactsA, maxContactsB,
                   transmitRate, newContRate, looseContRate, deathRate, birthRate, simulationNumber);

    }

}

int main(int argc, char* argv[])
{
    contactDynamics(argc, argv);
    //viralDynamics(argc, argv);

    /*std::random_device rDev;
    std::mt19937_64 generator = std::mt19937_64(rDev());
    std::uniform_int_distribution<> distrib(0, 200);
    std::uniform_int_distribution<> distrib2(0, 250);
    int N = 50;
    for (int i = 0; i < 1; i++)
    {
        int nAdd = 10;//distrib(generator);
        int nDel = 50;//distrib2(generator);

        std::cout << "nAdd = " << nAdd << ", nDel = " << nDel<<  ",  " <<nAdd - nDel << ", ";
        int a = splitRandomNumber( nDel,  nAdd,  N,  generator);
        std::cout << a << std::endl;

    }*/
    return 0;
}
