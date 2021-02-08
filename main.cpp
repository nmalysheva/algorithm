#include <iostream>
#include "ContactNetwork.h"
#include "SSA.h"
#include "NSA.h"
#include <chrono>
#include <unistd.h>
#include <fstream>
#include <string>
#include <vector>
#include "AndersonTauLeap.h"
#include "types.h"

using namespace lemon;
using namespace std;

void executeSSA(double tStart, double tEnd, size_t nInfected, size_t nSusceptible, size_t nEdges,int maxContactsA, int maxContactsB,
        double transmRate, double newContRate, double looseContRate, double diagnRate, double dRate, double bRate, size_t simulationNumber)
{
    ContactNetwork contNetwork(nInfected,
                               nSusceptible,
                               nEdges,
                               maxContactsA,
                               maxContactsB,
                               transmRate,
                               newContRate,
                               looseContRate,
                               diagnRate,
                               dRate,
                               bRate);

    //std::vector<BenStructure> benToFile = contNetwork.getBenStructure(tStart);
    size_t nPopulation = contNetwork.size();
    size_t startEdges = contNetwork.countEdges();
    std::cout <<"Nodes: "  << nPopulation <<std::endl;
    std::cout <<"Edges: "  << startEdges <<std::endl;


    std::vector<double> timeInfection; // vector stores times of infection
    timeInfection.reserve(1e4 + 1);

    std::vector<uint32_t> numberOfTransmitEdges;// vector stores num. of edges eligible for transmission
    numberOfTransmitEdges.reserve(1e+4 + 1);

    NetworkStorage nwStorage;
    nwStorage.reserve(1e6 + 1);

    SSA ssa;
    auto start_time = std::chrono::high_resolution_clock::now();
    ssa.execute(tStart, tEnd, contNetwork, nwStorage,timeInfection,numberOfTransmitEdges,"v"/*, benToFile*/);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    char buffer[50];

    std::ostringstream strs;
    sprintf( buffer, "%.1e", newContRate);
    std::string newContRateStr(buffer);

    sprintf( buffer, "%.1e", looseContRate);
    std::string looseContRateStr(buffer);

    sprintf( buffer, "%.1e", diagnRate);
    std::string diagnRateStr(buffer);

    sprintf( buffer, "%.1e", dRate);
    std::string dRateStr(buffer);

    sprintf( buffer, "%.1e", transmRate);
    std::string transmRateStr(buffer);

    string fileName = "SSA_"  + std::to_string(nPopulation) + "_nInf_" + std::to_string(nInfected) +
                      "_MaxCont_" + std::to_string(maxContactsA) + "-" + std::to_string(maxContactsB) +
                      "_addR_" + newContRateStr + "_delR_" + looseContRateStr +
                      "_diagnR_" + diagnRateStr  + "_deathR_" + dRateStr +
                      "_trR_" + transmRateStr + "_" +
                      std::to_string(simulationNumber) + ".txt";

    ofstream newFile;
    newFile.open(fileName);

    newFile << "{";
    newFile << "\"duration_in_milliseconds\" : " << chrono::duration <double, milli> (time).count() << ", ";
    newFile << "\"population_start\" : " << nPopulation << ",";
    newFile << "\"population_end\": " << contNetwork.size() << ",";
    newFile << "\"edges\": " << startEdges << ",";
    newFile << "\"infected_start\": " << nInfected << ",";
    newFile << "\"infected_end\": " << contNetwork.countByState(Specie::I) + contNetwork.countByState(Specie::D) << ",";


    newFile << "\"max_contact_range\" : [" << maxContactsA << ", "  << maxContactsB << "],";
    newFile << "\"transmission_rate\" : " << transmRate << ",";
    newFile << "\"rate_of_make_a_new_contact\" : " << newContRate << ",";
    newFile << "\"rate_of_loose_a_contact\" : " << looseContRate << ",";
    newFile << "\"death_rate\" : " << dRate << ",";
    newFile << "\"birth_rate\" : " << bRate << ",";
    newFile << "\"diagnosis_rate\" : " << diagnRate << ",";
    newFile << "\"transmission_rate\" : " << transmRate << ",";
    newFile << "\"infect_time\" : [";
    std::string delim = "";
    for (size_t i = 0; i < timeInfection.size(); i++)
    {
        newFile << delim << timeInfection.at(i);
        delim = ", ";
    }
    newFile << "], ";

    /*newFile << "number of edges with inf.: "<< std::endl;
    for (size_t i = 0; i < numberOfTransmitEdges.size(); i++)
    {
        newFile << numberOfTransmitEdges.at(i) << ' ';
    }*/


    newFile << "\"networkStates\" : " <<  "[";


    delim = "";
    for (const auto &item : nwStorage)
    {
        newFile << delim;
        delim = ",";
        newFile << "{";

        newFile << "\"time\" : " << item.first << ", ";
        newFile << "\"nw_states\" : " <<  "[";

        std::string delim2 = "";
        for (const auto &spcs : item.second)
        {
            newFile << delim2;
            delim2 = ", ";
            newFile << "{";
            newFile << "\"id\" : " << spcs.id << ", ";
            newFile << "\"state\" : " << spcs.state << ", ";
            newFile << "\"neighbors\" : " <<  "[";
            std::string delim3 = "";
            for (const auto &nghbr : spcs.contacts)
            {
                newFile  << delim3 << nghbr;
                delim3 = ",";
            }
            newFile << "]}";
        }

        newFile << "]}";
    }
    newFile << "]}";
    newFile << std::endl;


    newFile.close();

}

void executeNSA(double tStart, double tEnd, size_t nInfected, size_t nSusceptible, size_t nEdges,int maxContactsA, int maxContactsB,
                double transmRate, double newContRate, double looseContRate, double diagnRate, double dRate, double bRate, double epsilon, size_t simulationNumber)
{
    ///UniqueID().reset();

    ContactNetwork contNetwork(nInfected,
                               nSusceptible,
                               nEdges,
                               maxContactsA,
                               maxContactsB,
                               transmRate,
                               newContRate,
                               looseContRate,
                               diagnRate,
                               dRate,
                               bRate);

    //std::vector<BenStructure> benToFile = contNetwork.getBenStructure(tStart);
    size_t nPopulation = contNetwork.size();
    size_t startEdges = contNetwork.countEdges();


    std::vector<double> timeInfection;
    timeInfection.reserve(1e4 + 1);

    NetworkStorage nwStorage;
    nwStorage.reserve(1e6 + 1);

    size_t nRejections = 0;
    size_t nAcceptance = 0;
    size_t nThin = 0;


    NSA nsa;
    auto start_time = std::chrono::high_resolution_clock::now();
    nsa.execute(tStart, tEnd, contNetwork, nwStorage, timeInfection,  "v",  epsilon,nRejections, nAcceptance, nThin/*, benToFile*/);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    char buffer[50];

    std::ostringstream strs;
    sprintf( buffer, "%.1e", newContRate);
    std::string newContRateStr(buffer);

    sprintf( buffer, "%.1e", looseContRate);
    std::string looseContRateStr(buffer);

    sprintf( buffer, "%.1e", diagnRate);
    std::string diagnRateStr(buffer);

    sprintf( buffer, "%.1e", dRate);
    std::string dRateStr(buffer);

    sprintf( buffer, "%.1e", transmRate);
    std::string transmRateStr(buffer);

    string fileName = "NSA_"  + std::to_string(nPopulation) + "_nInf_" + std::to_string(nInfected) +
                      "_MaxCont_" + std::to_string(maxContactsA) + "-" + std::to_string(maxContactsB) +
                      "_addR_" + newContRateStr + "_delR_" + looseContRateStr +
                      "_diagnR_" + diagnRateStr + "_deathR_" + dRateStr +
                      "_trR_" + transmRateStr + "_" +
                      std::to_string(simulationNumber) + ".txt";

    ofstream newFile;
    newFile.open(fileName);

    newFile << "{";
    newFile << "\"duration_in_milliseconds\" : " << chrono::duration <double, milli> (time).count() << ", ";
    newFile << "\"population_start\" : " << nPopulation << ",";
    newFile << "\"population_end\": " << contNetwork.size() << ",";
    newFile << "\"edges\": " << startEdges << ",";
    newFile << "\"infected_start\": " << nInfected << ",";
    newFile << "\"infected_end\": " << contNetwork.countByState(Specie::I) + contNetwork.countByState(Specie::D) << ",";

    newFile << "\"max_contact_range\" : [" << maxContactsA << ", "  << maxContactsB << "],";
    newFile << "\"transmission_rate\" : " << transmRate << ",";
    newFile << "\"rate_of_make_a_new_contact\" : " << newContRate << ",";
    newFile << "\"rate_of_loose_a_contact\" : " << looseContRate << ",";
    newFile << "\"death_rate\" : " << dRate << ",";
    newFile << "\"birth_rate\" : " << bRate << ",";
    newFile << "\"diagnosis_rate\" : " << diagnRate << ",";
    newFile << "\"transmission_rate\" : " << transmRate << ",";
    newFile << "\"epsilon\" : " << epsilon <<  ",";
    newFile << "\"rejected\" : " << nRejections << ",";
    newFile << "\"accepted\" : " << nAcceptance << ",";
    newFile << "\"thined\" : "   << nThin << ",";

    newFile << "\"infect_time\" : [";
    std::string delim = "";
    for (size_t i = 0; i < timeInfection.size(); i++)
    {
        newFile << delim << timeInfection.at(i);
        delim = ", ";
    }
    newFile << "], ";

    newFile << "\"networkStates\" : " <<  "[";


    delim = "";
    for (const auto &item : nwStorage)
    {
        newFile << delim;
        delim = ",";
        newFile << "{";

        newFile << "\"time\" : " << item.first << ", ";
        newFile << "\"nw_states\" : " <<  "[";

        std::string delim2 = "";
        for (const auto &spcs : item.second)
        {
            newFile << delim2;
            delim2 = ", ";
            newFile << "{";
            newFile << "\"id\" : " << spcs.id << ", ";
            newFile << "\"state\" : " << spcs.state << ", ";
            newFile << "\"neighbors\" : " <<  "[";
            std::string delim3 = "";
            for (const auto &nghbr : spcs.contacts)
            {
                newFile  << delim3 << nghbr;
                delim3 = ",";
            }
            newFile << "]}";
        }

        newFile << "]}";
    }
    newFile << "]}";
    newFile << std::endl;

    newFile.close();

}

void executeSSAOnlyContactUpdate(double tStart, double tEnd, size_t nPopulation, size_t nEdges,int maxContactsA, int maxContactsB,
                 double newContRate, double looseContRate,  size_t simulationNumber)
{
    ContactNetwork contNetwork(0,
                               nPopulation,
                               nEdges,
                               maxContactsA,
                               maxContactsB,
                               0,
                               newContRate,
                               looseContRate,
                               0,
                               0,
                               0);

    //std::vector<BenStructure> benToFile = contNetwork.getBenStructure(tStart);


    size_t startEdges = contNetwork.countEdges();
    std::cout <<"Nodes: "  << nPopulation <<std::endl;
    std::cout <<"Edges: "  << startEdges <<std::endl;

    std::vector<double> timeInfection;
    timeInfection.reserve(1e4 + 1);

    std::vector<uint32_t> numberOfTransmitEdges;
    numberOfTransmitEdges.reserve(1e+4 + 1);

    NetworkStorage nwStorage;
    nwStorage.reserve(1e6 + 1);

    SSA ssa;
    auto start_time = std::chrono::high_resolution_clock::now();
    ssa.execute(tStart, tEnd, contNetwork, nwStorage, /*timeSteps,*/ timeInfection, /*populationState,*/ numberOfTransmitEdges, /*degreeDistr,*/ "c"/*, benToFile*/);
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


    std::cout << "start save SSA" << std::endl;
    ofstream newFile;
    newFile.open(fileName);


    newFile << "{";
    newFile << "\"duration_in_milliseconds\" : " << chrono::duration <double, milli> (time).count() << ", ";
    newFile << "\"population_start\" : " << nPopulation << ",";
    newFile << "\"population_end\": " << contNetwork.size() << ",";
    newFile << "\"edges\": " << startEdges << ",";

    newFile << "\"max_contact_range\" : [" << maxContactsA << ", "  << maxContactsB << "],";
    newFile << "\"rate_of_make_a_new_contact\" : " << newContRate << ",";
    newFile << "\"rate_of_loose_a_contact\" : " << looseContRate << ",";


    //newFile << "{";
    newFile << "\"networkStates\" : " <<  "[";


    std::string delim = "";
    for (const auto &item : nwStorage)
    {
        newFile << delim;
        delim = ",";
        newFile << "{";

        newFile << "\"time\" : " << item.first << ", ";
        newFile << "\"nw_states\" : " <<  "[";

        std::string delim2 = "";
        for (const auto &spcs : item.second)
        {
            newFile << delim2;
            delim2 = ", ";
            newFile << "{";
            newFile << "\"id\" : " << spcs.id << ", ";
            newFile << "\"state\" : " << spcs.state << ", ";
            newFile << "\"neighbors\" : " <<  "[";
            std::string delim3 = "";
            for (const auto &nghbr : spcs.contacts)
            {
                newFile  << delim3 << nghbr;
                delim3 = ",";
            }
            newFile << "]}";
        }

        newFile << "]}";
    }
    newFile << "]}";
    newFile << std::endl;


    newFile.close();
    std::cout << "end save SSA" << std::endl;
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
    //UniqueID().reset();
    ContactNetwork contNetwork(0,
                               nPopulation,
                               nEdges,
                               maxContactsA,
                               maxContactsB,
                               0,
                               newContRate,
                               looseContRate,
                               0,
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

    NetworkStorage nwStorage;
    //PoissonTauleap(tStart, tEnd, contNetwork, epsilon, timeSteps, degreeDistr, true/*, benToFile*/, generator);
    AndersonTauLeap(tStart, tEnd, contNetwork, epsilon, nwStorage, "c"/*, benToFile*/, generator);

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

    std::cout << "start save NSA" << std::endl;
    ofstream newFile;
    newFile.open(fileName);


    newFile << "{";
    newFile << "\"duration_in_milliseconds\" : " << chrono::duration <double, milli> (time).count() << ", ";
    newFile << "\"population_start\" : " << nPopulation << ", ";
    newFile << "\"population_end\": " << contNetwork.size() << ", ";
    newFile << "\"edges\": " << startEdges << ", ";

    newFile << "\"max_contact_range\" : [" << maxContactsA << ", "  << maxContactsB << "], ";
    newFile << "\"rate_of_make_a_new_contact\" : " << newContRate << ", ";
    newFile << "\"rate_of_loose_a_contact\" : " << looseContRate << ", ";


    newFile << "\"networkStates\" : " <<  "[";


    std::string delim = "";
    for (const auto &item : nwStorage)
    {
        newFile << delim;
        delim = ",";
        newFile << "{";

        newFile << "\"time\" : " << item.first << ", ";
        newFile << "\"nw_states\" : " <<  "[";

        std::string delim2 = "";
        for (const auto &spcs : item.second)
        {
            newFile << delim2;
            delim2 = ", ";
            newFile << "{";
            newFile << "\"id\" : " << spcs.id << ", ";
            newFile << "\"state\" : " << spcs.state << ", ";
            newFile << "\"neighbors\" : " <<  "[";
            std::string delim3 = "";
            for (const auto &nghbr : spcs.contacts)
            {
                newFile  << delim3 << nghbr;
                delim3 = ",";
            }
            newFile << "]}";
        }

        newFile << "]}";
    }
    newFile << "]}";
    newFile << std::endl;

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

    maxContactsA = nPopulation - 1;
    maxContactsB = nPopulation - 1;
    std::cout << "Max. contact boundaries =[" << maxContactsA << ", " << maxContactsB << "]" << std::endl;

    newContRate = std::strtod(arr[2], 0);
    std::cout << "New Contact Rate: " << newContRate << std::endl;

    looseContRate = std::strtod(arr[3], 0);
    std::cout << "Loose Contact Rate: " << looseContRate << std::endl;


}

void contactDynamics(int argc, char* argv[])
{
    if (argc != 9)
    {
        std::string msg = "Invalid parameters";
        throw std::domain_error(msg);
    }
    else
    {
        size_t simulationNumber = 0;
        char* arr[] = {argv[2], argv[3], argv[8]};
        std::pair<double, double> simulationTime = getAndPrintSimulationParameters(arr, simulationNumber);

        size_t nPopulation = 0;
        size_t nEdges = 0;
        int maxContactsA = 0;
        int maxContactsB = 0;
        double newContRate = 0;
        double looseContRate = 0;

        char* arr1[]  = {argv[4], argv[5], argv[6], argv[7]};
        getAndPrintNetworkParameters(nPopulation, nEdges, maxContactsA, maxContactsB,
                                      newContRate, looseContRate, arr1);
        double epsilon = 0.03;

        executeNSAOnlyContactUpdate(simulationTime.first, simulationTime.second, nPopulation, nEdges, maxContactsA, maxContactsB,
                         newContRate, looseContRate, epsilon, simulationNumber);
        executeSSAOnlyContactUpdate(simulationTime.first, simulationTime.second, nPopulation, nEdges, maxContactsA, maxContactsB,
                                   newContRate, looseContRate, simulationNumber);
        std::cout << "*********************************" << std::endl;
        //executeRKF45OnlyContactUpdate(tStart, tEnd, nPopulation, nEdges, maxContactsA, MaxContactsB,
        //                              newContRate, looseContRate, epsilon, simulationNumber);


    }

}


void viralDynamics(int argc, char* argv[])
{
    if (argc != 13)
    {
        std::string msg = "Invalid parameters 2";
        throw std::domain_error(msg);
    }
    else
    {

        char* arr[] = {argv[2], argv[3], argv[12]};
        size_t simulationNumber = 0;
        std::pair<double, double> simulationTime = getAndPrintSimulationParameters(arr, simulationNumber);

        size_t nPopulation = 0;
        size_t nEdges = 0;
        int maxContactsA = 0;
        int maxContactsB = 0;
        double newContRate = 0;
        double looseContRate = 0;

        char* arr1[]  = {argv[4], argv[6], argv[7], argv[8],};
        getAndPrintNetworkParameters(nPopulation, nEdges, maxContactsA, maxContactsB,
                                     newContRate, looseContRate, arr1);


        size_t nInfected = std::stoi(argv[5]);
        std::cout << "Number of infected: " << nInfected << std::endl;


        double diagnosisRate = std::strtod(argv[9], 0);
        std::cout << "Diagnosis rate: " << diagnosisRate << std::endl;

        double deathRate = std::strtod(argv[10], 0);
        std::cout << "Death rate: " << deathRate << std::endl;

        double transmitRate = std::strtod(argv[11], 0);
        std::cout << "Transmission rate: " << transmitRate << std::endl;

        double birthRate = 0;


        double epsilon = 0.03;


        executeNSA(simulationTime.first, simulationTime.second, nInfected, nPopulation - nInfected, nEdges, maxContactsA, maxContactsB,
                   transmitRate, newContRate, looseContRate, diagnosisRate, deathRate, birthRate, epsilon, simulationNumber);
        //executeSSA(simulationTime.first, simulationTime.second, nInfected, nPopulation - nInfected, nEdges, maxContactsA, maxContactsB,
          //         transmitRate, newContRate, looseContRate, diagnosisRate, deathRate, birthRate, simulationNumber);

    }

}

int main(int argc, char* argv[])
{
    if (std::string(argv[1]) == "-c")
    {
        /*parameters:
         * tStart
         * tEnd
         * nPopulation
         * nEdges_at_the_Begining
         * add_cont_rate
         * loose_cont_rate
         * simulation_Number
        */
        contactDynamics(argc, argv);
    }
    else if (std::string(argv[1])  == "-v")
    {
        /*parameters:
         * tStart
         * tEnd
         * nPopulation
         * nInfected
         * nEdges_at_the_Begining
         * add_cont_rate
         * loose_cont_rate
         * diagnosis_rate
         * death_rate
         * transmission_rate
         * simulation_Number
        */
        viralDynamics(argc, argv);
    }
    else
    {
        std::string msg = "Invalid parameters 1";
        throw std::domain_error(msg);
    }
    return 0;
}
