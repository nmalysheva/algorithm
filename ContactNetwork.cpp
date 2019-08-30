//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include "ContactNetwork.h"
#include "UniqueID.h"
#include <random>
#include <lemon/full_graph.h>
#include <map>


size_t  ContactNetwork::countByState(Specie::State st) const
{
    size_t result = 0;

    for (auto& it: population)
    {
        if (it.second.getState() == st)
        {
            result++;
        }
    }
    return result;
}


double  ContactNetwork::getTransmissionRateSum()const
{
    double result = 0;
    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        result += transmissionRates[eIt];
    }
    return result;
}

double  ContactNetwork::getEdgeDeletionRateSum()const
{
    double result = 0;
    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        result +=getEdgeDelitionRate(eIt);
    }
    return result;
}

double  ContactNetwork::getEdgeAdditionRateSum()const
{
    double result = 0;
    for (lemon::ListGraph::EdgeIt eIt(complement); eIt != lemon::INVALID; ++eIt)
    {
        result += getEdgeAdditionRate(eIt);

    }
    return result;
}

double  ContactNetwork::getDeathRateSum()const
{
    double result = 0;
    for (auto& it: population)
    {
        result += it.second.getDeathRate();

    }

    return result;
}

double ContactNetwork::getTransmissionRateLimit() const
{
    return transmissionRate;
}

size_t ContactNetwork::getMaxContactsLimitOfInfected()const
{
    size_t result = 0;

    for (auto& it: population)
    {
        if (it.second.getState() == Specie::State::I)
        {
            result += it.second.getMaxNumberOfContacts();
        }
    }
    return result;
}

double ContactNetwork::getBirthRateSum()const
{
    return birthRate;
}

size_t ContactNetwork::size() const
{
    return lemon::countNodes(network);
}

void ContactNetwork::addEdge(lemon::ListGraph::Edge &complementEdge, double trRate)
{
    lemon::ListGraph::Node u = complement.u(complementEdge);
    lemon::ListGraph::Node v = complement.v(complementEdge);

    // find nodes in network
    lemon::ListGraph::Node cu = network.nodeFromId(complement.id(u));
    lemon::ListGraph::Node cv = network.nodeFromId(complement.id(v));

    complement.erase(complementEdge);

    lemon::ListGraph::Edge e = network.addEdge(cu, cv);
    transmissionRates[e] = trRate;

    uint32_t sourceUID = nodeUIDs[cu];
    uint32_t targetUID = nodeUIDs[cv];
    population.at(sourceUID).incNumberOfContacts();
    population.at(targetUID).incNumberOfContacts();
}

void ContactNetwork::removeEdge(lemon::ListGraph::Edge &edge)
{
    //lemon::ListGraph::Edge e = network.edgeFromId(edgeID);

    lemon::ListGraph::Node u = network.u(edge);
    lemon::ListGraph::Node v = network.v(edge);

    int sourceID = network.id(u);
    int targetID = network.id(v);

    network.erase(edge); // erase from network

    u = complement.nodeFromId(sourceID);
    v = complement.nodeFromId(targetID);
    complement.addEdge(u, v);   // add to complement graph

    uint32_t sourceUID = nodeUIDs[u];
    uint32_t targetUID = nodeUIDs[v];
    population.at(sourceUID).decNumberOfContacts();
    population.at(targetUID).decNumberOfContacts();
}

void ContactNetwork::removeNode(lemon::ListGraph::Node &node)
{
    //std::cout <<nodeID <<std::endl;
    //std::cout <<network.valid(node) <<std::endl;

    lemon::ListGraph::IncEdgeIt ieIt(network, node);
    while (ieIt != lemon::INVALID)
    {
        lemon::ListGraph::IncEdgeIt tmpIt = ieIt;
        ++ieIt;
        network.erase(tmpIt);
    }


    uint32_t sourceUID = nodeUIDs[node];

    lemon::ListGraph::Node complementNode = complement.nodeFromId(network.id(node));

    lemon::ListGraph::IncEdgeIt cieIt(complement, complementNode);
    while (cieIt != lemon::INVALID)
    {
        lemon::ListGraph::IncEdgeIt tmpIt = cieIt;
        ++cieIt;
        complement.erase(tmpIt);
    }

    network.erase(node);
    complement.erase(complementNode);
    population.erase(sourceUID);
}

void ContactNetwork::executeEdgeDelition(double rStart, double rBound)
{
    double result = rStart;
    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        result += getEdgeDelitionRate(eIt);
        //result += disassembleRates[eIt];
        if (result >= rBound)
        {
            removeEdge(eIt);
            break;
        }
    }
}
void ContactNetwork::executeEdgeAddition(double rStart, double rBound)
{
    double result = rStart;
    //std::cout << result << std::endl;
    for (lemon::ListGraph::EdgeIt eIt(complement); eIt != lemon::INVALID; ++eIt)
    {

        result += getEdgeAdditionRate(eIt);
        //std::cout << result << std::endl;
        //result += assembleRates[eIt];
        if (result >= rBound)
        {
            double tRate = 0;
            lemon::ListGraph::Node u = complement.u(eIt);
            lemon::ListGraph::Node v = complement.v(eIt);

            // find nodes in network
            lemon::ListGraph::Node cu = network.nodeFromId(complement.id(u));
            lemon::ListGraph::Node cv = network.nodeFromId(complement.id(v));

            uint32_t sourceUID = nodeUIDs[cu];
            uint32_t targetUID = nodeUIDs[cv];

            if ((population.at(sourceUID).getState()  == Specie::S  && population.at(targetUID).getState() == Specie::I) ||
                (population.at(sourceUID).getState()  == Specie::I  && population.at(targetUID).getState() == Specie::S))
            {

                //tRate = transmitDistribution(generator);
                tRate = transmissionRate;
            }

            //std::cout<<"trate " << tRate << std::endl;
            addEdge(eIt, tRate);
            break;
        }

    }
}


void ContactNetwork::executeTransmission(double rStart, double rBound, double time)
{
    double result = rStart;

    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        result += transmissionRates[eIt];
        //std::cout << transmissionRates[eIt] << std::endl;
        if (result >= rBound)
        {
            lemon::ListGraph::Node u = network.u(eIt);
            lemon::ListGraph::Node v = network.v(eIt);
            uint32_t sourceUID = nodeUIDs[u];
            uint32_t targetUID = nodeUIDs[v];

            lemon::ListGraph::Node nd;
            if (population.at(sourceUID).getState()  == Specie::S)
            {
                population.at(sourceUID).infect(time);
                nd = u;

            }
            else if (population.at(targetUID).getState()  == Specie::S)
            {
                population.at(targetUID).infect(time);
                nd = v;
            }

            for(lemon::ListGraph::IncEdgeIt ieIt(network, nd); ieIt!=lemon::INVALID; ++ieIt)
            {
                transmissionRates[ieIt] = 0;
                lemon::ListGraph::Node incNode = network.v(ieIt);
                uint32_t nodeUID = nodeUIDs[incNode];
                if(population.at(nodeUID).getState() == Specie::S)
                {
                    transmissionRates[ieIt] = transmissionRate;//transmitDistribution(generator);

                }
            }
            break;

        }

    }
}

void ContactNetwork::executeDeath(double rStart, double rBound)
{
    double result = rStart;
    for(lemon::ListGraph::NodeIt nIt(network); nIt!=lemon::INVALID; ++nIt)
    {
        uint32_t nodeUID = nodeUIDs[nIt];
        result += population.at(nodeUID).getDeathRate();

        if (result >= rBound)
        {
            removeNode(nIt);
            break;
        }
    }
}

void ContactNetwork::executeBirth(double rStart, double rBound)
{
    double result = rStart;
    result += birthRate;
    if (result >= rBound)
    {
        unsigned char maxContacts = maxContactsDistribution(generator);
        double  dRate = deathRate;//deathRateRateDistribution(generator);
        double newContRate = newContactRate;//newContactRateDistribution(generator);
        double looseContRate = looseContactRate;//looseContactRateDistribution(generator);
        Specie sp = Specie(maxContacts, 0, dRate, newContRate, looseContRate);
        uint32_t  id = UniqueID().id;
        population.emplace(id, sp);
        lemon::ListGraph::Node n = network.addNode();
        nodeUIDs[n] = id;

        n = complement.addNode();

        //break;
    }
}

double  ContactNetwork::getEdgeAdditionRate(lemon::ListGraph::Edge complementEdge) const
{
    lemon::ListGraph::Node u = complement.u(complementEdge);
    lemon::ListGraph::Node v = complement.v(complementEdge);

    uint32_t sourceUID = nodeUIDs[network.nodeFromId(complement.id(u))];
    uint32_t targetUID = nodeUIDs[network.nodeFromId(complement.id(v))];

    double sourceRate = population.at(sourceUID).getNewContactRate();
    double targetRate = population.at(targetUID).getNewContactRate();

    double result = 0;
    if (sourceRate > 0 && targetRate > 0)
    {
        result = (sourceRate + targetRate) / 2;
    }
    return result;
}

double  ContactNetwork::getEdgeDelitionRate(lemon::ListGraph::EdgeIt networkEdgeIt) const
{
    lemon::ListGraph::Node u = network.u(networkEdgeIt);
    lemon::ListGraph::Node v = network.v(networkEdgeIt);

    uint32_t sourceUID = nodeUIDs[u];
    uint32_t targetUID = nodeUIDs[v];
    double result = (population.at(sourceUID).getLooseContactRate() + population.at(targetUID).getLooseContactRate()) / 2;

    return result;

}

void ContactNetwork::init(size_t nInfected, size_t nSusceptible, size_t nEdges, int maxContactsA, int MaxContactsB,
                          double transmRate, double newContRate, double looseContRate, double dRate, double bRate)
{

    std::random_device rDev;
    generator = std::mt19937_64(rDev());
    //generator.seed(::time(NULL)); //to change the seed for every run
    generator.seed(3);


    //transmitDistribution = std::exponential_distribution<double> (8);
    maxContactsDistribution = std::uniform_int_distribution<unsigned char>(maxContactsA, MaxContactsB);
    transmissionRate = transmRate;
    //newContactRateDistribution = std::exponential_distribution<double>(0.5);
    newContactRate = newContRate;
    //looseContactRateDistribution = std::exponential_distribution<double>(0.5);
    looseContactRate = looseContRate;
    //deathRateRateDistribution = std::exponential_distribution<double> (10);
    deathRate = dRate;
    birthRate = bRate;

    size_t  nPopulation = nInfected + nSusceptible;

    lemon::FullGraph fullG(nPopulation);
    lemon::GraphCopy<lemon::FullGraph, lemon::ListGraph> cg(fullG, complement);
    cg.run();


    for (size_t i = 0; i < nPopulation; i ++)
    {
        unsigned char maxContacts = maxContactsDistribution(generator);
        double  dRate = deathRate;//deathRateRateDistribution(generator);
        double newContRate = newContactRate;//newContactRateDistribution(generator);
        double looseContRate = looseContactRate;//looseContactRateDistribution(generator);
        //std::cout << "new: " << newContRate << " loose: " << looseContRate << std::endl;
        Specie::State st = Specie::S;
        Specie sp = Specie(maxContacts, 0, dRate, newContRate, looseContRate, st);
        uint32_t  id = UniqueID().id;
        population.emplace(id, sp);
        lemon::ListGraph::Node n = network.addNode();
        nodeUIDs[n] = id;
    }

    auto it = population.begin();
    for (size_t i = 0; i < nInfected; i ++)
    {
        it->second.infect(0);
        it++;
    }


    //size_t nComplEdges = nPopulation * (nPopulation - 1) / 2;
    int maxEdgeId = complement.maxEdgeId();

    for (size_t i = 0; i < nEdges; i ++)
    {
        std::uniform_int_distribution<int> dist(1, maxEdgeId);
        int edgeId = dist(generator);

        lemon::ListGraph::Edge cEdge = complement.edgeFromId(edgeId);
        if (complement.valid(cEdge) &&  getEdgeAdditionRate(cEdge) > 0)
        {

            double tRate = 0;
            lemon::ListGraph::Node u = complement.u(cEdge);
            lemon::ListGraph::Node v = complement.v(cEdge);

            // find nodes in network
            lemon::ListGraph::Node cu = network.nodeFromId(complement.id(u));
            lemon::ListGraph::Node cv = network.nodeFromId(complement.id(v));

            uint32_t sourceUID = nodeUIDs[cu];
            uint32_t targetUID = nodeUIDs[cv];

            if ((population.at(sourceUID).getState()  == Specie::S  && population.at(targetUID).getState() == Specie::I) ||
                (population.at(sourceUID).getState()  == Specie::I  && population.at(targetUID).getState() == Specie::S))
            {

                tRate = transmissionRate;//transmitDistribution(generator);//expDist(generator);
            }
            addEdge(cEdge, tRate);

            //std::cout<<edgeId<<" " << lemon::countEdges(complement) << " "<<complement.maxEdgeId()<<std::endl;
        }
    }
    generator.seed(::time(NULL));
    //std::cout << "infected: " << countByState(Specie::I) <<std::endl;
    //std::cout << "edges: " << lemon::countEdges(network) <<std::endl;


}



size_t  ContactNetwork::countEdges() const
{
    return lemon::countEdges(network);
}

ContactNetwork& ContactNetwork::operator=(const ContactNetwork& other)
{
    if (this != &other)
    { // self-assignment check expected
       this->population = other.population;
       lemon::GraphCopy<lemon::ListGraph, lemon::ListGraph> cg(other.network, this->network);
       cg.edgeMap(other.transmissionRates, this->transmissionRates);
       cg.nodeMap(other.nodeUIDs, this->nodeUIDs);
       cg.run();

       lemon::GraphCopy<lemon::ListGraph, lemon::ListGraph> cg2(other.complement, this->complement);
       cg2.run();

       this->birthRate = other.birthRate;

       this->generator = other.generator;
       //this->transmitDistribution = other.transmitDistribution;
       this->transmissionRate = other.transmissionRate;
       this->maxContactsDistribution = other.maxContactsDistribution;
       //this->newContactRateDistribution = other.newContactRateDistribution;
       this->newContactRate = other.newContactRate;
       //this->looseContactRateDistribution = other.looseContactRateDistribution;
        this->looseContactRate = other.looseContactRate;
       //this->deathRateRateDistribution = other.deathRateRateDistribution;
        this->deathRate = other.deathRate;
    }
    return *this;
}
void ContactNetwork::copyNetwork (const ContactNetwork& other)
{
    *this = other;

}

std::vector<size_t> ContactNetwork::getDegreeDistribution()
{
    std::vector<size_t> result;
    for(lemon::ListGraph::NodeIt nIt(network); nIt!=lemon::INVALID; ++nIt)
    {
        size_t degree =0 ;
        for(lemon::ListGraph::IncEdgeIt e(network, nIt); e!=lemon::INVALID; ++e)
        {
            ++ degree;
        }
        result.push_back(degree);
    }
    return result;
}