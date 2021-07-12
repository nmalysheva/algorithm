//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include "ContactNetwork.h"
#include "UniqueID.h"
#include <random>
#include <lemon/full_graph.h>
#include <unistd.h>
#include <algorithm>
#include <lemon/adaptors.h>
#include <stdexcept>
#include <iostream>
#include <fstream>

void ContactNetwork::init(size_t nInfected, size_t nSusceptible, size_t nEdges, int maxContactsL, int MaxContactsU,
                          double transmRate, double newContRate, double looseContRate, double diagnRate, double dRate, double bRate)
{
    initRandomGenerator();
    //generator.seed(::time(NULL)); //to change the seed for every run
    generator.seed(3);

    initRates(maxContactsL, MaxContactsU, transmRate, newContRate, looseContRate,  diagnRate, dRate, bRate);

    size_t  nPopulation = nInfected + nSusceptible;
    initComplementNetwork(nPopulation);

    size_t nInf = 0;

    /* rates of assemble and disassemble edges are sampled from distributions
     *
     */
    std::uniform_real_distribution<double> lcd(looseContactRate*0.2, looseContactRate);
    std::uniform_real_distribution<double> ncd(newContactRate*0.2, newContactRate);
    isContactsLimited = false;
    for (size_t i = 0; i < nPopulation; i ++)
    {
        size_t maxContacts = maxContactsL;
        if (maxContacts < nPopulation - 1)
        {
            isContactsLimited = true;
        }

        double dRate = deathRate;//deathRateRateDistribution(generator);

        double looseContRate = looseContactRate;//lcd(generator);

        double newContRate = newContactRate;//ncd(generator);

        Specie::State st = Specie::S;
        Specie sp = Specie(maxContacts, 0, 0, newContRate, looseContRate, st);
        lemon::ListGraph::Node newNode = network.addNode();
        diagnosisRates[newNode] = 0;

        if (nInf < nInfected)
        {
            sp.changeState(Specie::I, 0);
            diagnosisRates[newNode] = diagnosisRate;
            sp.setDeathRate(dRate);
            nInf ++;

        }

        population.emplace(nodeIdMap[newNode], sp);


    }

    int maxNumberOfEdges = complement.maxEdgeId();

    //for ben structure
    //std::ofstream out;
    //out.open("Ben_Cont_Dyn_SSA.txt");

    for (size_t i = 0; i < nEdges; i ++)
    {
        std::uniform_int_distribution<int> dist(0, maxNumberOfEdges);
        int edgeId = dist(generator);
        lemon::ListGraph::Edge cEdge = complement.edgeFromId(edgeId);
        if (complement.valid(cEdge) &&  getEdgeAdditionRate(cEdge) > 0)
        {
            std::pair<int, int> b = addEdge(cEdge);
            //for ben structure
            //out << "0 " << b.first << " " << b.second << " True" << std::endl ;
        }
    }
//    out.close();
    generator.seed(::time(nullptr) * getpid()); //reset generator
}

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


std::vector<std::pair<double, lemon::ListGraph::Edge>> ContactNetwork::getTransmissionRateSum() const
{
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propCumSum;
    propCumSum.reserve(1e6 + 1);

    //element <0, INVALID>
    std::pair<double, lemon::ListGraph::Edge> invalidElem {0, lemon::ListGraph::Edge (lemon::INVALID)};
    propCumSum.push_back(invalidElem);

    size_t index = 1;
    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        double rate = transmissionRates[eIt];
        if (rate > 0)
        {
            propCumSum.emplace_back(propCumSum.at(index - 1).first + rate, eIt);
            index ++;
        }
    }
    propCumSum.shrink_to_fit();
    return propCumSum;

}

std::vector<std::pair<double, lemon::ListGraph::Node>> ContactNetwork:: getDiagnosisRateSum()const
{
    std::vector<std::pair<double, lemon::ListGraph::Node>> propCumSum;
    propCumSum.reserve(1e6 + 1);

    //element <0, INVALID>
    std::pair<double, lemon::ListGraph::Node> invalidElem {0, lemon::ListGraph::Node (lemon::INVALID)};
    propCumSum.push_back(invalidElem);

    size_t index = 1;
    for (lemon::ListGraph::NodeIt nIt(network); nIt != lemon::INVALID; ++nIt)
    {
        double rate = diagnosisRates[nIt];
        if (rate > 0)
        {
            propCumSum.emplace_back(propCumSum.at(index - 1).first + rate, nIt);
            index ++;
        }
    }
    propCumSum.shrink_to_fit();
    return propCumSum;
}

std::vector<std::pair<double, lemon::ListGraph::Edge>> ContactNetwork::getEdgeDeletionRateSum()const
{
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propCumSum;
    propCumSum.reserve(population.size() * (population.size() - 1) / 2); //TODO ADJUST TO POPULATION SIZE

    //element <0, INVALID>
    std::pair<double, lemon::ListGraph::Edge> invalidElem {0, lemon::ListGraph::Edge (lemon::INVALID)};
    propCumSum.push_back(invalidElem);

    size_t index = 1;
    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        double rate = getEdgeDeletionRate(eIt);
        if (rate > 0)
        {
            propCumSum.emplace_back(propCumSum.at(index - 1).first + rate, eIt);
            //std::cout << "elem = " << elem.first << std::endl;
            index ++;
        }
    }
    propCumSum.shrink_to_fit();
    return propCumSum;

}


std::vector<std::pair<double, lemon::ListGraph::Edge>> ContactNetwork::getEdgeAdditionRateSum()const
{
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propCumSum;
    propCumSum.reserve(population.size() * (population.size() - 1) / 2);

    //element <0, INVALID>
    std::pair<double, lemon::ListGraph::Edge> invalidElem {0, lemon::ListGraph::Edge (lemon::INVALID)};
    propCumSum.push_back(invalidElem);

    size_t index = 1;
    for (lemon::ListGraph::EdgeIt eIt(complement); eIt != lemon::INVALID; ++eIt)
    {
        double rate = getEdgeAdditionRate(eIt);
        if (rate > 0)
        {
            propCumSum.emplace_back(propCumSum.at(index - 1).first + rate, eIt);
            index++;
        }
    }
    propCumSum.shrink_to_fit();
    return propCumSum;
}

std::vector<std::pair<double, lemon::ListGraph::Node>> ContactNetwork::getDeathRateSum()const
{
    std::vector<std::pair<double, lemon::ListGraph::Node>> propCumSum;
    propCumSum.reserve(1e6 + 1);

    //element <0, INVALID>
    std::pair<double, lemon::ListGraph::Node> invalidElem {0, lemon::ListGraph::Node (lemon::INVALID)};
    propCumSum.push_back(invalidElem);

    size_t index = 1;
    for(lemon::ListGraph::NodeIt nIt(network); nIt!=lemon::INVALID; ++nIt)
    {
        uint32_t nodeUID = nodeIdMap[nIt];
        double rate = population.at(nodeUID).getDeathRate();

        propCumSum.emplace_back(propCumSum.at(index - 1).first + rate, nIt);
        index++;
    }

    propCumSum.shrink_to_fit();
    return propCumSum;

}

double ContactNetwork::getTransmissionRateLimit() const
{
    return transmissionRate;
}


double ContactNetwork::getBirthRateSum()const
{
    return birthRate;
}

size_t ContactNetwork::size() const
{
    return lemon::countNodes(network);
}

std::pair<int, int> ContactNetwork::addEdge(lemon::ListGraph::Edge &complementEdge)
{
    if (!complement.valid(complementEdge))
    {
        std::string msg = "ERROR: INVALID EDGE TO ADD!";
        throw std::domain_error(msg);
    }
    //nodes of the given edge in a complement graph
    lemon::ListGraph::Node complU = complement.u(complementEdge);
    lemon::ListGraph::Node complV = complement.v(complementEdge);

    // find respective nodes in network
    lemon::ListGraph::Node networkU = network.nodeFromId(complement.id(complU));
    lemon::ListGraph::Node networkV = network.nodeFromId(complement.id(complV));


    //erase given edge from complement graph
    complement.erase(complementEdge);

    //add new edge to the network with given transmission rate
    lemon::ListGraph::Edge newEdge = network.addEdge(networkU, networkV);

    std::pair<int, int> result = std::make_pair(network.id(networkU), network.id(networkV));

    //get IDs of nodes. Have type int because it was originally in the definition of the IdMap class
    int sourceUID = nodeIdMap[networkU];
    int targetUID = nodeIdMap[networkV];

    //for these nodes increase number of contacts
    population.at(sourceUID).incNumberOfContacts();
    population.at(targetUID).incNumberOfContacts();

    //calculate transmission rate
    double trRate = 0;
    if ((population.at(sourceUID).getState()  == Specie::S  && population.at(targetUID).getState() == Specie::I) ||
        (population.at(sourceUID).getState()  == Specie::I  && population.at(targetUID).getState() == Specie::S))
    {
        trRate = transmissionRate;
    }

    if ((population.at(sourceUID).getState()  == Specie::S  && population.at(targetUID).getState() == Specie::D) ||
        (population.at(sourceUID).getState()  == Specie::D  && population.at(targetUID).getState() == Specie::S))
    {
        //for edges Diagnosed "D"- Susceptible "S" transmission rate is 50% lower than for edges "I" - "S"
        trRate = transmissionRate * 0.5;
    }

    transmissionRates[newEdge] = trRate;
    return result;
}

std::pair<int, int> ContactNetwork::removeEdge(lemon::ListGraph::Edge &edge)
{
    if (!network.valid(edge))
    {
        std::string msg = "ERROR: INVALID EDGE TO DEL!";
        throw std::domain_error(msg);
    }
    //nodes of the given edge in a graph
    lemon::ListGraph::Node networkU = network.u(edge);
    lemon::ListGraph::Node networkV = network.v(edge);

    int sourceUID = network.id(networkU);
    int targetUID = network.id(networkV);

    std::pair<int, int> result = std::make_pair(sourceUID, targetUID);

    network.erase(edge); // erase from network

    // find respective nodes in complement graph
    lemon::ListGraph::Node complU = complement.nodeFromId(sourceUID);
    lemon::ListGraph::Node complV = complement.nodeFromId(targetUID);
    complement.addEdge(complU, complV);   // add to complement graph

    // after removing an edge decrease num. of actual contacts for each of the incident nodes
    population.at(sourceUID).decNumberOfContacts();
    population.at(targetUID).decNumberOfContacts();

    return result;

}

void ContactNetwork::removeNode(lemon::ListGraph::Node &node)
{
    //delete all edges to given node in network
    lemon::ListGraph::IncEdgeIt ieIt(network, node);
    while (ieIt != lemon::INVALID)
    {
        lemon::ListGraph::IncEdgeIt tmpIt = ieIt;
        lemon::ListGraph::Node oppositeNode = network.oppositeNode(node, ieIt);
        population.at(nodeIdMap[oppositeNode]).decNumberOfContacts();
        ++ieIt;
        network.erase(tmpIt);

    }

    //find respective node in complement graph
    lemon::ListGraph::Node complementNode = complement.nodeFromId(network.id(node));

    //delete all edges to given node in complement graph
    lemon::ListGraph::IncEdgeIt cieIt(complement, complementNode);
    while (cieIt != lemon::INVALID)
    {
        lemon::ListGraph::IncEdgeIt tmpIt = cieIt;
        ++cieIt;
        complement.erase(tmpIt);
    }

    int sourceUID = nodeIdMap[node];
    population.erase(sourceUID);
    network.erase(node);
    complement.erase(complementNode);

}

void ContactNetwork::executeTransmission(lemon::ListGraph::Edge & edge, double time)
{
    lemon::ListGraph::Node networkU = network.u(edge);
    lemon::ListGraph::Node networkV = network.v(edge);

    int sourceUID = nodeIdMap[networkU];
    int targetUID = nodeIdMap[networkV];
    lemon::ListGraph::Node infectedNode;

    //TODO Fool check: only susceptible (S) can be infected (I), etc. Prevent all impossible transactions aka
    // D->I, etc. Is it really necessary here?
    if (population.at(sourceUID).getState()  == Specie::S)
    {
        population.at(sourceUID).changeState(Specie::I, time);
        population.at(sourceUID).setDeathRate(deathRate);
        infectedNode = networkU;
    }

    else if (population.at(targetUID).getState()  == Specie::S)
    {
        population.at(targetUID).changeState(Specie::I, time);
        population.at(targetUID).setDeathRate(deathRate);
        infectedNode = networkV;
    }

    diagnosisRates[infectedNode] = diagnosisRate;
    // no adaptivity
    for(lemon::ListGraph::IncEdgeIt ieIt(network, infectedNode); ieIt!=lemon::INVALID; ++ieIt)
    {
        transmissionRates[ieIt] = 0;
        lemon::ListGraph::Node neighbourNode = network.v(ieIt);
        uint32_t nodeUID = nodeIdMap[neighbourNode];
        if(population.at(nodeUID).getState() == Specie::S)
        {
            transmissionRates[ieIt] = transmissionRate;//transmitDistribution(generator);

        }
    }

}

void ContactNetwork::executeDiagnosis(lemon::ListGraph::Node & node, double time)
{
    int nodeUID = nodeIdMap[node];
    population.at(nodeUID).changeState(Specie::D, time);

    diagnosisRates[node] = 0;

    //std::cout << "DIAGNOSIS!!!!" << std::endl;
    //adaptivity: as soon as diagnosed, cut all contacts and reduce
    //new contact rate to 3%
    lemon::ListGraph::IncEdgeIt ieIt(network, node);
    while (ieIt != lemon::INVALID)
    {
        lemon::ListGraph::Edge tmpIt(ieIt);
        ++ieIt;
        removeEdge(tmpIt);
    }
    //std::cout << "num. of contacts: " << population.at(nodeUID).getNumberOfContacts() << std::endl;
    //std::cout << "old rate: " << population.at(nodeUID).getNewContactRate() << std::endl;
    population.at(nodeUID).setNewContactRate(
            population.at(nodeUID).getNewContactRate() * 0.3);
    //std::cout << "new rate: " << population.at(nodeUID).getNewContactRate() << std::endl;

    //population.at(nodeUID).setLooseContactRate(
      //      population.at(nodeUID).getLooseContactRate() * 4);

}

void ContactNetwork::executeDeath(lemon::ListGraph::Node & node)
{
    removeNode(node);

    // after removing node from the population decrease max. number of contacts for each specie.
    //TODO: may be not integrate max. num of contacts to the specie. Leave it in the network
    for (auto& it: population)
    {
        it.second.setMaxNumberOfContacts(it.second.getMaxNumberOfContacts() - 1);
    }

}


//TODO outdated
void ContactNetwork::executeBirth(double rStart, double rBound)
{
    double result = rStart;
    result += birthRate;
    if (result >= rBound)
    {
        size_t maxContacts = maxContactsDistribution(generator);
        double  dRate = deathRate;//deathRateRateDistribution(generator);
        double newContRate = newContactRate;//newContactRateDistribution(generator);
        double looseContRate = looseContactRate;//looseContactRateDistribution(generator);
        Specie sp = Specie(maxContacts, 0, dRate, newContRate, looseContRate);

        lemon::ListGraph::Node newNode = network.addNode();
        population.emplace(nodeIdMap[newNode], sp);

        newNode = complement.addNode();
    }

    if (! isContactsLimited)
    {
        isContactsLimited = false;
        for (auto& it: population)
        {
            if (it.second.getMaxNumberOfContacts() < population.size() - 1)
            {
                isContactsLimited = true;
                break;
            }

        }
    }
}

double  ContactNetwork::getEdgeAdditionRate(const lemon::ListGraph::Edge &complementEdge) const
{
    lemon::ListGraph::Node complU = complement.u(complementEdge);
    lemon::ListGraph::Node complV = complement.v(complementEdge);

    int sourceUID = nodeIdMap[network.nodeFromId(complement.id(complU))];
    double sourceRate = population.at(sourceUID).getNewContactRate();
    //sourceRate = sourceRate * (population.at(sourceUID).getMaxNumberOfContacts() - population.at(sourceUID).getNumberOfContacts());

    int targetUID = nodeIdMap[network.nodeFromId(complement.id(complV))];
    double targetRate = population.at(targetUID).getNewContactRate();

    // rate of adding an edge is a multiplication of rates of both incident nodes
    double result = sourceRate * targetRate;
    return result;
}

double  ContactNetwork::getEdgeDeletionRate(const lemon::ListGraph::Edge &networkEdge) const
{
    lemon::ListGraph::Node networkU = network.u(networkEdge);
    lemon::ListGraph::Node networkV = network.v(networkEdge);

    int sourceUID = nodeIdMap[networkU];
    int targetUID = nodeIdMap[networkV];

    double sourceRate = population.at(sourceUID).getLooseContactRate();
    double targetRate = population.at(targetUID).getLooseContactRate();

    // rate of deleting an edge is a multiplication of rates of both incident nodes
    double result = sourceRate * targetRate;
    return result;

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

       //cg.nodeMap(other.nodeIdMap, this->nodeIdMap);
       cg.run();

       lemon::GraphCopy<lemon::ListGraph, lemon::ListGraph> cg2(other.complement, this->complement);
       //cg2.nodeMap(other.complementAdjacentEdges, this->complementAdjacentEdges);
       cg2.run();

       this->nodeIdMap = lemon::IdMap<lemon::ListGraph, lemon::ListGraph::Node>(this->network);

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

        this ->isContactsLimited = other.isContactsLimited;
    }
    return *this;
}

std::vector<size_t> ContactNetwork::getDegreeDistribution() const
{
    std::vector<size_t> result;
    result.reserve(1e6); //reserving space for vector.

    for(lemon::ListGraph::NodeIt nIt(network); nIt!=lemon::INVALID; ++nIt)
    {
        size_t degree = 0 ;

        for(lemon::ListGraph::IncEdgeIt e(network, nIt); e!=lemon::INVALID; ++e)
        {
            degree++;
        }

        result.push_back(degree);
    }

    return result;
}

std::vector<specieState> ContactNetwork::getNetworkState() const
{
    std::vector<specieState> result;
    result.reserve(this->size()); //reserving space for vector.

    for(lemon::ListGraph::NodeIt nIt(network); nIt!=lemon::INVALID; ++nIt)
    {
        specieState spState;

        int nodeId = nodeIdMap[nIt];

        spState.id = nodeId;
        Specie::State st = population.at(nodeId).getState();

        spState.state = st;

        std::vector<int> neighbors;
        neighbors.reserve(maxContactsLimitU);
        for(lemon::ListGraph::IncEdgeIt e(network, nIt); e!=lemon::INVALID; ++e)
        {
            lemon::ListGraph::Node neighbor = network.oppositeNode(nIt, e);
            neighbors.push_back(nodeIdMap[neighbor]);
        }
        neighbors.shrink_to_fit();

        spState.contacts = neighbors;

        result.push_back(spState);
    }

    return result;
}

size_t ContactNetwork::getMaxContactsLimitOfInfected(double t)const
{
    size_t result = 0;

    for (auto& it: population)
    {
        if (it.second.getState() == Specie::State::I || it.second.getState() == Specie::State::D)
        {
            result += it.second.getNumberOfContactsLimit(t);
        }
    }
    return result;

}

size_t ContactNetwork::getMaxContactsLimitOfSusceptible(double t)const
{
    size_t result = 0;

    for (auto& it: population)
    {
        if (it.second.getState() == Specie::State::S)
        {
            result += it.second.getNumberOfContactsLimit(t);
            //std::cout << "max cont: " << it.second.getNumberOfContactsLimit(t) << std::endl;
        }
    }
    return result;

}

void ContactNetwork::initRandomGenerator()
{
    std::random_device rDev;
    generator = std::mt19937_64(rDev());
}

void ContactNetwork::initRates(int maxContactsL, int MaxContactsU, double transmRate, double newContRate, double looseContRate,
                               double diagnRate, double dRate, double bRate)
{
    maxContactsLimitL = maxContactsL;
    maxContactsLimitU = MaxContactsU;
    transmissionRate = transmRate;

    //std::uniform_real_distribution<> dis(0.5, newContRate);
    newContactRate = newContRate;//dis(generator);//newContRate;
    //std::uniform_real_distribution<> dis2(0.5, looseContRate);
    looseContactRate = looseContRate;
    deathRate = dRate;
    birthRate = bRate;
    diagnosisRate = diagnRate;
}

void ContactNetwork::initComplementNetwork(size_t nPopulation)
{
    lemon::FullGraph fullG(nPopulation);
    lemon::GraphCopy<lemon::FullGraph, lemon::ListGraph> cg(fullG, complement);
    cg.run();

}


std::vector<BenStructure> ContactNetwork::getBenStructure(double t)
{
    std::vector<BenStructure> result;
    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        lemon::ListGraph::Node networkU = network.u(eIt);
        lemon::ListGraph::Node networkV = network.v(eIt);

        int sourceUID = nodeIdMap[networkU];
        int targetUID = nodeIdMap[networkV];
        if (sourceUID < targetUID )
        {
            result.emplace_back(t, sourceUID, targetUID,true);
        }
        else
        {
            result.emplace_back(t, targetUID, sourceUID,true);
        }
    }
    return result;

}

BenStructure::BenStructure(double _t,int _u, int _v, bool _st)
{
    t = _t;
    u = _u;
    v = _v;
    state = _st;
}

bool BenStructure::isValid()
{
    return (u > -1) && (v > -1);
}

lemon::ListGraph::Edge ContactNetwork::getComplementEdge(int a, int b)
{
    lemon::ListGraph::Node complU = complement.nodeFromId(a);
    lemon::ListGraph::Node complV = complement.nodeFromId(b);

    lemon::ListGraph::Edge e = lemon::findEdge(complement, complU, complV);
    return e;
}
lemon::ListGraph::Edge ContactNetwork::getEdge(int a, int b)
{

    lemon::ListGraph::Node networkU = network.nodeFromId(a);
    lemon::ListGraph::Node networkV = network.nodeFromId(b);
    lemon::ListGraph::Edge e = lemon::findEdge(network, networkU, networkV);
    return e;

}
