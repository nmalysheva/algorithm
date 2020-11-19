//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include "ContactNetwork.h"
#include "UniqueID.h"
#include <random>
#include <lemon/full_graph.h>
#include <map>
#include <set>
#include <unistd.h>
#include <algorithm>
#include <lemon/adaptors.h>
#include <stdexcept>
#include <iostream>

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

    isContactsLimited = false;
    for (size_t i = 0; i < nPopulation; i ++)
    {
        size_t maxContacts = maxContactsDistribution(generator);
        if (maxContacts < nPopulation - 1)
        {
            isContactsLimited = true;
        }

        double dRate = deathRate;//deathRateRateDistribution(generator);
        double newContRate = newContactRate;//ncrDistribution(generator);//newContactRateDistribution(generator);

        //double looseContRate = looseContactRate;//looseContactRate;//looseContactRateDistribution(generator);

        std::uniform_real_distribution<double> lcd(looseContactRate*0.2, looseContactRate*3);
        double looseContRate = lcd(generator);


        Specie::State st = Specie::S;
        //Specie sp = Specie(maxContacts, 0, dRate, newContRate, looseContRate, st);
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
    /*for(lemon::ListGraph::NodeIt nIt(complement); nIt!=lemon::INVALID; ++nIt)
    {
        complementAdjacentEdges[nIt] = countAdjacentEdges(nIt);
    }*/

    int maxNumberOfEdges = complement.maxEdgeId();

    // TODO: change edge adding with regarding to the amount of edges possible to add safely.
    for (size_t i = 0; i < nEdges; i ++)
    {
        std::uniform_int_distribution<int> dist(0, maxNumberOfEdges);
        int edgeId = dist(generator);
        lemon::ListGraph::Edge cEdge = complement.edgeFromId(edgeId);
        if (complement.valid(cEdge) &&  getEdgeAdditionRate(cEdge) > 0)
        {
            addEdge(cEdge);
        }
    }

    generator.seed(::time(nullptr) * getpid());
}

size_t ContactNetwork::countAdjacentEdges(const lemon::ListGraph::Node &complementNode) const
{
    size_t result = 0;

    int sourceUID = nodeIdMap[network.nodeFromId(complement.id(complementNode))];
    double sourceRate = population.at(sourceUID).getNewContactRate();
    if (sourceRate > 0)
    {
        for(lemon::ListGraph::IncEdgeIt eIt(complement, complementNode); eIt!=lemon::INVALID; ++eIt)
        {
            lemon::ListGraph::Node oppositeNode = complement.oppositeNode(complementNode, eIt);
            int oppositeUID = nodeIdMap[network.nodeFromId(complement.id(oppositeNode))];

            double oppositeRate = population.at(oppositeUID).getNewContactRate();
            if (oppositeRate > 0)
            {
                result++;
            }
        }
    }
    else
    {
        std::cout << "cae else " <<std::endl;
    }


    return result;
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


/*double  ContactNetwork::getTransmissionRateSum()const
{
    double result = 0;
    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        result += transmissionRates[eIt];
    }
    return result;
}*/
std::vector<std::pair<double, lemon::ListGraph::Edge>> ContactNetwork::getTransmissionRateSum()
{
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propCumSum;
    propCumSum.reserve(1e6 + 1); //TODO ADJUST TO POPULATION SIZE

    std::pair<double, lemon::ListGraph::Edge> elem2 {0, lemon::ListGraph::Edge (lemon::INVALID)};
    propCumSum.push_back(elem2);

    size_t index = 1;
    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        double rate = transmissionRates[eIt];
        if (rate > 0)
        {
            std::pair<double, lemon::ListGraph::Edge> elem {propCumSum.at(index - 1).first + rate, eIt};
            propCumSum.push_back(elem);
            index ++;
        }
    }
    propCumSum.shrink_to_fit();
    return propCumSum;

}

std::vector<std::pair<double, lemon::ListGraph::Node>> ContactNetwork:: getDiagnosisRateSum()const
{
    std::vector<std::pair<double, lemon::ListGraph::Node>> propCumSum;
    propCumSum.reserve(1e6 + 1); //TODO ADJUST TO POPULATION SIZE

    std::pair<double, lemon::ListGraph::Node> elem2 {0, lemon::ListGraph::Node (lemon::INVALID)};
    propCumSum.push_back(elem2);

    size_t index = 1;
    for (lemon::ListGraph::NodeIt nIt(network); nIt != lemon::INVALID; ++nIt)
    {
        double rate = diagnosisRates[nIt];
        if (rate > 0)
        {
            std::pair<double, lemon::ListGraph::Node> elem {propCumSum.at(index - 1).first + rate, nIt};
            propCumSum.push_back(elem);
            index ++;
        }
    }
    propCumSum.shrink_to_fit();
    return propCumSum;
}

std::vector<std::pair<double, lemon::ListGraph::Edge>> ContactNetwork::getEdgeDeletionRateSum()const
{
    std::vector<std::pair<double, lemon::ListGraph::Edge>> propCumSum;
    propCumSum.reserve(1e6 + 1); //TODO ADJUST TO POPULATION SIZE

    std::pair<double, lemon::ListGraph::Edge> elem2 {0, lemon::ListGraph::Edge (lemon::INVALID)};
    propCumSum.push_back(elem2);

    size_t index = 1;
    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        double rate = getEdgeDeletionRate(eIt);
        if (rate > 0)
        {
            std::pair<double, lemon::ListGraph::Edge> elem {propCumSum.at(index - 1).first + rate, eIt};
            propCumSum.push_back(elem);
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
    propCumSum.reserve(1e6 + 1); //TODO ADJUST TO POPULATION SIZE

    std::pair<double, lemon::ListGraph::Edge> elem2 {0, lemon::ListGraph::Edge (lemon::INVALID)};
    propCumSum.push_back(elem2);

    size_t index = 1;
    for (lemon::ListGraph::EdgeIt eIt(complement); eIt != lemon::INVALID; ++eIt)
    {
        double rate = getEdgeAdditionRate(eIt);
        if (rate > 0)
        {
            std::pair<double, lemon::ListGraph::Edge> elem {propCumSum.at(index - 1).first + rate, eIt};
            propCumSum.push_back(elem);
            index++;
        }
    }
    propCumSum.shrink_to_fit();
    return propCumSum;
}

/*double  ContactNetwork::getDeathRateSum()const
{
    double result = 0;
    for (auto& it: population)
    {
        result += it.second.getDeathRate();

    }

    return result;
}*/

std::vector<std::pair<double, lemon::ListGraph::Node>> ContactNetwork::getDeathRateSum()const
{
    std::vector<std::pair<double, lemon::ListGraph::Node>> propCumSum;
    propCumSum.reserve(1e6 + 1); //TODO ADJUST TO POPULATION SIZE

    std::pair<double, lemon::ListGraph::Node> elem2 {0, lemon::ListGraph::Node (lemon::INVALID)};
    propCumSum.push_back(elem2);

    size_t index = 1;
    for(lemon::ListGraph::NodeIt nIt(network); nIt!=lemon::INVALID; ++nIt)
    {
        uint32_t nodeUID = nodeIdMap[nIt];
        double rate = population.at(nodeUID).getDeathRate();

        std::pair<double, lemon::ListGraph::Node> elem {propCumSum.at(index - 1).first + rate, nIt};
        propCumSum.push_back(elem);
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


    //size_t before = lemon::countEdges(complement);
    //erase given edge from complement graph
    complement.erase(complementEdge);

    /*size_t after = lemon::countEdges(complement);

    if (before == after)
    {
        std::cout << "err;" << complement.id(complU) << "; " << complement.id(complV) << std::endl;
    }*/
    //std::cout << "before: " << lemon::countEdges(complement) << "; ";
    //std::cout << "after: " << lemon::countEdges(complement) << std::endl;
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
        trRate = transmissionRate * 0.5;
    }

    transmissionRates[newEdge] = trRate;

    /*if (population.at(sourceUID).getNewContactRate() == 0)
    {
        for (lemon::ListGraph::IncEdgeIt eIt(complement, complU); eIt != lemon::INVALID; ++eIt)
        {
            lemon::ListGraph::Node oppositeNode = complement.oppositeNode(complU, eIt);
            if (complementAdjacentEdges[oppositeNode] > 0)
            {
                complementAdjacentEdges[oppositeNode]--;
            }


        }
    }

    if (population.at(targetUID).getNewContactRate() == 0)
    {
        for (lemon::ListGraph::IncEdgeIt eIt(complement, complV); eIt != lemon::INVALID; ++eIt)
        {
            lemon::ListGraph::Node oppositeNode = complement.oppositeNode(complV, eIt);
            if (complementAdjacentEdges[oppositeNode] > 0)
            {
                complementAdjacentEdges[oppositeNode]--;
            }

        }
    }

    complementAdjacentEdges[complU] = countAdjacentEdges(complU);
    complementAdjacentEdges[complV] = countAdjacentEdges(complV);*/

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

    //double oldSourceRate = population.at(sourceUID).getNewContactRate();
    //double oldTargetRate = population.at(targetUID).getNewContactRate();
    std::pair<int, int> result = std::make_pair(sourceUID, targetUID);

    network.erase(edge); // erase from network

    // find respective nodes in complement graph
    lemon::ListGraph::Node complU = complement.nodeFromId(sourceUID);
    lemon::ListGraph::Node complV = complement.nodeFromId(targetUID);
    complement.addEdge(complU, complV);   // add to complement graph

    population.at(sourceUID).decNumberOfContacts();
    population.at(targetUID).decNumberOfContacts();


    /*if (oldSourceRate == 0)
    {
        for (lemon::ListGraph::IncEdgeIt eIt(complement, complU); eIt != lemon::INVALID; ++eIt)
        {
            lemon::ListGraph::Node oppositeNode = complement.oppositeNode(complU, eIt);
            if (population.at(nodeIdMap[oppositeNode]).getNewContactRate() > 0)
            {
                complementAdjacentEdges[oppositeNode]++;
            }
        }
    }

    if (oldTargetRate == 0)
    {
        for (lemon::ListGraph::IncEdgeIt eIt(complement, complV); eIt != lemon::INVALID; ++eIt)
        {
            lemon::ListGraph::Node oppositeNode = complement.oppositeNode(complV, eIt);
            if (population.at(nodeIdMap[oppositeNode]).getNewContactRate() > 0)
            {
                complementAdjacentEdges[oppositeNode]++;
            }
        }
    }
    complementAdjacentEdges[complU] = countAdjacentEdges(complU);
    complementAdjacentEdges[complV] = countAdjacentEdges(complV);*/
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
        //lemon::ListGraph::Node oppositeNode = complement.oppositeNode(complementNode, tmpIt);
        ++cieIt;
        complement.erase(tmpIt);
        //complementAdjacentEdges[oppositeNode] = countAdjacentEdges(oppositeNode);
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
    // D->I, etc.
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

    //adaptivity: as soon as diagnosed, cut all contacts and decrease
    //new contact rate to 30%
    lemon::ListGraph::IncEdgeIt ieIt(network, node);
    while (ieIt != lemon::INVALID)
    {
        lemon::ListGraph::Edge tmpIt(ieIt);
        ++ieIt;
        removeEdge(tmpIt);
    }
    population.at(nodeIdMap[node]).setNewContactRate(
            population.at(nodeIdMap[node]).getNewContactRate() * 0.3);

}

void ContactNetwork::executeDeath(lemon::ListGraph::Node & node)
{
    removeNode(node);

    for (auto& it: population)
    {
        it.second.setMaxNumberOfContacts(it.second.getMaxNumberOfContacts() - 1);
    }

}

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
    double result = 0;

    lemon::ListGraph::Node complU = complement.u(complementEdge);
    lemon::ListGraph::Node complV = complement.v(complementEdge);

    int sourceUID = nodeIdMap[network.nodeFromId(complement.id(complU))];
    double sourceRate = population.at(sourceUID).getNewContactRate();
    sourceRate = sourceRate * (population.at(sourceUID).getMaxNumberOfContacts() - population.at(sourceUID).getNumberOfContacts());

    int targetUID = nodeIdMap[network.nodeFromId(complement.id(complV))];
    double targetRate = population.at(targetUID).getNewContactRate();
    targetRate = targetRate * (population.at(targetUID).getMaxNumberOfContacts() - population.at(targetUID).getNumberOfContacts());
    if ( (sourceRate > 0) && (targetRate > 0) )
    {

        double sourceAddEdges = population.at(sourceUID).getMaxNumberOfContacts() - population.at(sourceUID).getNumberOfContacts(); //complementAdjacentEdges[complU];
        double targetAddEdges = population.at(targetUID).getMaxNumberOfContacts() - population.at(targetUID).getNumberOfContacts();// complementAdjacentEdges[complV];

        if ( (sourceAddEdges == 0) || (targetAddEdges == 0) )
        {
            //TODO throw exception
            std::cout << "LOLLL" <<std::endl;
            //std::cout << "sourceRate = " << sourceRate << "; targetRate = " << targetRate << std::endl;
            std::cout << "sourceAddEdges = " << sourceAddEdges << "; targetAddEdges = " << targetAddEdges << std::endl;
            std::cout << "--------------------------!!!!-----------------------" << std::endl;
            std::string msg = "ERROR: IcomplementAdjacentEdges is 0!";
            throw std::domain_error(msg);


        }
        else
        {
            sourceRate = sourceRate / sourceAddEdges;
            targetRate = targetRate / targetAddEdges;
        }

        result = sourceRate * targetRate;
        if (result > 1e+10)
        {
            std::cout << "sourceAddEdges = " << sourceAddEdges << "; targetAddEdges = " << targetAddEdges << std::endl;
            std::cout << "sourcerate = " << sourceRate << "; trate = " << targetRate << std::endl;
            std::cout << "s.max = " << population.at(sourceUID).getMaxNumberOfContacts() << "; t.max = " << population.at(targetUID).getMaxNumberOfContacts() << std::endl;
            std::cout << "s.curr = " << population.at(sourceUID).getNumberOfContacts() << "; t.curr = " << population.at(targetUID).getNumberOfContacts() << std::endl;
        }
    }

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

    //double result = (sourceRate + targetRate) / 2;
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
    result.reserve(maxContactsLimitU); //reserving space for vector.

    for(lemon::ListGraph::NodeIt nIt(network); nIt!=lemon::INVALID; ++nIt)
    {
        size_t degree = 0 ;

        for(lemon::ListGraph::IncEdgeIt e(network, nIt); e!=lemon::INVALID; ++e)
        {
            ++ degree;
        }

        result.push_back(degree);
    }

    return result;
}

/*size_t ContactNetwork::getAmountOfEdgesToAddSafe()
{
    //subgraph();


    std::vector<size_t> capacitiesCount = countCapacities();
    size_t nEdgesToAdd = 0;
    for (size_t i = 0; i < capacitiesCount.size(); i ++)
    {
        //std::cout << capacitiesCount.at(i) << " ";
        if (capacitiesCount.at(i) >= i + 1 + 1)
        {
            size_t nFullBlock = (i + 1 + 1);
            size_t eEdgesProBlock = nFullBlock * (nFullBlock - 1) / 2;
            nEdgesToAdd += capacitiesCount.at(i) / nFullBlock * eEdgesProBlock;
            capacitiesCount.at(i) = capacitiesCount.at(i) % nFullBlock;
        }
    }
    //std::cout <<std::endl;
    //std::cout << "after blocks:" << nEdgesToAdd <<std::endl;
    size_t curNodes = 0;
    for (size_t i = 0; i < capacitiesCount.size(); i ++)
    {
        //std::cout <<capacitiesCount.at(i) << " ";
        curNodes +=capacitiesCount.at(i);
    }

    if (capacitiesCount.size() > 0)
    {
        size_t in = capacitiesCount.size();
        //std::cout << in <<std::endl;

        while (in > 0)
        {
            if (capacitiesCount.at(in - 1) > 0)
            {
                size_t nRequiredForFullBlock = in - 1 + 1 + 1;
                //std::cout <<"nRequiredForFullBlock " << nRequiredForFullBlock<<std::endl;
                //std::cout <<"curNodes " << curNodes<<std::endl;
                if (nRequiredForFullBlock > curNodes)
                {
                    nRequiredForFullBlock = curNodes;
                    size_t currCap = curNodes - 1;
                    size_t toAdd = curNodes * (curNodes - 1) / 2;
                    //std::cout << toAdd << std::endl;
                    for (size_t i = 0; i < in - 1 + 1; i++)
                    {
                        if (capacitiesCount.at(i) > 0)
                        {
                            //std::cout << "i=" <<i <<" " <<capacitiesCount.at(i)  <<std::endl;;
                            toAdd -= currCap - std::min(i, currCap) * capacitiesCount.at(i);
                        }
                    }
                    nEdgesToAdd += toAdd;
                    curNodes = 0;
                    //std::cout << toAdd << std::endl;
                    break;
                }
                else
                {
                    size_t currentBlock = capacitiesCount.at(in - 1);
                    capacitiesCount.at(in - 1) = 0;
                    curNodes -= currentBlock;

                    size_t currCap = in - 1;
                    size_t toAdd = nRequiredForFullBlock * (nRequiredForFullBlock - 1) / 2;

                    int ind = in - 1;
                    while (currentBlock < nRequiredForFullBlock)
                    {
                        ind--;
                        if (ind < 0)
                        {
                            break;
                        }

                        size_t nCuradd = std::min(nRequiredForFullBlock - currentBlock, capacitiesCount.at(ind));
                        toAdd -= std::min(nRequiredForFullBlock - currentBlock, capacitiesCount.at(ind)) *
                                 (currCap - ind);
                        capacitiesCount.at(ind) -= nCuradd;
                        curNodes -= nCuradd;
                        currentBlock += nCuradd;

                    }
                    //std::cout << "toAdd = " << toAdd <<std::endl;
                    nEdgesToAdd += toAdd;

                }

            }
            in--;
        }
    }
    return nEdgesToAdd;
} */

size_t ContactNetwork::getAmountOfEdgesToAdd() const
{
    size_t nAdd = 0;
    for (lemon::ListGraph::EdgeIt eIt(complement); eIt != lemon::INVALID; ++eIt)
    {
        if (getEdgeAdditionRate(eIt) > 0)
        {
            nAdd++;
        }

    }
    return nAdd;

    size_t result = 0;

    /*for (lemon::ListGraph::NodeIt nIt(complement); nIt != lemon::INVALID; ++nIt)
    {
        result += complementAdjacentEdges[nIt];
    }

    result = result / 2;*/
    return result;
}

/*std::vector<size_t> ContactNetwork::countCapacities() const
{
    std::vector<size_t> capacitiesCount (maxContactsLimitU, 0);

    size_t  nNodes = lemon::countNodes(network); //TODO may be nNodes can be provided to the function already and there is no need to calc. it again
    for (lemon::ListGraph::NodeIt nIt(network); nIt != lemon::INVALID; ++nIt)
    {
        int nodeUID = nodeIdMap[nIt];
        size_t cap = population.at(nodeUID).getMaxNumberOfContacts() - population.at(nodeUID).getNumberOfContacts();
        if (cap > maxContactsLimitU)
        {
            //TODO throw exception
            std::cout << "capacity error" <<std::endl;
            std::cout <<nodeUID <<std::endl;
            std::cout << population.at(nodeUID).getMaxNumberOfContacts() <<" " << population.at(nodeUID).getNumberOfContacts() <<" ";
            std::cout << cap << std::endl;
        }

        if (cap > 0)
        {
            if (cap > nNodes - 1)
            {
                cap = nNodes - 1;
            }
            capacitiesCount.at(cap - 1)++;
        }

    }
    return capacitiesCount;
}
*/


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

size_t ContactNetwork::getAmountOfEdgesToAddSafe()
{
    // TODO some fancy map like lessMap from Lemon
    //std::cout << "check subgraph: " << std::endl;
    //lemon::ListGraph::NodeMap<bool> filter(network, true);
    //lemon::ListGraph::NodeIt nIt(network);
    //filter[nIt] = false;
    //lemon::ListGraph::Node nn = network.addNode();
    //filter[nn] = true;
    //std::cout << filter[nn] << std::endl;
    //network.erase(nn);
    /* for (lemon::ListGraph::NodeIt nIt(network); nIt != lemon::INVALID; ++nIt)
    {
       int nodeUID = nodeIdMap[nIt];
        size_t cap = population.at(nodeUID).getMaxNumberOfContacts() - population.at(nodeUID).getNumberOfContacts();

        if (cap == 0)
        {
            filter[nIt] = false;
        }

       std::cout << filter[nIt] << " ";
    }
    std::cout << std::endl;

    lemon::FilterNodes<lemon::ListGraph> subgr(network, filter);
    std::cout << "nodes subgr: " << lemon::countNodes(subgr) << std::endl;
    lemon::ListGraph::Node nn = network.addNode();
    filter[nn] = true;
    std::cout << "nodes subgr1: " << lemon::countNodes(subgr) << std::endl;*/

    //lemon::ListGraph::EdgeMap<bool> edge_filter(network);
    //lemon::SubGraph<lemon::ListGraph> sg(network, filter, edge_filter);

    //std::cout << "netw. edges: "  << lemon::countEdges(network) << std::endl;
    //std::cout << "subgr. edges: " << lemon::countEdges(subgr) << std::endl;
    //std::cout << "subgr. nodes: " << lemon::countNodes(subgr) << std::endl;
    size_t num = 0;
    if (isContactsLimited)
    {
        using vertex = std::pair<lemon::ListGraph::Node, int>;
        struct CustomCompare
        {
            bool operator()(const vertex& lhs, const vertex& rhs)
            {
                bool result = false;

                if (lhs.second == rhs.second)
                {
                    if (lhs.first != rhs.first)
                    {
                        result = true;
                    }
                }
                else
                {
                    result = lhs.second < rhs.second;
                }

                return result;

            }

        };

        std::set<vertex, CustomCompare> s;
        //for (lemon::FilterNodes<lemon::ListGraph>::NodeIt snIt(subgr); snIt != lemon::INVALID; ++snIt)
        for (lemon::ListGraph::NodeIt snIt(network); snIt != lemon::INVALID; ++snIt)
        {
            int nodeUID = nodeIdMap[snIt];
            //std::cout << "nodeUID =  " << nodeUID << std::endl;
            size_t cap = population.at(nodeUID).getMaxNumberOfContacts() - population.at(nodeUID).getNumberOfContacts();
            //std::cout << cap << " ";
            if (cap > 0)
            {
                s.emplace(std::make_pair(snIt, cap));
            }
        }
        //std::cout << std::endl;
        while (true)
        {
            auto it = s.begin();
            if (it == s.end())
            {
                break;
            }
            int cap = it->second;
            lemon::ListGraph::Node nd = it->first;
            auto it2 = it;
            it2++;
            while (cap > 0 && it2 != s.end())
            {
                if (it2 == s.end())
                {
                    break;
                }
                lemon::ListGraph::Edge ed =  lemon::findEdge(network, nd, it2->first);
                if (ed == lemon::INVALID) // if nodes are not connected
                {
                    cap--;
                    num ++;
                    auto it3 = it2;
                    if (it2->second - 1 > 0)
                    {
                        s.emplace_hint(it2, it2->first, it2->second - 1);
                    }
                    it2++;
                    s.erase(it3);
                }
                else
                {
                    it2++;
                }
            }
            s.erase(it);
        }
    }

    else
    {
        size_t N = population.size();
        num = N * (N - 1) / 2 - countEdges();
    }

    return num;

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
    maxContactsDistribution = std::uniform_int_distribution<size_t>(maxContactsLimitL, maxContactsLimitU);

    //std::uniform_real_distribution<> dis(newContRate, newContRate *100);
    transmissionRate = transmRate;
    newContactRate = newContRate;//dis(generator);//newContRate;
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

double ContactNetwork::getExpectedEdgeAdditionRate(size_t &nAdd)const {
    double sum = 0;
    double sumSquare = 0;
    nAdd = 0;
    for (lemon::ListGraph::EdgeIt eIt(complement); eIt != lemon::INVALID; ++eIt)
    {
        double eRate = getEdgeAdditionRate(eIt);
        sum += eRate;
        sumSquare += eRate * eRate;
        if (eRate > 0)
        {
            nAdd++;
        }
    }

    double result = 0;
    if (sumSquare > 0)
    {
        result = sumSquare / sum;
    }
    return result;
}

double  ContactNetwork::getExpectedEdgeDeletionRate()const
{
    double sum = 0;
    double sumSquare = 0;
    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        double eRate = getEdgeDeletionRate(eIt);
        //std::cout << "eRate: " << eRate << std::endl;
        sum += eRate;
        sumSquare += eRate * eRate;
    }
    double result = 0;
    if (sumSquare > 0)
    {
        result = sumSquare / sum;
    }
    return result;
}

size_t ContactNetwork::updateSurvivalProbability(size_t nDeletions, size_t nAdditions,
        std::vector<BenStructure> &benToFile, double time)
{
    bool isPositive = (nAdditions >= nDeletions);
    double change = 0;
    if (isPositive)
    {
        change = static_cast<double>(nAdditions - nDeletions);
    }
    else
    {
        change = static_cast<double>(nDeletions - nAdditions);
    }

    double nEdges = static_cast<double>(countEdges());
    double xAvg = nEdges * 2;

    if (isPositive)
    {
        xAvg += change / 2;
    }

    else
    {
        xAvg -= change / 2;
    }

    double pSurvival = std::pow(1.0 - 1.0 / xAvg, static_cast<double>(nDeletions));
    std::cout << "pSurvival = " << pSurvival << std::endl;
    std::uniform_real_distribution<> randuni;

    lemon::ListGraph::EdgeIt eIt(network);
    size_t nSurvived = 0;
    while (eIt != lemon::INVALID)
    {
        double r = randuni(generator);
        while (r == 0)
        {
            r = randuni(generator);
        }
        if (r > pSurvival)
        {
            lemon::ListGraph::EdgeIt eTemp = eIt;
            ++eIt;
            std::pair<int, int> bP = removeEdge(eTemp);

            BenStructure b(time, -1, -1, false);
            if (bP.first < bP.second)
            {
                b.u = bP.first;
                b.v = bP.second;
            }
            else
            {
                b.v = bP.first;
                b.u = bP.second;
            }
            benToFile.push_back(b);
        }
        else
        {
            nSurvived++;
            ++eIt;
        }
    }

    return nSurvived;
}


/*bool ContactNetwork::isCapacityLimited()
{
    bool result = false;

    size_t popSize = population.size();
    for (auto& it: population)
    {
        if (it.second.getMaxNumberOfContacts() < popSize - 1)
        {
            result = true;
            break;
        }

    }
    return result;
}*/

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

    // find respective nodes in network
    lemon::ListGraph::Node networkU = network.nodeFromId(a);
    lemon::ListGraph::Node networkV = network.nodeFromId(b);
    lemon::ListGraph::Edge e = lemon::findEdge(network, networkU, networkV);
    return e;

}
