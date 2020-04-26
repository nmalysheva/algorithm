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

double  ContactNetwork::getEdgeDeletionRateSum(size_t &nDel)const
{
    double result = 0;
    nDel = 0;
    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        result += getEdgeDeletionRate(eIt);
        nDel++;
    }
    return result;
}

double  ContactNetwork::getEdgeAdditionRateSum(size_t &nAdd)const
{
    double result = 0;
    nAdd = 0;
    for (lemon::ListGraph::EdgeIt eIt(complement); eIt != lemon::INVALID; ++eIt)
    {
        result += getEdgeAdditionRate(eIt);
        if (getEdgeAdditionRate(eIt) > 0)
        {
            nAdd++;
        }

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

void ContactNetwork::addEdge(lemon::ListGraph::Edge &complementEdge /*, double trRate*/)
{
    //nodes of the given edge in a complement graph
    lemon::ListGraph::Node complU = complement.u(complementEdge);
    lemon::ListGraph::Node complV = complement.v(complementEdge);

    // find respective nodes in network
    lemon::ListGraph::Node networkU = network.nodeFromId(complement.id(complU));
    lemon::ListGraph::Node networkV = network.nodeFromId(complement.id(complV));

    //erase given edge from complement graph
    complement.erase(complementEdge);
    nEdgesPossibleAdd[complU] = 0;
    for(lemon::ListGraph::IncEdgeIt eIt(complement, complU); eIt!=lemon::INVALID; ++eIt)
    {
        lemon::ListGraph::Node oppositeNode = complement.oppositeNode(complU, eIt);
        int oppositeUID = nodeIdMap[network.nodeFromId(complement.id(oppositeNode))];
        double oppositeRate = population.at(oppositeUID).getNewContactRate();
        if (oppositeRate > 0)
        {
            nEdgesPossibleAdd[complU] = nEdgesPossibleAdd[complU] + 1;
        }
    }
    nEdgesPossibleAdd[complV] = 0;
    for(lemon::ListGraph::IncEdgeIt eIt(complement, complV); eIt!=lemon::INVALID; ++eIt)
    {
        lemon::ListGraph::Node oppositeNode = complement.oppositeNode(complV, eIt);
        int oppositeUID = nodeIdMap[network.nodeFromId(complement.id(oppositeNode))];
        double oppositeRate = population.at(oppositeUID).getNewContactRate();
        if (oppositeRate > 0)
        {
            nEdgesPossibleAdd[complV] = nEdgesPossibleAdd[complV] + 1;
        }
    }


    //add new edge to the network with given transmission rate
    lemon::ListGraph::Edge newEdge = network.addEdge(networkU, networkV);

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

        //tRate = transmitDistribution(generator);
        trRate = transmissionRate;
    }

    transmissionRates[newEdge] = trRate;


}

void ContactNetwork::removeEdge(lemon::ListGraph::Edge &edge)
{
    //lemon::ListGraph::Edge e = network.edgeFromId(edgeID);

    //nodes of the given edge in a graph
    lemon::ListGraph::Node networkU = network.u(edge);
    lemon::ListGraph::Node networkV = network.v(edge);

    int sourceUID = network.id(networkU);
    int targetUID = network.id(networkV);

    network.erase(edge); // erase from network

    // find respective nodes in complement graph
    lemon::ListGraph::Node complU = complement.nodeFromId(sourceUID);
    lemon::ListGraph::Node complV = complement.nodeFromId(targetUID);
    complement.addEdge(complU, complV);   // add to complement graph
    nEdgesPossibleAdd[complU] = 0;
    for(lemon::ListGraph::IncEdgeIt eIt(complement, complU); eIt!=lemon::INVALID; ++eIt)
    {
        lemon::ListGraph::Node oppositeNode = complement.oppositeNode(complU, eIt);
        int oppositeUID = nodeIdMap[network.nodeFromId(complement.id(oppositeNode))];
        double oppositeRate = population.at(oppositeUID).getNewContactRate();
        if (oppositeRate > 0)
        {
            nEdgesPossibleAdd[complU] = nEdgesPossibleAdd[complU] + 1;
        }
    }
    nEdgesPossibleAdd[complV] = 0;
    for(lemon::ListGraph::IncEdgeIt eIt(complement, complV); eIt!=lemon::INVALID; ++eIt)
    {
        lemon::ListGraph::Node oppositeNode = complement.oppositeNode(complV, eIt);
        int oppositeUID = nodeIdMap[network.nodeFromId(complement.id(oppositeNode))];
        double oppositeRate = population.at(oppositeUID).getNewContactRate();
        if (oppositeRate > 0)
        {
            nEdgesPossibleAdd[complV] = nEdgesPossibleAdd[complV] + 1;
        }
    }

    //decrtease number of contacts of the disconnected nodes
    population.at(sourceUID).decNumberOfContacts();
    population.at(targetUID).decNumberOfContacts();
}

void ContactNetwork::removeNode(lemon::ListGraph::Node &node)
{
    //delete all edges to given node in network
    lemon::ListGraph::IncEdgeIt ieIt(network, node);
    while (ieIt != lemon::INVALID)
    {
        lemon::ListGraph::IncEdgeIt tmpIt = ieIt;
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

void ContactNetwork::executeEdgeDeletion(double rStart, double rBound)
{
    double result = rStart;
    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        result += getEdgeDeletionRate(eIt);
        //result += disassembleRates[eIt];
        if (result >= rBound)
        {
            removeEdge(eIt);
            break;
        }
    }
}

void ContactNetwork::executeEdgeDeletion(size_t edgeNumber, size_t maxEdgesToDelete)
{

    lemon::ListGraph::Edge eIt = edgesIdMapInv[edgeNumber];

    if (! network.valid(eIt))
    {
        //TODO Throw exception
    }
    else
    {
        removeEdge(eIt);
    }
}


void ContactNetwork::executeEdgeAddition(double rStart, double rBound)
{
    double result = rStart;

    for (lemon::ListGraph::EdgeIt eIt(complement); eIt != lemon::INVALID; ++eIt)
    {
        result += getEdgeAdditionRate(eIt);

        if (result >= rBound)
        {
            addEdge(eIt);
            break;
        }

    }
}

void ContactNetwork::executeEdgeAddition(size_t edgeNumber,  size_t maxEdgesToAdd)
{
    //std::cout << "edge Num: " << edgeNumber << "; ";
    //std::cout << "max edges to add: " << getAmountOfEdgesToAdd() << std::endl;

    size_t currentEdge = 0;

    lemon::ListGraph::EdgeIt edgeToAddIt(lemon::INVALID);

    for(lemon::ListGraph::EdgeIt eIt(complement); eIt != lemon::INVALID; ++eIt)
    {

        double edgeAdditionRate = getEdgeAdditionRate(eIt);
        //std::cout  << edgeAdditionRate << ", ";
        if (edgeAdditionRate > 0)
        {
            //std::cout  << currentEdge << ", ";
            if (currentEdge == edgeNumber)
            {
               // std::cout << "currentEdge: " << currentEdge << std::endl;
                edgeToAddIt = eIt;
                break;

            }
            currentEdge ++;
        }
    }


    if (! complement.valid(edgeToAddIt))
    {
        std::cout << "INVALID" << std::endl;
    }
    else
    {
        addEdge(edgeToAddIt);
    }
}


/*void ContactNetwork::executeEdgeDeletionUniform()
{
    int maxEdgeId = network.maxEdgeId();

    //std::cout <<"maxEdgeId " << maxEdgeId << " " << lemon::countEdges(network) <<std::endl;
    std::uniform_int_distribution<int> dist(0, maxEdgeId);


    int edgeId = dist(generator);

    while (! network.valid(network.edgeFromId(edgeId)))
    {
        edgeId = dist(generator);
    }
    lemon::ListGraph::Edge e = network.edgeFromId(edgeId);
    removeEdge(e);


}*/

/*void ContactNetwork::executeEdgeAdditionUniform()
{
    //int maxEdgeId = complement.maxEdgeId();
    int maxEdges = lemon::countEdges(complement);
    std::uniform_int_distribution<int> dist(0, maxEdges - 1);
    double edgeAddRate = 0;
    int edgeId = dist(generator);

    /*lemon::ListGraph::EdgeIt eIt(complement);
    eIt = eIt + edgeId;

    //std::cout<<"trate " << tRate << std::endl;
    addEdge(e, tRate);
}*/

void ContactNetwork::executeTransmission(double rStart, double rBound, double time)
{
    double result = rStart;

    for (lemon::ListGraph::EdgeIt eIt(network); eIt != lemon::INVALID; ++eIt)
    {
        result += transmissionRates[eIt];
        //std::cout << transmissionRates[eIt] << std::endl;
        if (result >= rBound)
        {
            lemon::ListGraph::Node networkU = network.u(eIt);
            lemon::ListGraph::Node networkV = network.v(eIt);

            int sourceUID = nodeIdMap[networkU];
            int targetUID = nodeIdMap[networkV];

            lemon::ListGraph::Node infectedNode;
            if (population.at(sourceUID).getState()  == Specie::S)
            {
                population.at(sourceUID).infect(time);
                infectedNode = networkU;

            }
            else if (population.at(targetUID).getState()  == Specie::S)
            {
                population.at(targetUID).infect(time);
                infectedNode = networkV;
            }

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
            break;

        }

    }
}

void ContactNetwork::executeDeath(double rStart, double rBound)
{
    double result = rStart;
    for(lemon::ListGraph::NodeIt nIt(network); nIt!=lemon::INVALID; ++nIt)
    {
        uint32_t nodeUID = nodeIdMap[nIt];
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
        size_t maxContacts = maxContactsDistribution(generator);
        double  dRate = deathRate;//deathRateRateDistribution(generator);
        double newContRate = newContactRate;//newContactRateDistribution(generator);
        double looseContRate = looseContactRate;//looseContactRateDistribution(generator);
        Specie sp = Specie(maxContacts, 0, dRate, newContRate, looseContRate);

        lemon::ListGraph::Node newNode = network.addNode();
        population.emplace(nodeIdMap[newNode], sp);

        newNode = complement.addNode();

        //break;
    }
}

double  ContactNetwork::getEdgeAdditionRate(lemon::ListGraph::Edge complementEdge) const
{
    double result = 0;

    lemon::ListGraph::Node complU = complement.u(complementEdge);
    lemon::ListGraph::Node complV = complement.v(complementEdge);

    int sourceUID = nodeIdMap[network.nodeFromId(complement.id(complU))];
    double sourceRate = population.at(sourceUID).getNewContactRate();

    int targetUID = nodeIdMap[network.nodeFromId(complement.id(complV))];
    double targetRate = population.at(targetUID).getNewContactRate();
    if ( (sourceUID > 0) && (targetRate > 0) )
    {
        double sourceAddEdges = nEdgesPossibleAdd[complU];


        double targetAddEdges = nEdgesPossibleAdd[complV];

        //std::cout << "sourceAddEdges = " << sourceAddEdges << "; targetAddEdges = " << targetAddEdges << std::endl;
        sourceRate = sourceRate / sourceAddEdges;
        targetRate = targetRate / targetAddEdges;

        result = sourceRate * targetRate;
    }

    return result;
}

double  ContactNetwork::getEdgeDeletionRate(lemon::ListGraph::EdgeIt networkEdgeIt) const
{
    lemon::ListGraph::Node networkU = network.u(networkEdgeIt);
    lemon::ListGraph::Node networkV = network.v(networkEdgeIt);

    int sourceUID = nodeIdMap[networkU];
    int targetUID = nodeIdMap[networkV];

    double sourceRate = population.at(sourceUID).getLooseContactRate();
    double targetRate = population.at(targetUID).getLooseContactRate();

    //double result = (population.at(sourceUID).getLooseContactRate() + population.at(targetUID).getLooseContactRate()) / 2;
    double result = sourceRate * targetRate;
    return result;

}

void ContactNetwork::init(size_t nInfected, size_t nSusceptible, size_t nEdges, int maxContactsA, int MaxContactsB,
                          double transmRate, double newContRate, double looseContRate, double dRate, double bRate)
{
    std::random_device rDev;
    generator = std::mt19937_64(rDev());
    //generator.seed(::time(NULL)); //to change the seed for every run
    generator.seed(3);

    maxContactsLimitA = maxContactsA;
    maxContactsLimitB = MaxContactsB;
    //transmitDistribution = std::exponential_distribution<double> (8);
    maxContactsDistribution = std::uniform_int_distribution<size_t>(maxContactsA, MaxContactsB);
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

    size_t nInf = 0;
    std::cout << "maxContacts: " ;
    size_t sumTmp = 0;
    for (size_t i = 0; i < nPopulation; i ++)
    {
        size_t maxContacts = maxContactsDistribution(generator);
        std::cout << maxContacts << " ";
        double  dRate = deathRate;//deathRateRateDistribution(generator);
        double newContRate = newContactRate;//newContactRateDistribution(generator);
        double looseContRate = looseContactRate;//looseContactRateDistribution(generator);

        Specie::State st = Specie::S;
        Specie sp = Specie(maxContacts, 0, dRate, newContRate, looseContRate, st);

        if (nInf < nInfected)
        {
            sp.infect(0);
            nInf ++;

        }

        lemon::ListGraph::Node newNode = network.addNode();
        population.emplace(nodeIdMap[newNode], sp);

    }
    std::cout <<std::endl;

    //lemon::RangeIdMap<lemon::ListGraph, lemon::ListGraph::Edge> em (network);
    //lemon::RangeIdMap<lemon::ListGraph, lemon::ListGraph::Edge>::InverseMap inmap(lemon::RangeIdMap<lemon::ListGraph, lemon::ListGraph::Edge> (network));
    //std::cout << "ttt check:" << sizeof(ttt) << std::endl;
    //minNumOfEdges = getAmountOfEdgesToAddSafe();
    minNumOfEdges = subgraph();
    std::cout << "minNumOfEdges: " << minNumOfEdges <<std::endl;
    //int maxNumberOfEdges = lemon::countEdges(complement);
    int maxNumberOfEdges = complement.maxEdgeId();
    // TODO: change edge adding with regarding to the amount of edges possible to add safely.
    for (size_t i = 0; i < nEdges; i ++)
    {
        //std::cout << maxNumberOfEdges - 1 <<std::endl;
        std::uniform_int_distribution<int> dist(0, maxNumberOfEdges);
        int edgeId = dist(generator);
        //std::cout << edgeId <<std::endl;

        lemon::ListGraph::Edge cEdge = complement.edgeFromId(edgeId);
        if (complement.valid(cEdge) &&  getEdgeAdditionRate(cEdge) > 0)
        {
            addEdge(cEdge);
        }
    }


    std::cout <<lemon::countNodes(complement) <<std::endl;

    for(lemon::ListGraph::NodeIt nIt(complement); nIt!=lemon::INVALID; ++nIt)
    {
        nEdgesPossibleAdd[nIt] = 0;
        for(lemon::ListGraph::IncEdgeIt eIt(complement, nIt); eIt!=lemon::INVALID; ++eIt)
        {
            lemon::ListGraph::Node oppositeNode = complement.oppositeNode(nIt, eIt);
            int oppositeUID = nodeIdMap[network.nodeFromId(complement.id(oppositeNode))];
            double oppositeRate = population.at(oppositeUID).getNewContactRate();
            if (oppositeRate > 0)
            {
                nEdgesPossibleAdd[nIt] = nEdgesPossibleAdd[nIt] + 1;
            }
        }
    }

    for(lemon::ListGraph::NodeIt nIt(complement); nIt!=lemon::INVALID; ++nIt)
    {
        std::cout << nEdgesPossibleAdd[nIt] <<" ";
    }

    std::cout <<std::endl;
    //std::cout << "edgesIdMap check:" << network.valid(edgesIdMapInv[1]) << std::endl;
    generator.seed(::time(NULL) * getpid());
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

       //cg.nodeMap(other.nodeIdMap, this->nodeIdMap);
       cg.run();

       lemon::GraphCopy<lemon::ListGraph, lemon::ListGraph> cg2(other.complement, this->complement);
       cg2.nodeMap(other.nEdgesPossibleAdd, this->nEdgesPossibleAdd);
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
    }
    return *this;
}

/*void ContactNetwork::copyNetwork (const ContactNetwork& other)
{
    *this = other;

}*/

std::vector<size_t> ContactNetwork::getDegreeDistribution()
{
    std::vector<size_t> result; //TODO: reserve place in a begining

    //size_t sumTmp  = 0;

    for(lemon::ListGraph::NodeIt nIt(network); nIt!=lemon::INVALID; ++nIt)
    {
        size_t degree = 0 ;

        for(lemon::ListGraph::IncEdgeIt e(network, nIt); e!=lemon::INVALID; ++e)
        {
            ++ degree;
        }

        //sumTmp += degree;
        result.push_back(degree);
    }

    //std::cout << "edges: " << sumTmp / 2 << std::endl;
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
    size_t result = 0;

    for (lemon::ListGraph::EdgeIt eIt(complement); eIt != lemon::INVALID; ++eIt)
    {
        if (getEdgeAdditionRate(eIt) > 0)
        {
            result++;
        }
    }

    return result;
}

std::vector<size_t> ContactNetwork::countCapacities() const
{
    std::vector<size_t> capacitiesCount (maxContactsLimitB, 0);

    size_t  nNodes = lemon::countNodes(network); //TODO may be nNodes can be provided to the function already and there is no need to calc. it again
    for (lemon::ListGraph::NodeIt nIt(network); nIt != lemon::INVALID; ++nIt)
    {
        int nodeUID = nodeIdMap[nIt];
        size_t cap = population.at(nodeUID).getMaxNumberOfContacts() - population.at(nodeUID).getNumberOfContacts();
        if (cap > maxContactsLimitB)
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


size_t ContactNetwork::getNumberContactsOfInfected()const
{
    size_t result = 0;

    for (auto& it: population)
    {
        if (it.second.getState() == Specie::State::I)
        {
            result += it.second.getNumberOfContacts();
        }
    }
    return result;
}

size_t ContactNetwork::getMaxContactsLimitOfInfected(double t)const
{
    size_t result = 0;

    for (auto& it: population)
    {
        if (it.second.getState() == Specie::State::I)
        {
            result += it.second.getNumberOfContactsLimit(t);
            //std::cout << "max cont: " << it.second.getNumberOfContactsLimit(t) << std::endl;
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

size_t ContactNetwork::subgraph()
{
    // TODO some fancy map like lessMap from Lemon
    lemon::ListGraph::NodeMap<bool> filter(network, true);
    /*for (lemon::ListGraph::NodeIt nIt(network); nIt != lemon::INVALID; ++nIt)
    {
        int nodeUID = nodeIdMap[nIt];
        size_t cap = population.at(nodeUID).getMaxNumberOfContacts() - population.at(nodeUID).getNumberOfContacts();

        if (cap == 0)
        {
            filter[nIt] = false;
        }

    }*/

    //lemon::FilterNodes<lemon::ListGraph> subgr(network, filter);

    //lemon::ListGraph::EdgeMap<bool> edge_filter(network);
    //lemon::SubGraph<lemon::ListGraph> sg(network, filter, edge_filter);

    //std::cout << "netw. edges: "  << lemon::countEdges(network) << std::endl;
    //std::cout << "subgr. edges: " << lemon::countEdges(subgr) << std::endl;
    //std::cout << "subgr. nodes: " << lemon::countNodes(subgr) << std::endl;

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

    std::set<vertex, CustomCompare> s; //TODO WHY MULTISET? Set  doesn't work
    //std::cout << "cap: ";

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

    size_t num = 0;
    while (true)
    {
        auto it = s.begin();
        if (it == s.end())
        {
            //std::cout << "nummm: " << num << std::endl;
            break;
        }
        int cap = it->second;
        lemon::ListGraph::Node nd = it->first;
        auto it2 = it;
        it2++;
        while (cap > 0 && it2 != s.end())
        {
            //++it2;
            if (it2 == s.end())
            {
                break;
            }
            //lemon::ListGraph::Edge ed =  lemon::findEdge(subgr, nd, it2->first);
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
                //std::cout <<"connected"<<std::endl;
                it2++;
            }
        }
        s.erase(it);
    }
    //std::cout <<"lol"<<std::endl;
    return num;



}