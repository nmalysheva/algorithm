//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#ifndef ALGO_CONTACTNETWORK_H
#define ALGO_CONTACTNETWORK_H

#include <unordered_map>
#include <lemon/list_graph.h>
#include <lemon/maps.h>
#include "Specie.h"
#include <random>


class ContactNetwork {

public:

    ContactNetwork(size_t nInfected,
                   size_t nSusceptible,
                   size_t nEdges,
                   int maxContactsA,
                   int MaxContactsB,
                   double transmRate,
                   double newContRate,
                   double looseContRate,
                   double dRate,
                   double bRate) : network(lemon::ListGraph()),
                                   complement(lemon::ListGraph()),
                                   transmissionRates(network),
                                   nodeIdMap(network),
                                   edgesIdMap (network),
                                   edgesIdMapInv(edgesIdMap),
                                   nEdgesPossibleAdd(complement, 10)

                                   {
                                       init(nInfected, nSusceptible, nEdges, maxContactsA, MaxContactsB,
                                               transmRate, newContRate, looseContRate, dRate, bRate);
                                   };

    ContactNetwork& operator=(const ContactNetwork& other);

    /*ContactNetwork (const ContactNetwork& other):  network(lemon::ListGraph()),
                      transmissionRates(lemon::ListGraph::EdgeMap<double>(network)),
                      //nEdgesComplement(lemon::ListGraph::NodeMap<size_t>(complement)),
                      nodeIdMap(lemon::IdMap<lemon::ListGraph, lemon::ListGraph::Node>(network)) ,
                      //edgesIdMap (lemon::RangeIdMap<lemon::ListGraph, lemon::ListGraph::Edge>(network)){copyNetwork(other);};
                      edgesIdMap (network),
                      edgesIdMapInv(edgesIdMap)
                      {
                          copyNetwork(other);
                      };*/

    size_t  size() const; //amount of nodes
    size_t  countByState(Specie::State st) const;  //return amount of infected species in network
    size_t  countEdges() const;


    //sum of rates of particular reactions
    double  getTransmissionRateSum()const; //return sum of all transmition rates
    double  getEdgeDeletionRateSum(size_t &nDel)const;
    double  getEdgeAdditionRateSum(size_t &nAdd)const;
    double  getDeathRateSum()const;
    double  getBirthRateSum()const;
    double getTransmissionRateLimit() const;
    size_t getMaxContactsLimitOfInfected()const;
    size_t getMaxContactsLimitOfInfected(double t)const;
    size_t getMaxContactsLimitOfSusceptible(double t)const;
    size_t getNumberContactsOfInfected()const;
    size_t getAmountOfEdgesToAddSafe();
    size_t getAmountOfEdgesToAdd() const;

    void addEdge(lemon::ListGraph::Edge & complementEdge/*, double trRate*/);


    void removeEdge(lemon::ListGraph::Edge & edge);

    void removeNode(lemon::ListGraph::Node &node);

    void executeEdgeDeletion(double rStart, double rBound);
    void executeEdgeDeletion(size_t edgeNumber,  size_t maxEdgesToDelete);

    void executeEdgeAddition(double rStart, double rBound);
    void executeEdgeAddition(size_t edgeNumber,  size_t maxEdgesToAdd);

    void executeTransmission(double rStart, double rBound, double time);
    void executeDeath(double rStart, double rBound);
    void executeBirth(double rStart, double rBound);

    std::vector<size_t> getDegreeDistribution();
    std::vector<size_t> countCapacities() const;



    size_t minNumOfEdges;
    //!!!
    size_t subgraph();


private:

    double  getEdgeAdditionRate(lemon::ListGraph::Edge complementEdge) const;
    double  getEdgeDeletionRate(lemon::ListGraph::EdgeIt networkEdgeIt) const;
    void init(size_t nInfected, size_t nSusceptible, size_t nEdges, int maxContactsA, int MaxContactsB,
            double transmRate, double newContRate, double looseContRate, double dRate, double bRate);

    void copyNetwork (const ContactNetwork& other);


    lemon::ListGraph network;
    lemon::ListGraph complement;

    lemon::ListGraph::EdgeMap<double> transmissionRates; //TODO: calculate automatically instead of storing?

    lemon::IdMap<lemon::ListGraph, lemon::ListGraph::Node> nodeIdMap;

    lemon::RangeIdMap<lemon::ListGraph, lemon::ListGraph::Edge>::InverseMap edgesIdMapInv;
    lemon::RangeIdMap<lemon::ListGraph, lemon::ListGraph::Edge> edgesIdMap;


    // node map for complement graph; stores amount of adjucent edges with addition rate > 0
    lemon::ListGraph::NodeMap<size_t> nEdgesPossibleAdd;

    std::unordered_map<int, Specie> population;
    //std::random_device rDev;
    std::mt19937_64 generator;


    //std::exponential_distribution<double> transmitDistribution;
    double transmissionRate;
    std::uniform_int_distribution<size_t> maxContactsDistribution;
    //std::exponential_distribution<double> newContactRateDistribution;
    double newContactRate;
    //std::exponential_distribution<double> looseContactRateDistribution;
    double looseContactRate;
    //std::exponential_distribution<double> deathRateRateDistribution;
    double deathRate;

    double birthRate;



    size_t maxContactsLimitA;
    size_t maxContactsLimitB;





};


#endif //ALGO_CONTACTNETWORK_H
