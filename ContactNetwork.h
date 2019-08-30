//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#ifndef ALGO_CONTACTNETWORK_H
#define ALGO_CONTACTNETWORK_H

#include <unordered_map>
#include <lemon/list_graph.h>
#include "Specie.h"
#include <random>


class ContactNetwork {

public:
    /*ContactNetwork(): network(lemon::ListGraph()),
                      transmissionRates(lemon::ListGraph::EdgeMap<double>(network)),
                      nodeUIDs(lemon::ListGraph::NodeMap<double>(network)) {};*/


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
                                   transmissionRates(lemon::ListGraph::EdgeMap<double>(network)),
                                   nodeUIDs(lemon::ListGraph::NodeMap<double>(network))
                                   {
                                       init(nInfected, nSusceptible, nEdges, maxContactsA, MaxContactsB, transmRate, newContRate, looseContRate, dRate, bRate);
                                   };

    ContactNetwork& operator=(const ContactNetwork& other);

    ContactNetwork (const ContactNetwork& other):  network(lemon::ListGraph()),
                      transmissionRates(lemon::ListGraph::EdgeMap<double>(network)),
                      nodeUIDs(lemon::ListGraph::NodeMap<double>(network)) {copyNetwork(other);};

    size_t  size() const; //amount of nodes
    size_t  countByState(Specie::State st) const;  //return amount of infected species in network
    size_t  countEdges() const;
    //size_t const getNumberOfSusceptible() const; //return amount of susceptible species in network
    //size_t const getNumberOfRecovered() const;  //return amount of recovered species in network


    //sum of rates of particular reactions
    double  getTransmissionRateSum()const; //return sum of all transmition rates
    double  getEdgeDeletionRateSum()const;
    double  getEdgeAdditionRateSum()const;
    double  getDeathRateSum()const;
    double  getBirthRateSum()const;
    double getTransmissionRateLimit() const;
    size_t getMaxContactsLimitOfInfected()const;

    void addEdge(lemon::ListGraph::Edge & complementEdge, double trRate);


    void removeEdge(lemon::ListGraph::Edge & edge);

    void removeNode(lemon::ListGraph::Node &node);

    void executeEdgeDelition(double rStart, double rBound);
    void executeEdgeAddition(double rStart, double rBound);
    void executeTransmission(double rStart, double rBound, double time);
    void executeDeath(double rStart, double rBound);
    void executeBirth(double rStart, double rBound);

    std::vector<size_t> getDegreeDistribution();


private:

    double  getEdgeAdditionRate(lemon::ListGraph::Edge complementEdge) const;
    double  getEdgeDelitionRate(lemon::ListGraph::EdgeIt networkEdgeIt) const;
    void init(size_t nInfected, size_t nSusceptible, size_t nEdges, int maxContactsA, int MaxContactsB,
            double transmRate, double newContRate, double looseContRate, double dRate, double bRate);

    void copyNetwork (const ContactNetwork& other);


    lemon::ListGraph network;
    lemon::ListGraph complement;

    lemon::ListGraph::EdgeMap<double> transmissionRates;
    lemon::ListGraph::NodeMap<double> nodeUIDs;

    std::unordered_map<uint32_t, Specie> population;

    //std::random_device rDev;
    std::mt19937_64 generator;


    //std::exponential_distribution<double> transmitDistribution;
    double transmissionRate;
    std::uniform_int_distribution<unsigned char> maxContactsDistribution;
    //std::exponential_distribution<double> newContactRateDistribution;
    double newContactRate;
    //std::exponential_distribution<double> looseContactRateDistribution;
    double looseContactRate;
    //std::exponential_distribution<double> deathRateRateDistribution;
    double deathRate;

    double birthRate;



};


#endif //ALGO_CONTACTNETWORK_H
