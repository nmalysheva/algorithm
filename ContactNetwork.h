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


struct BenStructure
{
    double t;
    int u;
    int v;
    bool state;

    BenStructure(double _t,int _u, int _v, bool _st);
    bool isValid();

};

class ContactNetwork {

public:

    ContactNetwork(size_t nInfected,
                   size_t nSusceptible,
                   size_t nEdges,
                   int maxContactsL,
                   int MaxContactsU,
                   double transmRate,
                   double newContRate,
                   double looseContRate,
                   double diagnRate,
                   double dRate,
                   double bRate) : network(lemon::ListGraph()),
                                   complement(lemon::ListGraph()),
                                   transmissionRates(network),
                                   diagnosisRates(network),
                                   nodeIdMap(network)
                                   //complementAdjacentEdges(complement)

                                   {
                                       init(nInfected, nSusceptible, nEdges, maxContactsL, MaxContactsU,
                                               transmRate, newContRate, looseContRate, diagnRate, dRate, bRate);
                                   };

    ContactNetwork& operator=(const ContactNetwork& other);

    size_t  size() const; //amount of nodes
    size_t  countByState(Specie::State st) const;  //return amount of I/S/R etc. species in network
    size_t  countEdges() const;


    //sum of rates of particular reactions
    std::vector<std::pair<double, lemon::ListGraph::Edge>> getTransmissionRateSum();
    std::vector<std::pair<double, lemon::ListGraph::Edge>> getEdgeDeletionRateSum()const;
    double  getExpectedEdgeDeletionRate()const;

    std::vector<std::pair<double, lemon::ListGraph::Edge>> getEdgeAdditionRateSum()const;
    double  getExpectedEdgeAdditionRate(size_t &nAdd)const;

    std::vector<std::pair<double, lemon::ListGraph::Node>> getDeathRateSum()const;
    double  getBirthRateSum()const;

    std::vector<std::pair<double, lemon::ListGraph::Node>> getDiagnosisRateSum()const;

    double  getTransmissionRateLimit() const;

    size_t  getMaxContactsLimitOfInfected(double t)const;
    size_t  getMaxContactsLimitOfSusceptible(double t)const;

    size_t  getAmountOfEdgesToAddSafe();
    size_t  getAmountOfEdgesToAdd() const;

    std::pair<int, int> addEdge(lemon::ListGraph::Edge & complementEdge);
    std::pair<int, int> removeEdge(lemon::ListGraph::Edge & edge);
    void removeNode(lemon::ListGraph::Node & node);


    void executeTransmission(lemon::ListGraph::Edge & edge, double time);
    void executeDiagnosis(lemon::ListGraph::Node & node, double time);
    void executeDeath(lemon::ListGraph::Node & node);
    void executeBirth(double rStart, double rBound);

    size_t updateSurvivalProbability(size_t nDeletions, size_t nAdditions, std::vector<BenStructure> &benToFile, double time);

    std::vector<size_t> getDegreeDistribution() const;

    std::vector<BenStructure> getBenStructure(double t); //


    double  getEdgeAdditionRate(const lemon::ListGraph::Edge &complementEdge) const;
    double  getEdgeDeletionRate(const lemon::ListGraph::Edge &networkEdge) const;

    lemon::ListGraph::Edge getComplementEdge(int a, int b);
    lemon::ListGraph::Edge getEdge(int a, int b);


private:

    void initRandomGenerator();
    void initRates(int maxContactsL, int MaxContactsU, double transmRate, double newContRate, double looseContRate,
                   double diagnRate, double dRate, double bRate);
    void initComplementNetwork(size_t nPopulation);


    void init(size_t nInfected, size_t nSusceptible, size_t nEdges, int maxContactsL, int MaxContactsU,
            double transmRate, double newContRate, double looseContRate, double diagnRate,
            double dRate, double bRate);

    size_t countAdjacentEdges(const lemon::ListGraph::Node &complementNode) const;


    lemon::ListGraph network;
    lemon::ListGraph complement;

    lemon::ListGraph::EdgeMap<double> transmissionRates; //TODO: calculate automatically instead of storing?
    lemon::ListGraph::NodeMap<double> diagnosisRates; //TODO: calculate automatically instead of storing?

    lemon::IdMap<lemon::ListGraph, lemon::ListGraph::Node> nodeIdMap;

    // node map for complement graph; stores amount of adjucent edges with addition rate > 0
    //lemon::ListGraph::NodeMap<size_t> complementAdjacentEdges;

    std::unordered_map<int, Specie> population;

    //std::random_device rDev;
    std::mt19937_64 generator;
    
    // TODO: replace all rates with rate distributions to have different rates sampled for different species/edges
    std::uniform_int_distribution<size_t> maxContactsDistribution;
    double transmissionRate;
    double newContactRate;
    double looseContactRate;
    double deathRate;
    double diagnosisRate;

    double birthRate;

    bool isContactsLimited;



    size_t maxContactsLimitL;
    size_t maxContactsLimitU;





};


#endif //ALGO_CONTACTNETWORK_H
