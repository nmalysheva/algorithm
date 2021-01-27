/**
 * Created by Malysheva, Nadezhda on 2019-07-28.
 *
 * Class ContactNetwork describes a time evolving network.
 * Nodes - Species are connected to each other with edges. Graph structure change through time by
 * assembling /disassembling edges.
 * Graph is simple, i.e only one edg between each pair of nodes. Graph is undirected.
 * TODO: may be change to mixed graph? where edges between two species where
 * infection is not possible (for example, S-S, I-I edges) are undirected and ones with transmission rate > 0
 * for instance edge between Infected (I) and Susceprible (S) is directed (as transmission only possible one direction)
 * TODO: estimete pros & contras of this aproach
 * Uses LEMON library for graph representation.
 * Graph is represented by two graphs - actual network & complement network. Complement NW is used for
 * more convenient and direct addition of the edges. TODO:may be store all edges in one graph with different labels/colors?
*/

#ifndef ALGO_CONTACTNETWORK_H
#define ALGO_CONTACTNETWORK_H

#include <unordered_map>
#include <lemon/list_graph.h>
#include <lemon/maps.h>
#include "Specie.h"
#include <random>
#include "types.h"


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
                                   {
                                       init(nInfected, nSusceptible, nEdges, maxContactsL, MaxContactsU,
                                               transmRate, newContRate, looseContRate, diagnRate, dRate, bRate);
                                   };

    ContactNetwork& operator=(const ContactNetwork& other);

    size_t  size() const; //@return amount of nodes
    size_t  countByState(Specie::State st) const;  //@return amount of I/S/R etc. species in network
    size_t  countEdges() const;//@return amount of edges

/* calculates sum & cumsum of rates of particular reactions.
 * Used in SSA & NSA to find reaction being executed
 * @return a vector of pairs <cumulative sum, instance> for instances of question.
 * first element of the vector is always pair <0, INVALID> for convenience
 */
    std::vector<std::pair<double, lemon::ListGraph::Edge>> getTransmissionRateSum() const; //transmission
    std::vector<std::pair<double, lemon::ListGraph::Edge>> getEdgeDeletionRateSum()const;
    std::vector<std::pair<double, lemon::ListGraph::Edge>> getEdgeAdditionRateSum()const;
    std::vector<std::pair<double, lemon::ListGraph::Node>> getDeathRateSum()const;
    std::vector<std::pair<double, lemon::ListGraph::Node>> getDiagnosisRateSum()const;
    //double  getExpectedEdgeDeletionRate()const;
    //double  getExpectedEdgeAdditionRate(size_t &nAdd)const;

    double  getBirthRateSum()const; //TODO:NOT USED, OUTDATED

/*
 * @return highest possible transmission rate in network. Since now rate is static,
 * just return transmission rate given by initiating Contact Network.
 * TODO: for diagn. nodes new infect. rate is lower than for inf. split calc. according to the status.
 */
    double  getTransmissionRateLimit() const;

/*
 * @return  max. number of contacts infected AND DIAGNOSED nodes can have during time t.
 * TODO: for diagn. nodes new contact rate is lower than for inf. split calc. according to the status.
 */
    size_t  getMaxContactsLimitOfInfected(double t)const;
/*
 * @return  max. number of contacts susceptible nodes can have during time t.
 */
    size_t  getMaxContactsLimitOfSusceptible(double t)const;

    //size_t  getAmountOfEdgesToAddSafe();

 /*
 * Adding edge to the network. input - reference to the edge from complement network
 * @return pair of ids of nodes that was connected by given edge
 */
    std::pair<int, int> addEdge(lemon::ListGraph::Edge & complementEdge);

    /*
    * Removing edge from the network. input - reference to the edge from actual network
    * @return pair of ids of nodes that was disconnected
    */

    std::pair<int, int> removeEdge(lemon::ListGraph::Edge & edge);

    /*
    * Removing node from the network. input - reference to the node from actual network
    */
    void removeNode(lemon::ListGraph::Node & node);


    /*
    * Executing particular reactions. TODO: take in account reactions can be defined bu user.
     * Do as a  template?!
    */
    void executeTransmission(lemon::ListGraph::Edge & edge, double time);
    void executeDiagnosis(lemon::ListGraph::Node & node, double time);
    void executeDeath(lemon::ListGraph::Node & node);
    void executeBirth(double rStart, double rBound);

    //size_t updateSurvivalProbability(size_t nDeletions, size_t nAdditions, std::vector<BenStructure> &benToFile, double time);

    /*
     * Gets degree distribution of the network.
     * @return vector of the degrees of each node in the network
    */
    std::vector<size_t> getDegreeDistribution() const;
    std::vector<specieState> getNetworkState() const;

    std::vector<BenStructure> getBenStructure(double t); //


    double  getEdgeAdditionRate(const lemon::ListGraph::Edge &complementEdge) const;
    double  getEdgeDeletionRate(const lemon::ListGraph::Edge &networkEdge) const;

    lemon::ListGraph::Edge getComplementEdge(int a, int b); //@return complement edge by given nodes ids
    lemon::ListGraph::Edge getEdge(int a, int b);//@return edge of actua network by given nodes ids


private:

    void initRandomGenerator();
    void initRates(int maxContactsL, int MaxContactsU, double transmRate, double newContRate, double looseContRate,
                   double diagnRate, double dRate, double bRate);
    void initComplementNetwork(size_t nPopulation);


    void init(size_t nInfected, size_t nSusceptible, size_t nEdges, int maxContactsL, int MaxContactsU,
            double transmRate, double newContRate, double looseContRate, double diagnRate,
            double dRate, double bRate);

    //size_t countAdjacentEdges(const lemon::ListGraph::Node &complementNode) const;


    lemon::ListGraph network;
    lemon::ListGraph complement;

    lemon::ListGraph::EdgeMap<double> transmissionRates; //TODO: calculate automatically instead of storing or not?
    lemon::ListGraph::NodeMap<double> diagnosisRates; //TODO: calculate automatically instead of storing or not?

    lemon::IdMap<lemon::ListGraph, lemon::ListGraph::Node> nodeIdMap;

    std::unordered_map<int, Specie> population;

    //std::random_device rDev;
    std::mt19937_64 generator;
    
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
