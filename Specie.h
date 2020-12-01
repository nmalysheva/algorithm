/***
 * Created by Malysheva, Nadezhda on 2019-07-28.
 * Class describes a specie (aka node in network).
***/

#ifndef ALGO_SPECIE_H
#define ALGO_SPECIE_H


#include <cstddef>

class Specie {

/**
 * Epidemic state.
 * Depending on the model, can have variety of states,
 * for example in SIR model - states are "S" - susceptible,
 * "I" - infected, "R" - recovered.
 * TODO: state set can be defined by user, i.e. free set of statuses, not necessary described as above.
 * TODO: possibility to set a known models state sets: SIR, SI, etc.
 */
public:
    enum State {S, I, R, D}; //Susceptible, Infected, Recovered, Diagnosed

public:

/**
 * Default constructor
 * Create a specie with status "S" by default and all other params set to "0"
 */
    Specie();

    Specie(size_t maxNumberOfContacts,
           size_t numberOfContacts,
           double deathRate,
           double newContactRate,
           double looseContactRate,
           State st = S);

    size_t  getMaxNumberOfContacts() const; //@return max. number of contacts specie "allowed to have"
    size_t  getNumberOfContacts() const; //@return current number of contacts specie has

    double getNewContactRate() const;   //@return rate of establishing a new contact
    double getLooseContactRate() const; //@return rate of loosing already existing contacts

    double  getDeathRate() const;  //@return specie's death rate

/**
 * @return time of the last state change.
 * for instance, if infection occurred and state changed from "S" to "I",
 * this function will return the tome of infection.
 */
    double  getLastStateChangeTime() const;

    State  getState() const; // return epidemic state of the specie

/**
 * calculates analytical solution for max. number of contacts specie can have after time t,
 * given current specie parameters - rates of establishing and loosing contacts & number of current contacts.
 * Considering specie has:
 * a - rate of establish a contact
 * b - rate of loose a contact
 * Cmax - maximum number of contacts
 * c - current number of contacts.
 *
 * then analytical ODE for number of contacts of separated specie w/o external factors:
 * @see PAPER
 *
 */
    double getNumberOfContactsLimit(double t) const;

    //setters

    void setMaxNumberOfContacts(unsigned int maxNumOfCont);

    void incNumberOfContacts(); // increase number of current contacts to 1
    void decNumberOfContacts(); // decrease number of current contacts to 1


    void setNewContactRate  (double newContRate);
    void setLooseContactRate(double looseContRate);

    void setDeathRate(double dRate);


    //status changes
    void changeState (State st, double time);  //change status of the specie



    //comparison operator for map
    bool operator== (const Specie &sp) const;

    //destructor
    ~Specie() {};


private:

    void setNumberOfContacts(size_t nOfCont);

/**
 * calculates analytical expectation  for number of contacts specie can have after time t,
 * given current specie parameters - rates of establishing and loosing contacts & number of current contacts.
 * Considering specie has:
 * a - rate of establish a contact
 * b - rate of loose a contact
 * Cmax - maximum number of contacts
 * c - current number of contacts.
 *
 * then analytical ODE for expectation  of contacts of separated specie w/o external factors:
 * @see PAPER
 * @return expectation
 *
 */
    double ExpectationOfContacts(double a, double b, double t)const;

/**
 * calculates analytical expectation  for number of contacts specie can have after time t,
 * given current specie parameters - rates of establishing and loosing contacts & number of current contacts.
 * Considering specie has:
 * a - rate of establish a contact
 * b - rate of loose a contact
 * Cmax - maximum number of contacts
 * c - current number of contacts.
 *
 * then analytical ODE for expectation  of contacts of separated specie w/o external factors:
 * @see PAPER
 * @return variance
 */
    double VarianceOfContacts(double a, double b, double t)const;

private:

    size_t  maxNumberOfContacts; //max. number of contacts specie "allowed to have"
    size_t  numberOfContacts; //current number of contacts
    double deathRate;         // death rate of the specie
    double newContactRate;  //rate of establishing a new contact, i.e. new edge
    double looseContactRate;//rate to loose a contact, i.e. an edge

    double stateChangeTime;  //last time of state change (i.e. S->I, I->R, etc.)


    State state;

};


#endif //ALGO_SPECIE_H
