//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#ifndef ALGO_SPECIE_H
#define ALGO_SPECIE_H


#include <cstddef>

class Specie {
public:
    enum State {S, I, R, D}; //Susceptible, Infected, Recovered, Diagnosed

public:

    //constructors

    Specie();  //default constructor
    Specie(//unsigned char age,
           size_t maxNumberOfContacts,
           size_t numberOfContacts,
           double deathRate,
           double newContactRate,
           double looseContactRate,
           State st = S);

    //unsigned char const getAge() const; //return age of the specie
    size_t  getMaxNumberOfContacts() const; //return max. number of contacts per any time instance
    size_t  getNumberOfContacts() const; //return current number of contacts per any time instance

    double getNewContactRate() const;   //return rate for establishing a new contact
    double getLooseContactRate() const; //return rate for loosing of the one of already existing contacts

    double  getDeathRate() const;  //return specie's death rate

    double  getLastStateChangeTime() const; //return last infection time if specie is infected,  -infinity if not

    State  getState() const; // return epidemic state of the specie

    double getNumberOfContactsLimit(double t) const;

    //setters
    //void setAge(unsigned char age);
    void setMaxNumberOfContacts(unsigned int maxNumOfCont);

    void incNumberOfContacts();
    void decNumberOfContacts();


    void setNewContactRate  (double newContRate);
    void setLooseContactRate(double looseContRate);

    void setDeathRate(double dRate);


    //status changes

    void changeState (State st, double time);  //change status of the specie to INFECTED



    //comparison operator for map
    bool operator== (const Specie &sp) const;

    //destructor
    ~Specie() {};


private:

    void setNumberOfContacts(size_t nOfCont);

    double ExpectationOfContacts(double a, double b, double t)const;
    double VarianceOfContacts(double a, double b, double t)const;

private:

    size_t  maxNumberOfContacts; //max. number of contacts at any timepoint
    size_t  numberOfContacts; //current number of contacts
    double deathRate;
    double newContactRate;  //rate to establish a new contact
    double looseContactRate;//rate to loose a contact

    double stateChangeTime;  //last time of state change (i.e. S->I, I->R, etc.)


    State state;

};


#endif //ALGO_SPECIE_H
