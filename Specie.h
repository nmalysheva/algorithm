//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#ifndef ALGO_SPECIE_H
#define ALGO_SPECIE_H


class Specie {
public:
    enum State {S, I, R};

public:

    //constructors

    Specie();  //default constructor
    Specie(//unsigned char age,
           unsigned char maxNumberOfContacts,
           //unsigned int cumNumberOfContacts,
           unsigned char numberOfContacts,
           double deathRate,
           double newContactRate,
           double looseContactRate,
           State st = S);

    //getters
    //unsigned char const getAge() const; //return age of the specie
    //unsigned int  const getMaxNumberOfContacts() const; //return max. number of contacts per year
    //unsigned int  const getCumNumberOfContacts() const; //return cumulative number of contacts for the current year

    unsigned char  getMaxNumberOfContacts() const; //return max. number of contacts per any time instance
    unsigned char  getNumberOfContacts() const; //return max. number of contacts per any time instance

    double getNewContactRate() const; //return rate for establishing a new contact
    double getLooseContactRate() const; ////return rate for loosing of the one of already existing contacts

    double  getDeathRate() const;  //return specie's death rate

    double  getInfectionTime() const; //return last infection time if specie is infected,  -infinity if not
    double  getRecoveryTime() const;  //return last recovery time if specie was recovered, -infinity if not

    State  getState() const; // return epidemic state of the specie

    //setters
    //void setAge(unsigned char age);
    void setMaxNumberOfContacts(unsigned int maxNumOfCont);
    //void setCumNumberOfContacts(unsigned int cumNumOfCont);
    //void incCumNumberOfContacts();
    //void resetCumNumberOfContacts();
    void incNumberOfContacts();
    void decNumberOfContacts();


    void setNewContactRate  (double newContRate);
    void setLooseContactRate(double looseContRate);

    void setDeathRate(double dRate);


    //status changes

    bool infect (double infectTime);  //change status of the specie to INFECTED
    bool recover(double recTime);     //change status of the specie to RECOVERED



    //comparison operator for map
    bool operator== (const Specie &sp) const;

    //destructor
    ~Specie() {};


private:

    void setNumberOfContacts(unsigned char nOfCont);
    //unsigned char age;
    //unsigned char  maxNumberOfContacts; // per year
    //unsigned int  cumNumberOfContacts; //cumulative number of contacts. Reset every year.

    unsigned char  maxNumberOfContacts; //at any timepoint
    unsigned char  numberOfContacts; //current number of contacts
    double deathRate;
    double newContactRate;  //rate to establish a new contact
    double looseContactRate;//rate to loose a contact

    double infectionTime;  //last infection time
    double recoveryTime;   //last recovery time

    State state;

};


#endif //ALGO_SPECIE_H
