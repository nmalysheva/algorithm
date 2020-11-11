//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include "Specie.h"
#include <limits>
#include <iostream>
#include <cmath>
#include <stdexcept>



Specie::Specie()
{
    //age = static_cast<unsigned char> (-1);

    maxNumberOfContacts = 0;
    numberOfContacts = 0;

    deathRate        = 0;
    newContactRate   = 0;
    looseContactRate = 0;

    state = S;
    stateChangeTime =   0; // -infinity
}


Specie::Specie(/*unsigned char age,*/ size_t maxNumberOfContacts,
                                      size_t numberOfContacts, double deathRate,
               double newContactRate, double looseContactRate, State st)
{
    //setAge(age);
    setMaxNumberOfContacts(maxNumberOfContacts);
    setNumberOfContacts(numberOfContacts);
    setDeathRate(deathRate);
    setNewContactRate(newContactRate);
    setLooseContactRate(looseContactRate);

    state = st;
    stateChangeTime = 0;

}

size_t Specie::getMaxNumberOfContacts() const
{
    return maxNumberOfContacts;
}

size_t Specie::getNumberOfContacts() const
{
    return numberOfContacts;
}

double Specie::getNewContactRate() const
{
    double result = newContactRate * (maxNumberOfContacts - numberOfContacts);

    return result;
}

double Specie::getLooseContactRate() const
{
    double result = looseContactRate;// / numberOfContacts; //* numberOfContacts;
    /*if (numberOfContacts == 0)
    {
        result = 0;
    }*/
    return result;
}

double Specie::getDeathRate() const
{
    return deathRate;
}

Specie::State Specie::getState() const
{
    return state;
}

void Specie::setMaxNumberOfContacts(unsigned int maxNumOfCont)
{
    maxNumberOfContacts = maxNumOfCont;
}

void Specie::setNumberOfContacts(size_t nOfCont)
{
    numberOfContacts = nOfCont;
}

void Specie::incNumberOfContacts()
{
    numberOfContacts++;
}
void Specie::decNumberOfContacts()
{
    numberOfContacts--;
}

void Specie::setNewContactRate  (double newContRate)
{
    newContactRate = newContRate;
}

void Specie::setLooseContactRate(double looseContRate)
{
    looseContactRate = looseContRate;
}

void Specie::setDeathRate(double dRate)
{
    deathRate = dRate;
}


void Specie::changeState (State st, double time)
{
    state = st;
    stateChangeTime = time;
}

double Specie::getLastStateChangeTime() const
{
    return stateChangeTime;
}


bool Specie::operator== (const Specie &sp) const
{
    bool result = ( /*age == sp.getAge() &&*/
                   maxNumberOfContacts == sp.getMaxNumberOfContacts() &&
                   numberOfContacts == sp.getNumberOfContacts() &&
                   deathRate == sp.getDeathRate() &&
                   newContactRate == sp.getNewContactRate() &&
                   looseContactRate == sp.getLooseContactRate() &&
                   stateChangeTime == sp.getLastStateChangeTime() &&
                   state == sp.getState());
    return result;

}


double Specie::getNumberOfContactsLimit(double t) const
{
    double a = newContactRate * maxNumberOfContacts;
    double b = newContactRate + looseContactRate;
    double numConStart = numberOfContacts;
    double numConEnd   = ExpectationOfContacts(a, b, t) + 2 * sqrt(VarianceOfContacts(a, b, t));
    numConEnd = std::min(numConEnd, static_cast<double>(maxNumberOfContacts));
    double result = std::max(numConStart, numConEnd);

    return result;

}

double Specie::ExpectationOfContacts(double a, double b, double t)const
{

    return (a / b - (a / b - numberOfContacts) * exp(-b * t));
}

double Specie::VarianceOfContacts(double a, double b, double t)const
{
    return (ExpectationOfContacts(a, b, t) - numberOfContacts * exp(-2 * b * t));
}
