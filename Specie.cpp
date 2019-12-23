//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include "Specie.h"
#include <limits>
#include <iostream>

Specie::Specie()
{
    //age = static_cast<unsigned char> (-1);

    maxNumberOfContacts = static_cast<size_t> (-1);
    //cumNumberOfContacts = static_cast<unsigned char> (-1);
    numberOfContacts = static_cast<size_t> (-1);

    deathRate        = std::numeric_limits<double>::min();
    newContactRate   = std::numeric_limits<double>::min();
    looseContactRate = std::numeric_limits<double>::min();

    state = S;
    infectionTime =   std::numeric_limits<double>::infinity(); // -infinity
    recoveryTime  = - std::numeric_limits<double>::infinity();
}


Specie::Specie(/*unsigned char age, */size_t maxNumberOfContacts,
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
    infectionTime = std::numeric_limits<double>::min();
    recoveryTime  = std::numeric_limits<double>::min();

}


//getters
/*unsigned char const Specie::getAge() const
{
    return age;
}*/

size_t Specie::getMaxNumberOfContacts() const
{
    return maxNumberOfContacts;
}


/*unsigned int const Specie::getCumNumberOfContacts() const
{
    return cumNumberOfContacts;
}*/

size_t Specie::getNumberOfContacts() const
{
    return numberOfContacts;
}

double Specie::getNewContactRate() const
{
    double result = newContactRate;

    if (numberOfContacts == maxNumberOfContacts)
    {
        result = 0;
    }
    return result;
}

double Specie::getLooseContactRate() const
{
    double result = looseContactRate;
    if (numberOfContacts == 0)
    {
        result = 0;
    }
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

//setters
/*void Specie::setAge(unsigned char age)
{
    this->age = age;
}*/

void Specie::setMaxNumberOfContacts(unsigned int maxNumOfCont)
{
    maxNumberOfContacts = maxNumOfCont;
}

void Specie::setNumberOfContacts(size_t nOfCont)
{
    numberOfContacts = nOfCont;
}

/*void Specie::setCumNumberOfContacts(unsigned int cumNumOfCont)
{
    cumNumberOfContacts = cumNumOfCont;
}*/

/*void Specie::incCumNumberOfContacts()
{
    cumNumberOfContacts++;
    //if (cumNumberOfContacts >= maxNumberOfContacts)
    //{
    //    newContactRate = 0;
    //}
}*/

/*void Specie::resetCumNumberOfContacts()
{
    cumNumberOfContacts = 0;
}*/

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


bool Specie::infect (double infectTime)
{
    bool result = false;
    if (state == S)
    {
        state = I;
        infectionTime = infectTime;
        result = true;

    }
    else //throw exception
    {
        std::cout << "EXCEPTION. ONLY SUSCEPTIBLE CAN BE INFECTED";
    }
    return result;
}

bool Specie::recover(double recTime)
{
    bool result = false;
    if (state == I)
    {
        state = R;
        recoveryTime = recTime;
        result = true;

    }
    else //throw exception
    {
        std::cout << "EXCEPTION. ONLY Infected CAN RECOVER";
    }
    return result;
}

double Specie::getInfectionTime() const
{
    return infectionTime;
}
double Specie::getRecoveryTime() const
{
    return recoveryTime;
}


bool Specie::operator== (const Specie &sp) const
{
    //bool result = age == sp.getAge();

    bool result = ( /*age == sp.getAge() &&*/
                   maxNumberOfContacts == sp.getMaxNumberOfContacts() &&
                   //cumNumberOfContacts == sp.getCumNumberOfContacts() &&
                   numberOfContacts == sp.getNumberOfContacts() &&
                   deathRate == sp.getDeathRate() &&
                   newContactRate == sp.getNewContactRate() &&
                   looseContactRate == sp.getLooseContactRate() &&
                   infectionTime == sp.getInfectionTime() &&
                   recoveryTime == sp.getRecoveryTime() &&
                   state == sp.getState());
    return result;

}
