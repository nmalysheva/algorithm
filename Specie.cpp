//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include "Specie.h"
#include <limits>
#include <iostream>
#include <math.h>
#include <stdexcept>



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

    /*if (numberOfContacts == maxNumberOfContacts)
    {
        result = 0;
    }*/
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


bool Specie::infect (double infectTime)
{
    bool result = false;
    if (state == S)
    {
        state = I;
        infectionTime = infectTime;
        setNewContactRate(getNewContactRate() * 0.03);
        setLooseContactRate(getLooseContactRate() * 2.3);
        result = true;

    }
    else
    {
        std::string msg = "EXCEPTION. ONLY SUSCEPTIBLE CAN BE INFECTED";
        throw std::domain_error(msg);
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
    else
    {
        throw std::domain_error("Current status is .ONLY Infected CAN RECOVER");
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
    bool result = ( /*age == sp.getAge() &&*/
                   maxNumberOfContacts == sp.getMaxNumberOfContacts() &&
                   numberOfContacts == sp.getNumberOfContacts() &&
                   deathRate == sp.getDeathRate() &&
                   newContactRate == sp.getNewContactRate() &&
                   looseContactRate == sp.getLooseContactRate() &&
                   infectionTime == sp.getInfectionTime() &&
                   recoveryTime == sp.getRecoveryTime() &&
                   state == sp.getState());
    return result;

}


double Specie::getNumberOfContactsLimit(double t) const
{
    double a = newContactRate * maxNumberOfContacts;
    double b = newContactRate + looseContactRate;
    /*std::cout << "a: " << a << "; b: " << b << std::endl;

    std::cout << "numberOfContacts: " << numberOfContacts << std::endl;
    std::cout << "ExpectationOfContacts: " << ExpectationOfContacts(a, b, t) << std::endl;
    std::cout << "t: " << t << std::endl;*/
    double numConStart = numberOfContacts;
    double numConEnd   = ExpectationOfContacts(a, b, t) + 2 * sqrt(VarianceOfContacts(a, b, t));
    numConEnd = std::min(numConEnd, static_cast<double>(maxNumberOfContacts));
    /*std::cout << "numConStart:" << numConStart << std::endl;
    std::cout << "numConEnd:" << numConEnd << std::endl;*/
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
