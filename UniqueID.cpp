//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#include "UniqueID.h"

uint32_t UniqueID::nextID = 0;

UniqueID::UniqueID()
{
    id = ++nextID;
}

void UniqueID::reset()
{
    nextID = 0;
}

