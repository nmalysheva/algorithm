//
// Created by Malysheva, Nadezhda on 2019-07-28.
//

#ifndef ALGO_UNIQUEID_H
#define ALGO_UNIQUEID_H


#include <cstdint>

class UniqueID {
protected:
    static uint32_t nextID;

public:
    uint32_t id;
    UniqueID();
    void reset();
};



#endif //ALGO_UNIQUEID_H
