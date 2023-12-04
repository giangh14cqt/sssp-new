//
// Created by Truong Giang Do on 28/11/2023.
//

#ifndef SSSP_NEW_RANDD_H
#define SSSP_NEW_RANDD_H

#include "Graph.h"

class Random {
public:
private:
    // Make the default constructor private.
    Random() {}

public:
    // Delete the copy constructor function.
    Random(const Random &) = delete;

    // Delete the overloading of assignment operator
    Random &operator=(const Random &) = delete;

    static Random &Get() {
        static Random inst;
        return inst;
    }

    // Seed the random number generator.
    static void Seed() {
        srand(time(nullptr));
    }

    static int GenInt() {
        return rand();
    }

    static int GenInt(int min, int max) {
        return rand() % (max - min + 1) + min;
    }
};

#endif //SSSP_NEW_RANDD_H
