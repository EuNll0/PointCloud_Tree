#ifndef _MY_COMMON_H
#define _MY_COMMON_H
#include <iostream>
#include <memory>
#include <limits>
#include <vector>
#include <glog/logging.h>

#include "Aggrate/Baseshape.h"
#include "Aggrate/Memory.h"

static float MachineEpsilon = std::numeric_limits<float>::epsilon() * 0.5;
inline float gamma(int n) { return (n * MachineEpsilon) / (1 - n * MachineEpsilon); }

template <typename T>
inline T getgamma(int n)
{
    T MachineEpsilon = std::numeric_limits<float>::epsilon() * 0.5;
    return 1+ 2 * (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

#endif
