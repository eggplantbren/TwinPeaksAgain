// Includes
#include <ctime>
#include <iostream>
#include "RNG.h"
#include "Sampler.hpp"
#include "Examples/SimpleExample.h"

int main()
{
    // Create an RNG
    TwinPeaks::RNG rng(time(0));

    // Create and initialise a Sampler
    TwinPeaks::Sampler<TwinPeaks::SimpleExample> sampler(10, 1000);

    while(sampler.get_depth() < 500.0)
        sampler.do_iteration(rng);

    return 0;
}

