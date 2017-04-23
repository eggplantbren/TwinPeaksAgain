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

    for(int i=0; i<10*100; ++i)
        sampler.do_iteration(rng);

    return 0;
}

