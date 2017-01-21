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
    TwinPeaks::Sampler<TwinPeaks::SimpleExample> sampler(1000);
    sampler.initialise(rng);

    // Do some iterations.
    for(int i=0; i<1000; ++i)
        sampler.do_iteration(rng);

    return 0;
}

