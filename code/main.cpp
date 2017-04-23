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

    // Create a Sampler
    TwinPeaks::Sampler<TwinPeaks::SimpleExample> sampler(10, 1000);

    // Run it
    sampler.run_to_depth(rng, 500.0);

    return 0;
}

