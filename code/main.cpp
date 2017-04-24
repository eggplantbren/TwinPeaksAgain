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
    TwinPeaks::Sampler<TwinPeaks::SimpleExample> sampler(100, 1000, 10);

    // Ascend each scalar
    while(true)
    {
        sampler.run_to_depth(rng, 70.0);
        bool more_to_do = sampler.next_task();

        if(!more_to_do)
            break;
    }

    return 0;
}

