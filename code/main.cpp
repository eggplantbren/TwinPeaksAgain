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

    // Create a Sampler and ascend the two scalars
    TwinPeaks::Sampler<TwinPeaks::SimpleExample> sampler(1, 1000);

    while(true)
    {
        sampler.run_to_depth(rng, 500.0);
        bool more_to_do = sampler.next_task(rng);

        if(!more_to_do)
            break;
    }

    return 0;
}

