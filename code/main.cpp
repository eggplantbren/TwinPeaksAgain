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
    TwinPeaks::Sampler<TwinPeaks::SimpleExample> sampler(100);
    sampler.initialise(rng);

    // Print some information
    auto scalars = sampler.get_scalars();
    auto uccs    = sampler.get_uccs();
    for(size_t i=0; i<scalars.size(); ++i)
    {
        std::cout<<std::get<0>(scalars[i])<<' '<<std::get<1>(scalars[i])<<' ';
        std::cout<<uccs[i]<<'\n';
    }

    return 0;
}

