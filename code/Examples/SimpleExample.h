#ifndef TwinPeaks_SimpleExample
#define TwinPeaks_SimpleExample

#include <vector>

#include "RNG.h"

namespace TwinPeaks
{

/*
* A simple example where the two scalars become incompatible
* rather quickly. i.e. the old TwinPeaks algorithm breaks.
*/

class SimpleExample
{
    private:
        // The parameters
        std::vector<double> xs;

    public:
        // Constructor doesn't do much, just sets up xs with the right size.
        SimpleExample();

        // Generate from the prior
        void from_prior(RNG& rng);

};


} // namespace TwinPeaks

#endif

