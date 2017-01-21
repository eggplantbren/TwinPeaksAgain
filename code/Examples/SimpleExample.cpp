#include "SimpleExample.h"

namespace TwinPeaks
{

SimpleExample::SimpleExample()
:xs(100)
{

}

void SimpleExample::from_prior(RNG& rng)
{
    for(double& x: xs)
        x = rng.rand();
}





} // namespace TwinPeaks

