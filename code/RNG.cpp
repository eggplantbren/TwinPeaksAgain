#include "RNG.h"
#include <cmath>

namespace TwinPeaks
{

RNG::RNG()
:uniform(0., 1.)
,normal(0., 1.)
{

}

RNG::RNG(unsigned int seed)
:uniform(0., 1.)
,normal(0., 1.)
{
    set_seed(seed);
}

void RNG::set_seed(unsigned int seed)
{
    twister.seed(seed);
}

double RNG::rand()
{
    return uniform(twister);
}

double RNG::randn()
{
    return normal(twister);
}

double RNG::randt2()
{
    // t-distribution with 2 degrees of freedom (less ridiculous than Cauchy)
    return this->randn()/sqrt(-log(this->rand()));
}

double RNG::randh()
{
    return pow(10.0, 1.5 - 3*std::abs(this->randt2()))*this->randn();
}

int RNG::rand_int(int N)
{
    return static_cast<int>(floor(N*this->rand()));
}

} // namespace InfoNest

