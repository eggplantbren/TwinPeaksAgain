#include "SpikeSlab.h"
#include "../Utils.h"
#include <sstream>

namespace TwinPeaks
{

SpikeSlab::SpikeSlab()
:xs(size)
{

}

void SpikeSlab::from_prior(RNG& rng)
{
    for(double& x: xs)
        x = rng.rand();
}

double SpikeSlab::perturb(RNG& rng)
{
    // How many coordinates to perturb
    int reps = 1;
    if(rng.rand() <= 0.5)
        reps += (int)(pow((double)xs.size(), rng.rand()));

    // Perturb them
    int which;
    for(int i=0; i<reps; ++i)
    {
        which = rng.rand_int(xs.size());
        xs[which] += rng.randh();
        wrap(xs[which], 0.0, 1.0);
    }

    return 0.0;
}

double SpikeSlab::get_scalar(size_t which) const
{
    if(which >= num_scalars)
        throw std::invalid_argument("Invalid scalar requested.");

    // Scalar 0 is the SpikeSlab likelihood,
    // Scalar 1 is one that induces periodicities.

    double s = 0.0;
    if(which == 0)
    {
        static constexpr double u = 0.01;
        static constexpr double v = 0.1;
        static constexpr double C = log(1.0/sqrt(2*M_PI));
        static constexpr double log_half = log(0.5);

        double logl1 = size*(C - log(u));
        double logl2 = size*(C - log(v));

        for(double x: xs)
        {
            logl1 += -0.5*pow(x/u, 2);
            logl2 += -0.5*pow(x/v, 2);
        }

        s = logsumexp(logl1 + log_half, logl2 + log_half);
    }
    else
    {
        for(double x: xs)
            s += -pow(sin(4*M_PI*x), 2);
    }

    return s;
}

void SpikeSlab::print(std::ostream& out) const
{
    for(size_t i=0; i<xs.size(); ++i)
    {
        out << xs[i];
        if(i != xs.size()-1)
            out << ",";
    }
}

std::string SpikeSlab::description()
{
    std::stringstream s;

    for(size_t i=0; i<size; ++i)
    {
        s << "x[" << i << ']';
        if(i != size - 1)
            s << ",";
    }

    return s.str();
}

} // namespace TwinPeaks

