#include "SimpleExample.h"
#include "../Utils.h"

namespace TwinPeaks
{

SimpleExample::SimpleExample()
:xs(10)
{

}

void SimpleExample::from_prior(RNG& rng)
{
    for(double& x: xs)
        x = rng.rand();
}

double SimpleExample::perturb(RNG& rng)
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

double SimpleExample::get_scalar(size_t which) const
{
    if(which >= num_scalars)
        throw std::invalid_argument("Invalid scalar requested.");

    double s = 0.0;
    if(which == 0)
    {
        for(double x: xs)
            s += -pow(x - 0.5, 2);
    }
    else
    {
        for(double x: xs)
            s += -std::abs(x);
    }

    return s;
}

void SimpleExample::print(std::ostream& out) const
{
    for(size_t i=0; i<xs.size(); ++i)
    {
        out << xs[i];
        if(i != xs.size()-1)
            out << ",";
    }
}

} // namespace TwinPeaks

