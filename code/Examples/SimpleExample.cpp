#include "SimpleExample.h"
#include "Utils.h"

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
            s += -pow(sin(4.0*M_PI*x), 2);
    }

    return s;
}

std::vector<double> SimpleExample::get_scalars() const
{
    std::vector<double> result(num_scalars);
    for(size_t i=0; i<num_scalars; ++i)
        result[i] = get_scalar(i);
    return result;
}

} // namespace TwinPeaks

