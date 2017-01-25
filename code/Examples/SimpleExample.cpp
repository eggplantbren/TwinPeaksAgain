#include "SimpleExample.h"
#include "Utils.h"

namespace TwinPeaks
{

SimpleExample::SimpleExample()
:xs(2)
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

std::tuple<double, double> SimpleExample::get_scalars() const
{
    double s1 = 0.0;
    double s2 = 0.0;

    s1 = xs[0];
    s2 = xs[1];

//    for(double x: xs)
//    {
//        s1 += x[0];
//        s2 += x[1];
//    }

    return {s1, s2};
}

} // namespace TwinPeaks

