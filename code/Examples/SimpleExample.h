#ifndef TwinPeaks_SimpleExample
#define TwinPeaks_SimpleExample

// Includes
#include <tuple>
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
    public:
        // Number of scalars in this problem
        static constexpr size_t num_scalars = 2;

    private:
        // The parameters
        std::vector<double> xs;

    public:
        // Constructor doesn't do much,
        // just sets up xs with the right size.
        SimpleExample();

        // Generate from the prior
        void from_prior(RNG& rng);

        // Metropolis proposal
        double perturb(RNG& rng);

        // Evaluate and return one of the scalars
        double get_scalar(size_t which) const;

        // Get all of the scalars
        std::vector<double> get_scalars() const;
};

} // namespace TwinPeaks

#endif

