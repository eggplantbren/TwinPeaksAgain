#ifndef TwinPeaks_Atoms
#define TwinPeaks_Atoms

#include "../RNG.h"
#include <ostream>
#include <string>
#include <vector>

namespace TwinPeaks
{

class Atoms
{
    public:
        static constexpr size_t num_scalars = 2;

    private:
        static constexpr int num_atoms = 30;
        static constexpr double L = 100.0;

        // Positions
        std::vector<double> x, y, z;

    public:
        Atoms();

        // Generate the point from the prior
        void from_prior(RNG& rng);

        // Metropolis-Hastings proposals
        double perturb(RNG& rng);

        // Getter
        double get_scalar(size_t which) const;

        // Print to output stream (CSV-friendly format)
        void print(std::ostream& out) const;

        // For CSV header
        static std::string description();
};

} // namespace TwinPeaks

#endif

