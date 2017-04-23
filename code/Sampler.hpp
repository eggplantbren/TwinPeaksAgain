#ifndef TwinPeaks_Sampler
#define TwinPeaks_Sampler

// Includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>     // For size_t
#include <vector>
#include "RNG.h"
#include "Utils.h"

namespace TwinPeaks
{

/*
* A template class, an object of which represents a sampler
* (of a given particle type).
*/

template<class ParticleType>
class Sampler
{
    private:
        // The particles
        std::vector<ParticleType> particles;

        // Corresponding scalars
        std::vector<double> scalars;
        std::vector<double> tiebreakers;

        // Iteration counter
        unsigned int iteration;

        // Number of MCMC steps to use
        size_t mcmc_steps;

    public:
        // Constructor. You must specify the number of particles
        // and MCMC steps per iteration.
        Sampler(size_t num_particles, size_t mcmc_steps);

        // Do a single NS iteration
        void do_iteration(RNG& rng,
                          bool generate_new_particle=true);

    /***** Private helper functions *****/
    private:

        // Initialise all the particles from the prior
        void initialise(RNG& rng);

        // Find the worst particle
        size_t find_worst_particle() const;
};


/*************************************************/
/*           IMPLEMENTATIONS BELOW               */
/*************************************************/

template<class ParticleType>
Sampler<ParticleType>::Sampler(size_t num_particles,
                               size_t mcmc_steps)
:particles(num_particles)
,scalars(num_particles)
,tiebreakers(num_particles)
,iteration(0)
,mcmc_steps(mcmc_steps)
{

}

template<class ParticleType>
void Sampler<ParticleType>::initialise(RNG& rng)
{
    std::cout << "Generating " << particles.size() << " ";
    std::cout << "particles from the prior..." << std::flush;

    for(size_t i=0; i<particles.size(); ++i)
    {
        particles[i].from_prior(rng);
        scalars[i] = particles[i].get_scalar(0);
        tiebreakers[i] = rng.rand();
    }

    std::cout << "done.\n" << std::endl;
}

template<class ParticleType>
void Sampler<ParticleType>::do_iteration(RNG& rng,
                                         bool generate_new_particle)
{
    // On first iteration, generate particles from prior
    if(iteration == 0)
        initialise(rng);
    ++iteration;

    // Worst particle
    size_t worst = find_worst_particle();

    // Print some messages to stdout
    std::cout << std::setprecision(12);
    std::cout << "Iteration " << iteration << ". ";
    std::cout << "log(scalar) = " << scalars[worst] << '.';
    std::cout << std::endl;

}

template<class ParticleType>
size_t Sampler<ParticleType>::find_worst_particle() const
{
    // Zip log likelihoods with tiebreakers
    std::vector< std::pair<double, double> > lltbs(particles.size());
    for(size_t i=0; i<particles.size(); ++i)
    {
        lltbs[i] = std::pair<double, double>
                            {scalars[i], tiebreakers[i]};
    }

    // Argsort
    std::vector<size_t> indices = argsort(lltbs);
    return indices[0];
}

} // namespace TwinPeaks

#endif

