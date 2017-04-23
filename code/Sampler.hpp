#ifndef TwinPeaks_Sampler
#define TwinPeaks_Sampler

// Includes
#include <fstream>
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

    public:
        // Constructor. You must specify the number of particles.
        Sampler(size_t num_particles);

        // Do a single NS iteration
        void do_iteration(RNG& rng,
                          bool generate_new_particle=true);

    /***** Private helper functions *****/
    private:

        // Initialise all the particles from the prior
        void initialise(RNG& rng);
};


/*************************************************/
/*           IMPLEMENTATIONS BELOW               */
/*************************************************/

template<class ParticleType>
Sampler<ParticleType>::Sampler(size_t num_particles)
:particles(num_particles)
,scalars(num_particles)
,tiebreakers(num_particles)
,iteration(0)
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

    std::cout << "done." << std::endl;
}

template<class ParticleType>
void Sampler<ParticleType>::do_iteration(RNG& rng,
                                         bool generate_new_particle)
{
    if(iteration == 0)
        initialise(rng);


    ++iteration;
}

} // namespace TwinPeaks

#endif

