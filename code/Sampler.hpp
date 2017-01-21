#ifndef TwinPeaks_Sampler
#define TwinPeaks_Sampler

#include <stdexcept>
#include <stdlib.h>     // For size_t
#include <vector>

#include "RNG.h"


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
        // How many particles to use
        const size_t num_particles;

        // The particles
        std::vector<double> particles;

        // Flag for whether the Sampler is ready to go.
        bool initialised;

    public:
        // Constructor. You must specify the number of particles.
        Sampler(size_t num_particles);

        // Generate all particles from the prior, prepare the sampler
        void initialise(RNG& rng);



};



/* IMPLEMENTATIONS BELOW */
template<class ParticleType>
Sampler<ParticleType>::Sampler(size_t num_particles)
:num_particles(num_particles)
,particles(num_particles)
,initialised(false)
{
    if(num_particles == 0)
        throw std::invalid_argument("Invalid number of particles.");
}

template<class ParticleType>
void Sampler<ParticleType>::initialise(RNG& rng)
{
    for(ParticleType& particle: particles)
        particle.from_prior(rng);

    initialised = true;
}

} // namespace TwinPeaks

#endif

