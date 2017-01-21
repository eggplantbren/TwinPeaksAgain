#ifndef TwinPeaks_Sampler
#define TwinPeaks_Sampler

// Includes
#include <iostream>
#include <stdexcept>
#include <stdlib.h>     // For size_t
#include <tuple>
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
        std::vector<ParticleType> particles;

        // The scalars of the particles
        std::vector<std::tuple<double, double>> scalars;

        // Particles that define the context
        std::vector<ParticleType> context_particles;
        std::vector<std::tuple<double, double>> context_scalars;

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
,scalars(num_particles)
,context_particles(num_particles)
,context_scalars(num_particles)
,initialised(false)
{
    if(num_particles == 0)
        throw std::invalid_argument("Invalid number of particles.");
}

template<class ParticleType>
void Sampler<ParticleType>::initialise(RNG& rng)
{
    // Generate the particles
    std::cout<<"# Generating "<<num_particles<<" particles from the prior...";
    std::cout<<std::flush;
    for(size_t i=0; i<num_particles; ++i)
    {
        particles[i].from_prior(rng);
        scalars[i] = particles[i].get_scalars();
    }
    std::cout<<"done."<<std::endl;

    // Generate the *context* particles
    std::cout<<"# Generating "<<num_particles<<' ';
    std::cout<<"context particles from the prior...";
    std::cout<<std::flush;
    for(size_t i=0; i<num_particles; ++i)
    {
        context_particles[i].from_prior(rng);
        context_scalars[i] = context_particles[i].get_scalars();
    }
    std::cout<<"done."<<std::endl;

    initialised = true;
}

} // namespace TwinPeaks

#endif

