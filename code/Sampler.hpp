#ifndef TwinPeaks_Sampler
#define TwinPeaks_Sampler

// Includes
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>     // For size_t
#include <tuple>
#include <vector>
#include "Context.h"
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

        // The scalars of the particles
        std::vector<std::tuple<double, double>> scalars;

        // The context
        Context context;

        // Flag for whether the Sampler is ready to go.
        bool initialised;

        // Iteration counter
        unsigned int iteration;

    public:
        // Constructor. You must specify the number of particles.
        Sampler(size_t num_particles);

        // Generate all particles from the prior, prepare the sampler
        void initialise(RNG& rng);

        // Do a single NS iteration
        void do_iteration(RNG& rng);

        // Some getters
        const std::vector<std::tuple<double, double>> get_scalars()
        { return scalars; }

    /***** Private helper functions *****/
    private:

        // Find the worst two particles in terms of first scalar.
        std::tuple<size_t, size_t> find_worst() const;

        // Replace the given particle
        void replace_particle(size_t which_particle, RNG& rng);
};


/*************************************************/
/*           IMPLEMENTATIONS BELOW               */
/*************************************************/

template<class ParticleType>
Sampler<ParticleType>::Sampler(size_t num_particles)
:particles(num_particles)
,scalars(num_particles)
,context()
,initialised(false)
,iteration(0)
{
    if(num_particles == 0)
        throw std::invalid_argument("Invalid number of particles.");
}

template<class ParticleType>
void Sampler<ParticleType>::initialise(RNG& rng)
{
    // Generate the particles
    std::cout<<"# Generating "<<particles.size()<<' ';
    std::cout<<"particles from the prior...";
    std::cout<<std::flush;
    for(size_t i=0; i<particles.size(); ++i)
    {
        particles[i].from_prior(rng);
        scalars[i] = particles[i].get_scalars();
    }
    std::cout<<"done."<<std::endl;

    // Clear output file
    std::fstream fout("sample_info.txt", std::ios::out);
    fout.close();

    initialised = true;
    iteration = 0;
}

template<class ParticleType>
void Sampler<ParticleType>::do_iteration(RNG& rng)
{
    // Increment iteration
    ++iteration;

    // Print message
    std::cout<<"# Iteration "<<iteration<<". ";

    // Find the worst particle.
    std::tuple<size_t, size_t> worst = find_worst();
//    auto worst_scalars = scalars[worst];

    // Write out iteration and worst particle's scalars.
    std::fstream fout("sample_info.txt", std::ios::out | std::ios::app);
    fout<<iteration<<' ';
//    fout<<std::get<0>(worst_scalars)<<' ';
//    fout<<std::get<1>(worst_scalars)<<std::endl;
    fout.close();

    // Generate replacement particle
    std::cout<<"Generating replacement particle..."<<std::flush;
//    replace_particle(worst, rng);
    std::cout<<"done."<<std::endl;
}

template<class ParticleType>
void Sampler<ParticleType>::replace_particle(size_t which_particle, RNG& rng)
{
    // Copy a survivor
    if(particles.size() > 1)
    {
        size_t copy;
        do
        {
            copy = rng.rand_int(particles.size());
        }while(copy == which_particle);

        particles[which_particle] = particles[copy];
        scalars[which_particle] = scalars[copy];
    }

    // Do MCMC steps
    constexpr unsigned int mcmc_steps = 1000;
    unsigned int accepted = 0;
    for(unsigned int i=0; i<mcmc_steps; ++i)
    {
        // Generate proposal
        ParticleType proposal = particles[which_particle];
        double logH = proposal.perturb(rng);

        // Evaluate proposal
        auto proposal_scalars = proposal.get_scalars();

        // Accept?
        if(context.is_okay(proposal_scalars) && rng.rand() <= exp(logH))
        {
            particles[which_particle] = proposal;
            scalars[which_particle] = proposal_scalars;
            ++accepted;
        }
    }
    std::cout<<"accepted "<<accepted<<'/'<<mcmc_steps<<"...";
}

template<class ParticleType>
std::tuple<size_t, size_t> Sampler<ParticleType>::find_worst() const
{
    // Vector of values of scalar 1
    std::vector<double> s(scalars.size());
    for(size_t i=0; i<scalars.size(); ++i)
        s[i] = std::get<0>(scalars[i]);

    // Argsort
    auto ii = argsort(s);

    return {ii[0], ii[1]};
}

} // namespace TwinPeaks

#endif

