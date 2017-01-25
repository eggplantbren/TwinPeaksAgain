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

        // Keep track of remaining mass
        double log_mass;

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

        // Find the worst two particles in terms of the chosen scalar.
        std::tuple<size_t, size_t> find_worst_two(size_t scalar) const;

        // Replace the given particle
        unsigned int replace_particle(size_t which_particle, RNG& rng);
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
,log_mass(0.0)
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

    // Choose a scalar
    size_t scalar = rng.rand_int(2);

    // Find the worst two particles in terms of that scalar.
    size_t i1, i2;
    std::tie(i1, i2) = find_worst_two(scalar);

    // Extract scalars
    double x1, y1, x2, y2;
    std::tie(x1, y1) = scalars[i1];
    std::tie(x2, y2) = scalars[i2];

    // Write out particle information.
    std::fstream fout("sample_info.txt", std::ios::out | std::ios::app);
    fout<<iteration<<' '<<log_mass<<' ';
    fout<<x1<<' '<<y1<<std::endl;
    fout.close();

    // Add to Context
    if(scalar == 0)
        context.add_rectangle({x2, y1});
    else
        context.add_rectangle({x1, y2});

    // Decrement remaining mass
    log_mass -= 1.0/particles.size();

    // Generate replacement particles
    std::cout<<"Generating replacement particles..."<<std::flush;
    unsigned int accepts = 0;
    accepts += replace_particle(i1, rng);
    accepts += replace_particle(i2, rng);
    std::cout<<"accepted "<<accepts<<"/10000...";
    std::cout<<"done."<<std::endl;
}

template<class ParticleType>
unsigned int Sampler<ParticleType>::
                    replace_particle(size_t which_particle, RNG& rng)
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
    constexpr unsigned int mcmc_steps = 10000;
    unsigned int accepted = 0;
    for(unsigned int i=0; i<mcmc_steps/2; ++i)
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

    return accepted;
}

template<class ParticleType>
std::tuple<size_t, size_t> Sampler<ParticleType>::
                    find_worst_two(size_t scalar) const
{
    if(scalar != 0 && scalar != 1)
        throw std::invalid_argument("Invalid argument to find_worst_two.");

    // Vector of values of the chosen scalar
    std::vector<double> s(scalars.size());

    double x, y;
    for(size_t i=0; i<scalars.size(); ++i)
    {
        std::tie(x, y) = scalars[i];
        s[i] = (scalar == 0) ? (x) : (y);
    }

    // Argsort
    auto ii = argsort(s);

    return {ii[0], ii[1]};
}

} // namespace TwinPeaks

#endif

