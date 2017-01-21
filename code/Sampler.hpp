#ifndef TwinPeaks_Sampler
#define TwinPeaks_Sampler

// Includes
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>     // For size_t
#include <tuple>
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

        // The scalars of the particles
        std::vector<std::tuple<double, double>> scalars;

        // The uccs of the particles, and tiebreakers
        std::vector<size_t> uccs;
        std::vector<double> ucc_tiebreakers;

        // Particles that define the context
        std::vector<ParticleType> context_particles;
        std::vector<std::tuple<double, double>> context_scalars;

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
        const std::vector<size_t> get_uccs() const
        { return uccs; }


    /***** Private helper functions *****/
    private:

        // Calculate the ucc wrt the context scalars
        size_t calculate_ucc(const std::tuple<double, double>& s) const;

        // Calculate uccs of all particles
        void calculate_uccs();

        // Find the worst particle and return its index.
        size_t find_worst() const;

        // Replace the given particle with one above the given threshold
        void replace_particle(size_t which_particle,
                              const std::tuple<size_t, double>& threshold,
                              RNG& rng);

};


/*************************************************/
/*           IMPLEMENTATIONS BELOW               */
/*************************************************/

template<class ParticleType>
Sampler<ParticleType>::Sampler(size_t num_particles)
:particles(num_particles)
,scalars(num_particles)
,uccs(num_particles)
,ucc_tiebreakers(num_particles)
,context_particles(num_particles)
,context_scalars(num_particles)
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

    // Generate the *context* particles
    std::cout<<"# Generating "<<context_particles.size()<<' ';
    std::cout<<"context particles from the prior...";
    std::cout<<std::flush;
    for(size_t i=0; i<context_particles.size(); ++i)
    {
        context_particles[i].from_prior(rng);
        context_scalars[i] = context_particles[i].get_scalars();
    }
    std::cout<<"done."<<std::endl;

    // Calculate the uccs and generate tiebreakers
    calculate_uccs();
    for(double& tb: ucc_tiebreakers)
        tb = rng.rand();

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
    size_t worst = find_worst();
    auto worst_scalars = scalars[worst];

    // Write out iteration and worst particle's scalars.
    std::fstream fout("sample_info.txt", std::ios::out | std::ios::app);
    fout<<iteration<<' ';
    fout<<std::get<0>(worst_scalars)<<' ';
    fout<<std::get<1>(worst_scalars)<<std::endl;
    fout.close();

    // Generate replacement particle
    std::cout<<"Generating replacement particle..."<<std::flush;
    replace_particle(worst,
                     {uccs[worst], ucc_tiebreakers[worst]},
                     rng);
    std::cout<<"done."<<std::endl;
}

template<class ParticleType>
void Sampler<ParticleType>::replace_particle
                             (size_t which_particle,
                              const std::tuple<size_t, double>& threshold,
                              RNG& rng)
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
        uccs[which_particle] = uccs[copy];
        ucc_tiebreakers[which_particle] = ucc_tiebreakers[copy];
    }

    // Threshold, smooshed together
    double thresh = std::get<0>(threshold) + std::get<1>(threshold);

    // Do MCMC steps
    constexpr unsigned int mcmc_steps = 1000;
    unsigned int accepted = 0;
    for(unsigned int i=0; i<mcmc_steps; ++i)
    {
        // Generate proposal
        ParticleType proposal = particles[which_particle];
        double logH = proposal.perturb(rng);
        double proposal_ucc_tiebreaker = ucc_tiebreakers[which_particle]
                                                            + rng.randh();
        wrap(proposal_ucc_tiebreaker, 0.0, 1.0);

        // Evaluate proposal
        auto proposal_scalars = proposal.get_scalars();
        size_t proposal_ucc = calculate_ucc(proposal_scalars);
        double ucc_smooshed = proposal_ucc + proposal_ucc_tiebreaker;

        // Accept?
        if(ucc_smooshed < thresh && rng.rand() <= exp(logH))
        {
            particles[which_particle] = proposal;
            scalars[which_particle] = proposal_scalars;
            uccs[which_particle] = proposal_ucc;
            ucc_tiebreakers[which_particle] = proposal_ucc_tiebreaker;
            ++accepted;
        }
    }
    std::cout<<"accepted "<<accepted<<'/'<<mcmc_steps<<"...";
}

template<class ParticleType>
size_t Sampler<ParticleType>::calculate_ucc
                            (const std::tuple<double, double>& s) const
{
    size_t ucc = 0;
    for(size_t j=0; j<context_particles.size(); ++j)
    {
        if(both_above(context_scalars[j], s))
            ++ucc;
    }

    return ucc;
}

template<class ParticleType>
void Sampler<ParticleType>::calculate_uccs()
{
    for(size_t i=0; i<particles.size(); ++i)
        uccs[i] = calculate_ucc(scalars[i]);
}

template<class ParticleType>
size_t Sampler<ParticleType>::find_worst() const
{
    // Worst particle found so far
    size_t worst_index = 0;
    double worst_value = uccs[worst_index] + ucc_tiebreakers[worst_index];

    for(size_t i=1; i<particles.size(); ++i)
    {
        double current_value = uccs[i] + ucc_tiebreakers[i];
        if(current_value > worst_value)
        {
            worst_index = i;
            worst_value = current_value;
        }
    }

    return worst_index;
}

} // namespace TwinPeaks

#endif

