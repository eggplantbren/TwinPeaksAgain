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
        std::vector<double> scalars1;
        std::vector<double> scalars2;

        // The context
        Context context;

        // Flag for whether the Sampler is ready to go.
        bool initialised;

        // Iteration counter
        unsigned int iteration;

        // Keep track of remaining mass
        double log_mass;

        // Argsort results
        std::vector<size_t> indices1;
        std::vector<size_t> indices2;


    public:
        // Constructor. You must specify the number of particles.
        Sampler(size_t num_particles);

        // Generate all particles from the prior, prepare the sampler
        void initialise(RNG& rng);

        // Do a single NS iteration
        void do_iteration(RNG& rng);

    /***** Private helper functions *****/
    private:

        // compute ranks
        void compute_ranks();

        // compute rank map
        void compute_rank_map();

        // Argsorts
        void compute_sort_indices();

        // Replace the given particle
        unsigned int replace_particle(size_t which_particle, RNG& rng);
};


/*************************************************/
/*           IMPLEMENTATIONS BELOW               */
/*************************************************/

template<class ParticleType>
Sampler<ParticleType>::Sampler(size_t num_particles)
:particles(num_particles)
,scalars1(num_particles)
,scalars2(num_particles)
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
        std::tie(scalars1[i], scalars2[i]) = particles[i].get_scalars();
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
//    std::cout<<"# Iteration "<<iteration<<". ";

    // Ranks and rank map
    compute_sort_indices();

    // Print all the particles' scalars
    for(size_t i=0; i<particles.size(); ++i)
    {
        std::cout<<scalars1[i]<<' '<<scalars2[i]<<std::endl;
    }
    std::cout<<'\n';
    std::cout<<"Worst s1 = "<<scalars1[indices1[0]]<<'\n';
    std::cout<<"Worst s2 = "<<scalars2[indices2[0]]<<std::endl;

    exit(0);

//    // Write out particle information.
//    std::fstream fout("sample_info.txt", std::ios::out | std::ios::app);
//    fout<<iteration<<' '<<log_mass<<' ';
//    fout<<x1<<' '<<y1<<std::endl;
//    fout.close();

//    // Add to Context
//    if(scalar == 0)
//        context.add_rectangle({x2, y1});
//    else
//        context.add_rectangle({x1, y2});

//    // Decrement remaining mass
//    log_mass -= 1.0/particles.size();

//    // Generate replacement particles
//    std::cout<<"Generating replacement particles..."<<std::flush;
//    unsigned int accepts = 0;
//    accepts += replace_particle(i1, rng);
//    accepts += replace_particle(i2, rng);
//    std::cout<<"accepted "<<accepts<<"/10000...";
//    std::cout<<"done."<<std::endl;
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
        scalars1[which_particle] = scalars1[copy];
        scalars2[which_particle] = scalars2[copy];
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
            std::tie(scalars1[which_particle], scalars2[which_particle])
                            = proposal_scalars;
            ++accepted;
        }
    }

    return accepted;
}

template<class ParticleType>
void Sampler<ParticleType>::compute_sort_indices()
{
    // Argsort by each scalar
    indices1 = argsort(scalars1);
    indices2 = argsort(scalars2);
}

} // namespace TwinPeaks

#endif

