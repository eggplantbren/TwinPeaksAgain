#ifndef TwinPeaks_Sampler
#define TwinPeaks_Sampler

// Includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <stdlib.h>
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
        size_t more_particles;

        // Current scalar being done
        size_t which_scalar;

        // Results from the single-scalar runs
        std::vector<std::vector<double>> results;

        // Output files
        std::fstream sample_file;
        std::fstream sample_info_file;

        // Thinning for sample_file
        size_t thin;

    public:
        // Constructor. You must specify the number of particles
        // and MCMC steps per iteration. If you want to boost the
        // number of particles for the last run, use more_particles>1.
        Sampler(size_t num_particles, size_t mcmc_steps,
                size_t more_particles=1, size_t thin=1);

        // Run until the specified depth
        void run_to_depth(RNG& rng, double depth);

        // Move on to next task. Returns `true` if there
        // actually was a next task to move on to.
        bool next_task();

        // Print the results vectors to the given stream.
        void print_results(std::ostream& out) const;

    /***** Private helper functions *****/
    private:

        // Do a single NS iteration
        void do_iteration(RNG& rng,
                          bool generate_new_particle=true);

        // Get the current depth
        double get_depth() const { return (double)iteration/particles.size(); }

        // Initialise all the particles from the prior
        void initialise(RNG& rng);

        // Find the worst particle
        size_t find_worst_particle() const;

        // Replace worst particle
        void replace_particle(RNG& rng, size_t which);

        // Combined scalar wrt results
        double combined_scalar(const std::vector<double>& ss) const;
};


/*************************************************/
/*           IMPLEMENTATIONS BELOW               */
/*************************************************/

template<class ParticleType>
Sampler<ParticleType>::Sampler(size_t num_particles,
                               size_t _mcmc_steps,
                               size_t _more_particles,
                               size_t _thin)
:particles(num_particles)
,scalars(num_particles)
,tiebreakers(num_particles)
,iteration(0)
,mcmc_steps(_mcmc_steps)
,more_particles(_more_particles)
,which_scalar(0)
,results(ParticleType::num_scalars)
,thin(_thin)
{
    sample_file.open("sample.csv", std::ios::out);
    sample_file << "iteration,logX,";
    sample_file << ParticleType::description() << std::endl;

    sample_info_file.open("sample_info.csv", std::ios::out);

    sample_info_file << "which_scalar,iteration,logX,";
    for(size_t i=0; i<ParticleType::num_scalars; ++i)
    {
        sample_info_file << "scalars[" << i << "]";
        if(i != ParticleType::num_scalars - 1)
            sample_info_file << ",";
    }
    sample_info_file << std::endl;
}

template<class ParticleType>
void Sampler<ParticleType>::run_to_depth(RNG& rng, double depth)
{
    initialise(rng);
    while(get_depth() < depth)
        do_iteration(rng);
}

template<class ParticleType>
bool Sampler<ParticleType>::next_task()
{
    if(which_scalar == ParticleType::num_scalars)
        return false;

    ++which_scalar;

    if(which_scalar == ParticleType::num_scalars)
    {
        particles.resize(more_particles*particles.size());
        scalars.resize(particles.size());
        tiebreakers.resize(particles.size());
    }

    return true;
}

template<class ParticleType>
void Sampler<ParticleType>::initialise(RNG& rng)
{
    std::cout << "Generating " << particles.size() << " ";
    std::cout << "particles from the prior..." << std::flush;

    for(size_t i=0; i<particles.size(); ++i)
    {
        particles[i].from_prior(rng);

        if(which_scalar < ParticleType::num_scalars)
            scalars[i] = particles[i].get_scalar(which_scalar);
        else
            scalars[i] = combined_scalar(get_scalars(particles[i]));

        tiebreakers[i] = rng.rand();
    }

    std::cout << "done.\n" << std::endl;
    iteration = 0;
}

template<class ParticleType>
void Sampler<ParticleType>::do_iteration(RNG& rng,
                                         bool generate_new_particle)
{
    ++iteration;

    // Worst particle
    size_t worst = find_worst_particle();

    // Print some messages to stdout
    std::cout << std::setprecision(12);
    std::cout << "Iteration " << iteration << ". ";
    std::cout << "log(X) ~= " << -(double)iteration/particles.size() << ", ";
    std::cout << "log(scalar" << which_scalar << ") = ";
    std::cout << scalars[worst] << '.';
    std::cout << std::endl;

    // Accumulate results
    if(which_scalar < ParticleType::num_scalars)
        results[which_scalar].push_back(scalars[worst]);

    // Write to files
    if((which_scalar == ParticleType::num_scalars) &&
       (iteration % thin == 0))
    {
        sample_file << iteration << ',';
        sample_file << (-(double)iteration/particles.size()) << ',';
        particles[worst].print(sample_file);
        sample_file << std::endl;
    }

    sample_info_file << std::setprecision(12);
    sample_info_file << which_scalar << ',';
    sample_info_file << iteration << ',';
    sample_info_file << (-(double)iteration/particles.size()) << ',';
    std::vector<double> ss = get_scalars(particles[worst]);
    for(size_t i=0; i<ss.size(); ++i)
    {
        sample_info_file << ss[i];
        if(i != ss.size() - 1)
            sample_info_file << ',';
    }
    sample_info_file << std::endl;

    // Generate replacement
    if(generate_new_particle)
        replace_particle(rng, worst);
}

template<class ParticleType>
void Sampler<ParticleType>::replace_particle(RNG& rng, size_t which)
{
    std::pair<double, double> threshold{scalars[which],
                                        tiebreakers[which]};

    std::cout << "    Generating replacement particle...";
    std::cout << std::flush;

    // Copy a survivor
    if(particles.size() > 1)
    {
        int copy;
        do
        {
           copy = rng.rand_int(particles.size());
        }while(copy == (int)which);

        particles[which] = particles[copy];
        scalars[which] = scalars[copy];
        tiebreakers[which] = tiebreakers[copy];
    }

    // Proposal stuff
    ParticleType proposal;
    double s_proposal, tb_proposal;

    // Acceptance counter
    unsigned int accepted = 0;

    // Do the MCMC
    for(size_t i=0; i<mcmc_steps; ++i)
    {
        // Generate proposal
        proposal = particles[which];
        double logH = proposal.perturb(rng);
        if(rng.rand() <= exp(logH))
        {
            if(which_scalar < ParticleType::num_scalars)
                s_proposal = proposal.get_scalar(which_scalar);
            else
                s_proposal = combined_scalar(get_scalars(proposal));

            tb_proposal = tiebreakers[which] + rng.randh();
            wrap(tb_proposal, 0.0, 1.0);
            std::pair<double, double> lltb_proposal{s_proposal,
                                                    tb_proposal};

            // Accept?
            if(threshold < lltb_proposal)
            {
                particles[which] = proposal;
                scalars[which] = s_proposal;
                tiebreakers[which] = tb_proposal;

                ++accepted;
            }
        }
    }

    std::cout << "done. ";
    std::cout << "Accepted " << accepted << '/';
    std::cout << mcmc_steps << " proposals.\n" << std::endl;
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

template<class ParticleType>
void Sampler<ParticleType>::print_results(std::ostream& out) const
{
    for(size_t i=0; i<results[0].size(); ++i)
    {
        for(size_t j=0; j<results.size(); ++j)
            out << results[j][i] << ' ';
        out << '\n';
    }
}

template<class ParticleType>
double Sampler<ParticleType>::combined_scalar
                    (const std::vector<double>& ss) const
{
    unsigned int result = 0;
    for(size_t i=0; i<ParticleType::num_scalars; ++i)
    {
        for(size_t j=0; j<results[i].size(); ++j)
        {
            if(ss[i] > results[i][j])
                ++result;
            else
                break;
        }
    }
    return static_cast<double>(result)/(particles.size() / more_particles);
}

} // namespace TwinPeaks

#endif

