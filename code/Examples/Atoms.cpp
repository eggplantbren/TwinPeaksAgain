#include "Atoms.h"
#include "../Utils.h"
#include <cmath>

using namespace std;

namespace TwinPeaks
{

Atoms::Atoms()
:x(num_atoms), y(num_atoms), z(num_atoms)
{

}

void Atoms::from_prior(RNG& rng)
{
    for(size_t i=0; i<x.size(); i++)
    {
        x[i] = L*rng.rand();
        y[i] = L*rng.rand();
        z[i] = L*rng.rand();
    }
}

double Atoms::get_scalar(size_t which_scalar) const
{
    double s = 0.0;

    if(which_scalar == 0)
    {
        for(int i=0; i<num_atoms; ++i)
        {
            for(int j=(i+1); j<num_atoms; ++j)
            {
                double rsq = pow(x[i] - x[j], 2) + pow(z[i] - z[j], 2);
                s += -4*(pow(1./rsq, 6) - 2.*pow(1./rsq, 3));
            }
        }
    }
    else
    {

        // Polymer potential
        double k = 36*pow(2., 2./3);
        double c = pow(2., 1./6);
        double r;
        // First sum
        for(int i=0; i<(num_atoms-1); ++i)
        {
            r = sqrt(pow(x[i] - x[i+1], 2) + pow(y[i] - y[i+1], 2)
                            + pow(z[i] - z[i+1], 2));
            s += -0.5*k*pow(r - c, 2);
        }

        // Second sum
        double dotprod;
        for(int i=1; i<(num_atoms-1); ++i)
        {
            dotprod = (x[i] - x[i-1])*(x[i+1] - x[i])
                        + (y[i] - y[i-1])*(y[i+1] - y[i])
                        + (z[i] - z[i-1])*(z[i+1] - z[i]);
            s += -0.5*k*pow(dotprod - 1., 2);
        }
    }

    return s;
}

double Atoms::perturb(RNG& rng)
{
    int reps = pow((double)num_atoms*3, rng.rand());
    for(int rep=0; rep<reps; ++rep)
    {
        int i = rng.rand_int(num_atoms);

        int which_coord = rng.rand_int(3);
        double* coord = nullptr;
        if(which_coord == 0)
            coord = &(x[i]);
        else if(which_coord == 1)
            coord = &(y[i]);
        else
            coord = &(z[i]);

        *coord += L*rng.randh();
        wrap(*coord, 0.0, L);
    }
    return 0.0;
}


void Atoms::print(std::ostream& out) const
{
    for(size_t i=0; i<x.size(); ++i)
        out << x[i] << ",";

    for(size_t i=0; i<y.size(); ++i)
        out << y[i] << ",";

    for(size_t i=0; i<z.size(); ++i)
    {
        out << z[i];
        if(i != z.size() - 1)
            out << ",";
    }
}

} // namespace TwinPeaks

