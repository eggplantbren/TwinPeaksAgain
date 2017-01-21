#include "Context.h"

namespace TwinPeaks
{

Context::Context()
{

}

void Context::add_point(const std::tuple<double, double>& new_point)
{
    points.push_front(new_point);
}

size_t Context::calculate_ucc(const std::tuple<double, double>& x) const
{
    size_t ucc = 0;

    for(auto it = points.begin(); it != points.end(); ++it)
    {
        if(both_above(*it, x))
            ++ucc;
    }

    return ucc;
}

} // namespace TwinPeaks

