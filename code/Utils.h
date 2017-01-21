#ifndef InfoNest_TwinPeaks
#define InfoNest_TwinPeaks

#include <tuple>

// Useful functions, some of which were copied from DNest4

namespace TwinPeaks
{

// This is non-standard, gcc supports it but mingw64 doesn't,
// so putting it here
#ifndef M_PI
    constexpr double M_PI = 3.141592653589793;
#endif

double mod(double y, double x);
int    mod(int y, int x);
void   wrap(double& x, double min, double max);

// Are both values in the first tuple above both in the second?
bool both_above(const std::tuple<double, double>& x,
                const std::tuple<double, double>& y);

} // namespace TwinPeaks

#endif

