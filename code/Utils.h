#ifndef InfoNest_TwinPeaks
#define InfoNest_TwinPeaks

#include <algorithm>
#include <tuple>
#include <vector>

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

// Argsort from
// http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
std::vector<size_t> argsort(const std::vector<T>& v)
{
	// initialize original index locations
	std::vector<size_t> idx(v.size());
	for(size_t i=0; i<idx.size(); i++)
		idx[i] = i;

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

	return idx;
}

} // namespace TwinPeaks

#endif

