#include "Utils.h"
#include <algorithm>
#include <cassert>
#include <cmath>

namespace TwinPeaks
{

double mod(double y, double x)
{
	assert(x > 0.0);
	return (y/x - floor(y/x))*x;
}

void wrap(double& x, double min, double max)
{
	x = TwinPeaks::mod(x - min, max - min) + min;
}

int mod(int y, int x)
{
	assert(x > 0);
	if(y >= 0)
		return y - (y/x)*x;
	else
		return (x-1) - TwinPeaks::mod(-y-1, x);
}


double logsumexp(const std::vector<double>& values)
{
    double top = *max_element(values.begin(), values.end());
    double result = 0.0;
    for(double v: values)
        result += exp(v - top);
    return top + log(result);
}

double logsumexp(double x, double y)
{
    return logsumexp(std::vector<double>{x, y});
}

} // namespace TwinPeaks

