#ifndef TwinPeaks_Context
#define TwinPeaks_Context

// Includes
#include <list>
#include <tuple>
#include "Utils.h"

namespace TwinPeaks
{

/*
* A class that manages the context points.
*/
class Context
{
    private:
        // The context points
        // Using a list; recently added points stored at the front.
        std::list<std::tuple<double, double>> points;

    public:
        // Constructor; starts with zero points
        Context();

        // Add the given point to the context.
        void add_point(const std::tuple<double, double>& new_point);

        // Calculate an UCC
        size_t calculate_ucc(const std::tuple<double, double>& x) const;
};

} // namespace TwinPeaks

#endif

