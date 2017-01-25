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
        // The forbidden rectangles
        // Using a list; most recent ones stored at the front.
        std::list<std::tuple<double, double>> rectangles;

    public:
        // Constructor; starts with zero points
        Context();

        // Add the given point to the context.
        void add_rectangle(const std::tuple<double, double>& new_rectangle);

        // Is the given point okay?
        bool is_okay(const std::tuple<double, double>& point) const;
};

} // namespace TwinPeaks

#endif

