#include "Context.h"

namespace TwinPeaks
{

Context::Context()
{

}

void Context::add_rectangle(const std::tuple<double, double>& new_rectangle)
{
    rectangles.push_front(new_rectangle);
}

bool Context::is_okay(const std::tuple<double, double>& point) const
{
    double x, y, x_rect, y_rect;

    // Unpack the given point
    std::tie(x, y) = point;

    for(auto it = rectangles.begin(); it != rectangles.end(); ++it)
    {
        // Unpack rectangle coordinates
        std::tie(x_rect, y_rect) = *it;

        if(x < x_rect || y < y_rect)
            return false;
    }

    return true;
}

} // namespace TwinPeaks

