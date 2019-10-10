#ifndef KNITROEXAMPLES_REFERENCEPATH_H
#define KNITROEXAMPLES_REFERENCEPATH_H

#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include "ReferencePath/cubic_spline.h"
#include "Waypoint/Waypoints.h"

class ReferencePath
{
public:
    ReferencePath(const std::vector<double>& x, const std::vector<double>& y, const double ds) : spline_(x, y)
    {
        x_.reserve(spline_.s.size());
        y_.reserve(spline_.s.size());
        yaw_.reserve(spline_.s.size());

        for(float s=0; s<spline_.s.back(); s+=ds)
        {
            std::array<double, 2> point = spline_.calc_position(s);

            x_.push_back(point[0]);
            y_.push_back(point[1]);
            yaw_.push_back(spline_.calc_yaw(s));
        }
    }

    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> yaw_;
    Spline2D spline_;
};


#endif //KNITROEXAMPLES_REFERENCEPATH_H
