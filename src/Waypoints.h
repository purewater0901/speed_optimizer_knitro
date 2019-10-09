#ifndef KNITROEXAMPLES_WAYPOINTS_H
#define KNITROEXAMPLES_WAYPOINTS_H

#include <vector>
#include <iostream>

class Waypoints
{
public:
    Waypoints(const std::vector<double>& x,
              const std::vector<double>& y,
              const std::vector<double>& yaw)
              : x_(x), y_(y), yaw_(yaw){}

    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> yaw_;
};


#endif //KNITROEXAMPLES_WAYPOINTS_H
