#ifndef KNITROEXAMPLES_WAYPOINTS_H
#define KNITROEXAMPLES_WAYPOINTS_H

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

class Waypoints
{
public:
    Waypoints(const double predictDistance, const std::string& fileName);

    void readCsv(const std::string& fileName);
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> yaw_;

private:
    std::vector<std::string> split(std::string& input, char delimiter);
    double predictDistance_;

};


#endif //KNITROEXAMPLES_WAYPOINTS_H
