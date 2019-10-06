#include <iostream>
#include <vector>
#include <array>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "ExPositionOptimizer.h"

int main()
{
    double Smin = 0.0;
    double Smax = 50.0;
    double T = 10;
    double dt = 0.1;
    int N = int(T/dt);

    std::array<double, 3> weight = {0.6, 0.2, 0.2};

    //set speed
    std::vector<double> v;
    v.resize(N);
    /*
    for(size_t i=0; i<v.size()-1; ++i)
        v[i] = 5.0;
    v[v.size()-1] = 0.0;
     */

    for(size_t i=1; i<v.size()/3; ++i)
        v[i] = 3.0;

    for(size_t i=v.size()/3; i<v.size()*2/3; ++i)
        v[i] = 4.0;

    for(size_t i=v.size()*2/3; i<v.size()*3/4; ++i)
        v[i] = 2.0;

    for(size_t i=v.size()*3/4; i<v.size(); ++i)
        v[i] = 1.0;

    v[0] = 1.0;
    v[v.size()-1] = 0.0;

    PositionOptimizer instance(N, v, weight, dt, Smin, Smax);
    knitro::KTRSolver solver(&instance);
    int solveStatus = solver.solve();

    const std::vector<double>& point = solver.getXValues();
    for(size_t i=0; i<point.size(); ++i)
    {
        std::cout << point[i] << std::endl;
    }


    return 0;
}
