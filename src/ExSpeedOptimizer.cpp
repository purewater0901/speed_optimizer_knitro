#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "ExSpeedOptimizer.h"
#include "ExampleHelpers.h"

cv::Point2i cv_offset(float x, float y,int image_width=1000, int image_height=1000)
{
    cv::Point2i output;
    output.x = int(x * 20) + 10;
    output.y = int(y * 150) + 800;
    return output;
}

cv::Point2i cv_offset2(float x, float y)
{
    cv::Point2i output;
    output.x = int(x * 50) + 10;
    output.y = int(y * 50) + 400;
    return output;
}


int main()
{
    int n = 500;
    std::vector<double> v;
    v.resize(n);
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

    double ds = 0.1;
    double omega = 10;
    double a = 0.0;

    // Create a problem instance.
    ProblemQCQP instance(n, v, omega, ds, a);

    // Create a solver
    //knitro::KTRSolver solver(&instance, KTR_GRADOPT_FORWARD, KTR_HESSOPT_BFGS);
    knitro::KTRSolver solver(&instance);

    int solveStatus = solver.solve();

    printSolutionResults(solver, solveStatus);

    const std::vector<double>& point = solver.getXValues();

    cv::namedWindow("optimized_speed", cv::WINDOW_NORMAL);
    cv::Mat bg(1000, 2200, CV_8UC3, cv::Scalar(255, 255, 255));

    for(size_t i=0; i<v.size()-1; ++i)
    {
        cv::line(
                bg,
                cv_offset(i*ds, -v[i], bg.cols, bg.rows),
                cv_offset((i+1)*ds, -v[i+1], bg.cols, bg.rows),
                cv::Scalar(0, 0, 255),
                10);
    }
    for(size_t i=0; i<n-1; ++i)
    {
        cv::line(
                bg,
                cv_offset(i*ds, -std::sqrt(point[i]), bg.cols, bg.rows),
                cv_offset((i+1)*ds, -std::sqrt(point[i+1]), bg.cols, bg.rows),
                cv::Scalar(0, 0, 0),
                10);
    }

    cv::imshow("optimized_speed", bg);
    cv::waitKey();
    cv::destroyAllWindows();

    return 0;
}
