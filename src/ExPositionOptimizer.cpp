#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "ExPositionOptimizer.h"

cv::Point2i cv_offset(float x, float y,int image_width=1000, int image_height=1000)
{
    cv::Point2i output;
    output.x = int(x * 20) + 100;
    output.y = int(y * 20) + 300;
    return output;
}

int main()
{
    std::string filename = "../timed_result.csv";
    double Smin = 0.0;
    double Smax = 70.0;
    double T = 20;
    double dt = 0.1;
    int N = int(T/dt);

    std::array<double, 3> weight = {10, 0.1, 0.1};

    //set speed
    std::vector<double> Vr;  //restricted speed
    std::vector<double> Vd;  //desired speed
    Vr.resize(N);
    Vd.resize(N);
    /*
    for(size_t i=0; i<v.size(); ++i)
        v[i] = 2.0;
        */

    for(size_t i=1; i<Vr.size()/3; ++i)
    {
        Vr[i] = 3.0;
        Vd[i] = Vr[i]-1.0;
    }

    for(size_t i=Vr.size()/3; i<Vr.size()*3/4; ++i)
    {
        Vr[i] = 4.0;
        Vd[i] = Vr[i]-1.0;
    }

    for(size_t i=Vr.size()*3/4; i<Vr.size(); ++i)
    {
        Vr[i] = 2.0;
        Vd[i] = Vr[i]-1.0;
    }

    PositionOptimizer instance(N, Vr, Vd, weight, dt, Smin, Smax);
    knitro::KTRSolver solver(&instance, KTR_GRADOPT_FORWARD, KTR_HESSOPT_BFGS);
    //knitro::KTRSolver solver(&instance);
    int solveStatus = solver.solve();

    const std::vector<double>& point = solver.getXValues();
    for(size_t i=0; i<point.size(); ++i)
        std::cout << point[i] << std::endl;

    std::vector<double> vresult;
    vresult.resize(N);
    for(size_t i=0; i<N-1; ++i)
        vresult[i] =  (point[i+1] - point[i])/dt;
    vresult[vresult.size()-1] = vresult[vresult.size()-2];

    std::cout << "------------------ Speed -----------------------" << std::endl;
    for(int i=0; i<N; ++i)
        std::cout << vresult[i] << std::endl;

    /*  Write to file*/
    std::ofstream writing_file;
    writing_file.open(filename, std::ios::out);
    for(int i=0; i<N; ++i)
        writing_file << i*dt << "," << vresult[i] << std::endl;


    /* open cv output */
    cv::namedWindow("timed_optimized_speed", cv::WINDOW_NORMAL);
    cv::Mat bg(1000, 2000, CV_8UC3, cv::Scalar(255, 255, 255));

    for(size_t i=0; i<Vr.size()-1; ++i)
    {
        cv::line(
                bg,
                cv_offset(i*dt, -Vr[i], bg.cols, bg.rows),
                cv_offset((i+1)*dt, -Vr[i+1], bg.cols, bg.rows),
                cv::Scalar(0, 0, 255),
                10);
    }
    for(size_t i=0; i<N-1; ++i)
    {
        cv::line(
                bg,
                cv_offset(i*dt, -vresult[i], bg.cols, bg.rows),
                cv_offset((i+1)*dt, -vresult[i+1], bg.cols, bg.rows),
                cv::Scalar(0, 0, 0),
                10);
    }

    cv::imshow("timed_optimized_speed", bg);
    cv::waitKey();
    cv::destroyAllWindows();

    return 0;
}
