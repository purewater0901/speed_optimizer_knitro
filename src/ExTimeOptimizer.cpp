#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <fstream>
#include "ExTimeOptimizer.h"

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
    int N = 100;
    std::string filename = "../position_result.csv";

    std::vector<double> Vr;
    std::vector<double> Ar;  //acceleration restriction
    std::vector<double> Ac;  //comfort acceleration restriction
    Vr.resize(N);
    Ar.resize(N);
    Ac.resize(N);

    for(size_t i=1; i<Vr.size()/3; ++i)
        Vr[i] = 3.0;

    for(size_t i=Vr.size()/3; i<Vr.size(); ++i)
        Vr[i] = 0.0;

    Vr[0] = 2.0;

    for(int i=0; i<Ar.size(); ++i)
    {
        Ar[i] = 1.2;
        Ac[i] = 0.8;
    }

    double ds = 0.1;
    double a0 = 0.0;

    std::array<double, 4> weight = {1.0, 0.025, 0, 0.02};

    // Create a problem instance.
    TimeOptimizer instance(N, Vr, Ar, Ac, weight, ds, a0);

    // Create a solver
    knitro::KTRSolver solver(&instance, KTR_GRADOPT_FORWARD, KTR_HESSOPT_BFGS);
    //knitro::KTRSolver solver(&instance);
    //solver.setParam("outlev", 0);

    int solveStatus = solver.solve();

    const std::vector<double>& point = solver.getXValues();

    /*  Write to file*/
    std::ofstream writing_file;
    writing_file.open(filename, std::ios::out);
    writing_file << "position" << "," << "reference_speed" <<  "," << "result_speed" <<  ","
                 << "acceleration" << "," << "longitudinal_slack_variable"<<std::endl;
    for(int i=0; i<N; ++i)
        writing_file << i*ds<< "," << Vr[i] << "," << std::sqrt(point[i]) << "," << point[i+N] << ","
                     << point[i+2*N] << std::endl;
    writing_file.close();

    cv::namedWindow("optimized_speed", cv::WINDOW_NORMAL);
    cv::Mat bg(1000, 2200, CV_8UC3, cv::Scalar(255, 255, 255));

    for(size_t i=0; i<Vr.size()-1; ++i)
    {
        cv::line(
                bg,
                cv_offset(i*ds, -Vr[i], bg.cols, bg.rows),
                cv_offset((i+1)*ds, -Vr[i+1], bg.cols, bg.rows),
                cv::Scalar(0, 0, 255),
                10);
    }
    for(size_t i=0; i<N-1; ++i)
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

