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
    double predictDistance = 20.0;
    double ds = 0.2;
    std::string waypointFilename = "../data/saved_waypoints.csv";
    std::string filename = "../position_result.csv";

    /* reading file */
    Waypoints waypoints(predictDistance, waypointFilename);
    ReferencePath referencePath(waypoints.x_, waypoints.y_, ds);
    int N = referencePath.x_.size();

    // another information
    double m = 1500;
    double V0 = 1.0;
    double a0 = 0.0;
    double mu = 0.8;

    // speed restriction and acceleration restriction
    std::vector<double> Vr;
    std::vector<double> Vd;
    std::vector<double> Arlon;  //acceleration longitudinal restriction
    std::vector<double> Arlat;  //acceleration lateral restriction
    std::vector<double> Aclon;  //comfort longitudinal acceleration restriction
    std::vector<double> Aclat;  //comfort lateral acceleration restriction
    Vr.resize(N);
    Vd.resize(N);
    Arlon.resize(N);
    Arlat.resize(N);
    Aclon.resize(N);
    Aclat.resize(N);

    for(size_t i=0; i<Vr.size()/2; ++i)
    {
        Vr[i] = 2.0;
        Vd[i] = 1.5;
    }

    for(size_t i=Vr.size()/2; i<Vr.size(); ++i)
    {
        Vr[i] = 4.0;
        Vd[i] = 3.5;
    }

    for(size_t i=0; i<N; ++i)
    {
        Arlon[i] = 0.5*mu*9.83;
        Arlat[i] = 0.5*mu*9.83;
        Aclon[i] = 0.4*mu*9.83;
        Aclat[i] = 0.4*mu*9.83;
    }


    std::array<double, 5> weight = {1.0, 5.0, 0.0, 10, 10};

    // Create a problem instance.
    TimeOptimizer instance(N, Vr, Vd, Arlon, Arlat, Aclon, Aclat, weight, referencePath, m, ds, V0, a0, mu);

    // Create a solver
    knitro::KTRSolver solver(&instance, KTR_GRADOPT_FORWARD, KTR_HESSOPT_BFGS);
    //knitro::KTRSolver solver(&instance);
    //solver.setParam("outlev", 0);

    int solveStatus = solver.solve();

    const std::vector<double>& point = solver.getXValues();

    /*  Write to file*/
    std::ofstream writing_file;
    writing_file.open(filename, std::ios::out);
    writing_file << "position" << ","
                 << "reference_speed" <<  "," << "desired_speed" << ","
                 << "result_speed" <<  ","
                 << "longitudinal acceleration" << "," << "lateral acceleration" << ","
                 << "longitudinal_slack_variable" << "," << "lateral_slack_variable" <<std::endl;
    for(int i=0; i<N; ++i)
        writing_file << i*ds<< ","
                     << Vr[i] << "," << Vd[i] << ","
                     << std::sqrt(point[i]) << ","
                     << point[i+N] << "," << point[i+5*N]/m << ","
                     << point[i+2*N] << "," << point[i+3*N]<< std::endl;
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
        cv::line(
                bg,
                cv_offset(i*ds, -Vd[i], bg.cols, bg.rows),
                cv_offset((i+1)*ds, -Vd[i+1], bg.cols, bg.rows),
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

