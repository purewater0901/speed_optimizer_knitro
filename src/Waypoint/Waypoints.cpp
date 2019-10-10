#include "Waypoint/Waypoints.h"

Waypoints::Waypoints(const double predictDistance, const std::string &fileName)
    : predictDistance_(predictDistance)
{
    readCsv(fileName);
}

void Waypoints::readCsv(const std::string& fileName)
{
    std::ifstream ifs(fileName);
    if(!ifs)
        return;

    std::string line;
    int count = 0;
    double distance = 0.0;
    while(getline(ifs, line))
    {
        if(count==0)
        {
            count++;
            continue;
        }
        std::vector<std::string> strvec = split(line,',');

        x_.push_back(std::stod(strvec.at(0)));
        y_.push_back(std::stod(strvec.at(1)));
        yaw_.push_back(std::stod(strvec.at(3)));

        if(x_.size()>1)
            distance += std::sqrt(std::pow(x_[x_.size()-1]-x_[x_.size()-2],2)
                                + std::pow(y_[y_.size()-1]-y_[y_.size()-2],2));

        if(distance>predictDistance_)
            return;
        count++;
    }

}

std::vector<std::string> Waypoints::split(std::string& input, char delimiter)
{
    std::istringstream stream(input);
    std::string field;
    std::vector<std::string> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}