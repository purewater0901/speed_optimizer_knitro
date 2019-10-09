#ifndef KNITROEXAMPLES_EXTIMEOPTIMIZER_H
#define KNITROEXAMPLES_EXTIMEOPTIMIZER_H

#include <vector>
#include <array>
#include <Eigen/Eigen>
#include <cmath>
#include "KTRSolver.h"
#include "KTRProblem.h"

class TimeOptimizer : public knitro::KTRProblem
{
public:
    TimeOptimizer(const int N,
                  const std::vector<double>& Vr,
                  const std::vector<double>& Ar,
                  const std::vector<double>& Ac,
                  const std::array<double, 4>& weight,
                  const double ds,
                  const double a0)
                  : KTRProblem(3*N, 2*N), N_(N), weight_(weight), Vr_(Vr), Ar_(Ar), Ac_(Ac), epsilon_(1.0e-6), ds_(ds), a0_(a0)
    {
        setObjectiveProperties();
        setVariableProperties();
        setConstraintProperties();
    }

    double evaluateFC(const double* const x,
                      double* const c,
                      double* const objGrad,
                      double* const jac)
    {
        // Linear equality constraint.
        for(int i=0; i<N_-1; ++i)
            c[i] = (x[i+1]-x[i])/ds_ - 2*x[i+N_];
        c[N_-1] = 0.0;

        for(int i=0; i<N_; ++i)
            c[i+N_] = Ac_[i] + x[i+2*N_] - std::fabs(x[i+N_]);


        double Jt = 0.0;
        double Js = 0.0;
        double Jv = 0.0;
        double Ls = 0.0; //Longitudinal slack variable

        for(int i=0; i<N_-1; ++i)
        {
            Jt += (2*ds_)/(std::sqrt(x[i])+std::sqrt(x[i+1])+epsilon_);
            Js += std::pow((x[i+N_+1] - x[i+N_])/ds_, 2);
            Jv += std::pow(x[i] - Vr_[i]*Vr_[i], 2)*ds_;
            Ls += std::fabs(x[i+2*N_]);
        }
        Jv += std::pow(x[N_-1] - Vr_[N_-1]*Vr_[N_-1], 2)*ds_;
        Ls += std::fabs(x[3*N_-1]);

        return weight_[0]*Jt + weight_[1]*Js + weight_[2]*Jv + weight_[3]*Ls;
    }

private:
    void setObjectiveProperties()
    {
        //setObjType(knitro::KTREnums::ObjectiveType::ObjQuadratic);
        setObjType(knitro::KTREnums::ObjectiveType::ObjGeneral);
        setObjType(knitro::KTREnums::ObjectiveGoal::Minimize);
    }

    void setVariableProperties()
    {
        for(int i=0; i<N_; ++i)
        {
            if (i == 0)
            {
                setVarLoBnds(i, Vr_[i] * Vr_[i]);
                setVarUpBnds(i, Vr_[i] * Vr_[i]);

            } else
            {
                setVarLoBnds(i, 0);
                setVarUpBnds(i, Vr_[i] * Vr_[i]);
            }
        }

        for(int i=N_+1; i<2*N_-1; ++i)
        {
            setVarLoBnds(i, -Ar_[i-N_]);
            setVarUpBnds(i, Ar_[i-N_]);
        }
        setVarLoBnds(N_, a0_);
        setVarUpBnds(N_, a0_);
        setVarLoBnds(2*N_-1, 0.0);
        setVarUpBnds(2*N_-1, 0.0);

        // constraint for longitudinal slack variables
        for(int i=2*N_; i<3*N_; ++i)
        {
            setVarLoBnds(i, 0);
            setVarUpBnds(i, KTR_INFBOUND);
        }

        for(int i=0; i<3*N_; ++i)
            setXInitial(i, 0.0);
    }

    void setConstraintProperties()
    {
        for(int i=0; i<N_; ++i)
        {
            setConTypes(i, knitro::KTREnums::ConstraintType::ConLinear);
            setConLoBnds(i, 0.0);
            setConUpBnds(i, 0.0);

            setConTypes(i+N_, knitro::KTREnums::ConstraintType::ConGeneral);
            setConLoBnds(i+N_, 0.0);
            setConUpBnds(i+N_, KTR_INFBOUND);
        }
    }

    std::array<double,4> weight_;
    std::vector<double> Vr_;
    std::vector<double> Ar_;
    std::vector<double> Ac_;

    const double epsilon_;
    double ds_;
    double a0_;
    int N_;

};


#endif //KNITROEXAMPLES_EXTIMEOPTIMIZER_H
