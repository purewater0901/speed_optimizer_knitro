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
                    const std::array<double, 3>& weight,
                    const double ds,
                    const double a0)
                  : KTRProblem(2*N, N), N_(N), weight_(weight), Vr_(Vr), epsilon_(1.0e-6), ds_(ds), a0_(a0)
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

        double Jt = 0.0;
        double Js = 0.0;
        double Jv = 0.0;

        for(int i=0; i<N_-1; ++i)
        {
            Jt += (2*ds_)/(std::sqrt(x[i])+std::sqrt(x[i+1])+epsilon_);
            Js += std::pow((x[i+N_+1] - x[i+N_])/ds_, 2);
            Jv += std::pow(x[i] - Vr_[i]*Vr_[i], 2)*ds_;
        }
        Jv += std::pow(x[N_-1] - Vr_[N_-1]*Vr_[N_-1], 2)*ds_;

        return weight_[0]*Jt + weight_[1]*Js + weight_[2]*Jv;
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
            if (i == 0 || i == N_ - 1)
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
            setVarLoBnds(i, -KTR_INFBOUND);
            setVarUpBnds(i, KTR_INFBOUND);
        }
        setVarLoBnds(N_, a0_);
        setVarUpBnds(N_, a0_);
        setVarLoBnds(2*N_-1, 0.0);
        setVarUpBnds(2*N_-1, 0.0);

        for(int i=0; i<2*N_; ++i)
            setXInitial(i, 0.0);
    }

    void setConstraintProperties()
    {
        for(int i=0; i<N_; ++i)
        {
            setConTypes(i, knitro::KTREnums::ConstraintType::ConLinear);
            setConLoBnds(i, 0.0);
            setConUpBnds(i, 0.0);
        }
    }

    std::array<double,3> weight_;
    std::vector<double> Vr_;
    const double epsilon_;
    double ds_;
    double a0_;
    int N_;

};


#endif //KNITROEXAMPLES_EXTIMEOPTIMIZER_H
