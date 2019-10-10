#ifndef KNITROEXAMPLES_EXTIMEOPTIMIZER_H
#define KNITROEXAMPLES_EXTIMEOPTIMIZER_H

#include <vector>
#include <array>
#include <Eigen/Eigen>
#include <cmath>
#include "KTRSolver.h"
#include "KTRProblem.h"
#include "ReferencePath/ReferencePath.h"
#include "ReferencePath/ReferencePath.h"

class TimeOptimizer : public knitro::KTRProblem
{
public:
    TimeOptimizer(const int N,
                  const std::vector<double>& Vr,
                  const std::vector<double>& Arlon,
                  const std::vector<double>& Arlat,
                  const std::vector<double>& Aclon,
                  const std::vector<double>& Aclat,
                  const std::array<double, 5>& weight,
                  const ReferencePath& referencePath,
                  const double m,
                  const double ds,
                  const double a0,
                  const double mu)
                  : KTRProblem(6*N, 5*N), N_(N), weight_(weight), referencePath_(referencePath),
                    Vr_(Vr), Arlon_(Arlon), Arlat_(Arlat), Aclon_(Aclon), Aclat_(Aclat),
                    m_(m), epsilon_(1.0e-6), ds_(ds), a0_(a0), mu_(mu), g_(9.80665)
    {
        /*
         * @ variables
         * b --- speed 0~N
         * a --- acceleration N~2N
         * slack lon --- longitudinal slack variables 2N~3N
         * slack lat --- lateral slack variables 3N~4N
         * input lon --- longitudinal input 4N~5N
         * input lat --- lateral input 5N~6N
         */
        setObjectiveProperties();
        setVariableProperties();
        setConstraintProperties();
    }

    double evaluateFC(const double* const x,
                      double* const c,
                      double* const objGrad,
                      double* const jac)
    {
        /* constraint1 speed and acceleration constraints*/
        // Linear equality constraint.
        for(int i=0; i<N_-1; ++i)
            c[i] = (x[i+1]-x[i])/ds_ - 2*x[i+N_];
        c[N_-1] = 0.0;

        /* constraint2 longitudinal and lateral acceleration*/
        for(int i=0; i<N_; ++i)
        {
            c[i+N_]   = Aclon_[i] + x[i+2*N_] - std::fabs(x[i+N_]);      //longitudinal
            c[i+2*N_] = Aclat_[i] + x[1+3*N_] - std::fabs(x[i+5*N_]/m_); //lateral
        }

        /* constraint3 friction circle */
        for(int i=0; i<N_; ++i)
        {
            c[i+3*N_] = mu_*m_*g_ - std::sqrt(x[i+4*N_]*x[i+4*N_] + x[i+5*N_]*x[i+5*N_]);
            c[i+4*N_] = m_*Arlon_[i] - x[i+4*N_];
        }

        double Jt = 0.0;
        double Js = 0.0;
        double Jv = 0.0;
        double LonSlack = 0.0; //Longitudinal slack variable
        double LatSlack = 0.0; //Longitudinal slack variable

        for(int i=0; i<N_-1; ++i)
        {
            Jt += (2*ds_)/(std::sqrt(x[i])+std::sqrt(x[i+1])+epsilon_);
            Js += std::pow((x[i+N_+1] - x[i+N_])/ds_, 2);
            Jv += std::pow(x[i] - Vr_[i]*Vr_[i], 2)*ds_;
            LonSlack += std::fabs(x[i+2*N_]);
            LatSlack += std::fabs(x[i+3*N_]);
        }
        Jv += std::pow(x[N_-1] - Vr_[N_-1]*Vr_[N_-1], 2)*ds_;
        LonSlack += std::fabs(x[3*N_-1]);
        LatSlack += std::fabs(x[4*N_-1]);

        return weight_[0]*Jt + weight_[1]*Js + weight_[2]*Jv + weight_[3]*LonSlack + weight_[4]*LatSlack;
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
        // constraint for speed variables
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

        // constraint for acceleration variables
        for(int i=N_+1; i<2*N_-1; ++i)
        {
            setVarLoBnds(i, -Arlon_[i-N_]);
            setVarUpBnds(i, Arlon_[i-N_]);
        }
        setVarLoBnds(N_, a0_);
        setVarUpBnds(N_, a0_);
        setVarLoBnds(2*N_-1, 0.0);
        setVarUpBnds(2*N_-1, 0.0);

        // constraint for longitudinal  and lateral slack variables
        for(int i=2*N_; i<4*N_; ++i)
        {
            setVarLoBnds(i, 0);
            setVarUpBnds(i, KTR_INFBOUND);
        }

        for(int i=0; i<6*N_; ++i)
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
            setConTypes(i+2*N_, knitro::KTREnums::ConstraintType::ConGeneral);
            setConLoBnds(i+2*N_, 0.0);
            setConUpBnds(i+2*N_, KTR_INFBOUND);

            setConTypes(i+3*N_, knitro::KTREnums::ConstraintType::ConGeneral);
            setConLoBnds(i+3*N_, 0.0);
            setConUpBnds(i+3*N_, KTR_INFBOUND);
            setConTypes(i+4*N_, knitro::KTREnums::ConstraintType::ConGeneral);
            setConLoBnds(i+4*N_, 0.0);
            setConUpBnds(i+4*N_, KTR_INFBOUND);
        }
    }

    std::array<double,5> weight_;
    std::vector<double> Vr_;
    std::vector<double> Arlon_;
    std::vector<double> Arlat_;
    std::vector<double> Aclon_;
    std::vector<double> Aclat_;
    ReferencePath referencePath_;

    const double epsilon_;
    double m_;
    double ds_;
    double a0_;
    int N_;
    double mu_;
    const double g_;
};


#endif //KNITROEXAMPLES_EXTIMEOPTIMIZER_H
