#pragma once
#include "KTRSolver.h"
#include "KTRProblem.h"
#include <vector>

class SpeedOptimizer : public knitro::KTRProblem
{
public:
    SpeedOptimizer(const int N,
                   const std::vector<double>& v,
                   const double omega,
                   const double ds,
                   const double a0)
                   : KTRProblem(2*N, N, 3*(N-1), 3*N-1), N_(N), v_(v), omega_(omega), ds_(ds), a0_(a0)
    {
        setObjectiveProperties();
        setVariableProperties();
        setConstraintProperties();
        setDerivativeProperties();
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

        double Js = 0.0;
        double Jv = 0.0;

        for(int i=0; i<N_-1; ++i)
        {
            Js += std::pow((x[i+N_+1] - x[i+N_]), 2)/ds_;
            Jv += std::pow(x[i] - v_[i]*v_[i], 2)*ds_;
        }
        Jv += std::pow(x[N_-1] - v_[N_-1]*v_[N_-1], 2)*ds_;

        return omega_*Js + Jv;
    }

    int evaluateGA(const double* const x, double* const objGrad, double* const jac)
    {
        for(int i=0; i<N_; ++i)
        {
            objGrad[i] =  (2*ds_)*(x[i]-v_[i]*v_[i]);

            if(i==0)
                objGrad[N_] = (-2*omega_/ds_)*x[N_];
            else if(i==N_-1)
                objGrad[2*N_-1] = 0.0;
            else
                objGrad[i+N_] = (2*omega_/ds_)*(2*x[i+N_]-x[i+N_-1]-x[i+N_+1]);
        }

        for(int i=0; i<N_-1; ++i)
        {
            jac[i]      = -1/ds_;
            jac[i+(N_-1)]   = 1/ds_;
            jac[i+2*(N_-1)] = -2;
        }

        return 0;
    }

    int evaluateHess(const double* const x,
                     double objScaler,
                     const double* const lambda,
                     double* const hess)
    {
        for(int i=0; i<N_; ++i)
        {
            hess[i] = 2 * ds_ * objScaler;

            if(i==0 || i==N_-1)
                hess[i+N_] = (2*omega_*objScaler)/ds_;
            else
                hess[i+N_] = (4*omega_*objScaler)/ds_;

            if(i<N_-1)
                hess[i+2*N_] = (-2*omega_*objScaler)/ds_;
        }

        return 0;
    }

private:

    void setObjectiveProperties()
    {
        setObjType(knitro::KTREnums::ObjectiveType::ObjQuadratic);
        setObjGoal(knitro::KTREnums::ObjectiveGoal::Minimize);
    }

    void setVariableProperties()
    {
        for(int i=0; i<N_; ++i)
        {
            if(i==0 || i==N_-1)
            {
                setVarLoBnds(i, v_[i]*v_[i]);
                setVarUpBnds(i, v_[i]*v_[i]);

            }
            else
            {
                setVarLoBnds(i, 0);
                setVarUpBnds(i, v_[i]*v_[i]);
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

    void setDerivativeProperties()
    {
        for(int i=0; i<N_-1; ++i)
        {
            setJacIndexCons(i, i);
            setJacIndexVars(i, i);

            setJacIndexCons(i+N_-1, i);
            setJacIndexVars(i+N_-1, i+1);

            setJacIndexCons(i+2*(N_-1), i);
            setJacIndexVars(i+2*(N_-1), i+N_);
        }

        for(int i=0; i<N_; ++i)
        {
            setHessIndexRows(i, i);
            setHessIndexCols(i, i);

            setHessIndexRows(i+N_, i+N_);
            setHessIndexCols(i+N_, i+N_);

            if(i<N_-1)
            {
                setHessIndexRows(i+2*N_, i+N_);
                setHessIndexCols(i+2*N_, i+N_+1);
            }
        }
    }

    int N_;
    std::vector<double> v_;
    double omega_;
    double ds_;
    double a0_;

};
