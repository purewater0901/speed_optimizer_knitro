#ifndef KNITROEXAMPLES_EXPOSITIONOPTIMIZER_H
#define KNITROEXAMPLES_EXPOSITIONOPTIMIZER_H

#include <vector>
#include <array>
#include <Eigen/Eigen>
#include "KTRSolver.h"
#include "KTRProblem.h"

class PositionOptimizer : public knitro::KTRProblem
{
public:
    PositionOptimizer(const int N,
                      const std::vector<double>& v,
                      const std::array<double, 3>& weight,
                      const double dt,
                      const double Smin,
                      const double Smax)
            //: KTRProblem(N+3, N, 0, 4*N+6), N_(N), v_(v), weight_(weight), dt_(dt), Smin_(Smin), Smax_(Smax), p_(0.0)
    : KTRProblem(N+3, N), N_(N), v_(v), weight_(weight), dt_(dt), Smin_(Smin), Smax_(Smax), p_(0.0)
    {
        setObjectiveProperties();
        setVariableProperties();
        setConstraintProperties();
        //setDerivativeProperties();

        calcHess();
        for(int i=0; i<N_; ++i)
            p_ += v_[i]*v_[i];
        p_ = p_ * weight_[0]*dt_;
    }

    double evaluateFC(const double* const x,
                      double* const c,
                      double* const objGrad,
                      double* const jac)
    {

        /* Constraints */
        //1. Linear Constraints
        for(int i=0; i<N_; ++i)
            c[i] = (x[i+1]-x[i])/dt_;


        /* Objective Function */
        /*
        Eigen::VectorXd s = Eigen::VectorXd::Zero(N_+3);
        for(int i=0; i<N_+3; ++i)
            s(i) = x[i];

        return (1.0/2.0)*s.dot(H_*s) + q_.dot(s) + p_;
         */


        std::vector<double> v;
        std::vector<double> acc;
        std::vector<double> jerk;
        v.resize(N_);
        acc.resize(N_);
        jerk.resize(N_);
        double Jv = 0.0;
        double Ja = 0.0;
        double Jj = 0.0;

        for(size_t i=0; i<N_; ++i)
        {
            v[i] = (x[i+1]-x[i])/dt_;
            acc[i] = (x[i+2]-2*x[i+1]+x[i])/std::pow(dt_,2);
            jerk[i] = (x[i+3]-3*x[i+2]+3*x[i+1]-x[i])/std::pow(dt_,3);
            Jv += std::pow((v[i]-v_[i]),2);
            Ja += std::pow(acc[i],2);
            Jj += std::pow(jerk[i],2);
        }

        return weight_[0]*Jv + weight_[1]*Ja + weight_[2]*Jj;
    }

    int evaluateGA(const double* const x, double* const objGrad, double* const jac)
    {
        Eigen::VectorXd s = Eigen::VectorXd::Zero(N_+3);
        for(int i=0; i<N_+3; ++i)
            s(i) = x[i];

        Eigen::VectorXd grad = H_*s + q_;

        for(int i=0; i<N_+3; ++i)
            objGrad[i] = grad(i);

        return 0;
    }

    int evaluateHess(const double* const x,
                     double objScaler,
                     const double* const lambda,
                     double* const hess)
    {
        for(int i=0; i<N_+3; ++i)
        {
            hess[i] = objScaler*H_(i,i);

            if(i<N_+2)
                hess[i+N_+3] = objScaler*H_(i,i+1);

            if(i<N_+1)
                hess[i+2*N_+5] = objScaler*H_(i,i+2);

            if(i<N_)
                hess[i+3*N_+6] = objScaler*H_(i, i+3);
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
        setVarLoBnds(0, Smin_);
        setVarUpBnds(0, Smin_);

        for(int i=1; i<N_+3; ++i)
        {
            setVarLoBnds(i, Smin_);
            setVarUpBnds(i, Smax_);
        }

        for(int i=0; i<N_+3; ++i)
            setXInitial(i, Smin_);
    }

    void setConstraintProperties()
    {
        for(size_t i=0; i<N_; ++i)
        {
            setConTypes(i, knitro::KTREnums::ConstraintType::ConLinear);
            setConLoBnds(i, 0.0);
            setConUpBnds(i, v_[i]);
        }
    }

    void setDerivativeProperties()
    {
        for(int i=0; i<N_+3; ++i)
        {
            setHessIndexRows(i, i);
            setHessIndexCols(i, i);

            if(i<N_+2)
            {
                setHessIndexRows(i+N_+3, i);
                setHessIndexCols(i+N_+3, i+1);
            }

            if(i<N_+1)
            {
                setHessIndexRows(i+2*N_+5, i);
                setHessIndexCols(i+2*N_+5, i+2);
            }

            if(i<N_)
            {
                setHessIndexRows(i+3*N_+6, i);
                setHessIndexCols(i+3*N_+6, i+3);
            }
        }
    }

    void calcHess()
    {
        Eigen::MatrixXd Hvel  = Eigen::MatrixXd::Zero(N_+3, N_+3);
        Eigen::MatrixXd Hacc  = Eigen::MatrixXd::Zero(N_+3, N_+3);
        Eigen::MatrixXd Hjerk = Eigen::MatrixXd::Zero(N_+3, N_+3);

        Eigen::VectorXd q = Eigen::VectorXd::Zero(N_+3);

        for(int i=0; i<N_+1; ++i)
        {
            if(i==0)
            {
                Hvel(0,0) = 1;
                Hvel(0,1) = -1;
                Hvel(1,0) = -1;

                q[0] = -v_[0];
            }
            else if(i==N_)
            {
                Hvel(N_, N_) = 1;

                q[N_] = v_[N_-1];
            }
            else
            {
                Hvel(i, i)   = 2;
                Hvel(i, i+1) = -1;
                Hvel(i+1, i) = -1;

                q[i] = v_[i-1] - v_[i];
            }
        }

        for(int i=0; i<N_; ++i)
        {
            if(i==0)
            {
                Hacc(0, 0)  = 1.0;
                Hacc(0, 1) = -2.0;
                Hacc(0, 2) = 1.0;
                Hacc(1, 0) = -2.0;
                Hacc(1, 1) = 1.0;
            }
            else if(i==1)
            {
                Hacc(1, 1) = 5.0;
                Hacc(1, 2) = -4.0;
                Hacc(1, 3) = 1.0;
                Hacc(2, 1) = -4.0;
                Hacc(3, 1) = 1.0;
            }
            else
            {
                Hacc(i,i) = 6;
                Hacc(i, i+1) = -4;
                Hacc(i, i+2) = 1;
                Hacc(i+1, i) = -4;
                Hacc(i+2, i) = 1;
            }
        }
        Hacc(N_,N_)     = 5;
        Hacc(N_,N_+1)   = -2;
        Hacc(N_+1,N_)   = -2;
        Hacc(N_+1,N_+1) = 1;

        for(int i=0; i<N_; ++i)
        {
            if(i==0)
            {
                Hjerk(i,i)=1;
                Hjerk(i,i+1)=-3;
                Hjerk(i,i+2)=3;
                Hjerk(i,i+3)=-1;
                Hjerk(i+1,i)=-3;
                Hjerk(i+2,i)=3;
                Hjerk(i+3,i)=-1;
            }
            else if(i==1)
            {
                Hjerk(i,i)=10;
                Hjerk(i,i+1)=-12;
                Hjerk(i,i+2)=6;
                Hjerk(i,i+3)=-1;
                Hjerk(i+1,i)=-12;
                Hjerk(i+2,i)=6;
                Hjerk(i+3,i)=-1;
            }
            else if(i==2)
            {
                Hjerk(i,i)=19;
                Hjerk(i,i+1)=-15;
                Hjerk(i,i+2)=6;
                Hjerk(i,i+3)=-1;
                Hjerk(i+1,i)=-15;
                Hjerk(i+2,i)=6;
                Hjerk(i+3,i)=-1;
            }
            else
            {
                Hjerk(i,i)=20;
                Hjerk(i,i+1)=-15;
                Hjerk(i,i+2)=6;
                Hjerk(i,i+3)=-1;
                Hjerk(i+1,i)=-15;
                Hjerk(i+2,i)=6;
                Hjerk(i+3,i)=-1;
            }
        }
        Hjerk(N_,N_)   = 19;
        Hjerk(N_,N_+1) = -12;
        Hjerk(N_,N_+2) = 3;
        Hjerk(N_+1,N_) = -12;
        Hjerk(N_+2,N_) = 3;
        Hjerk(N_+1,N_+1) = 10;
        Hjerk(N_+1,N_+2) = -3;
        Hjerk(N_+2,N_+1) = -3;
        Hjerk(N_+2,N_+2) = 1;

        H_ = (2*weight_[0]/dt_)*Hvel + (2*weight_[1]/std::pow(dt_,3))*Hacc + (2*weight_[2]/std::pow(dt_,5))*Hjerk;
        q_ = -2 * weight_[0] * q;
    }

    int N_;
    std::vector<double> v_;
    std::array<double, 3> weight_;
    double dt_;
    double Smax_;
    double Smin_;
    Eigen::MatrixXd H_;
    Eigen::VectorXd q_;
    double p_;
};


#endif //KNITROEXAMPLES_EXPOSITIONOPTIMIZER_H
