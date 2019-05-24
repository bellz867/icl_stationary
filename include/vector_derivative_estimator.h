#ifndef VECTORDERIVATIVEESTIMATOR_H
#define VECTORDERIVATIVEESTIMATOR_H

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <ros/ros.h>

// LS estimator for a first order approximatoion of the derivative of a state vector wrt time, thanks Anup
struct VectorDerivativeEstimator
{
    bool firstUpdate;
    ros::Time tLast;
    Eigen::Matrix<float,6,1>  xHat;
    Eigen::Matrix<float,6,6> P;
    Eigen::Matrix<float,6,6> Q;
    Eigen::Matrix3f R;
    Eigen::Matrix<float,6,6> F;
    Eigen::Matrix<float,3,6> H;
    Eigen::Matrix<float,6,3> HT;

    VectorDerivativeEstimator();

    void initialize();

    Eigen::Vector3f update(Eigen::Vector3f newMeasure, ros::Time newTime);

    Eigen::Matrix<float,6,1> xDot(Eigen::Matrix<float,6,1> x);
};

#endif
