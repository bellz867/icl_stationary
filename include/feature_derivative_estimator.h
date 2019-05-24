#ifndef FEATUREDERIVATIVEESTIMATOR_H
#define FEATUREDERIVATIVEESTIMATOR_H

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <deque>
#include <ros/ros.h>

// LS estimator for a first order approximatoion of the derivative of a state vector wrt time, thanks Anup
struct FeatureDerivativeEstimator
{
    bool firstUpdate;
    ros::Time tLast;
    Eigen::Vector4f xHat;
    Eigen::Matrix4f P;
    Eigen::Matrix4f Q;
    Eigen::Matrix2f R;
    Eigen::Matrix4f F;
    Eigen::Matrix<float,2,4> H;
    Eigen::Matrix<float,4,2> HT;

    FeatureDerivativeEstimator();

    void initialize(float maxdtInit, int bufferSize);

    Eigen::Vector2f update(Eigen::Vector2f newMeasure, ros::Time newTime);

    Eigen::Vector4f xDot(Eigen::Vector4f x);
};

#endif
