#ifndef FEATUREDERIVATIVEESTIMATORLS_H
#define FEATUREDERIVATIVEESTIMATORLS_H

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <deque>
#include <ros/ros.h>

// LS estimator for a first order approximatoion of the derivative of a state vector wrt time, thanks Anup
struct FeatureDerivativeEstimator
{
    bool bufferFull; //Estimation will start after buffer is full for first time
    std::deque<ros::Time> timeBuff; //ring buffer for time data
    std::deque<Eigen::Vector2f> stateBuff; //ring buffer for position data
    int bufferSize;
    float maxdt;
    Eigen::MatrixXf A;

    FeatureDerivativeEstimator();

    void initialize(float maxdtInit, int bufferSize);

    Eigen::Vector2f update(Eigen::Vector2f newMeasure, ros::Time newTime);
};

#endif
