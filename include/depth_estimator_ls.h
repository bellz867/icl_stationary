#ifndef DEPTHESTIMATORLS_H
#define DEPTHESTIMATORLS_H

#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>
#include <feature_derivative_estimator.h>

// Kalman Filter that uses a LS estimator for input estimation
struct DepthEstimatorLS
{
	float zmin,zmax;
	float mzHat;//estimate of the states
	FeatureDerivativeEstimator mDotEsimator;//estimates for the derivatives for input into the kalman filter

	DepthEstimatorLS();

	void initialize(Eigen::Vector2f m, float zminInit, float zmaxInit, ros::Time t);

	float update(Eigen::Vector2f m, Eigen::Vector3f v, Eigen::Vector3f w, ros::Time t, float dt);

	float getmzDot(float mx, float my, float mz, float vz, float wx, float wy);

	Eigen::Vector2f getzeta(float mx, float my, float wx, float wy, float wz);

	Eigen::Vector2f getrho(float mx, float my, float vx, float vy, float vz);
};

#endif
