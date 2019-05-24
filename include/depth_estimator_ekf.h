#ifndef DEPTHESTIMATOREKF_H
#define DEPTHESTIMATOREKF_H

#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>

// Kalman Filter that uses a LS estimator for input estimation
struct DepthEstimatorEKF
{
	Eigen::Matrix3f Q;
	Eigen::Matrix3f P;
	Eigen::Matrix<float,2,3> H;
	Eigen::Matrix<float,3,2> HT;
	Eigen::Matrix2f R;
	float zmin,zmax;
	Eigen::Vector3f xHat;//estimate of the states

	DepthEstimatorEKF();

	void initialize(Eigen::Vector2f m, float zminInit, float zmaxInit);

	float update(Eigen::Vector2f m, Eigen::Vector3f v, Eigen::Vector3f w, float dt);

	Eigen::Vector3f getxDot(float mx, float my, float mz, float vx, float vy, float vz, float wx, float wy, float wz);

	Eigen::Matrix3f getF(float mx, float my, float mz, float vx, float vy, float vz, float wx, float wy, float wz);
};

#endif
