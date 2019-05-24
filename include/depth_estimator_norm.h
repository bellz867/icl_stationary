#ifndef DEPTHESTIMATORNORM_H
#define DEPTHESTIMATORNORM_H

#define _USE_MATH_DEFINES
#include <ros/ros.h>
#include <std_msgs/Time.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <deque>

#include <helper_functions.h>

struct DepthEstimatorNorm
{
  std::deque<Eigen::Vector3f> mzetaBuff;
  std::deque<Eigen::Vector3f> psizetaBuff;
  std::deque<Eigen::Vector3f> vBuff;
  Eigen::Vector3f vInt;
  Eigen::Vector3f psizetaInt;
  std::deque<ros::Time> tBuff;
  std::deque<float> dtBuff;
  std::vector<float> zkBuff;
  Eigen::Vector3f mk,mc,pik;
  float dkcHat;
  float zkHat;
  float zcHat;
  float zmin;
  float zmax;
  float tau;
  bool firstzk;
  bool zkKnown;

	DepthEstimatorNorm();

	void initialize(Eigen::Vector3f mInit, float zminInit, float zmaxInit, float tauInit);

	Eigen::Vector3f update(Eigen::Vector3f mc, Eigen::Vector3f ukc, Eigen::Matrix3f Rkc, Eigen::Vector3f v, Eigen::Vector3f w, Eigen::Vector3f pkc, ros::Time t, float dt);
};

#endif
