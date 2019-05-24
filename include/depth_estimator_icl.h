#ifndef DEPTHESTIMATORICL_H
#define DEPTHESTIMATORICL_H

#define _USE_MATH_DEFINES
#include <ros/ros.h>
#include <std_msgs/Time.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <deque>

#include <helper_functions.h>
#include <depth_estimator_ekf.h>
#include <depth_estimator_ls.h>
#include <depth_estimator_norm.h>

//depth estimator
struct DepthEstimatorICL
{
  std::deque<Eigen::Vector2f> zetaBuff;
  std::deque<Eigen::Vector2f> uvBuff;
  Eigen::Vector2f uvInt;
  std::deque<ros::Time> tBuff;
  std::deque<float> dtBuff;
  float yysum,yusum;
  Eigen::Vector3f uk;
  float dkHat;
  float dcHat;
  float dkcHat;
  float zmin;
  float zmax;
  float tau;
  bool firstzk;
  bool dkKnown;
  int numSaved;

  DepthEstimatorICL();
  void initialize(Eigen::Vector3f uInit, float zminInit, float zmaxInit, float tauInit);
  Eigen::Vector3f update(Eigen::Vector3f mc, Eigen::Vector3f tkc, Eigen::Matrix3f Rkc, Eigen::Vector3f v, Eigen::Vector3f pkc, ros::Time t, float dt, bool patternLost);
};

#endif
