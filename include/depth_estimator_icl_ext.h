#ifndef DEPTHESTIMATORICLExt_H
#define DEPTHESTIMATORICLExt_H

#define _USE_MATH_DEFINES
#include <ros/ros.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <deque>
#include <helper_functions.h>
#include <vector_derivative_estimator.h>

//depth estimator
struct DepthEstimatorICLExt
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
  VectorDerivativeEstimator uDotEstimator;

  DepthEstimatorICLExt();
  void initialize(Eigen::Vector3f uInit, float zminInit, float zmaxInit, float tauInit, ros::Time t);
  Eigen::Vector3f update(Eigen::Vector3f mc, Eigen::Vector3f tkc, Eigen::Matrix3f Rkc, Eigen::Vector3f v, Eigen::Vector3f w, Eigen::Vector3f pkc, ros::Time t, float dt, bool patternLost);
};

#endif
