#ifndef DEPTHESTIMATOR_H
#define DEPTHESTIMATOR_H

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
#include <depth_estimator_icl.h>
#include <depth_estimator_icl_ext.h>

//depth estimator
struct DepthEstimator
{
  int depthInd;
  std::deque<Eigen::Vector2f> zetaBuff;
  std::deque<Eigen::Vector2f> uvBuff;
  Eigen::Vector2f uvInt;
  std::deque<ros::Time> tBuff;
  std::deque<float> dtBuff;
  std::vector<float> dkBuff;
  Eigen::Vector3f mk,mc,pik,uk,uc;
  float dkHat;
  float dcHat;
  float dkcHat;
  float zcHatEKF;
  float zcHatLS;
  float zkHatNorm;
  float zcHatNorm;
  float dkcHatNorm;
  float dkHatICL;
  float dcHatICL;
  float dkcHatICL;
  float dkHatICLExt;
  float dcHatICLExt;
  float dkcHatICLExt;
  ros::Time lastt;
  ros::Time startt;
  float zmin;
  float zmax;
  float tau;
  bool firstzk;
  bool dkKnown;
  DepthEstimatorEKF depthEstimatorEKF;
  DepthEstimatorLS depthEstimatorLS;
  DepthEstimatorNorm depthEstimatorNorm;
  DepthEstimatorICL depthEstimatorICL;
  DepthEstimatorICLExt depthEstimatorICLExt;

  DepthEstimator();
  DepthEstimator(int depthIndInit, Eigen::Vector3f mInit, ros::Time t, float zminInit, float zmaxInit, Eigen::Vector3f pic, float tauInit);
  float update(Eigen::Matrix3f H, Eigen::Vector3f mcMeas, Eigen::RowVector3f nkT, Eigen::Vector3f tkc, Eigen::Matrix3f Rkc, Eigen::Vector3f v, Eigen::Vector3f w, ros::Time t, Eigen::Vector3f pkc, Eigen::Vector4f qkc, Eigen::Vector3f pic, bool patternLost);
  void predict(Eigen::Matrix3f Rkc, Eigen::Vector3f v, Eigen::Vector3f w, ros::Time t, Eigen::Vector3f pkc, Eigen::Vector4f qkc);
};

#endif
