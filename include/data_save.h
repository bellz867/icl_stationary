#ifndef DATASAVE_H
#define DATASAVE_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

struct DataSave
{
  float time;
  Eigen::Vector3f pkc;
  Eigen::Vector4f qkc;
  Eigen::Vector3f pkcICL;
  Eigen::Vector4f qkcICL;
  Eigen::Vector3f pik;
  Eigen::Vector3f pikICL;
  Eigen::Vector3f pic;
  Eigen::Vector3f picICL;
  Eigen::Vector3f picEKF;
  Eigen::Vector3f picLS;
  Eigen::Vector3f picICLExt;
  Eigen::Vector3f v;
  Eigen::Vector3f w;

  DataSave();
  DataSave(float timeInit, Eigen::Vector3f pkcInit, Eigen::Vector4f qkcInit, Eigen::Vector3f pkcICLInit, Eigen::Vector4f qkcICLInit, Eigen::Vector3f pikInit,
           Eigen::Vector3f pikICLInit, Eigen::Vector3f picInit, Eigen::Vector3f picICLInit, Eigen::Vector3f picEKFInit, Eigen::Vector3f picLSInit, Eigen::Vector3f picICLExtInit,
           Eigen::Vector3f vInit, Eigen::Vector3f wInit);
};

#endif
