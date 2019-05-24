#include <data_save.h>

DataSave::DataSave()
{}

DataSave::DataSave(float timeInit, Eigen::Vector3f pkcInit, Eigen::Vector4f qkcInit, Eigen::Vector3f pkcICLInit, Eigen::Vector4f qkcICLInit, Eigen::Vector3f pikInit,
                  Eigen::Vector3f pikICLInit, Eigen::Vector3f picInit, Eigen::Vector3f picICLInit, Eigen::Vector3f picEKFInit, Eigen::Vector3f picLSInit, Eigen::Vector3f picICLExtInit,
                  Eigen::Vector3f vInit, Eigen::Vector3f wInit)
{
  time = timeInit;
  pkc = pkcInit;
  qkc = qkcInit;
  pkcICL = pkcICLInit;
  qkcICL = qkcICLInit;
  pik = pikInit;
  pikICL = pikICLInit;
  pic = picInit;
  picICL = picICLInit;
  picEKF = picEKFInit;
  picLS = picLSInit;
  picICLExt = picICLExtInit;
  v = vInit;
  w = wInit;

  std::cout << "\n picICLExt \n" << picICLExt << std::endl;
}
