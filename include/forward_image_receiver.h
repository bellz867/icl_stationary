#ifndef FORWARDIMAGERECEIVER_H
#define FORWARDIMAGERECEIVER_H

#include <ctime>
#include <vector>
#include <mutex>

#include <ros/ros.h>
#include <ros/console.h>
#include <cv_bridge/cv_bridge.h>
#include <image_transport/image_transport.h>
#include <image_geometry/pinhole_camera_model.h>
#include <sensor_msgs/CameraInfo.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/video/tracking.hpp>

#include <patch_estimator.h>

struct ForwardImageReceiver
{
  //ros
  ros::NodeHandle nh;
  ros::Subscriber camInfoSub,odomSub;
  image_transport::ImageTransport it;
  image_transport::Subscriber imageSub;
  image_transport::Publisher undistortPub;

  //camera parameters
  cv::Mat camMat;//camera intrinsic matrix
  cv::Mat distCoeffs;//camera distortion coefficients
  cv::Mat map1,map2;//undistort maps
  bool gotCamParam;//indicate the camera intrinsic parameters are received
  std::string cameraName;//body name and camera name

  std::string markerName;//body name and camera name

  int imageWidth;
  int imageHeight;
  int minFeaturesBad;

  float fq,fp,ft,fn,fd;
  float zmin,zmax;

  bool saveExp;
  std::string expName;

  PatchEstimator* patch;

  ForwardImageReceiver();

  ~ForwardImageReceiver();

  //image callback
  void imageCB(const sensor_msgs::Image::ConstPtr& msg);

  // callback for getting camera intrinsic parameters
  void camInfoCB(const sensor_msgs::CameraInfo::ConstPtr& camInfoMsg);
};

#endif
