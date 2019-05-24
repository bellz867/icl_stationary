#ifndef PATCHESTIMATOR_H
#define PATCHESTIMATOR_H

#include <vector>
#include <deque>
#include <mutex>

#include <ros/ros.h>
#include <ros/console.h>
#include <cv_bridge/cv_bridge.h>
#include <image_transport/image_transport.h>
#include <geometry_msgs/Point32.h>
#include <nav_msgs/Odometry.h>

#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/video/tracking.hpp>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <depth_estimator.h>
#include <helper_functions.h>
#include <data_save.h>

typedef pcl::PointCloud<pcl::PointXYZRGB> PointCloudRGB;
typedef pcl::PointCloud<pcl::PointXYZ> PointCloud;

struct PatchEstimator
{
	ros::NodeHandle nh;
	std::string cameraName,markerName;
	image_transport::ImageTransport it;
	image_transport::Subscriber imageSub;
	image_transport::Publisher imagePub;
	cv::Mat kimage,pimage;
  ros::Subscriber odomSub,roiSub,markerOdomSub;
	ros::Publisher wallPub,poseDeltaPub,roiPub,odomPub,pointCloudPub,odomDelayedPub;
	Eigen::Vector3f tkcHat;
	Eigen::Vector3f nkHat;
	Eigen::Vector4f qkcHat;
	Eigen::Vector3f pkcHat;
	Eigen::Vector3f pckHat;
	Eigen::Vector4f qckHat;
	float dkcHat,tau,dkHat;
	float fx,fy,cx,cy,zmin,zmax;
	std::vector<DepthEstimator*> depthEstimators;
	std::deque<nav_msgs::Odometry> odomSync;
	std::deque<nav_msgs::Odometry> markerOdomSync;
	std::mutex odomMutex,roiMutex,pubMutex,markerOdomMutex;
	ros::Time tLast;
	float pTau,qTau,tTau,nTau,dTau;
	nav_msgs::Odometry keyOdom,imageOdom,markerKeyOdom,markerOdom;
	cv::Mat camMat;
	Eigen::Matrix3f camMatf,camMatIf;
	bool firstOdomImageCB,firstImageCB,firstMarkerOdomImageCB;
	bool dkEstimated;
	int imageWidth,imageHeight;
	int minFeaturesBad;
	bool patchShutdown;
	bool saveExp;
	std::string expName;
	std::vector<DataSave*> data;
	ros::Time tStart;
	bool alwaysSearch;
	bool patternLost;

	~PatchEstimator();

	PatchEstimator();

	PatchEstimator(int imageWidth, int imageHeight, int minFeaturesBad, float fxInit, float fyInit, float cxInit, float cyInit, float zminInit, float zmaxInit, float fq, float fp, float ft, float fn, float fd, std::string cameraNameInit, std::string markerNameInit, float tauInit, bool saveExpInit, std::string expNameInit, bool alwaysSearchInit);

	void markerOdomCB(const nav_msgs::Odometry::ConstPtr& msg);

	void odomCB(const nav_msgs::Odometry::ConstPtr& msg);

	void imageCB(const sensor_msgs::Image::ConstPtr& msg);

	bool match(cv::Mat& image, float dt, Eigen::Vector3f vc, Eigen::Vector3f wc, ros::Time t, std::vector<Eigen::Vector3f>& pics);

	void update(std::vector<cv::Point2f>& pPts, std::vector<cv::Point2f>& kPts, std::vector<cv::Point2f>& cPts, Eigen::Vector3f vc, Eigen::Vector3f wc, ros::Time t, float dt, std::vector<Eigen::Vector3f>& pics, bool patternFound);
};

#endif
