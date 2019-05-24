#include <mutex>

#include <ros/ros.h>
#include <ros/console.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/TwistStamped.h>
#include <geometry_msgs/PoseStamped.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <helper_functions.h>

struct OdomEstimator
{
	ros::NodeHandle nh;
  ros::Publisher odomPub,markerOdomPub,camPosePub;
  ros::Subscriber poseSub,velSub,markerPoseSub;
  std::string bodyName;
	std::string cameraName;
	std::string markerName;
	std::mutex poseMutex;
	std::mutex poseMarkerMutex;
	ros::Time tPoseLast,tVelLast,tMarkerPoseLast;//last time

	bool firstPose,firstVel,firstMarker;
	Eigen::Vector3f pcwHat;
	Eigen::Vector4f qcwHat;
	Eigen::Vector3f vcHat;
	Eigen::Vector3f wcHat;
	Eigen::Vector3f pfi;
	Eigen::Vector4f qfi;
	Eigen::Vector3f pmwHat;
	Eigen::Vector4f qmwHat;

	float piTau;
	float qiTau;
	float viTau;
	float wiTau;

  OdomEstimator();

	void markerPoseCB(const nav_msgs::Odometry::ConstPtr& msg);

	void poseCB(const nav_msgs::Odometry::ConstPtr& msg);

	void velCB(const nav_msgs::Odometry::ConstPtr& msg);
};
