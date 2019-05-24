#include <odom_estimator.h>

OdomEstimator::OdomEstimator()
{
	// Parameters
	ros::NodeHandle nhp("~");
	nhp.param<std::string>("bodyName", bodyName, "turtlebot3");
	nhp.param<std::string>("cameraName", cameraName, "camera");
	nhp.param<std::string>("markerName", markerName, "board");

	float pfix,pfiy,pfiz,qfiw,qfix,qfiy,qfiz,qmew,qmex,qmey,qmez;
	nhp.param<float>("pfix", pfix, 0.0);
	nhp.param<float>("pfiy", pfiy, 0.0);
	nhp.param<float>("pfiz", pfiz, 0.0);
	nhp.param<float>("qfiw", qfiw, 1.0);
	nhp.param<float>("qfix", qfix, 0.0);
	nhp.param<float>("qfiy", qfiy, 0.0);
	nhp.param<float>("qfiz", qfiz, 0.0);

	pfi = Eigen::Vector3f(pfix,pfiy,pfiz);

	qfi << qfiw,qfix,qfiy,qfiz;
	qfi /= qfi.norm();

	float fpi,fqi,fvi,fwi;
	nhp.param<float>("fpi", fpi, 1.0);
	nhp.param<float>("fqi", fqi, 1.0);
	nhp.param<float>("fvi", fvi, 1.0);
	nhp.param<float>("fwi", fwi, 1.0);
	piTau = 1.0/(2.0*M_PI*fqi);
	qiTau = 1.0/(2.0*M_PI*fqi);
	viTau = 1.0/(2.0*M_PI*fvi);
	wiTau = 1.0/(2.0*M_PI*fwi);

	firstPose = true;

	//publisher
	camPosePub = nh.advertise<geometry_msgs::PoseStamped>(cameraName+"/mocapPose",1);
  odomPub = nh.advertise<nav_msgs::Odometry>(cameraName+"/odom",1);
	markerOdomPub = nh.advertise<nav_msgs::Odometry>(markerName+"/odom",1);

	// Subscriber
	velSub = nh.subscribe(bodyName+"/odom",5,&OdomEstimator::velCB,this);
	poseSub = nh.subscribe(bodyName+"/odomEKF",5,&OdomEstimator::poseCB,this);
	markerPoseSub = nh.subscribe(markerName+"/odomEKF",5,&OdomEstimator::markerPoseCB,this);

	pcwHat = Eigen::Vector3f::Zero();
	qcwHat = Eigen::Vector4f::Zero();
	qcwHat(0) = 1.0;
	vcHat = Eigen::Vector3f::Zero();
	wcHat = Eigen::Vector3f::Zero();

	pmwHat = Eigen::Vector3f::Zero();
	qmwHat = Eigen::Vector4f::Zero();
	qmwHat(0) = 1.0;
}

void OdomEstimator::markerPoseCB(const nav_msgs::Odometry::ConstPtr& msg)
{
	//pose is of turtlebot in room frame x west, y south, z up
	ros::Time t = msg->header.stamp;
	Eigen::Vector3f pmw(msg->pose.pose.position.x,msg->pose.pose.position.y,msg->pose.pose.position.z);
	Eigen::Vector4f qmw(msg->pose.pose.orientation.w,msg->pose.pose.orientation.x,msg->pose.pose.orientation.y,msg->pose.pose.orientation.z);
	qmw /= qmw.norm();

	{
		std::lock_guard<std::mutex> poseMutexGuard(poseMarkerMutex);
		if (firstMarker)
		{
			tMarkerPoseLast = t;
			pmwHat = pmw;
			qmwHat = qmw;
			firstMarker = false;
			return;
		}

		if ((qmwHat + qmw).norm() < (qmwHat - qmw).norm())
		{
			qmw *= -1.0;
		}

		float dt = (t-tMarkerPoseLast).toSec();
		tMarkerPoseLast = t;

		// get the low pass gains
		float kpi = dt/(piTau+dt);
		float kqi = dt/(qiTau+dt);

		pmwHat += kpi*(pmw-pmwHat);
		qmwHat += kqi*(qmw-qmwHat);
		qmwHat /= qmwHat.norm();
	}
}

void OdomEstimator::velCB(const nav_msgs::Odometry::ConstPtr& msg)
{
	//check to make sure mocap recieved
	if (firstPose || firstMarker)
	{
		ROS_ERROR("NO MOCAP");
		return;
	}

	//velocity is of turtlebot in body frame of turtlebot x forward, y to left, z up
	// odom pose is in ENU, odom twist is in body
	ros::Time t = msg->header.stamp;
	Eigen::Vector3f vi(msg->twist.twist.linear.x,msg->twist.twist.linear.y,msg->twist.twist.linear.z);
	Eigen::Vector3f wi(msg->twist.twist.angular.x,msg->twist.twist.angular.y,msg->twist.twist.angular.z);

	{
		std::lock_guard<std::mutex> poseMutexGuard(poseMutex);
		Eigen::Vector3f vc = rotatevec((vi+getss(wi)*pfi),getqInv(qfi));
		Eigen::Vector3f wc = rotatevec(wi,getqInv(qfi));

		if (firstVel)
		{
			tVelLast = t;
			vcHat = vc;
			wcHat = wc;
			firstVel = false;
		}

		float dt = (t-tVelLast).toSec();
		tVelLast = t;

		// get the low pass gains
		float kvi = dt/(viTau+dt);
		float kwi = dt/(wiTau+dt);

		vcHat += kvi*(vc-vcHat);
		wcHat += kwi*(wc-wcHat);

		// build and publish odom message
		nav_msgs::Odometry odomMsg;
		odomMsg.header.stamp = t;
		odomMsg.header.frame_id = "world";
		odomMsg.child_frame_id = "camera";
		odomMsg.pose.pose.position.x = pcwHat(0);
		odomMsg.pose.pose.position.y = pcwHat(1);
		odomMsg.pose.pose.position.z = pcwHat(2);
		odomMsg.pose.pose.orientation.w = qcwHat(0);
		odomMsg.pose.pose.orientation.x = qcwHat(1);
		odomMsg.pose.pose.orientation.y = qcwHat(2);
		odomMsg.pose.pose.orientation.z = qcwHat(3);
		odomMsg.twist.twist.linear.x = vcHat(0);
		odomMsg.twist.twist.linear.y = vcHat(1);
		odomMsg.twist.twist.linear.z = vcHat(2);
		odomMsg.twist.twist.angular.x = wcHat(0);
		odomMsg.twist.twist.angular.y = wcHat(1);
		odomMsg.twist.twist.angular.z = wcHat(2);
		odomPub.publish(odomMsg);


		// build and publish odom message
		geometry_msgs::PoseStamped camPoseMsg;
		camPoseMsg.header.stamp = t;
		camPoseMsg.header.frame_id = "world";
		camPoseMsg.pose.position.x = pcwHat(0);
		camPoseMsg.pose.position.y = pcwHat(1);
		camPoseMsg.pose.position.z = pcwHat(2);
		camPoseMsg.pose.orientation.w = qcwHat(0);
		camPoseMsg.pose.orientation.x = qcwHat(1);
		camPoseMsg.pose.orientation.y = qcwHat(2);
		camPoseMsg.pose.orientation.z = qcwHat(3);
		camPosePub.publish(camPoseMsg);
	}

	{
		std::lock_guard<std::mutex> poseMutexGuard(poseMarkerMutex);
		// build and publish odom message
		nav_msgs::Odometry odomMsg;
		odomMsg.header.stamp = t;
		odomMsg.header.frame_id = "world";
		odomMsg.child_frame_id = "marker";
		odomMsg.pose.pose.position.x = pmwHat(0);
		odomMsg.pose.pose.position.y = pmwHat(1);
		odomMsg.pose.pose.position.z = pmwHat(2);
		odomMsg.pose.pose.orientation.w = qmwHat(0);
		odomMsg.pose.pose.orientation.x = qmwHat(1);
		odomMsg.pose.pose.orientation.y = qmwHat(2);
		odomMsg.pose.pose.orientation.z = qmwHat(3);
		odomMsg.twist.twist.linear.x = 0.0;
		odomMsg.twist.twist.linear.y = 0.0;
		odomMsg.twist.twist.linear.z = 0.0;
		odomMsg.twist.twist.angular.x = 0.0;
		odomMsg.twist.twist.angular.y = 0.0;
		odomMsg.twist.twist.angular.z = 0.0;
		markerOdomPub.publish(odomMsg);
	}
}

void OdomEstimator::poseCB(const nav_msgs::Odometry::ConstPtr& msg)
{
	//pose is of turtlebot in room frame x west, y south, z up
	ros::Time t = msg->header.stamp;
	Eigen::Vector3f piw(msg->pose.pose.position.x,msg->pose.pose.position.y,msg->pose.pose.position.z);
	Eigen::Vector4f qiw(msg->pose.pose.orientation.w,msg->pose.pose.orientation.x,msg->pose.pose.orientation.y,msg->pose.pose.orientation.z);
	qiw /= qiw.norm();

	{
		std::lock_guard<std::mutex> poseMutexGuard(poseMutex);
		Eigen::Vector3f pcw = piw + rotatevec(pfi,qiw);
		Eigen::Vector4f qcw = getqMat(qiw)*qfi;
		qcw /= qcw.norm();

		if (firstPose)
		{
			tPoseLast = t;
			pcwHat = pcw;
			qcwHat = qcw;
			firstPose = false;
		}

		if ((qcwHat + qcw).norm() < (qcwHat - qcw).norm())
		{
			qcw *= -1.0;
		}

		float dt = (t-tPoseLast).toSec();
		tPoseLast = t;

		// get the low pass gains
		float kpi = dt/(piTau+dt);
		float kqi = dt/(qiTau+dt);

		pcwHat += kpi*(pcw-pcwHat);
		qcwHat += kqi*(qcw-qcwHat);
		qcwHat /= qcwHat.norm();
	}
}
