#include <forward_image_receiver.h>

ForwardImageReceiver::~ForwardImageReceiver()
{
	delete patch;
}

ForwardImageReceiver::ForwardImageReceiver() : it(nh)
{
	// Parameters
	ros::NodeHandle nhp("~");
	nhp.param<std::string>("cameraName", cameraName, "camera");
	nhp.param<std::string>("markerName", markerName, "board");
	nhp.param<float>("fq", fq, 1.0);
	nhp.param<float>("fp", fp, 1.0);
	nhp.param<float>("ft", ft, 1.0);
	nhp.param<float>("fn", fn, 1.0);
	nhp.param<float>("fd", fd, 1.0);
	nhp.param<float>("zmin", zmin, 1.0);
	nhp.param<float>("zmax", zmax, 1.0);
	nhp.param<int>("minFeaturesBad", minFeaturesBad, 8);
	float tau;
	nhp.param<float>("tau", tau, 1.0);

	bool alwaysSearch;
	nhp.param<bool>("saveExp", saveExp, false);
	nhp.param<std::string>("expName", expName, "exp");
	nhp.param<bool>("alwaysSearch", alwaysSearch, false);

	// Get camera parameters
	gotCamParam = false;
	camInfoSub = nh.subscribe(cameraName+"/camera_info",1,&ForwardImageReceiver::camInfoCB,this);
	ROS_INFO("Waiting for camera parameters on topic %s/camera_info",cameraName.c_str());
	do
	{
		ros::spinOnce();
		ros::Duration(0.3).sleep();
	} while (!(ros::isShuttingDown()) && !gotCamParam);
	ROS_INFO("Got forward camera parameters");

	patch = new PatchEstimator(imageWidth, imageHeight, minFeaturesBad, camMat.at<float>(0,0),camMat.at<float>(1,1),camMat.at<float>(0,2),camMat.at<float>(1,2),zmin,zmax,fq,fp,ft,fn,fd,cameraName,markerName,tau,saveExp,expName,alwaysSearch);

	// Publishers
	undistortPub = it.advertise(cameraName+"/image_undistort",1);

	// Subscribers
	imageSub = it.subscribe(cameraName+"/image_raw", 60, &ForwardImageReceiver::imageCB,this);
}

//image callback
void ForwardImageReceiver::imageCB(const sensor_msgs::Image::ConstPtr& msg)
{
	clock_t processTime = clock();

	// convert to opencv image
	cv_bridge::CvImagePtr cv_ptr;
	ros::Time t = msg->header.stamp;
	// std::cout << "\n t forward " << t << std::endl;

	cv::Mat image; // current image
	try
	{
		cv_ptr = cv_bridge::toCvCopy(msg, sensor_msgs::image_encodings::MONO8);
		image = cv_ptr->image;
	}
	catch (cv_bridge::Exception& e)
	{
		ROS_ERROR("cv_bridge exception: %s", e.what());
		return;
	}

	// Prepare image for processing
	cv::Mat imageUndistort(image.size(),CV_8UC1);
	cv::remap(image,imageUndistort,map1,map2,CV_INTER_LINEAR,cv::BORDER_CONSTANT,cv::Scalar(0,0,0));

	//publish image
	cv_bridge::CvImage undistort_msg;
	undistort_msg.header = msg->header; // Same timestamp and tf frame as input image
	undistort_msg.encoding = sensor_msgs::image_encodings::MONO8; // Or whatever
	undistort_msg.image = imageUndistort; // Your cv::Mat
	undistortPub.publish(undistort_msg.toImageMsg());
}

// callback for getting camera intrinsic parameters
void ForwardImageReceiver::camInfoCB(const sensor_msgs::CameraInfo::ConstPtr& camInfoMsg)
{
	//get camera info
	image_geometry::PinholeCameraModel cam_model;
	cam_model.fromCameraInfo(camInfoMsg);
	cv::Mat cam_calib_matTemp = cv::Mat(cam_model.fullIntrinsicMatrix());
	cam_calib_matTemp.convertTo(camMat,CV_32F);
	cam_model.distortionCoeffs().convertTo(distCoeffs,CV_32F);

	imageWidth = camInfoMsg->width;
	imageHeight = camInfoMsg->height;

	cv::Size imageSize(imageWidth,imageHeight);

	cv::initUndistortRectifyMap(camMat,distCoeffs,cv::Mat::eye(3,3,CV_32F),camMat,imageSize,CV_16SC2,map1,map2);

	//unregister subscriber
	camInfoSub.shutdown();
	gotCamParam = true;
}
