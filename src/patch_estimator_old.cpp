#include <patch_estimator.h>

PatchEstimator::~PatchEstimator()
{
	if (saveExp)
	{
		std::cout << "\nsaving\n";
		std::ofstream saveFile("/home/ncr/ncr_ws/src/icl_stationary/experiment/"+expName+".txt");
		if (saveFile.is_open())
		{
			std::cout << "\nopen\n";
			saveFile << "time,";
			saveFile << "pkcx," << "pkcy," << "pkcz,";
			saveFile << "qkcw," << "qkcx," << "qkcy," << "qkcz,";
			saveFile << "pkcHatx," << "pkcHaty," << "pkcHatz,";
			saveFile << "qkcHatw," << "qkcHatx," << "qkcHaty," << "qkcHatz,";
			saveFile << "pikx," << "piky," << "pikz,";
			saveFile << "pikHatx," << "pikHaty," << "pikHatz,";
			saveFile << "picx," << "picy," << "picz,";
			saveFile << "picHatx," << "picHaty," << "picHatz,";
			saveFile << "picEKFx," << "picEKFy," << "picEKFz,";
			saveFile << "picLSx," << "picLSy," << "picLSz,";
			saveFile << "vx," << "vy," << "vz,";
			saveFile << "wx," << "wy," << "wz,";
			saveFile << "\n";

			for (int jj = 0; jj < data.size(); jj++)
			{
				float timej = data.at(jj)->time;
				Eigen::Vector3f pkcj = data.at(jj)->pkc;
				Eigen::Vector4f qkcj = data.at(jj)->qkc;
				Eigen::Vector3f pkcHatj = data.at(jj)->pkcHat;
				Eigen::Vector4f qkcHatj = data.at(jj)->qkcHat;
				Eigen::Vector3f pikj = data.at(jj)->pik;
				Eigen::Vector3f pikHatj = data.at(jj)->pikHat;
				Eigen::Vector3f picj = data.at(jj)->pic;
				Eigen::Vector3f picHatj = data.at(jj)->picHat;
				Eigen::Vector3f picEKFj = data.at(jj)->picEKF;
				Eigen::Vector3f picLSj = data.at(jj)->picLS;
				Eigen::Vector3f vj = data.at(jj)->v;
				Eigen::Vector3f wj = data.at(jj)->w;

				saveFile << timej << ",";
				saveFile << pkcj(0) << "," << pkcj(1) << "," << pkcj(2) << ",";
				saveFile << qkcj(0) << "," << qkcj(1) << "," << qkcj(2) << "," << qkcj(3) << ",";
				saveFile << pkcHatj(0) << "," << pkcHatj(1) << "," << pkcHatj(2) << ",";
				saveFile << qkcHatj(0) << "," << qkcHatj(1) << "," << qkcHatj(2) << "," << qkcHatj(3) << ",";
				saveFile << pikj(0) << "," << pikj(1) << "," << pikj(2) << ",";
				saveFile << pikHatj(0) << "," << pikHatj(1) << "," << pikHatj(2) << ",";
				saveFile << picj(0) << "," << picj(1) << "," << picj(2) << ",";
				saveFile << picHatj(0) << "," << picHatj(1) << "," << picHatj(2) << ",";
				saveFile << picEKFj(0) << "," << picEKFj(1) << "," << picEKFj(2) << ",";
				saveFile << picLSj(0) << "," << picLSj(1) << "," << picLSj(2) << ",";
				saveFile << vj(0) << "," << vj(1) << "," << vj(2) << ",";
				saveFile << wj(0) << "," << wj(1) << "," << wj(2) << ",";
				saveFile << "\n";

				delete data.at(jj);
			}
			saveFile.close();
			std::cout << "\nclose\n";
		}
		std::cout << "\nsaved\n";
	}
}

PatchEstimator::PatchEstimator() : it(nh)
{}

PatchEstimator::PatchEstimator(int imageWidthInit, int imageHeightInit, int minFeaturesBadInit, float fxInit, float fyInit, float cxInit, float cyInit, float zminInit, float zmaxInit, float fq, float fp, float ft, float fn, float fd, std::string cameraNameInit, std::string markerNameInit, float tauInit, bool saveExpInit, std::string expNameInit, bool alwaysSearchInit) : it(nh)
{
	imageWidth = imageWidthInit;
	imageHeight = imageHeightInit;
	minFeaturesBad = minFeaturesBadInit;
	patchShutdown = false;
	firstOdomImageCB = true;
	firstImageCB = true;
	patternLost = false;
	tkcHat = Eigen::Vector3f::Zero();
	nkHat = Eigen::Vector3f(0.0,0.0,1.0);
	qkcHat = Eigen::Vector4f(1.0,0.0,0.0,0.0);
	pkcHat = Eigen::Vector3f::Zero();
	pckHat = Eigen::Vector3f::Zero();
 	qckHat = Eigen::Vector4f(1.0,0.0,0.0,0.0);
	fx = fxInit;
	fy = fyInit;
	cx = cxInit;
	cy = cyInit;
	camMat = cv::Mat::zeros(3,3,CV_32F);
	camMat.at<float>(0,0) = fx;
	camMat.at<float>(0,2) = cx;
	camMat.at<float>(1,1) = fy;
	camMat.at<float>(1,2) = cy;
	camMat.at<float>(2,2) = 1.0;

	camMatf = Eigen::Matrix3f::Zero();
	camMatf(0,0) = fx;
	camMatf(0,2) = cx;
	camMatf(1,1) = fy;
	camMatf(1,2) = cy;
	camMatf(2,2) = 1.0;

	Eigen::JacobiSVD<Eigen::MatrixXf> svdcamMatf(camMatf, Eigen::ComputeThinU | Eigen::ComputeThinV);
	camMatIf = svdcamMatf.solve(Eigen::Matrix3f::Identity());

	zmin = zminInit;
	zmax = zmaxInit;
	dkcHat = zmin;
	dkHat = zmax;

	pTau = 1.0/(2.0*M_PI*fp);
	qTau = 1.0/(2.0*M_PI*fq);
	tTau = 1.0/(2.0*M_PI*ft);
	nTau = 1.0/(2.0*M_PI*fn);
	dTau = 1.0/(2.0*M_PI*fd);

	tau = tauInit;

	cameraName = cameraNameInit;
	markerName = markerNameInit;

	saveExp = saveExpInit;
	expName = expNameInit;
	alwaysSearch = alwaysSearchInit;

	imagePub = it.advertise(cameraName+"/image_output",1);
	imageSub = it.subscribe(cameraName+"/image_undistort", 1000, &PatchEstimator::imageCB,this);
	odomSub = nh.subscribe(cameraName+"/odom", 1000, &PatchEstimator::odomCB,this);
	markerOdomSub = nh.subscribe(markerName+"/odom", 1000, &PatchEstimator::markerOdomCB,this);
	odomPub = nh.advertise<nav_msgs::Odometry>(cameraName+"/odomHat",1);
	pointCloudPub = nh.advertise<PointCloud> ("wall_map", 1);
}

void PatchEstimator::markerOdomCB(const nav_msgs::Odometry::ConstPtr& msg)
{
	if (patchShutdown)
	{
		return;
	}
	std::lock_guard<std::mutex> odomMutexGuard(markerOdomMutex);
	markerOdomSync.push_back(*msg);
}

void PatchEstimator::odomCB(const nav_msgs::Odometry::ConstPtr& msg)
{
	if (patchShutdown)
	{
		return;
	}
	std::lock_guard<std::mutex> odomMutexGuard(odomMutex);
	odomSync.push_back(*msg);
}

//image callback
void PatchEstimator::imageCB(const sensor_msgs::Image::ConstPtr& msg)
{
	ROS_WARN("IMAGE START");
	if (patchShutdown)
	{
		return;
	}

	clock_t callbackTime = clock();
	clock_t processTime = clock();
	std::cout << std::endl;

	// sync the times to find the closest
	ros::Time t = msg->header.stamp;

	{
		std::lock_guard<std::mutex> odomMutexGuard(odomMutex);
		std::vector<float> timeDiff;
		int numMatch = odomSync.size();

		// std::cout << "\n t key " << t << std::endl;

		std::cout << "\n numMatch " << numMatch << std::endl;
		if (numMatch > 0)
		{
			// std::cout << "\n recent dt \n" << (odomSync.at(numMatch-1).header.stamp - t).toSec() << std::endl;
			for (int ii = 0; ii < numMatch; ii++)
			{
				float dtii = fabsf((odomSync.at(ii).header.stamp - t).toSec());
				timeDiff.push_back(fabsf(dtii));
			}

			int minTimeInd = std::distance(timeDiff.begin(),std::min_element(timeDiff.begin(),timeDiff.end()));
			imageOdom = odomSync.at(minTimeInd);

			std::cout << "\n odomSync size " << odomSync.size() << std::endl;
			for (int ii = 0; ii <= minTimeInd; ii++)
			{
				odomSync.pop_front();
			}
			std::cout << "\n odomSync size " << odomSync.size() << std::endl;
			if (firstOdomImageCB)
			{
				firstOdomImageCB = false;
			}
		}
		else
		{
			ROS_ERROR("NO NEW IMAGE ODOM");
			if (firstOdomImageCB)
			{
				return;
			}
		}
	}

	Eigen::Vector3f pcw(imageOdom.pose.pose.position.x,imageOdom.pose.pose.position.y,imageOdom.pose.pose.position.z);
	Eigen::Vector4f qcw(imageOdom.pose.pose.orientation.w,imageOdom.pose.pose.orientation.x,imageOdom.pose.pose.orientation.y,imageOdom.pose.pose.orientation.z);

	{
		std::lock_guard<std::mutex> odomMutexGuard(markerOdomMutex);
		std::vector<float> timeDiff;
		int numMatch = markerOdomSync.size();

		// std::cout << "\n t key " << t << std::endl;

		std::cout << "\n numMatch marker " << numMatch << std::endl;
		if (numMatch > 0)
		{
			// std::cout << "\n recent dt \n" << (odomSync.at(numMatch-1).header.stamp - t).toSec() << std::endl;
			for (int ii = 0; ii < numMatch; ii++)
			{
				float dtii = fabsf((markerOdomSync.at(ii).header.stamp - t).toSec());
				timeDiff.push_back(fabsf(dtii));
			}

			int minTimeInd = std::distance(timeDiff.begin(),std::min_element(timeDiff.begin(),timeDiff.end()));
			markerOdom = markerOdomSync.at(minTimeInd);

			std::cout << "\n markerSync size " << markerOdomSync.size() << std::endl;
			for (int ii = 0; ii <= minTimeInd; ii++)
			{
				markerOdomSync.pop_front();
			}
			std::cout << "\n markerSync size " << markerOdomSync.size() << std::endl;
			if (firstMarkerOdomImageCB)
			{
				firstMarkerOdomImageCB = false;
			}
		}
		else
		{
			ROS_ERROR("NO NEW MARKER ODOM");
			if (firstMarkerOdomImageCB)
			{
				return;
			}
		}
	}

	ROS_WARN("get odom time %2.4f",float(clock()-processTime)/CLOCKS_PER_SEC);
	processTime = clock();

	// std::cout << "\n fabsf(imageOdom.pose.pose.position.y) " << fabsf(imageOdom.pose.pose.position.y) << std::endl;
	// convert to opencv image
	cv_bridge::CvImagePtr cv_ptr;
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

	ROS_WARN("get image time %2.4f",float(clock()-processTime)/CLOCKS_PER_SEC);
	processTime = clock();

	Eigen::Vector3f pmw(markerOdom.pose.pose.position.x,markerOdom.pose.pose.position.y,markerOdom.pose.pose.position.z);
	Eigen::Vector4f qmw(markerOdom.pose.pose.orientation.w,markerOdom.pose.pose.orientation.x,markerOdom.pose.pose.orientation.y,markerOdom.pose.pose.orientation.z);
	std::vector<Eigen::Vector3f> pics;
	for (int ii = 0; ii < 48; ii++)
	{
		int rowi = ii/8;
		int coli = ii%8;
		Eigen::Vector3f pim(float(rowi+1)*0.06,float(coli+1)*0.06,0.007);
		Eigen::Vector3f piw = pmw + rotatevec(pim,qmw);
		Eigen::Vector3f pic = rotatevec((piw - pcw),getqInv(qcw));
		pics.push_back(pic);
	}

	if (firstImageCB)
	{
		// Check if chessboard is visible
		cv::Size patternSize(8,6);
		std::vector<cv::Point2f> pts;
		bool patternFound = cv::findChessboardCorners(image,patternSize,pts,cv::CALIB_CB_ADAPTIVE_THRESH+cv::CALIB_CB_NORMALIZE_IMAGE+cv::CALIB_CB_FAST_CHECK);

		//if found initialize the patch
		if (patternFound)
		{
				cv::cornerSubPix(image,pts,cv::Size(11,11),cv::Size(-1,-1),cv::TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1));

				keyOdom = imageOdom;
				markerKeyOdom = markerOdom;
				kimage = image.clone();
				pimage = image.clone();
				tLast = t;
				firstImageCB = false;
				tStart = t;

				for (int ii = 0; ii < pts.size(); ii++)
				{
					int rowi = ii/8;
					int coli = ii%8;
					Eigen::Vector3f pim(float(rowi+1)*0.06,float(coli+1)*0.06,0.007);
					Eigen::Vector3f piw = pmw + rotatevec(pim,qmw);
					Eigen::Vector3f pic = rotatevec((piw - pcw),getqInv(qcw));
					DepthEstimator* newDepthEstimator = new DepthEstimator(ii,Eigen::Vector3f((pts.at(ii).x-cx)/fx,(pts.at(ii).y-cy)/fy,1.0),tLast,zmin,zmax,pic,tau);
					// DepthEstimator* newDepthEstimator = new DepthEstimator(ii,Eigen::Vector3f((pts.at(ii).x-cx)/fx,(pts.at(ii).y-cy)/fy,1.0),tLast,zmin,zmax,R2I,R4I,Q,P);
					depthEstimators.push_back(newDepthEstimator);
				}

				return;
		}
		else
		{
			return;
		}
	}

	ROS_WARN("get checkerboard time %2.4f",float(clock()-processTime)/CLOCKS_PER_SEC);
	processTime = clock();

	float dt = (t-tLast).toSec();
	float timeFromStart = (t-tStart).toSec();
	tLast = t;

	Eigen::Vector3f vc(imageOdom.twist.twist.linear.x,imageOdom.twist.twist.linear.y,imageOdom.twist.twist.linear.z);
	Eigen::Vector3f wc(imageOdom.twist.twist.angular.x,imageOdom.twist.twist.angular.y,imageOdom.twist.twist.angular.z);

	// find the features
	bool patternFound = match(image,dt,vc,wc,t,pics);

	// pckHat += (rotatevec(vc,qckHat)*dt);
	// qckHat += (0.5*B(qckHat)*wc*dt);
	// qckHat /= qckHat.norm();
	// pkcHat = rotatevec(-pckHat,getqInv(qckHat));
	// qkcHat = getqInv(qckHat);
	// qkcHat /= qkcHat.norm();
	pkcHat += ((-vc-getss(wc)*pkcHat)*dt);
	qkcHat += (-0.5*B(qkcHat)*wc*dt);
	qkcHat /= qkcHat.norm();
	pckHat = rotatevec(-pkcHat,getqInv(qkcHat));
	qckHat = getqInv(qkcHat);
	qckHat /= qckHat.norm();

	std::cout << "\n dt \n" << dt << std::endl;
	std::cout << "\n vc \n" << vc << std::endl;
	std::cout << "\n wc \n" << wc << std::endl;
	std::cout << "\n pkcHat \n" << pkcHat << std::endl;
	std::cout << "\n qkcHat \n" << qkcHat << std::endl;

	// ros::shutdown();

	ROS_WARN("get match time %2.4f",float(clock()-processTime)/CLOCKS_PER_SEC);
	processTime = clock();

	pimage = image.clone();

	//send the pose estimate to the odom
	// Eigen::Vector4f qckHat = getqInv(qkcHat);
	// qckHat /= qckHat.norm();
	// Eigen::Vector3f pckHat = rotatevec(-pkcHat,qckHat);

	std::cout << "\n pckHat \n" << pckHat << std::endl;
	std::cout << "\n qckHat \n" << qckHat << std::endl;

	Eigen::Vector3f pkw(keyOdom.pose.pose.position.x,keyOdom.pose.pose.position.y,keyOdom.pose.pose.position.z);
	Eigen::Vector4f qkw(keyOdom.pose.pose.orientation.w,keyOdom.pose.pose.orientation.x,keyOdom.pose.pose.orientation.y,keyOdom.pose.pose.orientation.z);

	Eigen::Vector4f qcwHat = getqMat(qkw)*qckHat;
	qcwHat /= qcwHat.norm();

	Eigen::Vector3f pcwHat = rotatevec(pckHat,qkw)+pkw;

	Eigen::Vector3f pkc = rotatevec(pkw-pcw,getqInv(qcw));
	Eigen::Vector4f qkc = getqMat(getqInv(qcw))*qkw;
	qkc /= qkc.norm();

	std::cout << "\n pkc \n" << pkc << std::endl;
	std::cout << "\n qkc \n" << qkc << std::endl;

	// build and publish odom message
	nav_msgs::Odometry odomMsg;
	odomMsg.header.stamp = t;
	odomMsg.header.frame_id = "world";
	odomMsg.child_frame_id = "cameraHat";
	odomMsg.pose.pose.position.x = pcwHat(0);
	odomMsg.pose.pose.position.y = pcwHat(1);
	odomMsg.pose.pose.position.z = pcwHat(2);
	odomMsg.pose.pose.orientation.w = qcwHat(0);
	odomMsg.pose.pose.orientation.x = qcwHat(1);
	odomMsg.pose.pose.orientation.y = qcwHat(2);
	odomMsg.pose.pose.orientation.z = qcwHat(3);
	odomMsg.twist.twist.linear.x = vc(0);
	odomMsg.twist.twist.linear.y = vc(1);
	odomMsg.twist.twist.linear.z = vc(2);
	odomMsg.twist.twist.angular.x = wc(0);
	odomMsg.twist.twist.angular.y = wc(1);
	odomMsg.twist.twist.angular.z = wc(2);
	odomPub.publish(odomMsg);


	if (depthEstimators.size() > 0)
	{
		//plot the points
		// get all the points and draw them on the image
		// std::cout << "\n wall 1 \n";
		PointCloudRGB::Ptr map(new PointCloudRGB);
		pcl_conversions::toPCL(msg->header.stamp,map->header.stamp);

		// std::cout << "\n wall 1 1 \n";
		map->header.frame_id = "world";

		// std::cout << "\n wall 1 2 \n";
		map->height = 1;
		map->is_dense = true;
		map->points.clear();

		//publish the image
		pcl::PointXYZRGB pt,ptHat;
		for (uint16_t ii = 0; ii < depthEstimators.size(); ii++)
		{
			//get the points after update
			Eigen::Vector3f mci = depthEstimators.at(ii)->mc;
			if (patternFound)
			{
				cv::Point2f cPti(fx*mci(0)+cx,fy*mci(1)+cy);
				cv::circle(image, cPti, 10, cv::Scalar(150, 150, 150), -1);
			}

			Eigen::Vector3f mkiHat = depthEstimators.at(ii)->mk;
			cv::Point2i kPti(int(fx*mkiHat(0)+cx),int(fy*mkiHat(1)+cy));
			uint8_t colori = kimage.at<uint8_t>(kPti.y,kPti.x);

			Eigen::Vector3f pkiHat = depthEstimators.at(ii)->uk;
			Eigen::Vector3f pciHat = depthEstimators.at(ii)->uc;
			Eigen::Vector3f pkiHatICL = pkiHat;
			Eigen::Vector3f pciHatICL = pciHat;
			Eigen::Vector3f pkiHatNorm = depthEstimators.at(ii)->mk;
			Eigen::Vector3f pciHatEKF = depthEstimators.at(ii)->mc;
			Eigen::Vector3f pciHatLS = depthEstimators.at(ii)->mc;
			Eigen::Vector3f pciHatNorm = depthEstimators.at(ii)->mc;
			pkiHat *= depthEstimators.at(ii)->dkHat;
			pciHat *= depthEstimators.at(ii)->dcHat;
			pkiHatICL *= depthEstimators.at(ii)->dkHatICL;
			pciHatICL *= depthEstimators.at(ii)->dcHatICL;
			pciHatEKF *= depthEstimators.at(ii)->zcHatEKF;
			pciHatLS *= depthEstimators.at(ii)->zcHatLS;
			pkiHatNorm *= depthEstimators.at(ii)->zkHatNorm;
			pciHatNorm *= depthEstimators.at(ii)->zcHatNorm;
			Eigen::Vector3f pki = rotatevec(pics.at(ii)-pkc,getqInv(qkc));

			if (saveExp)
			{
				DataSave* dataSave = new DataSave(timeFromStart,pkc,qkc,pkcHat,qkcHat,pki,pkiHat,pics.at(ii),pciHat,pciHatEKF,pciHatLS,vc,wc);
				data.push_back(dataSave);
			}

			pkiHat = pkw + rotatevec(pkiHat,qkw);
			pkiHatICL = pkw + rotatevec(pkiHatICL,qkw);
			pkiHatNorm = pkw + rotatevec(pkiHatNorm,qkw);
			pciHat = pcw + rotatevec(pciHat,qcw);
			pciHatICL = pcw + rotatevec(pciHatICL,qcw);
			pciHatEKF = pcw + rotatevec(pciHatEKF,qcw);
			pciHatLS = pcw + rotatevec(pciHatLS,qcw);
			pciHatNorm = pcw + rotatevec(pciHatNorm,qcw);
			Eigen::Vector3f pci = pcw + rotatevec(pics.at(ii),qcw);

			pt.x = pci(0);
			pt.y = pci(1);
			pt.z = pci(2);
			pt.r = colori;
			pt.g = colori;
			pt.b = colori;
			map->points.push_back(pt);

			ptHat.x = pkiHat(0);
			ptHat.y = pkiHat(1);
			ptHat.z = pkiHat(2);
			ptHat.r = 255;
			ptHat.g = 0;
			ptHat.b = 0;
			map->points.push_back(ptHat);

			ptHat.x = pkiHatICL(0);
			ptHat.y = pkiHatICL(1);
			ptHat.z = pkiHatICL(2);
			ptHat.r = 175;
			ptHat.g = 0;
			ptHat.b = 0;
			map->points.push_back(ptHat);

			// ptHat.x = pkiHatNorm(0);
			// ptHat.y = pkiHatNorm(1);
			// ptHat.z = pkiHatNorm(2);
			// ptHat.r = 75;
			// ptHat.g = 0;
			// ptHat.b = 0;
			// map->points.push_back(ptHat);

			ptHat.x = pciHat(0);
			ptHat.y = pciHat(1);
			ptHat.z = pciHat(2);
			ptHat.r = 0;
			ptHat.g = 255;
			ptHat.b = 0;
			map->points.push_back(ptHat);

			ptHat.x = pciHatICL(0);
			ptHat.y = pciHatICL(1);
			ptHat.z = pciHatICL(2);
			ptHat.r = 0;
			ptHat.g = 175;
			ptHat.b = 0;
			map->points.push_back(ptHat);

			// ptHat.x = pciHatNorm(0);
			// ptHat.y = pciHatNorm(1);
			// ptHat.z = pciHatNorm(2);
			// ptHat.r = 0;
			// ptHat.g = 75;
			// ptHat.b = 0;
			// map->points.push_back(ptHat);

			ptHat.x = pciHatEKF(0);
			ptHat.y = pciHatEKF(1);
			ptHat.z = pciHatEKF(2);
			ptHat.r = 0;
			ptHat.g = 0;
			ptHat.b = 255;
			map->points.push_back(ptHat);

			ptHat.x = pciHatLS(0);
			ptHat.y = pciHatLS(1);
			ptHat.z = pciHatLS(2);
			ptHat.r = 0;
			ptHat.g = 255;
			ptHat.b = 255;
			map->points.push_back(ptHat);

			// std::cout << "\n index " << ii << std::endl;
			// std::cout << "\n dkHat " << depthEstimators.at(ii)->dkHat << std::endl;
			// std::cout << "\n dcHat " << depthEstimators.at(ii)->dcHat << std::endl;
			// std::cout << "\n zcHatEKF " << depthEstimators.at(ii)->zcHatEKF << std::endl;
			// std::cout << "\n dk " << depthEstimators.at(ii)->pik.norm() << std::endl;
		}
		map->width = map->points.size();
		pointCloudPub.publish(map);
	}

	ROS_WARN("get circle time %2.4f",float(clock()-processTime)/CLOCKS_PER_SEC);
	processTime = clock();

	// publish key image
	cv_bridge::CvImage out_msg;
	out_msg.header = msg->header; // Same timestamp and tf frame as input image
	out_msg.encoding = sensor_msgs::image_encodings::MONO8; // Or whatever
	out_msg.image = image; // Your cv::Mat

	{
		std::lock_guard<std::mutex> pubMutexGuard(pubMutex);
		imagePub.publish(out_msg.toImageMsg());
	}

	ROS_WARN("cb time %2.4f",float(clock()-callbackTime)/CLOCKS_PER_SEC);

	if (depthEstimators.size() < minFeaturesBad)
	{
		patchShutdown = true;
		ROS_WARN("shutdown after imagesub");
	}
}

//finds the features in the previous image in the new image and matches the features
bool PatchEstimator::match(cv::Mat& image, float dt, Eigen::Vector3f vc, Eigen::Vector3f wc, ros::Time t, std::vector<Eigen::Vector3f>& pics)
{
	clock_t estimatorUpdateTime = clock();

	//get the points from the previous image
	std::vector<cv::Point2f> pPts,kPts;
	for (int ii = 0; ii < depthEstimators.size(); ii++)
	{
		// Eigen::Vector4f xHati = depthEstimators.at(ii)->xHat;
		Eigen::Vector3f mpi = depthEstimators.at(ii)->mc;
		Eigen::Vector3f mki = depthEstimators.at(ii)->mk;
		pPts.push_back(cv::Point2f(fx*mpi(0)+cx,fy*mpi(1)+cy));
		kPts.push_back(cv::Point2f(fx*mki(0)+cx,fy*mki(1)+cy));
	}

	// ROS_WARN("keyframe %d patch %d pPts size before flow %d",keyInd,patchInd,int(pPts.size()));
	ROS_WARN("time for getting points %2.4f",float(clock()-estimatorUpdateTime)/CLOCKS_PER_SEC);
	estimatorUpdateTime = clock();

	// find the points using either approximated flow or looking for the board
	std::vector<cv::Point2f> cPts;//holds the estimated points
	bool patternFound = false;
	if (!alwaysSearch)
	{
		//use optical flow to find the features
		try
		{
			std::vector<uchar> status;//holds the status of a feature
			std::vector<float> err;//holds error of a feature
			cv::calcOpticalFlowPyrLK(pimage,image,pPts,cPts,status,err,cv::Size(21,21),1,cv::TermCriteria(cv::TermCriteria::COUNT|cv::TermCriteria::EPS,20,0.03),cv::OPTFLOW_LK_GET_MIN_EIGENVALS,0.0001);//optical flow to get new measurement
			cv::cornerSubPix(image,cPts,cv::Size(11,11),cv::Size(-1,-1),cv::TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1));
			ROS_WARN("time for optical flow %2.4f",float(clock()-estimatorUpdateTime)/CLOCKS_PER_SEC);
			estimatorUpdateTime = clock();
		}
		catch (cv::Exception e)
		{
			ROS_ERROR("optical flow failed");
		}
	}
	else
	{
			// Check if chessboard is visible
			cv::Size patternSize(8,6);
			// patternFound = cv::findChessboardCorners(image,patternSize,cPts,cv::CALIB_CB_ADAPTIVE_THRESH+cv::CALIB_CB_NORMALIZE_IMAGE+cv::CALIB_CB_FAST_CHECK);
			patternFound = cv::findChessboardCorners(image,patternSize,cPts,cv::CALIB_CB_FAST_CHECK);
			std::cout << "\n looking for pattern \n";
			//if found use subpix to correct to best corner individually
			if (patternFound)
			{
					cv::cornerSubPix(image,cPts,cv::Size(11,11),cv::Size(-1,-1),cv::TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1));
					std::cout << "\n pattern found \n";
			}
			else
			{
				patternLost = true;
				std::cout << "\n pattern lost \n";
			}
	}

	//update the estimtators using the estimted points
	update(pPts,kPts,cPts,vc,wc,t,dt,pics,patternFound);
	ROS_WARN("time for update call %2.4f",float(clock()-estimatorUpdateTime)/CLOCKS_PER_SEC);
	estimatorUpdateTime = clock();
	return patternFound;
}

void PatchEstimator::update(std::vector<cv::Point2f>& pPts, std::vector<cv::Point2f>& kPts, std::vector<cv::Point2f>& cPts, Eigen::Vector3f vc, Eigen::Vector3f wc, ros::Time t, float dt, std::vector<Eigen::Vector3f>& pics, bool patternFound)
{
	// //predict the orientation and position
	// pkcHat += (-vc*dt);
	// qkcHat += (-0.5*B(qkcHat)*wc*dt);
	// qkcHat /= qkcHat.norm();


	// std::cout << "\n dt \n" << dt << std::endl;
	// std::cout << "\n vc \n" << vc << std::endl;
	// std::cout << "\n wc \n" << wc << std::endl;
	// std::cout << "\n pkcHat \n" << pkcHat << std::endl;
	// std::cout << "\n qkcHat \n" << qkcHat << std::endl;
	float kp = 0.03/(pTau + 0.03);
	float kq = 0.03/(qTau + 0.03);
	float kt = 0.03/(tTau + 0.03);
	float kn = 0.03/(nTau + 0.03 + 100.0*fabsf(qkcHat(2)));
	float kd = 0.03/(dTau + 0.03);

	// return;

	// add in if the the pattern is found and always searching hjaslfahdahdfhasdhasdas
	if (!alwaysSearch)
	{
		clock_t updateClock = clock();
		try
		{
			assert(kPts.size()>0);
			assert(cPts.size()>0);

			cv::Mat inliersG,inliersGp;
			// cv::Mat G = cv::findHomography(kPts, cPts, cv::RANSAC, 2.0, inliersG, 2000, 0.99);//calculate homography using RANSAC
			cv::Mat G = cv::findHomography(kPts, cPts, 0);//calculate homography using RANSAC
			cv::Mat Gp = cv::findHomography(pPts, cPts, cv::RANSAC, 4.0, inliersGp, 2000, 0.99);//calculate homography using RANSAC

			ROS_WARN("time for homog %2.4f",float(clock()-updateClock)/CLOCKS_PER_SEC);
			updateClock = clock();

			// estimate the homography
			Eigen::Vector4f qkc = qkcHat;
			Eigen::Vector3f nk = nkHat;
			Eigen::Vector3f tkc = tkcHat;

			// adjust points based on homography
			// std::vector<cv::Point2f> cPtsAdj = cPts;
			// std::vector<cv::Point2f> cPtsAdj;
			// cv::perspectiveTransform(pPts,cPtsAdj,Gp);
			// cv::perspectiveTransform(pPts,cPtsAdj,G);

			if (!G.empty())
			{
				try
				{
					//find the solutions
					std::vector<cv::Mat> RkcH,tkcH,nkH;//holds the rotations, translations and normal vetors from decompose homography
					int numberHomogSols = cv::decomposeHomographyMat(G, camMat, RkcH, tkcH, nkH);// Decompose homography

					//check positive depth constraint on the inliers to find the solutions
					std::vector<Eigen::Vector3f> nks;
					std::vector<Eigen::Vector3f> tkcs;
					std::vector<Eigen::Vector4f> qkcs;
					std::vector<float> errors;

					// ROS_WARN("keyframe %d patch %d numberHomogSols %d",keyInd,patchInd,numberHomogSols);

					assert(numberHomogSols>0);
					assert(RkcH.size()>0);
					assert(tkcH.size()>0);
					assert(nkH.size()>0);

					for (int jj = 0; jj < numberHomogSols; jj++)
					{
						//convert normal to eigen
						Eigen::Vector3f nkj(nkH.at(jj).at<double>(0,0),nkH.at(jj).at<double>(1,0),nkH.at(jj).at<double>(2,0));

						// std::cout << "\n\n sol " << jj << std::endl;
						// std::cout << "\n nkjx " << nkj(0) << " nkjy " << nkj(1) << " nkjz " << nkj(2) << std::endl;

						if (nkj.norm() < 0.1)
						{
							nkj(0) = 0.0;
							nkj(1) = 0.0;
							nkj(2) = 1.0;
						}

						//if n^T*[0;0;1] > then solution in front of camera
						Eigen::Vector3f tkcj(tkcH.at(jj).at<double>(0,0),tkcH.at(jj).at<double>(1,0),tkcH.at(jj).at<double>(2,0));
						// std::cout << "\n tkcjx " << tkcj(0) << " tkcjy " << tkcj(1) << " tkcjz " << tkcj(2) << std::endl;
						if ((nkj(2) >= 0.0))
						{
							Eigen::Matrix3f Rkcj;
							for (int hh = 0; hh < 9; hh++)
							{
								Rkcj(hh/3,hh%3) = RkcH.at(jj).at<double>(hh/3,hh%3);
							}
							Eigen::Quaternionf qkcjq(Rkcj);// convert to quaternion
							Eigen::Vector4f qkcj(qkcjq.w(),qkcjq.x(),qkcjq.y(),qkcjq.z());
							qkcj /= qkcj.norm();

							if ((qkcHat + qkcj).norm() < (qkcHat - qkcj).norm())
							{
								qkcj *= -1.0;
							}

							// std::cout << "\n qkcjw " << qkcj(0) << " qkcjx " << qkcj(1) << " qkcjy " << qkcj(2) << " qkcjz " << qkcj(3) << std::endl;
							// std::cout << "\n tkcjx " << tkcj(0) << " tkcjy " << tkcj(1) << " tkcjz " << tkcj(2) << std::endl;
							// std::cout << "\n nkjx " << nkj(0) << " nkjy " << nkj(1) << " nkjz " << nkj(2) << std::endl;

							nks.push_back(nkj);
							tkcs.push_back(tkcj);
							qkcs.push_back(qkcj);
							if (tkcj.norm() > 0.001)
							{
								errors.push_back((qkcHat - qkcj).norm()+(pkcHat-tkcj).norm());
							}
							else
							{
								errors.push_back((qkcHat - qkcj).norm());
							}
						}
					}

					// ROS_WARN("keyframe %d patch %d errors size %d",keyInd,patchInd,int(errors.size()));

					int minqkcsErrorInd = std::distance(errors.begin(),std::min_element(errors.begin(),errors.end()));
					qkc = qkcs.at(minqkcsErrorInd);
					nk = nks.at(minqkcsErrorInd);
					tkc = tkcs.at(minqkcsErrorInd);

					// ROS_WARN("keyframe %d patch %d selected best",keyInd,patchInd);
				}
				catch (cv::Exception e)
				{
					ROS_ERROR("G failed");
				}
			}

			float kp = dt/(pTau + dt);
			float kq = dt/(qTau + dt);
			float kt = dt/(tTau + dt);
			float kn = dt/(nTau + dt + 100.0*fabsf(qkc(2)));
			float kd = dt/(dTau + dt);

			qkcHat += kq*(qkc - qkcHat);
			qkcHat /= qkcHat.norm();

			Eigen::Matrix3f RkcHat = getqRot(qkcHat);
			tkcHat += kt*(tkc - tkcHat);

			float alphank = 0.75;
			nk = (alphank*nk + (1.0-alphank)*(Eigen::Vector3f(0.0,0.0,1.0)));
			nk /= nk.norm();
			nkHat += kn*(nk - nkHat);
			nkHat /= nkHat.norm();

			ROS_WARN("time for update tnq %2.4f",float(clock()-updateClock)/CLOCKS_PER_SEC);
			updateClock = clock();

			Eigen::RowVector3f nkHatT = nkHat.transpose();
			Eigen::Matrix3f H = RkcHat+tkcHat*nkHatT;
			Eigen::RowVector3f H3 = H.block(2,0,1,3);
			Eigen::Matrix3f Gpf = Eigen::Matrix3f::Zero();
			for (int ii = 0; ii < 9; ii++)
			{
				Gpf(ii/3,ii%3) = Gp.at<double>(ii/3,ii%3);
			}

			Eigen::Matrix3f Hp = camMatIf*Gpf*camMatf;
			Eigen::Matrix<float,2,3> Hp12 = Hp.block(0,0,2,3);
			Eigen::RowVector3f Hp3 = Hp.block(2,0,1,3);
			// Eigen::RowVector3f Gp3 = Gpf.block(2,0,1,3);

			//remove the outliers
			std::vector<DepthEstimator*> depthEstimatorsIn;
			std::vector<float> dkcs;
			assert(depthEstimators.size() > 0);
			ROS_WARN("depthEstimators size before %d",int(depthEstimators.size()));
			// ROS_WARN("keyframe %d patch %d kPts size %d",keyInd,patchInd,int(kPts.size()));
			// ROS_WARN("keyframe %d patch %d cPts size %d",keyInd,patchInd,int(cPts.size()));
			Eigen::Vector3f mpi((pPts.at(0).x-cx)/fx,(pPts.at(0).y-cy)/fy,1.0);
			float H3mpi = Hp3*mpi;
			float alphapi = 1.0/H3mpi;
			Eigen::Vector3f mci1((cPts.at(0).x-cx)/fx,(cPts.at(0).y-cy)/fy,1.0);
			// Eigen::Vector3f mci = 0.1*alphapi*Hp*mpi+0.9*mci1;
			// mci /= mci(2);
			Eigen::Vector3f mci = mci1;

			float tlcx = (0.0-cx)/fx;
			float tlcy = (0.0-cy)/fy;
			float brcx = (imageWidth-cx)/fx;
			float brcy = (imageHeight-cy)/fy;

			// std::cout << "\n imageWidth " << imageWidth << std::endl;
			// std::cout << "\n imageHeight " << imageHeight << std::endl;
			// std::cout << "\n tlcx " << tlcx << std::endl;
			// std::cout << "\n tlcy " << tlcy << std::endl;
			// std::cout << "\n brcx " << brcx << std::endl;
			// std::cout << "\n brcy " << brcy << std::endl;
			// std::cout << "\n fx " << fx << std::endl;
			// std::cout << "\n fy " << fy << std::endl;
			// std::cout << "\n cx " << cx << std::endl;
			// std::cout << "\n cy " << cy << std::endl;

			for (int ii = 0; ii < kPts.size(); ii++)
			{
				mpi(0) = (pPts.at(ii).x-cx)/fx;
				mpi(1) = (pPts.at(ii).y-cy)/fy;
				H3mpi = Hp3*mpi;
				alphapi = 1.0/H3mpi;
				mci1(0) = (cPts.at(ii).x-cx)/fx;
				mci1(1) = (cPts.at(ii).y-cy)/fy;
				// mci.segment(0,2) = 0.25*alphapi*Hp12*mpi+0.75*mci1.segment(0,2);
				mci.segment(0,2) = mci1.segment(0,2);

				// std::cout << "\n cPts.at(ii).x " << cPts.at(ii).x << " cPts.at(ii).y " << cPts.at(ii).y << std::endl;
				// std::cout << "\n pPts.at(ii).x " << pPts.at(ii).x << " pPts.at(ii).y " << pPts.at(ii).y << std::endl;
				// std::cout << "\n mci1x " << mci1(0) << " mci1y " << mci1(1) << std::endl;
				// std::cout << "\n mcix " << mci(0) << " mciy " << mci(1) << std::endl;
				// std::cout << "\n inliersGp.at<uchar>(ii) " << int(inliersGp.at<uchar>(ii)) << std::endl;
				// std::cout << "\n (mci(0) >= tlcx) && (mci(0) < brcx) && (mci(1) >= tlcy) && (mci(1) < brcy) " << int((mci(0) >= tlcx) && (mci(0) < brcx) && (mci(1) >= tlcy) && (mci(1) < brcy)) << std::endl;

				if ((inliersGp.at<uchar>(ii)) && (mci(0) >= tlcx) && (mci(0) < brcx) && (mci(1) >= tlcy) && (mci(1) < brcy))
				{
					float dkcHati = depthEstimators.at(ii)->update(H,mci,nkHatT,tkcHat,RkcHat,vc,wc,t,pkcHat,qkcHat,pics.at(ii),false);
					depthEstimatorsIn.push_back(depthEstimators.at(ii));
					if (depthEstimators.at(ii)->dkKnown)
					{
						dkcs.push_back(dkcHati);
					}
				}
				else
				{
					delete depthEstimators.at(ii);

					// if (((mci(0)*fx+cx) >= 0) && ((mci(0)*fx+cx) < 640) && ((mci(1)*fy+cy) >= 0) && ((mci(1)*fy+cy) < 480))
					// {
					// 	// std::cout << "\n outlier pruned \n";
					// }
					// else
					// {
					// 	// std::cout << "\n out of image pruned \n";
					// }

				}
			}
			depthEstimators = depthEstimatorsIn;
			depthEstimatorsIn.clear();

			ROS_WARN("depthEstimators size after %d",int(depthEstimators.size()));
			ROS_WARN("time for estimators %2.4f",float(clock()-updateClock)/CLOCKS_PER_SEC);
			updateClock = clock();

			// assert(depthEstimators.size() > 0);

			if(dkcs.size() > 0)
			{
				std::sort(dkcs.begin(),dkcs.end());
				float dkcMed = dkcs.at(0);
				if (dkcs.size()%2 == 0)
				{
					dkcMed = (dkcs.at(dkcs.size()/2-1)+dkcs.at(dkcs.size()/2))/2.0;
				}
				else
				{
					dkcMed = dkcs.at(dkcs.size()/2);
				}

				dkcHat += kd*(dkcMed - dkcHat);
				pkcHat += kp*(tkcHat*(dkcHat/tkcHat.norm()) - pkcHat);
				dkEstimated = true;
			}

			ROS_WARN("dkcHat %2.1f",dkcHat);
		}
		catch (cv::Exception e)
		{
			ROS_ERROR("update failed");
		}
	}
	else
	{
		//if always searching check if features found, if so then use, if not then just predict using estimated position, orientation, and key position estimate
		clock_t updateClock = clock();
		try
		{
			if (patternFound)
			{
				assert(kPts.size()>0);
				assert(cPts.size()>0);

				cv::Mat inliersG,inliersGp;
				// cv::Mat G = cv::findHomography(kPts, cPts, cv::RANSAC, 2.0, inliersG, 2000, 0.99);//calculate homography using RANSAC
				cv::Mat G = cv::findHomography(kPts, cPts, 0);//calculate homography using RANSAC
				cv::Mat Gp = cv::findHomography(pPts, cPts, cv::RANSAC, 4.0, inliersGp, 2000, 0.99);//calculate homography using RANSAC

				ROS_WARN("time for homog %2.4f",float(clock()-updateClock)/CLOCKS_PER_SEC);
				updateClock = clock();

				// estimate the homography
				Eigen::Vector4f qkc = qkcHat;
				Eigen::Vector3f nk = nkHat;
				Eigen::Vector3f tkc = tkcHat;


				// adjust points based on homography
				// std::vector<cv::Point2f> cPtsAdj = cPts;
				// std::vector<cv::Point2f> cPtsAdj;
				// cv::perspectiveTransform(pPts,cPtsAdj,Gp);
				// cv::perspectiveTransform(pPts,cPtsAdj,G);

				std::cout << "\n start update \n";
				std::cout << "\n qkcHat \n" << qkcHat << std::endl;
				std::cout << "\n nkHat \n" << nkHat << std::endl;
				std::cout << "\n tkcHat \n" << tkcHat << std::endl;

				if (!G.empty())
				{
					try
					{
						//find the solutions
						std::vector<cv::Mat> RkcH,tkcH,nkH;//holds the rotations, translations and normal vetors from decompose homography
						int numberHomogSols = cv::decomposeHomographyMat(G, camMat, RkcH, tkcH, nkH);// Decompose homography

						//check positive depth constraint on the inliers to find the solutions
						std::vector<Eigen::Vector3f> nks;
						std::vector<Eigen::Vector3f> tkcs;
						std::vector<Eigen::Vector4f> qkcs;
						std::vector<float> errors;

						// ROS_WARN("keyframe %d patch %d numberHomogSols %d",keyInd,patchInd,numberHomogSols);

						assert(numberHomogSols>0);
						assert(RkcH.size()>0);
						assert(tkcH.size()>0);
						assert(nkH.size()>0);

						for (int jj = 0; jj < numberHomogSols; jj++)
						{
							//convert normal to eigen
							Eigen::Vector3f nkj(nkH.at(jj).at<double>(0,0),nkH.at(jj).at<double>(1,0),nkH.at(jj).at<double>(2,0));

							// std::cout << "\n\n sol " << jj << std::endl;
							// std::cout << "\n nkjx " << nkj(0) << " nkjy " << nkj(1) << " nkjz " << nkj(2) << std::endl;

							if (nkj.norm() < 0.1)
							{
								nkj(0) = 0.0;
								nkj(1) = 0.0;
								nkj(2) = 1.0;
							}

							//if n^T*[0;0;1] > then solution in front of camera
							Eigen::Vector3f tkcj(tkcH.at(jj).at<double>(0,0),tkcH.at(jj).at<double>(1,0),tkcH.at(jj).at<double>(2,0));
							// std::cout << "\n tkcjx " << tkcj(0) << " tkcjy " << tkcj(1) << " tkcjz " << tkcj(2) << std::endl;
							if ((nkj(2) >= 0.0))
							{
								Eigen::Matrix3f Rkcj;
								for (int hh = 0; hh < 9; hh++)
								{
									Rkcj(hh/3,hh%3) = RkcH.at(jj).at<double>(hh/3,hh%3);
								}
								Eigen::Quaternionf qkcjq(Rkcj);// convert to quaternion
								Eigen::Vector4f qkcj(qkcjq.w(),qkcjq.x(),qkcjq.y(),qkcjq.z());
								qkcj /= qkcj.norm();

								if ((qkcHat + qkcj).norm() < (qkcHat - qkcj).norm())
								{
									qkcj *= -1.0;
								}

								// std::cout << "\n qkcjw " << qkcj(0) << " qkcjx " << qkcj(1) << " qkcjy " << qkcj(2) << " qkcjz " << qkcj(3) << std::endl;
								// std::cout << "\n tkcjx " << tkcj(0) << " tkcjy " << tkcj(1) << " tkcjz " << tkcj(2) << std::endl;
								// std::cout << "\n nkjx " << nkj(0) << " nkjy " << nkj(1) << " nkjz " << nkj(2) << std::endl;

								nks.push_back(nkj);
								tkcs.push_back(tkcj);
								qkcs.push_back(qkcj);
								if (tkcj.norm() > 0.001)
								{
									errors.push_back((qkcHat - qkcj).norm()+(pkcHat-tkcj).norm());
								}
								else
								{
									errors.push_back((qkcHat - qkcj).norm());
								}
							}
						}

						// ROS_WARN("keyframe %d patch %d errors size %d",keyInd,patchInd,int(errors.size()));

						int minqkcsErrorInd = std::distance(errors.begin(),std::min_element(errors.begin(),errors.end()));
						qkc = qkcs.at(minqkcsErrorInd);
						nk = nks.at(minqkcsErrorInd);
						tkc = tkcs.at(minqkcsErrorInd);

						// ROS_WARN("keyframe %d patch %d selected best",keyInd,patchInd);
					}
					catch (cv::Exception e)
					{
						ROS_ERROR("G failed");
					}
				}

				qkcHat += kq*(qkc - qkcHat);
				qkcHat /= qkcHat.norm();
				qckHat = getqInv(qkcHat);

				Eigen::Matrix3f RkcHat = getqRot(qkcHat);
				tkcHat += kt*(tkc - tkcHat);

				float alphank = 0.75;
				nk = (alphank*nk + (1.0-alphank)*(Eigen::Vector3f(0.0,0.0,1.0)));
				nk /= nk.norm();
				nkHat += kn*(nk - nkHat);
				nkHat /= nkHat.norm();

				ROS_WARN("time for update tnq %2.4f",float(clock()-updateClock)/CLOCKS_PER_SEC);
				updateClock = clock();

				Eigen::RowVector3f nkHatT = nkHat.transpose();
				Eigen::Matrix3f H = RkcHat+tkcHat*nkHatT;
				Eigen::RowVector3f H3 = H.block(2,0,1,3);
				Eigen::Matrix3f Gpf = Eigen::Matrix3f::Zero();
				for (int ii = 0; ii < 9; ii++)
				{
					Gpf(ii/3,ii%3) = Gp.at<double>(ii/3,ii%3);
				}

				Eigen::Matrix3f Hp = camMatIf*Gpf*camMatf;
				Eigen::Matrix<float,2,3> Hp12 = Hp.block(0,0,2,3);
				Eigen::RowVector3f Hp3 = Hp.block(2,0,1,3);
				// Eigen::RowVector3f Gp3 = Gpf.block(2,0,1,3);

				//remove the outliers
				std::vector<DepthEstimator*> depthEstimatorsIn;
				std::vector<float> dkcs;
				assert(depthEstimators.size() > 0);
				ROS_WARN("depthEstimators size before %d",int(depthEstimators.size()));
				// ROS_WARN("keyframe %d patch %d kPts size %d",keyInd,patchInd,int(kPts.size()));
				// ROS_WARN("keyframe %d patch %d cPts size %d",keyInd,patchInd,int(cPts.size()));
				Eigen::Vector3f mpi((pPts.at(0).x-cx)/fx,(pPts.at(0).y-cy)/fy,1.0);
				float H3mpi = Hp3*mpi;
				float alphapi = 1.0/H3mpi;
				Eigen::Vector3f mci1((cPts.at(0).x-cx)/fx,(cPts.at(0).y-cy)/fy,1.0);
				// Eigen::Vector3f mci = 0.1*alphapi*Hp*mpi+0.9*mci1;
				// mci /= mci(2);
				Eigen::Vector3f mci = mci1;

				float tlcx = (0.0-cx)/fx;
				float tlcy = (0.0-cy)/fy;
				float brcx = (imageWidth-cx)/fx;
				float brcy = (imageHeight-cy)/fy;

				// std::cout << "\n imageWidth " << imageWidth << std::endl;
				// std::cout << "\n imageHeight " << imageHeight << std::endl;
				// std::cout << "\n tlcx " << tlcx << std::endl;
				// std::cout << "\n tlcy " << tlcy << std::endl;
				// std::cout << "\n brcx " << brcx << std::endl;
				// std::cout << "\n brcy " << brcy << std::endl;
				// std::cout << "\n fx " << fx << std::endl;
				// std::cout << "\n fy " << fy << std::endl;
				// std::cout << "\n cx " << cx << std::endl;
				// std::cout << "\n cy " << cy << std::endl;

				for (int ii = 0; ii < kPts.size(); ii++)
				{
					mpi(0) = (pPts.at(ii).x-cx)/fx;
					mpi(1) = (pPts.at(ii).y-cy)/fy;
					H3mpi = Hp3*mpi;
					alphapi = 1.0/H3mpi;
					mci1(0) = (cPts.at(ii).x-cx)/fx;
					mci1(1) = (cPts.at(ii).y-cy)/fy;
					// mci.segment(0,2) = 0.25*alphapi*Hp12*mpi+0.75*mci1.segment(0,2);
					mci.segment(0,2) = mci1.segment(0,2);

					// std::cout << "\n cPts.at(ii).x " << cPts.at(ii).x << " cPts.at(ii).y " << cPts.at(ii).y << std::endl;
					// std::cout << "\n pPts.at(ii).x " << pPts.at(ii).x << " pPts.at(ii).y " << pPts.at(ii).y << std::endl;
					// std::cout << "\n mci1x " << mci1(0) << " mci1y " << mci1(1) << std::endl;
					// std::cout << "\n mcix " << mci(0) << " mciy " << mci(1) << std::endl;
					// std::cout << "\n inliersGp.at<uchar>(ii) " << int(inliersGp.at<uchar>(ii)) << std::endl;
					// std::cout << "\n (mci(0) >= tlcx) && (mci(0) < brcx) && (mci(1) >= tlcy) && (mci(1) < brcy) " << int((mci(0) >= tlcx) && (mci(0) < brcx) && (mci(1) >= tlcy) && (mci(1) < brcy)) << std::endl;
					float dkcHati = depthEstimators.at(ii)->update(H,mci,nkHatT,tkcHat,RkcHat,vc,wc,t,pkcHat,qkcHat,pics.at(ii),patternLost);
					depthEstimatorsIn.push_back(depthEstimators.at(ii));
					if (depthEstimators.at(ii)->dkKnown)
					{
						dkcs.push_back(dkcHati);
					}
				}

				ROS_WARN("depthEstimators size after %d",int(depthEstimators.size()));
				ROS_WARN("time for estimators %2.4f",float(clock()-updateClock)/CLOCKS_PER_SEC);
				updateClock = clock();

				// assert(depthEstimators.size() > 0);

				if(dkcs.size() > 0)
				{
					std::sort(dkcs.begin(),dkcs.end());
					float dkcMed = dkcs.at(0);
					if (dkcs.size()%2 == 0)
					{
						dkcMed = (dkcs.at(dkcs.size()/2-1)+dkcs.at(dkcs.size()/2))/2.0;
					}
					else
					{
						dkcMed = dkcs.at(dkcs.size()/2);
					}

					dkcHat += kd*(dkcMed - dkcHat);
					pkcHat += kp*(tkcHat*(dkcHat/tkcHat.norm()) - pkcHat);
					pckHat = -rotatevec(pkcHat,qckHat);

					Eigen::RowVector3f tkcHatT = tkcHat.transpose();
					float tktk = tkcHatT*tkcHat;
					float tkpk = tkcHatT*pkcHat;
					if (tktk > 0.001)
					{
						dkHat += kd*(tkpk/tktk - dkHat);
					}
					dkEstimated = true;
				}

				patternLost = false;
				// ROS_WARN("dkcHat %2.2f",dkcHat);
				// ROS_WARN("dkHat %2.2f",dkHat);
			}
			else
			{
					Eigen::Matrix3f RkcHat = getqRot(qkcHat);
					dkcHat += kd*(pkcHat.norm() - dkcHat);
					tkcHat += (dt*(-vc/dkHat - getss(wc)*tkcHat));

					//predict all the estimators forward
					for (int ii = 0; ii < depthEstimators.size(); ii++)
					{
						depthEstimators.at(ii)->predict(RkcHat,vc,wc,t,pkcHat,qkcHat);
					}
			}

		}
		catch (cv::Exception e)
		{
			ROS_ERROR("update failed");
		}
	}

	std::cout << "\n after update \n";
	std::cout << "\n tkcHat \n" << tkcHat << std::endl;
	std::cout << "\n pkcHat \n" << pkcHat << std::endl;
	std::cout << "\n dkHat " << dkHat << std::endl;
	std::cout << "\n dkcHat " << dkcHat << std::endl;


}
