#include <feature_derivative_estimator.h>

FeatureDerivativeEstimator::FeatureDerivativeEstimator()
{
}

void FeatureDerivativeEstimator::initialize(float maxdtInit, int bufferSizeInit)
{
	xHat = Eigen::Matrix<float,4,1>::Zero();
	firstUpdate = true;
	P = Eigen::Matrix<float,4,4>::Zero();
	Q = Eigen::Matrix<float,4,4>::Zero();

	// feature variance
	R = 0.0001*Eigen::Matrix2f::Identity();//measurment covariance
	P.block(0,0,2,2) = 1.0*R;//covariance
	Q.block(0,0,2,2) = 1.0*R;//process covariance

	// flow variance
	P.block(2,2,2,2) = 1000.0*R;//covariance
	Q.block(2,2,2,2) = 1000.0*R;//process covariance


	//process jacobian
	F = Eigen::Matrix<float,4,4>::Identity();

	//measruement jacobian
	H = Eigen::Matrix<float,2,4>::Zero();
	H.block(0,0,2,2) = Eigen::Matrix2f::Identity();
	HT = H.transpose();
}

Eigen::Vector2f FeatureDerivativeEstimator::update(Eigen::Vector2f newMeasure, ros::Time newTime)
{
	ros::Time t = newTime;
	Eigen::Vector2f z = newMeasure;

	if (firstUpdate)
	{
		xHat.segment(0,2) = z;
		tLast = t;
		firstUpdate = false;
	}

	float dt = (t-tLast).toSec();
	tLast = t;

	//predict
	F.block(0,2,2,2) = dt*Eigen::Matrix2f::Identity();
	// std::cout << std::endl << "+++++++++++++++++" << std::endl;
	// std::cout << "\n xHat \n" << xHat <<std::endl;
	// std::cout << "\n xDot(xHat) \n" << xDot(xHat) <<std::endl;
	xHat += (xDot(xHat)*dt);
	// xHat.segment(3,4) /= xHat.segment(3,4).norm();
	// // P += ((F*P + P*F.transpose() + Q)*dt);
	P = (F*P*F.transpose() + Q);
	//
	// std::cout << "\n F \n" << F <<std::endl;
	// std::cout << "\n P \n" << P <<std::endl;

	// Eigen::Matrix2f argK = H*P*HT+R;
	// std::cout << "\n argK \n" << argK <<std::endl;
	// Eigen::JacobiSVD<Eigen::MatrixXf> svdargK(argK, Eigen::ComputeThinU | Eigen::ComputeThinV);
	// Eigen::Matrix2f argKI = svdargK.solve(Eigen::Matrix2f::Identity());
	Eigen::Matrix2f argK = H*P*HT + R;
	float argKdet = argK(0,0)*argK(1,1) - argK(0,1)*argK(1,0);
	Eigen::Matrix2f argKI;
	argKI(0,0) = argK(1,1)/argKdet;
	argKI(0,1) = -argK(0,1)/argKdet;
	argKI(1,0) = -argK(1,0)/argKdet;
	argKI(1,1) = argK(0,0)/argKdet;
	// std::cout << "\n argKI \n" << argKI <<std::endl;
	Eigen::Matrix<float,4,2> K = P*HT*argKI;
	// std::cout << "\n K \n" << K <<std::endl;

	// std::cout << "\n z \n" << z <<std::endl;
	// std::cout << "\n 88888 xHat ls \n" << xHat <<std::endl;
	// std::cout << "\n z-xHat \n" << (z-xHat.segment(0,7)) <<std::endl;
	// std::cout << "\n K*(z-xHat) \n" << (K*(z-xHat.segment(0,7))).segment(0,7) <<std::endl;

	// std::cout << std::endl << xHat.segment(0,7) << std::endl;
	// std::cout << std::endl << P << std::endl;
	xHat += (K*(z-xHat.segment(0,2)));

	// std::cout << std::endl << "----------------" << std::endl;
	// xHat.segment(3,4) /= xHat.segment(3,4).norm();
	P = (Eigen::Matrix4f::Identity() - K*H)*P;

	return xHat.segment(2,2);
}

Eigen::Vector4f FeatureDerivativeEstimator::xDot(Eigen::Vector4f x)
{
	Eigen::Vector4f xDot = Eigen::Vector4f::Zero();
	xDot.segment(0,2) = x.segment(2,2);
	return xDot;
}
