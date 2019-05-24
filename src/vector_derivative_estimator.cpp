#include <vector_derivative_estimator.h>

VectorDerivativeEstimator::VectorDerivativeEstimator()
{
}

void VectorDerivativeEstimator::initialize()
{
	xHat = Eigen::Matrix<float,6,1>::Zero();
	firstUpdate = true;
	P = Eigen::Matrix<float,6,6>::Zero();
	Q = Eigen::Matrix<float,6,6>::Zero();

	// feature variance
	R = 0.0001*Eigen::Matrix3f::Identity();//measurment covariance
	P.block(0,0,3,3) = 1.0*R;//covariance
	Q.block(0,0,3,3) = 1.0*R;//process covariance

	// flow variance
	P.block(3,3,3,3) = 1000.0*R;//covariance
	Q.block(3,3,3,3) = 1000.0*R;//process covariance

	//process jacobian
	F = Eigen::Matrix<float,6,6>::Identity();

	//measruement jacobian
	H = Eigen::Matrix<float,3,6>::Zero();
	H.block(0,0,3,3) = Eigen::Matrix3f::Identity();
	HT = H.transpose();
}

Eigen::Vector3f VectorDerivativeEstimator::update(Eigen::Vector3f newMeasure, ros::Time newTime)
{
	ros::Time t = newTime;
	Eigen::Vector3f z = newMeasure;

	if (firstUpdate)
	{
		xHat.segment(0,3) = z;
		tLast = t;
		firstUpdate = false;
	}

	float dt = (t-tLast).toSec();
	tLast = t;

	//predict
	F.block(0,3,3,3) = dt*Eigen::Matrix3f::Identity();
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
	Eigen::Matrix3f argK = H*P*HT + R;
	Eigen::JacobiSVD<Eigen::MatrixXf> svdargK(argK, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::Matrix3f argKI = svdargK.solve(Eigen::Matrix3f::Identity());
	// std::cout << "\n argKI \n" << argKI <<std::endl;
	Eigen::Matrix<float,6,3> K = P*HT*argKI;
	// std::cout << "\n K \n" << K <<std::endl;

	// std::cout << "\n z \n" << z <<std::endl;
	// std::cout << "\n 88888 xHat ls \n" << xHat <<std::endl;
	// std::cout << "\n z-xHat \n" << (z-xHat.segment(0,7)) <<std::endl;
	// std::cout << "\n K*(z-xHat) \n" << (K*(z-xHat.segment(0,7))).segment(0,7) <<std::endl;

	// std::cout << std::endl << xHat.segment(0,7) << std::endl;
	// std::cout << std::endl << P << std::endl;
	xHat += (K*(z-xHat.segment(0,3)));

	// std::cout << std::endl << "----------------" << std::endl;
	// xHat.segment(3,4) /= xHat.segment(3,4).norm();
	P = (Eigen::Matrix<float,6,6>::Identity() - K*H)*P;

	return xHat.segment(3,3);
}

Eigen::Matrix<float,6,1> VectorDerivativeEstimator::xDot(Eigen::Matrix<float,6,1> x)
{
	Eigen::Matrix<float,6,1> xDot = Eigen::Matrix<float,6,1>::Zero();
	xDot.segment(0,3) = x.segment(3,3);
	return xDot;
}
