#include <depth_estimator_ekf.h>

DepthEstimatorEKF::DepthEstimatorEKF()
{
	H = Eigen::Matrix<float,2,3>::Zero();
	H(0,0) = 1.0;
	H(1,1) = 1.0;
	HT = H.transpose();

	Q = Eigen::Matrix3f::Zero();
	// Q(0,0) = 0.005;
	// Q(1,1) = 0.005;
	// Q(2,2) = 0.005;
	//
	P = Eigen::Matrix3f::Zero();
	// P(0,0) = 0.000005;
	// P(1,1) = 0.000005;
	// P(2,2) = 1.5;
	//
	// R = Eigen::Matrix2f::Zero();
	// R(0,0) = 0.000005;
	// R(1,1) = 0.000005;

	// feature variance
	float rr = 0.00001;
	P(0,0) = 1.0*rr;//covariance
	P(1,1) = 1.0*rr;//covariance
	P(2,2) = 1.5;//covariance
	Q(0,0) = 100.0*rr;//process covariance
	Q(1,1) = 100.0*rr;//process covariance
	Q(2,2) = 100000.0*rr;//process covariance
	R = rr*Eigen::Matrix2f::Identity();//measurment covariancece
}

void DepthEstimatorEKF::initialize(Eigen::Vector2f m, float zminInit, float zmaxInit)
{
	zmin = zminInit;
	zmax = zmaxInit;
	xHat(0) = m(0);
	xHat(1) = m(1);
	xHat(2) = 1.0/zmin;
}

//update the kalman
float DepthEstimatorEKF::update(Eigen::Vector2f m, Eigen::Vector3f v, Eigen::Vector3f w, float dt)
{
	float mx = m(0);
	float my = m(1);
	float vx = v(0);
	float vy = v(1);
	float vz = v(2);
	float wx = w(0);
	float wy = w(1);
	float wz = w(2);

	//predict
	Eigen::Vector3f xHatDot = getxDot(xHat(0),xHat(1),xHat(2),vx,vy,vz,wx,wy,wz);
	Eigen::Matrix3f F = getF(xHat(0),xHat(1),xHat(2),vx,vy,vz,wx,wy,wz);
	// Eigen::Matrix3f PDot = A*P + P*A.transpose() + Q;

	// std::cout << "\n dt " << dt << std::endl;
	// std::cout << "\n xHat \n" << xHat << std::endl;
	// std::cout << "\n xHatDot \n" << xHatDot << std::endl;
	// std::cout << "\n PDot \n" << PDot << std::endl;
	// std::cout << "\n F \n" << F << std::endl;
	// std::cout << "\n F.transpose() \n" << F.transpose() << std::endl;
	// std::cout << "\n Q \n" << Q << std::endl;

	xHat += (xHatDot*dt);
	// P += (PDot*dt);
	P = (F*P*F.transpose() + Q);
	// std::cout << "\n P \n" << P << std::endl;

	// std::cout << "\n xHat \n" << xHat << std::endl;
	// std::cout << "\n P \n" << P << std::endl;

	//measurement update
	Eigen::Vector2f meas(mx,my);
	Eigen::Vector2f measHat(xHat(0),xHat(1));
	Eigen::Matrix2f argK = H*P*HT + R;
	float argKdet = argK(0,0)*argK(1,1) - argK(0,1)*argK(1,0);
	Eigen::Matrix2f argKI;
	argKI(0,0) = argK(1,1)/argKdet;
	argKI(0,1) = -argK(0,1)/argKdet;
	argKI(1,0) = -argK(1,0)/argKdet;
	argKI(1,1) = argK(0,0)/argKdet;
	Eigen::Matrix<float,3,2> K = P*HT*argKI;
	xHat += K*(meas - measHat);
	P = (Eigen::Matrix3f::Identity() - K*H)*P;

	if (xHat(2) > 1.0/zmin)
	{
		xHat(2) = 1.0/zmin;
	}

	if (xHat(2) < 1.0/zmax)
	{
		xHat(2) = 1.0/zmax;
	}

	return 1.0/xHat(2);
}

Eigen::Vector3f DepthEstimatorEKF::getxDot(float mx, float my, float mz, float vx, float vy, float vz, float wx, float wy, float wz)
{
	Eigen::Vector3f xDot;
  xDot(0) = -mz*vx + mx*mz*vz + my*wz - wy - mx*mx*wy + mx*my*wx;
  xDot(1) = -mz*vy + my*mz*vz - mx*wz + wx - mx*my*wy + my*my*wx;
  xDot(2) = mz*mz*vz - mx*mz*wy + my*mz*wx;
  return xDot;
}

Eigen::Matrix3f DepthEstimatorEKF::getF(float mx, float my, float mz, float vx, float vy, float vz, float wx, float wy, float wz)
{
	Eigen::Matrix3f F;
  // F(0,0) = mz*vz - 2.0*mx*wy + my*wx;
  // F(0,1) = wz + mx*wx;
  // F(0,2) = -vx + mx*vz;
  // F(1,0) = -wz - my*wy;
  // F(1,1) = mz*vz - mx*wy + 2.0*my*wx;
  // F(1,2) = -vy + my*vz;
  // F(2,0) = -mz*wy;
  // F(2,1) = mz*wx;
  // F(2,2) = 2.0*mz*vz - mx*wy + my*wx;
	F(0,0) = 1.0 + mz*vz - 2.0*mx*wy + my*wx;
	F(0,1) = wz + mx*wx;
	F(0,2) = -vx + mx*vz;
	F(1,0) = -wz - my*wy;
	F(1,1) = 1.0 + mz*vz - mx*wy + 2.0*my*wx;
	F(1,2) = -vy + my*vz;
	F(2,0) = -mz*wy;
	F(2,1) = mz*wx;
	F(2,2) = 1.0 + 2.0*mz*vz - mx*wy + my*wx;
  return F;
}
