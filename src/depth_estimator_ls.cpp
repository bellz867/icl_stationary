#include <depth_estimator_ls.h>

DepthEstimatorLS::DepthEstimatorLS()
{}

void DepthEstimatorLS::initialize(Eigen::Vector2f m, float zminInit, float zmaxInit, ros::Time t)
{
	zmin = zminInit;
	zmax = zmaxInit;
	mzHat = 1/zmin;
	mDotEsimator.initialize(0.05,5);
	Eigen::Vector2f mDot = mDotEsimator.update(m,t);
}

float DepthEstimatorLS::update(Eigen::Vector2f m, Eigen::Vector3f v, Eigen::Vector3f w, ros::Time t, float dt)
{
	float mx = m(0);
	float my = m(1);
	float vx = v(0);
	float vy = v(1);
	float vz = v(2);
	float wx = w(0);
	float wy = w(1);
	float wz = w(2);

	Eigen::Vector2f mDot = mDotEsimator.update(m,t);
	Eigen::Vector2f zeta = getzeta(mx,my,wx,wy,wz);
	Eigen::Vector2f rho = getrho(mx,my,vx,vy,vz);
	Eigen::RowVector2f rhoT = rho.transpose();

	float kmz = 100.0;

	float mzHatDot = getmzDot(mx,my,mzHat,vz,wx,wy)+kmz*(float(rhoT*(mDot-zeta)) - float(rhoT*rho)*mzHat);
	mzHat += (mzHatDot*dt);

	if (mzHat > 1.0/zmin)
	{
		mzHat = 1.0/zmin;
	}

	if (mzHat < 1.0/zmax)
	{
		mzHat = 1.0/zmax;
	}
	return (1.0/mzHat);
}

float DepthEstimatorLS::getmzDot(float mx, float my, float mz, float vz, float wx, float wy)
{
  return (mz*mz*vz - mx*mz*wy + my*mz*wx);
}

Eigen::Vector2f DepthEstimatorLS::getzeta(float mx, float my, float wx, float wy, float wz)
{
  float zetax = my*wz - wy - mx*mx*wy + mx*my*wx;
  float zetay = -mx*wz + wx - mx*my*wy + my*my*wx;
  return Eigen::Vector2f(zetax,zetay);
}

Eigen::Vector2f DepthEstimatorLS::getrho(float mx, float my, float vx, float vy, float vz)
{
	float rhox = -vx + mx*vz;
	float rhoy = -vy + my*vz;
  return Eigen::Vector2f(rhox,rhoy);
}
