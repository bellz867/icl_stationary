#include <depth_estimator_norm.h>

DepthEstimatorNorm::DepthEstimatorNorm()
{}

void DepthEstimatorNorm::initialize(Eigen::Vector3f mInit, float zminInit, float zmaxInit, float tauInit)
{
  mk = mInit;
  mc = mk;
  zmin = zminInit;
  zmax = zmaxInit;
  dkcHat = 0.0;
  zkHat = zmin;
  zkHat = zmin;
  tau = tauInit;
  firstzk = true;
  zkKnown = false;
  vInt = Eigen::Vector3f::Zero();
  psizetaInt = Eigen::Vector3f::Zero();
}

Eigen::Vector3f DepthEstimatorNorm::update(Eigen::Vector3f mc, Eigen::Vector3f ukc, Eigen::Matrix3f Rkc, Eigen::Vector3f v, Eigen::Vector3f w, Eigen::Vector3f pkc, ros::Time t, float dt)
{
  float kzk = 25.0;

  Eigen::RowVector3f mcT = mc.transpose();
  Eigen::RowVector3f ukcT = ukc.transpose();
  Eigen::Vector3f psi = getss(w)*mc;

  float zcDot = -v(2)-psi(2)*zcHat;
  float dkcDot = -ukcT*v;

  zcHat += (zcDot*dt);
  dkcHat += (dkcDot*dt);

  Eigen::Matrix<float,3,2> Yc;
  Eigen::Matrix<float,2,3> YcT;
  Yc.block(0,0,3,1) = mc;
  Yc.block(0,1,3,1) = -ukc;
  YcT.block(0,0,1,3) = mcT;
  YcT.block(1,0,1,3) = -ukcT;
  Eigen::Matrix2f YcYc = YcT*Yc;
  float YcYcdet = YcYc(0,0)*YcYc(1,1) - YcYc(0,1)*YcYc(1,0);
  Eigen::Matrix2f YcYcI;
  YcYcI(0,0) = YcYc(1,1)/YcYcdet;
  YcYcI(0,1) = -YcYc(0,1)/YcYcdet;
  YcYcI(1,0) = -YcYc(1,0)/YcYcdet;
  YcYcI(1,1) = YcYc(0,0)/YcYcdet;

  Eigen::Vector2f zeta = YcYcI*YcT*Rkc*mk;

  mzetaBuff.push_back(mc*zeta(0));
  psizetaBuff.push_back(psi*zeta(0));
  vBuff.push_back(v);
  tBuff.push_back(t);
  dtBuff.push_back(dt);
  vInt += (v*dt);
  psizetaInt += (psi*zeta(0)*dt);

  if (v.norm() < 0.05)
  {
    mzetaBuff.clear();
    psizetaBuff.clear();
    vBuff.clear();
    tBuff.clear();
    dtBuff.clear();
    vInt = Eigen::Vector3f();
    psizetaInt = Eigen::Vector3f();
  }

  if (pkc.norm() < 0.03)
  {
    mzetaBuff.clear();
    psizetaBuff.clear();
    vBuff.clear();
    tBuff.clear();
    dtBuff.clear();
    vInt = Eigen::Vector3f();
    psizetaInt = Eigen::Vector3f();
  }

  if (tBuff.size() > 3)
  {
    while ((tBuff.at(tBuff.size()-1) - tBuff.at(1)).toSec() > tau)
    {
      vInt -= (vBuff.at(0)*dtBuff.at(0));
      psizetaInt -= (psizetaBuff.at(0)*dtBuff.at(0));
      mzetaBuff.pop_front();
      psizetaBuff.pop_front();
      vBuff.pop_front();
      tBuff.pop_front();
      dtBuff.pop_front();
    }

    Eigen::Vector3f Y = mzetaBuff.at(mzetaBuff.size()-1) - mzetaBuff.at(0) + psizetaInt;
    Eigen::Vector3f U = -vInt;

    float Yx = Y(0);
    float Ux = U(0);

    float Yy = Y(1);
    float Uy = U(1);

    // std::cout << "\n mzetaBuffInt \n" << mzetaBuffInt << std::endl;
    // std::cout << "\n psizetaBuffInt \n" << psizetaBuffInt << std::endl;

    // std::cout << "\n (Ux/Yx) " << (Ux/Yx) << std::endl;
    // std::cout << "\n Yx " << Yx << std::endl;
    // std::cout << "\n Ux " << Ux << std::endl;
    //
    // std::cout << "\n (Uy/Yy)  " << (Uy/Yy) << std::endl;
    // std::cout << "\n Yy " << Yy << std::endl;
    // std::cout << "\n Uy " << Uy << std::endl;

    //check which estimates are good
    bool xGood = (fabsf(Ux) > 0.1) && (fabsf(Yx) > 0.1);
    bool yGood = (fabsf(Uy) > 0.1) && (fabsf(Yy) > 0.1);

    bool dirGoodzBad = false;

    //check the xs
    if (xGood)
    {
      float zk = (Yx*Ux)/(Yx*Yx);
      if ((zk > zmin) && (zk < zmax))
      {
        zkBuff.push_back(zk);
      }
      else
      {
        dirGoodzBad = true;
      }
    }

    //check the ys
    if (yGood)
    {
      float zk = (Yy*Uy)/(Yy*Yy);
      if ((zk > zmin) && (zk < zmax))
      {
        zkBuff.push_back(zk);
      }
      else
      {
        dirGoodzBad = true;
      }
    }

    if (dirGoodzBad)
    {
      mzetaBuff.clear();
      psizetaBuff.clear();
      vBuff.clear();
      tBuff.clear();
      dtBuff.clear();
      vInt = Eigen::Vector3f();
      psizetaInt = Eigen::Vector3f();
    }
  }

  if (zkBuff.size() > 1)
  {
    zkKnown = true;
    float zkMed = zkBuff.at(0);

    std::sort(zkBuff.begin(),zkBuff.end());
    if (zkBuff.size()%2 == 0)
    {
      zkMed = (zkBuff.at(zkBuff.size()/2-1)+zkBuff.at(zkBuff.size()/2))/2.0;
    }
    else
    {
      zkMed = zkBuff.at(zkBuff.size()/2);
    }

    // std::cout << "\n zkMed " << zkMed << std::endl;
    // std::cout << "\n zkHat " << zkHat << std::endl;

    if (zkMed < zmin)
    {
      zkMed = zmin;
    }

    if (zkMed > zmax)
    {
      zkMed = zmax;
    }

    zkHat += (kzk*(zkMed-zkHat)*dt);
    zcHat += (kzk*(zeta(0)*zkMed-zcHat)*dt);
    dkcHat += (kzk*(zeta(1)*zkMed-dkcHat)*dt);
  }

  return Eigen::Vector3f(zkHat,zcHat,dkcHat);
}
