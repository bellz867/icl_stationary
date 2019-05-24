#include <depth_estimator.h>

DepthEstimator::DepthEstimator()
{}

DepthEstimator::DepthEstimator(int depthIndInit, Eigen::Vector3f mInit, ros::Time t, float zminInit, float zmaxInit, Eigen::Vector3f pic, float tauInit)
{
  pik = pic;
  depthInd = depthIndInit;
  mk = mInit;
  mc = mk;
  uk = mk/mk.norm();
  uc = uk;
  lastt = t;
  zmin = zminInit;
  zmax = zmaxInit;
  dkcHat = 0.0;
  dkHat = zmin;
  dcHat = zmin;
  zcHatEKF = zmin;
  zcHatLS = zmin;
  zkHatNorm = zmin;
  zcHatNorm = zmin;
  dkHatICL = zmin;
  dkHatICLExt = zmin;
  dcHatICL = zmin;
  dcHatICLExt = zmin;
  dkcHatICL = 0.0;
  dkcHatICLExt = 0.0;
  tau = tauInit;
  startt = t;
  firstzk = true;
  dkKnown = false;
  uvInt = Eigen::Vector2f::Zero();

  depthEstimatorEKF.initialize(mk.segment(0,2),zmin,zmax);
  depthEstimatorLS.initialize(mk.segment(0,2),zmin,zmax,t);
  depthEstimatorNorm.initialize(mk,zmin,zmax,tau);
  depthEstimatorICL.initialize(uk,zmin,zmax,tau);
  depthEstimatorICLExt.initialize(uk,zmin,zmax,tau,t);
}

float DepthEstimator::update(Eigen::Matrix3f H, Eigen::Vector3f mcMeas, Eigen::RowVector3f nkT, Eigen::Vector3f tkc, Eigen::Matrix3f Rkc, Eigen::Vector3f v, Eigen::Vector3f w, ros::Time t, Eigen::Vector3f pkc, Eigen::Vector4f qkc, Eigen::Vector3f pic, bool patternLost)
{
  float dt = (t - lastt).toSec();
  lastt = t;

  zcHatEKF = depthEstimatorEKF.update(mcMeas.segment(0,2),v,w,dt);

  if (patternLost)
  {
    mc = mcMeas;
  }
  float mcTau = 1.0/(2.0*M_PI*120.0);
  float kmc = dt/(mcTau + dt);
  mc += kmc*(mcMeas - mc);
  mc(2) = 1.0;

  zcHatLS = depthEstimatorLS.update(mc.segment(0,2),v,w,t,dt);

  if (tkc.norm() < 0.001)
  {
    mk = mc;
    pik += kmc*(pic - pik);
    return dkcHat;
  }

  uc = mc/mc.norm();
  Eigen::Vector3f ukc = tkc/tkc.norm();

  // float zkTau = 1.0/(2.0*M_PI*50.0);
  // float kzk = dt/(zkTau + dt);
  float kzk = 20.0;

  Eigen::Vector3f zkzcdkcNorm = depthEstimatorNorm.update(mc,ukc,Rkc,v,w,pkc,t,dt);
  zkHatNorm = zkzcdkcNorm(0);
  zcHatNorm = zkzcdkcNorm(1);
  dkcHatNorm = zkzcdkcNorm(2);

  Eigen::Vector3f dkdcdkcICL = depthEstimatorICL.update(uc,ukc,Rkc,v,pkc,t,dt,patternLost);
  dkHatICL = dkdcdkcICL(0);
  dcHatICL = dkdcdkcICL(1);
  dkcHatICL = dkdcdkcICL(2);

  Eigen::Vector3f dkdcdkcICLExt = depthEstimatorICLExt.update(uc,ukc,Rkc,v,w,pkc,t,dt,patternLost);
  dkHatICLExt = dkdcdkcICLExt(0);
  dcHatICLExt = dkdcdkcICLExt(1);
  dkcHatICLExt = dkdcdkcICLExt(2);

  Eigen::RowVector3f ucT = uc.transpose();
  Eigen::RowVector3f ukcT = ukc.transpose();

  Eigen::Matrix<float,3,2> Yc;
  Eigen::Matrix<float,2,3> YcT;
  Yc.block(0,0,3,1) = uc;
  Yc.block(0,1,3,1) = -ukc;
  YcT.block(0,0,1,3) = ucT;
  YcT.block(1,0,1,3) = -ukcT;
  Eigen::Matrix2f YcYc = YcT*Yc;
  float YcYcdet = YcYc(0,0)*YcYc(1,1) - YcYc(0,1)*YcYc(1,0);
  Eigen::Matrix2f YcYcI;
  YcYcI(0,0) = YcYc(1,1)/YcYcdet;
  YcYcI(0,1) = -YcYc(0,1)/YcYcdet;
  YcYcI(1,0) = -YcYc(1,0)/YcYcdet;
  YcYcI(1,1) = YcYc(0,0)/YcYcdet;

  Eigen::Vector2f zeta = YcYcI*YcT*Rkc*uk;
  float dcDot = -ucT*v;
  float dkcDot = -ukcT*v;
  Eigen::Vector2f uv(dcDot,dkcDot);

  // std::cout << "\n vx " << v(0) << " vy " << v(1) << " vz " << v(2) << std::endl;
  // std::cout << "\n tx " << tkc(0) << " ty " << tkc(1) << " tz " << tkc(2) << std::endl;
  // std::cout << "\n ux " << ukc(0) << " uy " << ukc(1) << " uz " << ukc(2) << std::endl;

  dcHat += (dcDot*dt);
  dkcHat += (dkcDot*dt);

  zetaBuff.push_back(zeta);
  uvBuff.push_back(uv);
  tBuff.push_back(t);
  dtBuff.push_back(dt);
  uvInt += (uv*dt);

  if (v.norm() < 0.05)
  {
    zetaBuff.clear();
    uvBuff.clear();
    tBuff.clear();
    dtBuff.clear();
    uvInt = Eigen::Vector2f::Zero();
  }

  if (pkc.norm() < 0.03)
  {
    zetaBuff.clear();
    uvBuff.clear();
    tBuff.clear();
    dtBuff.clear();
    uvInt = Eigen::Vector2f::Zero();
  }

  if (tBuff.size() > 3)
  {
    while ((tBuff.at(tBuff.size()-1) - tBuff.at(1)).toSec() > tau)
    {
      uvInt -= (uvBuff.at(0)*dtBuff.at(0));
      zetaBuff.pop_front();
      uvBuff.pop_front();
      tBuff.pop_front();
      dtBuff.pop_front();
    }

    Eigen::Vector2f Y = zetaBuff.at(zetaBuff.size()-1) - zetaBuff.at(0);
    Eigen::Vector2f U = uvInt;

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
    bool measgood = (U.norm() > 0.1) && (Y.norm() > 0.1);
    // bool xGood = (fabsf(Ux) > 0.1) && (fabsf(Yx) > 0.1);
    // bool yGood = (fabsf(Uy) > 0.1) && (fabsf(Yy) > 0.1);

    bool dirGoodzBad = false;

    //check the ys
    if (measgood)
    {
      float dk = (Yx*Ux+Yy*Uy)/(Yx*Yx + Yy*Yy);
      if ((dk > zmin) && (dk < zmax))
      {
        dkBuff.push_back(dk);
      }
      else
      {
        dirGoodzBad = true;
      }
    }

    // // //check the xs
    // // if (xGood)
    // // {
    // //   float dk = (Yx*Ux)/(Yx*Yx);
    // //   if ((dk > zmin) && (dk < zmax))
    // //   {
    // //     dkBuff.push_back(dk);
    // //   }
    // //   else
    // //   {
    // //     dirGoodzBad = true;
    // //   }
    // // }
    //
    // //check the ys
    // if (yGood)
    // {
    //   float dk = (Yy*Uy)/(Yy*Yy);
    //   if ((dk > zmin) && (dk < zmax))
    //   {
    //     dkBuff.push_back(dk);
    //   }
    //   else
    //   {
    //     dirGoodzBad = true;
    //   }
    // }

    if (dirGoodzBad)
    {
      zetaBuff.clear();
      uvBuff.clear();
      tBuff.clear();
      dtBuff.clear();
      uvInt = Eigen::Vector2f::Zero();
    }
  }

  if (dkBuff.size() > 1)
  {
    dkKnown = true;
    float dkMed = dkBuff.at(0);

    std::sort(dkBuff.begin(),dkBuff.end());
    if (dkBuff.size()%2 == 0)
    {
      dkMed = (dkBuff.at(dkBuff.size()/2-1)+dkBuff.at(dkBuff.size()/2))/2.0;
    }
    else
    {
      dkMed = dkBuff.at(dkBuff.size()/2);
    }

    // std::cout << "\n zkMed " << zkMed << std::endl;
    // std::cout << "\n zkHat " << zkHat << std::endl;

    if (dkMed < zmin)
    {
      dkMed = zmin;
    }

    if (dkMed > zmax)
    {
      dkMed = zmax;
    }

    // dcHat += kzk*(zeta(0)*dkMed-dcHat);
    // dkcHat += kzk*(zeta(1)*dkMed-dkcHat);
    // dkHat += kzk*(dkMed-dkHat);
    dcHat += (kzk*(zeta(0)*dkMed-dcHat)*dt);
    dkcHat += (kzk*(zeta(1)*dkMed-dkcHat)*dt);
    dkHat += (kzk*(dkMed-dkHat)*dt);
  }

  return dkcHat;
}

void DepthEstimator::predict(Eigen::Matrix3f Rkc, Eigen::Vector3f v, Eigen::Vector3f w, ros::Time t, Eigen::Vector3f pkc, Eigen::Vector4f qkc)
{
  float dt = (t - lastt).toSec();
  lastt = t;

  Eigen::Vector3f mcMeas = pkc + rotatevec(mk*dkHatICL,qkc);
  if (fabsf(mcMeas(2)) > 0.001)
  {
    mc = mcMeas/mcMeas(2);
  }
  else
  {
    mc = mcMeas/(fabsf(mcMeas(2))+0.001);
  }
  mc(2) = 1.0;

  zcHatEKF = depthEstimatorEKF.update(mc.segment(0,2),v,w,dt);

  zcHatLS = depthEstimatorLS.update(mc.segment(0,2),v,w,t,dt);

  uc = mcMeas/mcMeas.norm();
  Eigen::Vector3f ukc = pkc/pkc.norm();

  // float zkTau = 1.0/(2.0*M_PI*50.0);
  // float kzk = dt/(zkTau + dt);
  float kzk = 25.0;

  Eigen::Vector3f zkzcdkcNorm = depthEstimatorNorm.update(mc,ukc,Rkc,v,w,pkc,t,dt);
  zkHatNorm = zkzcdkcNorm(0);
  zcHatNorm = zkzcdkcNorm(1);
  dkcHatNorm = zkzcdkcNorm(2);

  Eigen::Vector3f dkdcdkcICL = depthEstimatorICL.update(uc,ukc,Rkc,v,pkc,t,dt,true);
  dkHatICL = dkdcdkcICL(0);
  dcHatICL = dkdcdkcICL(1);
  dkcHatICL = dkdcdkcICL(2);

  Eigen::RowVector3f ucT = uc.transpose();
  Eigen::RowVector3f ukcT = ukc.transpose();

  Eigen::Matrix<float,3,2> Yc;
  Eigen::Matrix<float,2,3> YcT;
  Yc.block(0,0,3,1) = uc;
  Yc.block(0,1,3,1) = -ukc;
  YcT.block(0,0,1,3) = ucT;
  YcT.block(1,0,1,3) = -ukcT;
  Eigen::Matrix2f YcYc = YcT*Yc;
  float YcYcdet = YcYc(0,0)*YcYc(1,1) - YcYc(0,1)*YcYc(1,0);
  Eigen::Matrix2f YcYcI;
  YcYcI(0,0) = YcYc(1,1)/YcYcdet;
  YcYcI(0,1) = -YcYc(0,1)/YcYcdet;
  YcYcI(1,0) = -YcYc(1,0)/YcYcdet;
  YcYcI(1,1) = YcYc(0,0)/YcYcdet;

  Eigen::Vector2f zeta = YcYcI*YcT*Rkc*uk;
  float dcDot = -ucT*v;
  float dkcDot = -ukcT*v;
  Eigen::Vector2f uv(dcDot,dkcDot);

  // std::cout << "\n vx " << v(0) << " vy " << v(1) << " vz " << v(2) << std::endl;
  // std::cout << "\n tx " << tkc(0) << " ty " << tkc(1) << " tz " << tkc(2) << std::endl;
  // std::cout << "\n ux " << ukc(0) << " uy " << ukc(1) << " uz " << ukc(2) << std::endl;

  dcHat += (dcDot*dt);
  dkcHat += (dkcDot*dt);

  zetaBuff.clear();
  uvBuff.clear();
  tBuff.clear();
  dtBuff.clear();
  uvInt = Eigen::Vector2f::Zero();

  if (tBuff.size() > 3)
  {
    while ((tBuff.at(tBuff.size()-1) - tBuff.at(1)).toSec() > tau)
    {
      uvInt -= (uvBuff.at(0)*dtBuff.at(0));
      zetaBuff.pop_front();
      uvBuff.pop_front();
      tBuff.pop_front();
      dtBuff.pop_front();
    }

    Eigen::Vector2f Y = zetaBuff.at(zetaBuff.size()-1) - zetaBuff.at(0);
    Eigen::Vector2f U = uvInt;

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
    bool measgood = (U.norm() > 0.1) && (Y.norm() > 0.1);
    // bool xGood = (fabsf(Ux) > 0.1) && (fabsf(Yx) > 0.1);
    // bool yGood = (fabsf(Uy) > 0.1) && (fabsf(Yy) > 0.1);

    bool dirGoodzBad = false;

    //check the ys
    if (measgood)
    {
      float dk = (Yx*Ux+Yy*Uy)/(Yx*Yx + Yy*Yy);
      if ((dk > zmin) && (dk < zmax))
      {
        dkBuff.push_back(dk);
      }
      else
      {
        dirGoodzBad = true;
      }
    }

    // // //check the xs
    // // if (xGood)
    // // {
    // //   float dk = (Yx*Ux)/(Yx*Yx);
    // //   if ((dk > zmin) && (dk < zmax))
    // //   {
    // //     dkBuff.push_back(dk);
    // //   }
    // //   else
    // //   {
    // //     dirGoodzBad = true;
    // //   }
    // // }
    //
    // //check the ys
    // if (yGood)
    // {
    //   float dk = (Yy*Uy)/(Yy*Yy);
    //   if ((dk > zmin) && (dk < zmax))
    //   {
    //     dkBuff.push_back(dk);
    //   }
    //   else
    //   {
    //     dirGoodzBad = true;
    //   }
    // }

    if (dirGoodzBad)
    {
      zetaBuff.clear();
      uvBuff.clear();
      tBuff.clear();
      dtBuff.clear();
      uvInt = Eigen::Vector2f::Zero();
    }
  }

  if (dkBuff.size() > 1)
  {
    dkKnown = true;
    float dkMed = dkBuff.at(0);

    std::sort(dkBuff.begin(),dkBuff.end());
    if (dkBuff.size()%2 == 0)
    {
      dkMed = (dkBuff.at(dkBuff.size()/2-1)+dkBuff.at(dkBuff.size()/2))/2.0;
    }
    else
    {
      dkMed = dkBuff.at(dkBuff.size()/2);
    }

    // std::cout << "\n zkMed " << zkMed << std::endl;
    // std::cout << "\n zkHat " << zkHat << std::endl;

    if (dkMed < zmin)
    {
      dkMed = zmin;
    }

    if (dkMed > zmax)
    {
      dkMed = zmax;
    }

    // dcHat += kzk*(zeta(0)*dkMed-dcHat);
    // dkcHat += kzk*(zeta(1)*dkMed-dkcHat);
    // dkHat += kzk*(dkMed-dkHat);
    dcHat += (kzk*(zeta(0)*dkMed-dcHat)*dt);
    dkcHat += (kzk*(zeta(1)*dkMed-dkcHat)*dt);
    dkHat += (kzk*(dkMed-dkHat)*dt);
  }
}
