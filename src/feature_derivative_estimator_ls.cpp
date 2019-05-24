#include <feature_derivative_estimator.h>

FeatureDerivativeEstimator::FeatureDerivativeEstimator()
{
}

void FeatureDerivativeEstimator::initialize(float maxdtInit, int bufferSizeInit)
{
	bufferFull = false;
	bufferSize = bufferSizeInit;
	A = Eigen::MatrixXf::Zero(2*bufferSize,4);
	Eigen::MatrixXf I2 = Eigen::MatrixXf::Zero(2,4);
	I2(0,0) = 1.0;
	I2(1,2) = 1.0;
	for (int ii = 0; ii < bufferSize; ii++)
	{
		A.block(2*ii,0,2,4) = I2;
	}
	maxdt = maxdtInit;
}

Eigen::Vector2f FeatureDerivativeEstimator::update(Eigen::Vector2f newMeasure, ros::Time newTime)
{

	// Picture courtesy of Anup
	// Setting up least squares problem A*theta = P. theta is made up of the coefficients for the best fit line,
	// e.g., X = Mx*T + Bx, Y = My*t + By, Z = Mz*t + Bz. Velocity is estimated as the slope of the best fit line, i.e., Vx = Mx, Vy = My, Vz = Mz.
	// Each block of data is arranged like this:
	// [Xi]     [1, Ti,  0,  0,  0,  0] * [Bx]
	// [Yi]  =  [0,  0,  1, Ti,  0,  0]   [Mx]
	// [Zi]     [0,  0,  0,  0,  1, Ti]   [By]
	//  \/      \_____________________/   [My]
	//  Pi                 \/             [Bz]
	//                     Ai             [Mz]
	//                                     \/
	//                                   theta
	//
	// and then data is all stacked like this, where n is the buffer size:
	// [P1]     [A1] * [Bx]
	// [P2]  =  [A2]   [Mx]
	//  :        :     [By]
	// [Pn]     [An]   [My]
	//                 [Bz]
	//                 [Mz]

	//Fill buffers
	timeBuff.push_back(newTime);
	stateBuff.push_back(newMeasure);

	//If the index has rolled over once, the buffer is full
	if (!bufferFull && (timeBuff.size() >= bufferSize))
	{
		bufferFull = true;
	}

	while (timeBuff.size() > bufferSize)
	{
		timeBuff.pop_front();
		stateBuff.pop_front();
	}

	Eigen::Vector2f stateDerivative = Eigen::Vector2f::Zero();//initialize state derivative
	if (bufferFull)
	{
		Eigen::VectorXf b = Eigen::VectorXf::Zero(2*bufferSize);
		// Solve LLS for best fit line parameters
		for (int ii = 0; ii < bufferSize; ii++)
		{
			float dtii = (timeBuff.at(ii) - timeBuff.at(0)).toSec();
			A(2*ii,1) = dtii;
			A(2*ii+1,3) = dtii;
			b.segment(2*ii,2) = stateBuff.at(ii);
		}

		Eigen::JacobiSVD<Eigen::MatrixXf> svdA(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::VectorXf x = svdA.solve(b);
		stateDerivative(0) = x(1);
		stateDerivative(1) = x(3);
	}

	return stateDerivative;//return state derivative
}
