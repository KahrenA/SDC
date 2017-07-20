//**********************************************
// kalman_filter.cpp
//**********************************************
#include "kalman_filter.h"
#include <iostream>
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() 
{}

KalmanFilter::~KalmanFilter()
{}

//********************************************************************
//
void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) 
{
	x_ = x_in;	// changes constantly 
  	P_ = P_in;	// same
  	F_ = F_in;	// only dt impacts this 
  	H_ = H_in;	// not same for Lidar and Radar
  	R_ = R_in;  // not same ... see FEKF.cpp
  	Q_ = Q_in;	// same 
}

//**************************
//
void KalmanFilter::Predict() 
{
	//std::cout << "In KF:Predict\n";
	//std::cout << "x_ = " << x_ << "\n\n";
	//std::cout << "P_ = " << P_ << "\n\n";

  	// predict the state
  	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
  	
}

//******************************************
//
void KalmanFilter::Update(const VectorXd &z) 
{
	//std::cout << "IN KF:Update\n";

  	//**TODO:update the state by using Kalman Filter equations
  	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.
size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
  	
}


//*********************************************
//
void KalmanFilter::UpdateEKF(const VectorXd &z) 
{
	//std::cout << "In KF:UpdateEKF\n";

  	//**TODO:update the state by using Extended Kalman Filter equations
  	
  	// Lesson 5.20 indicates To calculate y, we use the equations that map 
	// the predicted location x​′​​ from Cartesian coordinates to polar coordinates:
	// y = z(of radar) - h(x')
	
	if ( (fabs(x_[0]) < 0.0001) && (fabs(x_[1]) < 0.0001) )
	{
		std::cout << "problem in converting to Polar\n\n";
		return;
	}	
	VectorXd h_polar(3);
	h_polar << 0, 0, 0;

	double px = x_[0];
	double py = x_[1];
	double vx = x_[2];
	double vy = x_[3];

	h_polar(0) = sqrt( (px * px) + (py * py));
	if(fabs(px) < 0.001)
		px = 0.001;	

	h_polar(1) = atan2(py, px);
	
	if(fabs(h_polar(1)) < 0.001)
		h_polar(1) = 0.001;

	h_polar(2) = (px * vx + py * vy) / h_polar(0);

	
	//---------------------------------------------
	VectorXd y = z - h_polar;

	// Let's normalize the y which is polar itself
	if(y(1) > M_PI)
	{	
		while(y(1) > M_PI)
		{
			std::cout << "y(1) > M_PI\n";
			y(1) -= 2*M_PI;
		}	
	}
	else if (y(1) < -M_PI)
	{
		while (y(1) < -M_PI)
		{
			std::cout << "y(1) < -M_PI\n";
			y(1) += 2 * M_PI;
		}
	}


	MatrixXd Ht = H_.transpose();
	
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
  	
}
