//**************************
// FusionEKF.cpp
//***************************
#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//********************
// Constructor 
//
FusionEKF::FusionEKF()
{
	is_initialized_ = false;

  	previous_timestamp_ = 0;

  	// initializing matrices
 	R_laser_ = MatrixXd(2, 2);
  	R_radar_ = MatrixXd(3, 3);
  	H_laser_ = MatrixXd(2, 4);
  	Hj_ = MatrixXd(3, 4);

	H_laser_ << 1, 0, 0, 0,
		 		0, 1, 0, 0;

  	//measurement covariance matrix - laser
  	R_laser_ << 0.0225, 0,
  		        0, 		0.0225;

  	//measurement covariance matrix - radar
  	R_radar_ << 0.09,  0,      0,
  	      	    0,     0.0009, 0,
  	      	    0,     0,      0.09;

  	//**	TODO:Finish initializing the FusionEKF.
    //  Set the process and measurement noises
  	noise_ax = 9;
	noise_ay = 9; 
  	
  	// Process noise/covariance matrix
	Q_matrix = MatrixXd(4, 4);
	Q_matrix << 0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0;

	P_matrix = MatrixXd(4, 4);
	P_matrix = MatrixXd::Zero(4, 4);
	P_matrix << 1, 0, 0, 0,
				0 ,1, 0, 0,
				0, 0, 1000, 0,
				0, 0, 0, 1000;

	F_matrix = MatrixXd(4, 4);
	F_matrix << 1, 0, 1, 0,
				0, 1, 0, 1,
				0, 0, 1, 0,
				0, 0, 0, 1;
  		 
  	// Measurement noise/covariance matrix
  	

}

//*****************
// Destructor
//
FusionEKF::~FusionEKF() 
{
}

//****************************************************************************
//     		        
void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) 
{
		
	//-------------------
	//  Initialization
	//-------------------
  	if (!is_initialized_) 
  	{
  		//	TODO:
  	    // Initialize the state ekf_.x_ with the first measurement.
  	  	// Create the covariance matrix.
	
		ekf_.H_ = H_laser_;
      
		// first measurement
	    cout << "EKF: Init" << endl;
	    ekf_.x_ = VectorXd(4);
	    ekf_.x_ << 1, 1, 1, 1;

	    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
	    {
			//	Convert radar from polar to cartesian coordinates and initialize state
			// px = py = phi; vx = vy = 0 per forum notes 
			ekf_.x_ << measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[1], 0, 0;  	
	
			previous_timestamp_ = measurement_pack.timestamp_;
 		
    	}
    	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
 		{
 			// Initialize state.
 			// set the state with the initial location and zero velocity
			ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

			previous_timestamp_ = measurement_pack.timestamp_;
 		}

		// Set ekf_ members initial values 
		ekf_.F_ = F_matrix;
		ekf_.P_ = P_matrix; 
		ekf_.Q_ = Q_matrix;

    	// done initializing, no need to predict or update
    	is_initialized_ = true;
		cout << "Init Done \n\n";
    	return;
  	}


  	//---------------
  	// Prediction
  	//---------------
	// cout << "EKF Predict \n\n";
	
  	//**TODO:
  	// Update the state transition matrix F according to the new elapsed time.
  	// Time is measured in seconds.
  	// Update the process noise covariance matrix.
  	// Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
  	//
	double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	
	previous_timestamp_ = measurement_pack.timestamp_;

	double dt_2 = dt * dt;
	double dt_3 = dt_2 * dt;
	double dt_4 = dt_3 * dt;

	//Modify the F matrix so that the time is integrated
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;
		
	
	//set the process covariance matrix Q
	Q_matrix << dt_4/4*noise_ax, 	0, 					dt_3/2*noise_ax, 	0,
				0, 					dt_4/4*noise_ay,	0,					dt_3/2*noise_ay,
			   	dt_3/2*noise_ax, 	0, 					dt_2*noise_ax, 		0,
			   	0, 					dt_3/2*noise_ay, 	0, 					dt_2*noise_ay;

	// cout << "Q_matrix = " << Q_matrix << "\n\n";

	// --------------------Predict ----------------------------------
   	ekf_.Init( ekf_.x_, ekf_.P_, ekf_.F_, ekf_.H_, ekf_.R_, Q_matrix);   // ekf_. means leave value as is
 	ekf_.Predict();

	//cout << "Proc_meas - Predict - x = " << ekf_.x_ << "\n\n";
	//cout << "Proc_meas - Predict - P = " << ekf_.P_ << "\n\n";

  	//=================
  	// Update
  	//=================

  	// TODO:
  	// Use the sensor type to perform the update step.
  	// Update the state and covariance matrices.
  	

	//VectorXd x_in = VectorXd(4);
  	//-------------------------------------------------------------
  	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
	{
	   	
	    // Get Hj Jacobian 
  		Hj_ = tools.CalculateJacobian( ekf_.x_ );	
  		
		//------------- Radar updates ---------------------------
	   	ekf_.Init( ekf_.x_, ekf_.P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_); 
		ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  	} 
	else
	{
		// -----------------Laser updates-------------------------
 		ekf_.Init( ekf_.x_, ekf_.P_, ekf_.F_, H_laser_, R_laser_, ekf_.Q_);
 		ekf_.Update(measurement_pack.raw_measurements_);

  	}

  	// print the output
  	//cout << "Proc_meas - Update x_ = " << ekf_.x_ << endl;
	//cout << "Proc_meas - Update P_ = " << ekf_.P_ << endl;
}
