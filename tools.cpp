//***************************
//
// Tools.cpp
//
//***************************
#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() 
{}

Tools::~Tools() 
{}

//*****************************************************************
//
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{

	//cout << "In Tools:CalculateRMSE\n\n";

	/**TODO: Calculate the RMSE here. */
 	 VectorXd rmse(4);
	rmse << 0,0,0,0;
	
    // check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0)
	{
		cout << "ground_truth vector size is 0 or not equal! \n";
		return rmse;
	}

	//cout << "Tools CalculateRMSE : estimations.size() = " << estimations.size() << "\n\n";
	
	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i)
	{
        VectorXd diff = estimations[i] - ground_truth[i];
		//cout << "diff1 = " << diff << "\n\n\n";

		// coefficient-wise multiplication
		diff = diff.array() * diff.array();
		//cout << "diff2 = " << diff << "\n\n\n";

		rmse += diff;
	}

	//calculate the mean
	rmse = rmse / estimations.size();
	
	//calculate the squared root
	rmse = rmse.array().sqrt();
	
	cout << "RMSE = " << rmse(0) << "\t" << rmse(1) << "\t" << rmse(2) << "\t" << rmse(3) << "\n\n\n\n";
	
	//return the result
	return rmse;

}

//********************************************************
//
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) 
{
	
	//cout << "In Tools::CalculateJacobian\n\n";

	/** TODO:  Calculate a Jacobian here. */
	
	MatrixXd Hj(3,4);

	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE 

	//check division by zero
	if( px && py == 0) 
	{
		cout << "CalculateJaobian() - Error - Division by Zero\n ";
		return Hj;
	}
	
	//compute the Jacobian matrix
	float c1 = (px * px) + (py * py);
	float c2 = sqrt( c1 ) ;
	float c3 = c1 * c2;

	Hj << px/c2,				py/c2,					0,		0,
		  -py/c1,				px/c1,					0,		0,
		  py*(vx*py-vy*px)/c3,	px*(vy*px-vx*py)/c3,	px/c2,	py/c2;

	// cout << "Jacobian Hj = " << Hj << "\n\n";
	return Hj;

 	
	
}
