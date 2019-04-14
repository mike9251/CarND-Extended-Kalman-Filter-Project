#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
    
    VectorXd rmse(4);//ground_truth[0].size());
    rmse << 0, 0, 0, 0;
    /*for(int i = 0; i < ground_truth[0](0).size(); i ++)
    {
    	rmse(i) = 0;
    }*/

    if(estimations.size() != ground_truth.size() || estimations.size() == 0){
    	return rmse;
    }

    for(int i = 0; i < ground_truth.size(); i++){
    	VectorXd diff = estimations[i] - ground_truth[i];

    	diff = diff.array() * diff.array();

    	rmse += diff;// / ground_truth.size();
    }

    rmse = rmse / ground_truth.size();

    return rmse.array().sqrt();


    /*for( int i = 0; i < estimations[0].size(); i++)
    {
    	float rmse = 0;
    	for(int j = 0; i < estimations.size(); j++)
    	{
    		rmse += (estimations[j](i) - ground_truth[j](i)) * (estimations[j](i) - ground_truth[j](i));
    	}
    	rmse_v(i) = rmse / estimations.size();
    }*/
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
	MatrixXd Hj(3, 4);

	// Unroll state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);
	
	// Pre-compute some term which recur in the Jacobian
	float c1 = px * px + py * py;
	float c2 = sqrt(c1);
	float c3 = c1 * c2;

	// Sanity check to avoid division by zero
	if (std::abs(c1) < 0.0001) {
		std::cout << "Error in CalculateJacobian. Division by zero." << std::endl;
		return Hj;
	}

	// Actually compute Jacobian matrix
	Hj << (px / c2),				(py / c2),					0,			0,
		-(py / c1),					(px / c1),					0,			0,
		py * (vx*py - vy*px) / c3,	px * (vy*px - vx*py) / c3,	px / c2,	py / c2;

	return Hj;
}
