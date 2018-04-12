#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	unsigned num_states = 4;
	VectorXd rmse(num_states);
	rmse.fill(0);

	if(estimations.size()==0 or ground_truth.size()==0 or estimations.size() != ground_truth.size()){
		std::cerr<<"Invalid input, returning 0"<<std::endl;
		return rmse;
	}
	for (unsigned idx=0; idx< estimations.size(); ++idx){
		VectorXd diff = estimations[idx] - ground_truth[idx];
		diff = diff.array().square();
		rmse += diff;
	}
	rmse /= estimations.size();

	rmse = rmse.array().sqrt();

	return rmse;
}
