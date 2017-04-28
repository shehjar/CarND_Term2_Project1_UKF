#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if ((estimations.size() < 0) &&(estimations.size() != ground_truth.size())){
    std::cout<<"RMSE cannot be calculated, please check vector sizes\n";
    return rmse;
  }

  for (size_t i = 0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // Average rmse
  rmse = rmse/estimations.size();

  // Taking square root of individual elements
  rmse = rmse.array().sqrt();

  return rmse;
}
