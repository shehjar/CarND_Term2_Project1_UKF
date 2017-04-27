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
  int dataCount = estimations.size();
  int chkCount = ground_truth.size();
  //std::cout<< "Data count = "<<dataCount<<std::endl;
  //std::cout<<"ground truth count = "<< chkCount<< std::endl;
  if((dataCount != chkCount)||(dataCount == 0)||(chkCount == 0)) return VectorXd::Zero(4);
  VectorXd rmse = VectorXd::Zero(estimations[0].size());

  //std::cout << "The size of rmse is "<<rmse.size();
  //std::cout << "The data is "<<dataCount;
  for (int i =0; i<dataCount; i++){
    VectorXd diff = estimations[i] - ground_truth[i];
    Eigen::ArrayXd diff_sq = diff.array();
    rmse += diff_sq.square().matrix();
  }

  // Averaging and taking a square root of the result
  rmse = rmse/dataCount;
  rmse = rmse.array().sqrt();

  return rmse;
}
