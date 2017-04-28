#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  is_initialized_ = false;
  time_us_ = 0;   // initial time for calculating change in time

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.0);
  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ = MatrixXd::Identity(5,5);     // setting covariance as Identity
  //P_.bottomRightCorner(3,3) << 100,0,0,
  //                              0,100,0,
  //                              0,0,100;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 7;
  double kappa = 0;
  double alpha = 0.001;
  double beta = 2;
  //lambda_ = 3 - n_aug_;
  lambda_ = alpha*alpha*(n_aug_+kappa) - n_aug_;

  // Initializing measurement covariance matrices
  R_radar = MatrixXd(3,3);
  R_radar << std_radr_*std_radr_, 0,0,
              0, std_radphi_*std_radphi_, 0,
              0, 0, std_radrd_*std_radrd_;

  R_lidar = MatrixXd(2,2);
  R_lidar << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

  // Initializing weights for sigma points
  weights_m_ = VectorXd(2*n_aug_ + 1);
  weights_c_ = VectorXd(2*n_aug_ + 1);
  double w_mc = 0.5/(lambda_+n_aug_);
  weights_m_.fill(w_mc);
  weights_c_.fill(w_mc);
  weights_m_(0) = lambda_/(lambda_+n_aug_);
  weights_c_(0) = weights_m_(0) + (1 - alpha*alpha + beta);

  //Initializing NIS
  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;
  // Initializing predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  Xsig_pred_.fill(0.0);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  // initialization
  if (!is_initialized_){
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_dot = meas_package.raw_measurements_(2);
      double px = rho*cos(phi);
      double py = rho*sin(phi);
      double vx = rho_dot*cos(phi);
      double vy = rho_dot*sin(phi);
      double v = sqrt(vx*vx + vy*vy);
      x_ << px, py, v, phi, 0.0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);
      double phi = atan2(py,px);
      x_ << px, py, 0.0, phi, 0.0;
    }
    // adjusting initial values close to zero:
    float eps = 0.001;
    if ((fabs(x_(0)) < eps) || (fabs(x_(1)) < eps)) {
      x_(0) = eps;
      x_(1) = eps;
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  //Prediction
  // Breaking the time step into an array of critical time steps for integration
  //double crit_t = 0.005;
  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
  std::cout<<"delta_t = "<<delta_t<<std::endl;
  //std::cout<<"Current timestamp = "<<meas_package.timestamp_<<std::endl;
  //std::cout<<"Previous timestamp = "<<time_us_<<std::endl;
  Prediction(delta_t);

  // Measurement Update
  if (use_laser_){
    if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      UpdateLidar(meas_package);
      time_us_ = meas_package.timestamp_;
    }
  }
  if (use_radar_){
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      UpdateRadar(meas_package);
      time_us_ = meas_package.timestamp_;
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // Generate Sigma Points
  MatrixXd Xsig_ = GenerateSigmaPoints();

  // Predict Sigma Points
  PredictSigmaPoints(Xsig_, delta_t);

  // Predict the next state through weighted mean
  x_ = Xsig_pred_*weights_m_;
  //Predict the covariance matrix through weighted mean
  MatrixXd mean = x_.replicate(1,2*n_aug_+1);
  MatrixXd diff = Xsig_pred_ - mean;
  /*Normalizing the angles*/
  for (int i=0; i< 2*n_aug_+1; i++){
    diff(3,i) = atan2(sin(Xsig_pred_(3,i)),cos(Xsig_pred_(3,i))) - atan2(sin(mean(3,i)),cos(mean(3,i)));
  }
  MatrixXd wt_rep = weights_c_.transpose().replicate(n_x_,1);
  MatrixXd wt_diff = wt_rep.array()*diff.array();
  P_ = wt_diff*diff.transpose();

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  MatrixXd Zsig = Xsig_pred_.topRows(2);

  //Predicted Mean
  VectorXd Zpred = Zsig*weights_m_;
  //Measurement Covariance Matrix S
  MatrixXd S = MatrixXd(2,2);
  S.fill(0.0);
  MatrixXd mean_z = Zpred.replicate(1,2*n_aug_+1);
  MatrixXd diff = Zsig - mean_z;
  MatrixXd wt_rep = weights_c_.transpose().replicate(2,1);
  MatrixXd wt_diff_z = wt_rep.array()*diff.array();
  S = wt_diff_z*diff.transpose();
  S += R_lidar;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, 2);
  MatrixXd mean_x = x_.replicate(1, 2*n_aug_+1);
  MatrixXd diff_x = Xsig_pred_ - mean_x;
  /*Normalizing the angles*/
  for (int i=0; i<2*n_aug_+1;i++){
    diff_x(3,i) = atan2(sin(Xsig_pred_(3,i)),cos(Xsig_pred_(3,i))) - atan2(sin(mean_x(3,i)),cos(mean_x(3,i)));
  }
  Tc = diff_x*wt_diff_z.transpose();

  // Kalman gain Matrix
  MatrixXd Si = S.inverse();
  MatrixXd K = Tc*Si;
  // Updating state and covariance matrix
  VectorXd z = meas_package.raw_measurements_;
  VectorXd meas_err = z - Zpred;
  x_ += K*meas_err;
  P_ -= K*S*K.transpose();

  // Calculating NIS for LIDAR
  NIS_laser_ = meas_err.transpose()*Si*meas_err;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //n_z = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_+1);
  for (int i = 0; i<2*n_aug_+1; i++){
    VectorXd sig_point = Xsig_pred_.col(i);
    double rho = sqrt(sig_point(0)*sig_point(0) + sig_point(1)*sig_point(1));
    double theta = atan2(sig_point(1),sig_point(0));
    double theta_dot = (sig_point(0)*cos(sig_point(3)) + sig_point(1)*sin(sig_point(3)))*sig_point(2)/rho;
    Zsig.col(i) << rho, theta, theta_dot;
  }

  // Calculate predicted mean measurement
  VectorXd Zpred = Zsig*weights_m_;
  // Calculate measurement covariance matrix S
  MatrixXd S = MatrixXd(3,3);
  S.fill(0.0);
  MatrixXd mean_z = Zpred.replicate(1,2*n_aug_+1);
  MatrixXd diff = Zsig - mean_z;
  //double d_ang =0.0;
  for (int i=0; i<2*n_aug_+1; i++){
    diff(1,i) = atan2(sin(Zsig(1,i)),cos(Zsig(1,i))) - atan2(sin(mean_z(1,i)),cos(mean_z(1,i)));
  }
  MatrixXd wt_rep = weights_c_.transpose().replicate(3,1);
  MatrixXd wt_diff_z = wt_rep.array()*diff.array();
  S = wt_diff_z*diff.transpose();
  S += R_radar;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, 3);
  Tc.fill(0.0);
  MatrixXd mean_x = x_.replicate(1, 2*n_aug_+1);
  MatrixXd diff_x = Xsig_pred_ - mean_x;
  for (int i=0; i<2*n_aug_+1; i++){
    diff_x(3,i) = atan2(sin(Xsig_pred_(3,i)),cos(Xsig_pred_(3,i))) - atan2(sin(mean_x(3,i)),cos(mean_x(3,i)));
  }
  Tc = diff_x*wt_diff_z.transpose();

  // Kalman gain Matrix
  MatrixXd Si = S.inverse();
  MatrixXd K = Tc*Si;
  // Updating state and covariance matrix
  VectorXd z = meas_package.raw_measurements_;
  VectorXd meas_err = z - Zpred;
  x_ += K*meas_err;
  P_ -= K*S*K.transpose();

  // Calculating NIS for RADAR
  NIS_radar_ = meas_err.transpose()*Si*meas_err;
}

MatrixXd UKF::GenerateSigmaPoints(){
  // Augmented state vector
  VectorXd X_aug = VectorXd(n_aug_);
  X_aug.fill(0.0);
  X_aug.head(n_x_) = x_;

  //Creating Augmented Covariance matrix P_aug
  MatrixXd P_aug = MatrixXd::Zero(n_aug_,n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  // Creating Sigma points
  MatrixXd A = P_aug.llt().matrixL();
  MatrixXd Xsig_ = MatrixXd(n_aug_, 2*n_aug_+1);
  Xsig_.fill(0.0);
  Xsig_.col(0)= X_aug;
  double coeff = sqrt(lambda_+n_aug_);
  for (int i=0; i< n_aug_; i++){
    Xsig_.col(i+1) = X_aug + coeff*A.col(i);
    Xsig_.col(i+1+n_aug_) = X_aug - coeff*A.col(i);
  }
  return(Xsig_);
}

void UKF::PredictSigmaPoints(MatrixXd Xsig_, double delta_t){
  // Evaluating sigma points through the non-linear process model
  VectorXd delta_x_ = VectorXd(n_x_);
  delta_x_.fill(0.0);
  VectorXd err = VectorXd(n_x_);
  err.fill(0.0);
  double delta_t_2 = delta_t*delta_t;
  double tol = 0.001;
  for(int i=0; i<2*n_aug_+1; ++i){
    VectorXd point = Xsig_.col(i);
    //double px = point(0);
    //double py = point(1);
    double v = point(2);
    double phi = point(3);
    double phi_dot = point(4);
    double v_err = point(5);
    double phi_err = point(6);
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);
    double phi_next = phi + phi_dot*delta_t;
    double sin_phi_next = sin(phi_next);
    double cos_phi_next = cos(phi_next);

    // Evaluate Error and delta_x terms
    err << 0.5*delta_t_2*cos_phi*v_err,
          0.5*delta_t_2*sin_phi*v_err,
          delta_t*v_err,
          0.5*delta_t_2*phi_err,
          delta_t*phi_err;

    if(fabs(phi_dot < tol)){
      delta_x_(0) = v*cos_phi*delta_t;
      delta_x_(1) = v*sin_phi*delta_t;
    }
    else{
      double v_by_phidot = v/phi_dot;
      delta_x_(0) = v_by_phidot*(sin_phi_next - sin_phi);
      delta_x_(1) = v_by_phidot*(-cos_phi_next + cos_phi);
    }
    delta_x_(3) = delta_t*phi_dot;
    // Evaluate predicted sigma points
    Xsig_pred_.col(i) = point.head(n_x_) + delta_x_ + err;
  }
}
