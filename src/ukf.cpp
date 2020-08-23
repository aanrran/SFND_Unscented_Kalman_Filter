#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 5.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;
  
  /**
   * These measurement noise values below are provided by the sensor manufacturer.
   */

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
  
  this->is_initialized_ =false;
  // State dimension
  this->n_x_ = 5;
  // Augmented state dimension
  this->n_aug_ = 7;
  // initial state vector
  this->x_ = VectorXd(n_x_);
  // initial covariance matrix
  this->P_ = MatrixXd(n_x_, n_x_);
  //the variance should be at the same scale of the other standard deviation
  this->P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
  				0, std_laspx_*std_laspx_, 0, 0, 0,
  				0, 0, std_radrd_*std_radrd_, 0, 0,
  				0, 0, 0, std_radphi_*std_radphi_, 0,
  				0, 0, 0, 0, std_yawdd_*std_yawdd_;
  // Sigma point spreading parameter
  this->lambda_ = 3 - n_aug_;
  
  // Weights of sigma points
  this->weights_ = VectorXd(n_aug_*2+1);
  // set weights
  this->weights_.fill(0.5/(this->n_aug_ + this->lambda_));
  this->weights_(0) = this->lambda_/(this->lambda_ + this->n_aug_);
  
  // predicted sigma points matrix
  this->Xsig_pred_ = MatrixXd(n_x_, n_aug_*2+1);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * step 1. Check first Measurement
   * yes -> Initialize State x and Covariance P
   * No  -> go to Prediction
   */

  //if the initialization is not done, start initialzation
  if(this->is_initialized_ == false) {
    // Convert radar from polar to cartesian coordinates and initialize state
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      //initialize the state vector x
      this->x_.fill(0.0);
      auto ro = static_cast<float>(meas_package.raw_measurements_(0));     
      auto phi = static_cast<float>(meas_package.raw_measurements_(1));

      this->x_(0) = ro * cos(phi);
      this->x_(1) = ro * sin(phi);


    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      //initialize the state vector x
      this->x_.fill(0.0);
      this->x_(0) = static_cast<float>(meas_package.raw_measurements_(0));
      this->x_(1) = static_cast<float>(meas_package.raw_measurements_(1));
    }
    
    //crossup the initialization
    this->is_initialized_ = true; 
    //intialize the starting time
    this->time_us_ = meas_package.timestamp_;
    
    /**
   * step 2. Prediction
   * a. Compute elapsed time dt
   * b. Compute the augmented Mean(x_aug) and Covariance(P_aug)
   * c. Compute the x_k|k, P_k|k
   */

  } else { //if the initialization is done, continue the prediction
    //update the delta t and the current time
    double dt = (meas_package.timestamp_ - this->time_us_) / 1000000.0;
    this->time_us_ = meas_package.timestamp_;
    //start the prediction
    Prediction(dt);
    
    /**
   * step 3. Check if the received data is Radar or Lidar
   * Radar-> set the measurement vector z = new Radar reading from the MeasurementPackage
   * Lidar-> set the meansurement vector z = new Lidar reading from the MeasurementPackage
   * a. calculate the pridicted measurement mean and cov z_k|k, S from the x_k|k, P_k|k
   * b. calculate the Kalman gain K, and the difference between real and predicted meansurement
   *	z_diff = z - z_k|k
   * c. calculate the new estimate x_k+1|k, P_k+1|k through Uncented Kalman filter
   */
    
    //update the state base on the data type
    if(this->use_radar_ == true && meas_package.sensor_type_ == MeasurementPackage::RADAR) UpdateRadar(meas_package);
    if(this->use_laser_ == true && meas_package.sensor_type_ == MeasurementPackage::LASER) UpdateLidar(meas_package);
  }

}

void UKF::Prediction(double delta_t) {
  /**
   * Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // create augmented mean state
  Eigen::VectorXd x_aug = VectorXd(this->n_aug_);
  x_aug.fill(0.0);
  x_aug.head(this->n_x_) = this->x_;

  // create augmented covariance matrix
  Eigen::MatrixXd P_aug = MatrixXd(this->n_aug_ , this->n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(this->n_x_,this->n_x_) = this->P_;

  // add Q matrice to the bottom right corner of P augment
  P_aug(this->n_x_, this->n_x_) = pow(this->std_a_, 2);
  P_aug(this->n_x_ +1, this->n_x_ +1) = pow(this->std_yawdd_, 2);

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Eigen::MatrixXd Xsig_aug = MatrixXd(this->n_aug_ , 2*this->n_aug_ +1);
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< this->n_aug_; ++i) {
    Xsig_aug.col(i+1)       = x_aug + sqrt(this->lambda_ + this->n_aug_) * L.col(i);
    Xsig_aug.col(i+1+ this->n_aug_) = x_aug - sqrt(this->lambda_ + this->n_aug_) * L.col(i);
  }

  // predict sigma points
  for (int i = 0; i< 2*this->n_aug_ +1; ++i) {
    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into the related column
    this->Xsig_pred_(0,i) = px_p;
    this->Xsig_pred_(1,i) = py_p;
    this->Xsig_pred_(2,i) = v_p;
    this->Xsig_pred_(3,i) = yaw_p;
    this->Xsig_pred_(4,i) = yawd_p;
  }

  // predicted state mean
  this->x_ = this->Xsig_pred_ * this->weights_;

  // predicted state covariance matrix
  this->P_.fill(0.0);
  for (int i = 0; i < 2 * this->n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = this->Xsig_pred_.col(i) - this->x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    this->P_ = this->P_ + this->weights_(i) * x_diff * x_diff.transpose() ;
  }

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   */

  //initialize the radar data size
  int n_z = 2;
  // create a vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  //add the incoming radar reading to z vector
  double px = meas_package.raw_measurements_(0);
  double py = meas_package.raw_measurements_(1);
  z << px, py;
  
  // transform sigma points into measurement space

  // create a vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z); 
  // mean predicted measurement
  z_pred = this->x_.head(n_z);

  // innovation covariance matrix S
  // create a matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z,n_z);

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  pow(this->std_laspx_, 2), 0,
  		0, pow(this->std_laspy_, 2);
  S = this->P_.topLeftCorner(n_z,n_z) + R;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(this->n_x_ , n_z);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * this->n_aug_ + 1; ++i) {  // 2n+1 simga points

    // state difference
    VectorXd x_diff = this->Xsig_pred_.col(i) - this->x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    
    // residual
    VectorXd z_diff = x_diff.head(n_z);

    Tc = Tc + this->weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // update state mean and covariance matrix
  this->x_ = this->x_ + K * z_diff;
  this->P_ = this->P_ - K*S*K.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   */
  
  //initialize the radar data size
  int n_z = 3;
  // create a vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  //add the incoming radar reading to z vector
  double rho = meas_package.raw_measurements_(0);
  double phi = meas_package.raw_measurements_(1);
  double rho_dot = meas_package.raw_measurements_(2);
  z << rho, phi, rho_dot;
  
  // transform sigma points into measurement space
  // create a matrix with sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * this->n_aug_ + 1);
  // create a vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  for (int i = 0; i < 2 * this->n_aug_ + 1; ++i) {  // 2n+1 simga points
    // extract values for better readability
    double p_x = this->Xsig_pred_(0,i);
    double p_y = this->Xsig_pred_(1,i);
    double v  = this->Xsig_pred_(2,i);
    double yaw = this->Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
    Zsig(1,i) = atan2(p_y,p_x);                                // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   // r_dot
  }

  // mean predicted measurement
  z_pred = Zsig * this->weights_;

  // innovation covariance matrix S
  // create a matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * this->n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + this->weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  pow(this->std_radr_, 2), 0, 0,
  		0, pow(this->std_radphi_, 2), 0,
  		0, 0, pow(this->std_radrd_, 2);
  S = S + R;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(this->n_x_ , n_z);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * this->n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = this->Xsig_pred_.col(i) - this->x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + this->weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  this->x_ = this->x_ + K * z_diff;
  this->P_ = this->P_ - K*S*K.transpose();

}