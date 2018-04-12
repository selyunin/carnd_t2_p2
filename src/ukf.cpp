#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  is_initialized_ = false;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initialize state vector size
  n_x_ = 5;

  // initialize augmented state vector size
  n_aug_ = n_x_ + 2;

  // number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix.square()
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  //define spreading parameter
  lambda_ = 3.0 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug_.fill(0.0);

  P_  << 1,    0,    0,    0,    0,
		 0,    1,    0,    0,    0,
		 0,    0,    1,    0,    0,
		 0,    0,    0,    1,    0,
		 0,    0,    0,    0,    1;

  sig_weights_ = VectorXd(2*n_aug_+1);
  CalculateWeights();
  R_radar_ = MatrixXd(3,3);
  R_radar_ << pow(std_radr_, 2),  0, 0,
		      0, pow(std_radphi_,2), 0,
		      0, 0, pow(std_radrd_, 2);

  R_lidar_ = MatrixXd(2,2);
  R_lidar_ << pow(std_laspx_, 2),  0,
		  	  0,   pow(std_laspy_,2);

  NIS_radar_ = 0.0;
  NIS_lidar_ = 0.0;
  NIS_radar_file_.open("NIS_radar.txt");
  NIS_lidar_file_.open("NIS_lidar.txt");

}

UKF::~UKF() {

	NIS_radar_file_.close();
	NIS_lidar_file_.close();

}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	if (!is_initialized_) {
		Initialize(meas_package);
		is_initialized_ = true;
		std::cout<<"Initialization done!"<<std::endl;
		return;
	}
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
	float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;
	Prediction(dt);
	/*****************************************************************************
	 *  Update
	 ****************************************************************************/
	//std::cout<<"Update step"<<endl;
	  /**
	     * Use the sensor type to perform the update step.
	     * Update the state and covariance matrices.
	   */
	  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
	    // Radar updates
		   	std::cout<<"RADAR:";
		   	for(int idx=0; idx< meas_package.raw_measurements_.size(); idx++){
		   		std::cout<<meas_package.raw_measurements_[idx]<<",";
		   	}
		   	std::cout<<endl;
		  UpdateRadar(meas_package);
	  } else {
		   	std::cout<<"LASER:";
		      	for(int idx=0; idx< meas_package.raw_measurements_.size(); idx++){
		       		std::cout<<meas_package.raw_measurements_[idx]<<",";
		       	}
		   	std::cout<<endl;
		   	// Laser updates
		  UpdateLidar(meas_package);
	  }
	  // print the output
	  // cout << "x_ = " << x_ << endl;
	  // cout << "P_ = " << P_ << endl;
}


void UKF::Prediction(double delta_t) {
	/**
	Prediction step:
	1. Calculating sigma points;
	2. Passing sigma points through process model;
	3. Predicting state mean and covariance.
	*/
	AugmentedSigmaPoints();
	SigmaPointPrediction(delta_t);
	PredictMeanAndCovariance();

	std::cout<<"Prediction step done"<<std::endl;
	std::cout<<"x = "<<std::endl;
	std::cout<<x_<<std::endl;
	std::cout<<"P = "<<std::endl;
	std::cout<<P_<<std::endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& meas_package) {
  /**
  Update step for lidar:
  1. predict lidar measurements
  2. compute innovation and Kalman gain
  3. update the belief about the object's position.
  */
	int n_z = 2;
	VectorXd z = meas_package.raw_measurements_;
	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z,n_z);
	S.fill(0.0);
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	std::cout<<"PredictLidarMeasurement:"<<std::endl;
	PredictLidarMeasurement(z_pred, S, Zsig);
	std::cout<<"PredictLidarMeasurement done"<<std::endl;
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	//calculate cross correlation matrix
	for (unsigned idx=0; idx < n_sig_; ++idx){
	  //measurement difference
	  VectorXd z_diff = Zsig.col(idx) - z_pred;
	  // state difference
	  VectorXd x_diff = Xsig_pred_.col(idx) - x_;
      //angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
	  Tc = Tc + sig_weights_(idx) * x_diff * z_diff.transpose();
	}
	//calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();
	//residual
	VectorXd z_diff = z - z_pred;
	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();
	//print result
	std::cout << "Updated state x: " << std::endl << x_ << std::endl;
	std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

	NIS_lidar_ = z.transpose() * S.inverse() * z;
//    NIS_lidar_file_ << meas_package.timestamp_<<","<<NIS_lidar_;
}
///**
// * Updates the state and the state covariance matrix using a radar measurement.
// * @param {MeasurementPackage} meas_package
// */
void UKF::UpdateRadar(const MeasurementPackage& meas_package) {
  /**
	Update step for radar:
	1. calculate predicted measurement for radar
	2. compute innovation and Kalman gain
	3. update Uthe belief about the object's position
  */
	  int n_z = 3;
	  VectorXd z = meas_package.raw_measurements_;
	  //mean predicted measurement
	  VectorXd z_pred = VectorXd(n_z);
	  z_pred.fill(0.0);
	  //measurement covariance matrix S
	  MatrixXd S = MatrixXd(n_z,n_z);
	  S.fill(0.0);
	  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	  //std::cout<<"PredictRadarMeas"<<std::endl;
	  PredictRadarMeasurement(z_pred, S, Zsig);
	  std::cout<<"PredictRadarMeas done"<<std::endl;
	  //create matrix for cross correlation Tc
	  MatrixXd Tc = MatrixXd(n_x_, n_z);
	  Tc.fill(0.0);
	  //calculate cross correlation matrix
	  for (unsigned idx=0; idx < Xsig_pred_.cols(); ++idx){
	      //measurement difference
	      VectorXd z_diff = Zsig.col(idx) - z_pred;
	      //angle normalization
	      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
	      // state difference
	      VectorXd x_diff = Xsig_pred_.col(idx) - x_;
	      //angle normalization
	      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
	      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

	      Tc = Tc + sig_weights_(idx) * x_diff * z_diff.transpose();
	  }
//	  std::cout<<"Kalman Gain:"<<std::endl;
	  //calculate Kalman gain K;
	  MatrixXd K = Tc * S.inverse();
	  //residual
	  VectorXd z_delta = z - z_pred;
	  //angle normalization
	  while (z_delta(1)> M_PI) z_delta(1)-=2.*M_PI;
	  while (z_delta(1)<-M_PI) z_delta(1)+=2.*M_PI;
	  //update state mean and covariance matrix
	  x_ = x_ + K * z_delta;
	  P_ = P_ - K * S * K.transpose();
	  //print result
	  std::cout << "Updated state x: " << std::endl << x_ << std::endl;
	  std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
	  //write result
	  //*x_out = x_;
	  //*P_out = P_;
	  NIS_radar_ = z.transpose() * S.inverse() * z;
//	  NIS_radar_file_ << meas_package.timestamp_<<","<<NIS_radar_;
}

void UKF::Initialize(const MeasurementPackage& measurement_pack){
  /*****************************************************************************
  *  Initialization
  ****************************************************************************/
   // first measurement
   x_.fill(0.0);
   if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
     /**
     Convert radar from polar to cartesian coordinates and initialize state.
     */
   	std::cout<<"RADAR:";
   	for(int idx=0; idx< measurement_pack.raw_measurements_.size(); idx++){
   		std::cout<<measurement_pack.raw_measurements_[idx]<<",";
   	}
   	std::cout<<endl;
   	double rho = measurement_pack.raw_measurements_(0);
   	double phi = measurement_pack.raw_measurements_(1);
   	double rho_dot = measurement_pack.raw_measurements_(2);
   	x_ << rho * cos(phi),
   		  rho * sin(phi),
   		  sqrt(pow(rho_dot * cos(phi),2) + pow(rho_dot * sin(phi),2)),
   		  0,
		  0;
   }
   else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
     /**
     Initialize state.
     */
   	std::cout<<"LASER:";
      	for(int idx=0; idx< measurement_pack.raw_measurements_.size(); idx++){
       		std::cout<<measurement_pack.raw_measurements_[idx]<<",";
       	}
   	std::cout<<endl;
   	x_ << measurement_pack.raw_measurements_(0),
		  measurement_pack.raw_measurements_(1),
		  0,
		  0,
		  0;
   }
   previous_timestamp_ = measurement_pack.timestamp_;
   // done initializing, no need to predict or update
   return;
 }


void UKF::AugmentedSigmaPoints(){
	//create augmented state covariance
	VectorXd x_aug(n_aug_);
	x_aug.fill(0.0);
	x_aug.head(n_x_) = x_;
	//std::cout<<"X_aug:\n"<<x_aug<<std::endl;
	//create augmented covariance matrix
	MatrixXd Q(2,2);
	Q << pow(std_a_,2), 0,
		 0, pow(std_yawdd_,2);
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_,n_x_) = P_;
	P_aug.block(P_.rows(), P_.cols(), Q.rows(), Q.cols()) = Q;
	//std::cout<<"P_aug:\n"<<P_aug<<std::endl;
	//create square root matrix
	MatrixXd A = P_aug.llt().matrixL();
	//create augmented sigma points
	//set first column of sigma point matrix
	Xsig_aug_.fill(0.0);
	Xsig_aug_.col(0)  = x_aug;
	//set remaining sigma points
	for (int i = 0; i < n_aug_; i++)
	{
	  Xsig_aug_.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * A.col(i);
	  Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A.col(i);
	}
}

void UKF::SigmaPointPrediction(double delta_t){
	VectorXd x_f(n_x_);
	VectorXd nu_f(n_x_);
	for(unsigned i=0; i < Xsig_aug_.cols(); ++i){
		double p_x = Xsig_aug_(0,i);
		double p_y = Xsig_aug_(1,i);
		double v = Xsig_aug_(2,i);
		double phi = Xsig_aug_(3,i);
		double phi_dot = Xsig_aug_(4,i);
		double nu_a = Xsig_aug_(5,i);
		double nu_phi_ddot = Xsig_aug_(6,i);
		//avoid division by zero
		if (fabs(phi_dot) < 1e-4){
		  x_f << v * cos(phi) * delta_t,
				 v * sin(phi) * delta_t,
				 0,
				 0,
				 0;

		}else{
		  x_f << v / phi_dot * (sin(phi + phi_dot * delta_t) - sin(phi)),
				 v / phi_dot * (-cos(phi + phi_dot * delta_t) + cos(phi)),
				 0,
				 phi_dot * delta_t,
				 0;
		}
		nu_f<<  0.5 * delta_t * delta_t * cos(phi) * nu_a,
			    0.5 * delta_t * delta_t * sin(phi) * nu_a,
			    delta_t * nu_a,
			    0.5 * delta_t * delta_t * nu_phi_ddot,
			    delta_t * nu_phi_ddot;
		//predict sigma points
		VectorXd x_k(n_x_);
		x_k << p_x,
			   p_y,
			   v,
			   phi,
			   phi_dot;
		VectorXd x_pred(n_x_);
		x_pred = x_k + x_f + nu_f;
		//write predicted sigma points into right column
		Xsig_pred_.col(i) = x_pred;
	}
	//print result
	//std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;
}

void UKF::PredictMeanAndCovariance(){
	  //predict state mean
	x_ = Xsig_pred_ * sig_weights_;
//	for (unsigned i=0; i<Xsig_pred_.cols(); ++i){
//		x_ = x_ + sig_weights_(i) * Xsig_pred_.col(i);
//		while (x_(3)> M_PI) x_(3)-=2.*M_PI;
//		while (x_(3)<-M_PI) x_(3)+=2.*M_PI;
//	}
	  //predict state covariance matrix
	  P_.fill(0.0);
	  for (unsigned i=0; i<Xsig_pred_.cols(); ++i){
		  VectorXd x_diff = Xsig_pred_.col(i) - x_;
		  //angle normalization
		  while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		  while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
		  P_ = P_ + sig_weights_(i) * x_diff * x_diff.transpose();
	  }
	  //print result
	  std::cout << "Predicted state" << std::endl;
	  std::cout << x_ << std::endl;
	  std::cout << "Predicted covariance matrix" << std::endl;
	  std::cout << P_ << std::endl;
	  //write result
}

void UKF::CalculateWeights(){
	sig_weights_(0) = lambda_ / (lambda_ + n_aug_);
	for (int i=1; i<2*n_aug_+1; i++) {
		sig_weights_(i) = 0.5/(n_aug_ + lambda_);
	}
}

void UKF::PredictLidarMeasurement(VectorXd& z_pred, MatrixXd& S, MatrixXd& Zsig){
	  //set measurement dimension, radar can measure r, phi, and r_dot
	  int n_z = 2;
	  //create matrix for sigma points in measurement space
//	  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	  //mean predicted measurement
//	  VectorXd z_pred = VectorXd(n_z);

	  //measurement covariance matrix S
//	  MatrixXd S = MatrixXd(n_z,n_z);

	  //transform sigma points into measurement space
	  Zsig = Xsig_pred_.block(0,0,n_z, 2*n_aug_ + 1);
//	  for (unsigned idx=0; idx < Xsig_pred_.cols(); ++idx){
//		  VectorXd x_sig = Xsig_pred_.col(idx);
//		  double px = x_sig(0);
//		  double py = x_sig(1);
//		  double v = x_sig(2);
//		  double phi = x_sig(3);
//		  double phi_dot = x_sig(4);
//
//		  VectorXd h(n_z);
//		  h(0) = px;
//		  h(1) = py;
//		  Zsig.col(idx) = h;
//	    }
	    //calculate mean predicted measurement
	  	z_pred = Zsig * sig_weights_;
//		for (unsigned idx=0; idx < Xsig_pred_.cols(); ++idx){
		    // state difference
//		    z_pred = z_pred + sig_weights_(idx) * Zsig.col(idx) ;
//		}
	    //calculate innovation covariance matrix S
		for (unsigned idx=0; idx < n_sig_; ++idx){
		    VectorXd z_diff = Zsig.col(idx) - z_pred;

			S = S + sig_weights_(idx) * z_diff * z_diff.transpose() ;
		}
		MatrixXd R(n_z, n_z);
//		R << pow(std_laspx_, 2),  0,
//			 0,   pow(std_laspy_,2);
		S = S + R_lidar_;
		//print result
		std::cout << "z_pred: " << std::endl << z_pred << std::endl;
		std::cout << "S: " << std::endl << S << std::endl;
}



void UKF::PredictRadarMeasurement(VectorXd& z_pred, MatrixXd& S, MatrixXd& Zsig){
	  //set measurement dimension, radar can measure r, phi, and r_dot
	  int n_z = 3;

	  //transform sigma points into measurement space
	  for (unsigned idx=0; idx < Xsig_pred_.cols(); ++idx){
		  double px = Xsig_pred_(0,idx);
		  double py = Xsig_pred_(1, idx);
		  double v = Xsig_pred_(2, idx);
		  double phi = Xsig_pred_(3, idx);
		  double phi_dot = Xsig_pred_(4, idx);

//		  VectorXd h(3);
//		  h(0) = sqrt(pow(px,2)+pow(py,2));
//		    if( fabs(h(0)) < 1.0e-4) // avoid infinitesimal
//			  h(0) = 1.0e-4;
//		    if( fabs(px) < 1.0e-4) // avoid infinitesimal
//			  px = 1.0e-4;
//		  h(1) = atan2(py, px);
//		  h(2) = ( px * v * cos(phi) + py * v * sin(phi) ) / h(0) ;
		  Zsig(0, idx) = sqrt(pow(px,2)+pow(py,2));
		  Zsig(1, idx) = atan2(py, px);
		  Zsig(2, idx) = ( px * v * cos(phi) + py * v * sin(phi) ) / Zsig(0, idx) ;
	    }
	    //calculate mean predicted measurement
//		for (unsigned idx=0; idx < Xsig_pred_.cols(); ++idx){
		    // state difference
//		    z_pred = z_pred + sig_weights_(idx) * Zsig.col(idx) ;
//			while (z_pred(1)> M_PI) z_pred(1)-=2.*M_PI;
//			while (z_pred(1)<-M_PI) z_pred(1)+=2.*M_PI;
//		}
		z_pred = Zsig * sig_weights_;
		S.fill(0.0);
	    //calculate innovation covariance matrix S
		for (unsigned idx=0; idx < Xsig_pred_.cols(); ++idx){
		    //angle normalization
		    VectorXd z_diff = Zsig.col(idx) - z_pred;
		    //angle normalization
		    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
			S = S + sig_weights_(idx) * z_diff * z_diff.transpose() ;
		}
		S = S + R_radar_;
		//print result
		std::cout << "z_pred: " << std::endl << z_pred << std::endl;
		std::cout << "S: " << std::endl << S << std::endl;
}
