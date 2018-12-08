#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

 H_laser_ << 1, 0, 0, 0,
            0, 1, 0, 0;

 ekf_.Q_ = MatrixXd(4,4);


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    // first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ << 1000.0, 0, 0, 0, // Initialize P with high uncertainty for all estimates.
              0, 1000.0, 0, 0,
              0, 0, 1000.0, 0,
              0, 0, 0, 1000.0;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
     double rho = measurement_pack.raw_measurements_[0];
     double phi = measurement_pack.raw_measurements_[1];
     double rho_dot = measurement_pack.raw_measurements_[2];

     double px = rho * cos(phi);
     double py = rho * sin(phi);
     double vx = 0;
     double vy = 0;
     ekf_.x_ << px, py, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
     double px = measurement_pack.raw_measurements_[0];
     double py = measurement_pack.raw_measurements_[1];
     double vx = 0;
     double vy = 0;
     ekf_.x_ << px, py, vx, vy;
  
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  
  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  // Find F 
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, dt, 0,
            0, 1, 0, dt,
            0, 0, 1, 0,
            0, 0, 0, 1;
  double noise_ax = 9;
  double noise_ay = 9;
  // Find covariance matrix
  ekf_.Q_ << pow(dt,4)/4*noise_ax, 0, pow(dt,3)/2*noise_ax, 0,
      0, pow(dt,4)/4*noise_ay, 0, pow(dt,3)/2*noise_ay,
      pow(dt,3)/2*noise_ax, 0, pow(dt,2)*noise_ax, 0,
      0, pow(dt,3)/2*noise_ay, 0, pow(dt,2)*noise_ay;
  // Predict 
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    
    ekf_.R_ = R_radar_;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates

    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
