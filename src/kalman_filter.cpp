#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;

  P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) { // LIDAR
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) { // RADAR

  // Transforn x_ to Polar coordinates
  double rho = sqrt(pow(x_[0],2) + pow(x_[1], 2));
  if(fabs(rho) < 0.0001){
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    rho = 0.0001;
  }
  double phi = atan2(x_[1], x_[0]);
  double rho_dot = (x_[0] * x_[2] + x_[1]*x_[3]) / rho;

  VectorXd hx = VectorXd(3);
  hx << rho, phi, rho_dot;
  
  VectorXd y = z - hx;
  // Normalize y[1] between -pi and pi
  while(y[1] > M_PI){
    y[1] -= M_PI * 2;
  }
  while(y[1] < -M_PI){
    y[1] += M_PI*2;
  }
  // Standard kalman filter updates
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
}
