#include "kalman_filter.h"

#include <iostream>
#include "tools.h"
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
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}


void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  //calculate Hj
  Tools tools; 
  MatrixXd Hj = tools.CalculateJacobian(x_);  
  cout << "Hj :" << endl << Hj << endl;
  //calculate h(x)
  VectorXd h(3);
  //recover state parameters
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  //pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = py/px;
  float c4 = px*vx + py*vy;
  float rho = c2;
  //float varphi = atan(c3);
  float varphi = atan2(py,px);
  //cout << "varphi  :" << endl << varphi << endl;
  float dotrho = c4/c2;
 	

  h << rho,varphi,dotrho;
  cout << "h  :" << endl << h << endl;
  

  //measurement covariance
  MatrixXd R;
  R = MatrixXd(3, 3);
  R<< 0.02, 0,0,
      0, 0.0005,0,
      0, 0, 0.02;	
   
  //VectorXd z_pred = Hj * x_;
  VectorXd y = z - h;
  if(y(1)>5){
    y(1) = y(1)-6.283;
  }else if(y(1)<-5){
    y(1) = y(1)+6.283;
  }
  cout << "z  :" << endl << z << endl;
  cout << "h  :" << endl << h << endl;
  cout << "y  :" << endl << y << endl;
  MatrixXd Ht = Hj.transpose();
  MatrixXd S = Hj * P_ * Ht + R;
  //MatrixXd S = Hj * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;


  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;


}
