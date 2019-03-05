#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // TODO: YOUR CODE HERE
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0)
  {
      std::cout << "Invalid estimation or ground_truth data" << std::endl;
      return rmse;
  }    
 
  // TODO: accumulate squared residuals
  for (int i=0; i < estimations.size(); ++i)
  {
    // ... your code here
    VectorXd residual = (estimations[i] - ground_truth[i]);
    residual = residual.array() * residual.array();
    rmse += residual;

  }

  // TODO: calculate the mean
  rmse = rmse / estimations.size();

  // TODO: calculate the squared root
  rmse = rmse.array().sqrt();
  cout << rmse;
  
  if( rmse(0) > .11 ||
      rmse(1) > .11 ||
      rmse(2) > .52 ||
      rmse(3) > .52 )
    cout << "Warning  rmse = " 
         << rmse(0) << "  " << rmse(1) << "  " 
         << rmse(2) << "  " << rmse(3) << endl
         << " currently exceeds tolerances of "
         << ".11, .11, .52, .52" << endl;
  
  // return the result
  return rmse;
  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) 
{
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // TODO: YOUR CODE HERE 
  float px2py2 = px*px + py *py;
  float sqrt_px2py2 = sqrt(px2py2);
  float vp = vx*py - vy*px;
  float px2py2pow3_2 = pow(px2py2, 3/2);
  // check division by zero
  if (px2py2 == 0)
    cout<< "Error!";

  else
    // compute the Jacobian matrix
    Hj << px/sqrt_px2py2, py/sqrt_px2py2, 0, 0,
        -py/px2py2, px/px2py2, 0, 0, 
        py * (vp)/(px2py2pow3_2), px * (vp)/(px2py2pow3_2), px/sqrt_px2py2, py/sqrt_px2py2;

  return Hj;
}
