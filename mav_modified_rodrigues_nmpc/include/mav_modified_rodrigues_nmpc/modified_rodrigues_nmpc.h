/*
 Copyright (c) 2015, Mina Kamel, ASL, ETH Zurich, Switzerland

 You can contact the author at <mina.kamel@mavt.ethz.ch>

 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of ETHZ-ASL nor the
 names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL ETHZ-ASL BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef INCLUDE_MAV_MODIFIED_RODRIGUES_NMPC_MODIFIED_RODRIGUES_NMPC_H_
#define INCLUDE_MAV_MODIFIED_RODRIGUES_NMPC_MODIFIED_RODRIGUES_NMPC_H_


#include <ros/ros.h>
#include <Eigen/Eigen>
#include <mav_msgs/conversions.h>
#include <mav_msgs/eigen_mav_msgs.h>
#include <stdio.h>
#include <mav_control_interface/mpc_queue.h>
#include "acado_common.h"
#include "acado_auxiliary_functions.h"
#include <mav_disturbance_observer/KF_disturbance_observer.h>
#include <std_srvs/Empty.h>
#include <lapacke.h>
#include <iostream>
#include <fstream>
#include <string>

ACADOvariables acadoVariables;
ACADOworkspace acadoWorkspace;

namespace mav_control {

lapack_logical select_lhp(const double *real, const double *imag)
{
  return *real < 0.0;
}

class NonlinearModelPredictiveControl
{
 public:
  NonlinearModelPredictiveControl(const ros::NodeHandle& nh, const ros::NodeHandle& private_nh);
  ~NonlinearModelPredictiveControl();


  // Dynamic parameters
  void setPositionPenality(const Eigen::Vector3d& q_position)
  {
    q_position_ = q_position;
  }
  void setVelocityPenality(const Eigen::Vector3d& q_velocity)
  {
    q_velocity_ = q_velocity;
  }
  void setAttitudePenality(const Eigen::Vector3d& q_attitude)
  {
    q_attitude_ = q_attitude;
  }
  void setAngVelPenality(const Eigen::Vector3d& q_angVel)
  {
	  q_angVel_ = q_angVel;
  }

  void setAngVelRef(double angVel_ref)
   {
 	  angVel_ref_ = angVel_ref;
   }

  void setCommandPenality(const Eigen::VectorXd& r_command)
  {
    r_command_ = r_command;
  }
 
  void setAltitudeIntratorGain(double Ki_altitude)
  {
    Ki_altitude_ = Ki_altitude;
  }

  void setXYIntratorGain(double Ki_xy)
  {
    Ki_xy_ = Ki_xy;
  }

  void setEnableOffsetFree(bool enable_offset_free)
  {
    enable_offset_free_ = enable_offset_free;
  }

  void setEnableShadowMRP(bool enable_shadow_mrp)
    {
      enable_shadow_mrp_ = enable_shadow_mrp;
    }

  void setEnableIntegrator(bool enable_integrator)
  {
    enable_integrator_ = enable_integrator;
  }

  void setControlLimits(const Eigen::VectorXd& control_limits) // constraints on f1,f2,f3,f4...
  {
	  force_min_ = control_limits(0);
	  force_max_ = control_limits(1);
  }

  void applyParameters();

  double getMass() const
  {
    return mass_;
  }

  double getForceConstant() const
  {
    return force_constant_;
  }

  // get reference and predicted state
  bool getCurrentReference(mav_msgs::EigenTrajectoryPoint* reference) const;
  bool getCurrentReference(mav_msgs::EigenTrajectoryPointDeque* reference) const;
  bool getPredictedState(mav_msgs::EigenTrajectoryPointDeque* predicted_state) const;

  // set odometry and commands
  void setOdometry(const mav_msgs::EigenOdometry& odometry);
  void setCommandTrajectoryPoint(const mav_msgs::EigenTrajectoryPoint& command_trajectory);
  void setCommandTrajectory(const mav_msgs::EigenTrajectoryPointDeque& command_trajectory);

  // compute control input
  void calculateForcesCommand(Eigen::VectorXd* ref_normforces, Eigen::VectorXd* ref_forces);
  double normalizeForce(double force);


  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
 private:


  // constants
  static constexpr double kGravity = 9.8066;
  static constexpr int kDisturbanceSize = 3;

  // ros node handles
  ros::NodeHandle nh_, private_nh_;

  // reset integrator service 
  ros::ServiceServer reset_integrator_service_server_;
  bool resetIntegratorServiceCallback(std_srvs::Empty::Request &req,
                                      std_srvs::Empty::Response &res);
  
  // sampling time parameters
  void initializeParameters();
  bool initialized_parameters_;

  void removeCSVfiles();

  // sampling time parameters
  double sampling_time_;
  double prediction_sampling_time_;

  // system model parameters
  double mass_;
  Eigen::Vector3d inertia_;
  double drag_coefficient_;
  double armlength_;
  double force_constant_;

  // controller parameters
  // state penalty
  Eigen::Vector3d q_position_;
  Eigen::Vector3d q_velocity_;
  Eigen::Vector3d q_attitude_;
  Eigen::Vector3d q_angVel_;

  // control penalty
  Eigen::VectorXd r_command_;

  // error integrator
  bool enable_integrator_;
  double Ki_altitude_;
  double Ki_xy_;
  double antiwindup_ball_;
  Eigen::Vector3d position_error_integration_;
  double position_error_integration_limit_;
 // double logistic_gain_;
  bool enable_shadow_mrp_;


  double force_min_;
  double force_max_;

  // reference queue
  MPCQueue mpc_queue_; 
  Vector3dDeque position_ref_, velocity_ref_, acceleration_ref_;
  std::deque<double> yaw_ref_, yaw_rate_ref_;
  double angVel_ref_;
  Eigen::Vector3d angVel_ref_vector_;


  // solver matrices
  Eigen::Matrix<double, ACADO_NY, ACADO_NY> W_;
  Eigen::Matrix<double, ACADO_NYN, ACADO_NYN> WN_;
  Eigen::Matrix<double, ACADO_N + 1, ACADO_NX> state_;
  Eigen::Matrix<double, ACADO_N, ACADO_NU> input_;
  Eigen::Matrix<double, ACADO_N, ACADO_NY> reference_;
  Eigen::Matrix<double, 1, ACADO_NYN> referenceN_;
  Eigen::Matrix<double, ACADO_N + 1, ACADO_NOD> acado_online_data_;

  // disturbance observer
  bool enable_offset_free_;

  // commands
  Eigen::VectorXd command_forces_;
  double u_ref_;

  // debug info
  bool verbose_;
  double solve_time_average_;
  ros::WallTime loop_start_time_;
  bool record_to_csv_;
  bool use_error_dynamics_;

  // most recent odometry information
  mav_msgs::EigenOdometry odometry_;
  Eigen::Vector3d position_error_ ;
  Eigen::Vector3d velocity_error_ ;
  Eigen::Quaterniond attitude_error_;
  Eigen::Vector3d mrp_error_;
  double  mrp_error_magnitude_;
  Eigen::Vector3d angVel_error_;
  bool received_first_odometry_;

  // initilize solver
  void initializeAcadoSolver(Eigen::VectorXd x0);

};

}

#endif /* INCLUDE_MAV_MODIFIED_RODRIGUES_NMPC_MODIFIED_RODRIGUES_NMPC_H_ */
