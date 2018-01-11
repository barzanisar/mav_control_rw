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

#include <mav_modified_rodrigues_nmpc/modified_rodrigues_nmpc.h>
#include <math.h>


#define PI 3.14159265

namespace mav_control {

constexpr double NonlinearModelPredictiveControl::kGravity;
constexpr int NonlinearModelPredictiveControl::kDisturbanceSize;

NonlinearModelPredictiveControl::NonlinearModelPredictiveControl(const ros::NodeHandle& nh,
                                                                 const ros::NodeHandle& private_nh)
    : nh_(nh),
      private_nh_(private_nh),
      initialized_parameters_(false),
      position_error_integration_(0, 0, 0),
      mpc_queue_(nh, private_nh, ACADO_N+1),
      verbose_(false),
      solve_time_average_(0),
      received_first_odometry_(false)
{


  acado_initializeSolver();

  W_.setZero();
  WN_.setZero();
  command_forces_.setConstant(1, ACADO_NU, 0);
  input_.setZero();
  state_.setZero();
  reference_.setZero();
  referenceN_.setZero();

  reset_integrator_service_server_ = nh_.advertiseService(
      "reset_integrator", &NonlinearModelPredictiveControl::resetIntegratorServiceCallback, this);

  initializeParameters();

  mpc_queue_.initializeQueue(sampling_time_, prediction_sampling_time_);

  removeCSVfiles();


}

NonlinearModelPredictiveControl::~NonlinearModelPredictiveControl()
{

}

 
bool NonlinearModelPredictiveControl::resetIntegratorServiceCallback(std_srvs::Empty::Request &req,
                                                                     std_srvs::Empty::Response &res)
{
  position_error_integration_.setZero();
  return true;
} 

void NonlinearModelPredictiveControl::removeCSVfiles()
{
	std::remove("/home/barza/bagfiles/CSVfiles/forces.csv");
	std::remove("/home/barza/bagfiles/CSVfiles/position.csv");
	std::remove("/home/barza/bagfiles/CSVfiles/velocity.csv");
	std::remove("/home/barza/bagfiles/CSVfiles/angVel.csv");
	std::remove("/home/barza/bagfiles/CSVfiles/orientation.csv");
	std::remove("/home/barza/bagfiles/CSVfiles/qnorm.csv");
	std::remove("/home/barza/bagfiles/CSVfiles/time.csv");
	std::remove("/home/barza/bagfiles/CSVfiles/fail_Status.csv");
	std::remove("/home/barza/bagfiles/CSVfiles/pred_State.csv");
	std::remove("/home/barza/bagfiles/CSVfiles/States.csv");
	std::remove("/home/barza/bagfiles/CSVfiles/cost.csv");
	std::remove("/home/barza/bagfiles/CSVfiles/avgSolveTime.csv");
	std::remove("/home/barza/bagfiles/CSVfiles/solveTime.csv");
}

//initializing onlineData parameters
void NonlinearModelPredictiveControl::initializeParameters()
{

	std::vector<double> inertia;

  //Get parameters from RosParam server
  private_nh_.param<bool>("verbose", verbose_, false);


  if (!private_nh_.getParam("mass", mass_)) {
    ROS_ERROR("mass in nonlinear MPC controller is not loaded from ros parameter "
              "server");
    abort();
  }


  if (!private_nh_.getParam("inertia", inertia)) {
	  ROS_ERROR(
		  "inertia in nonlinear MPC controller is not loaded from ros parameter server");
	  abort();
  }

  inertia_ << inertia.at(0), inertia.at(1), inertia.at(2);

  if (!private_nh_.getParam("drag_coefficient", drag_coefficient_)) {
    ROS_ERROR(
        "drag_coefficient in nonlinear MPC controller is not loaded from ros parameter server");
    abort();
  }

   if (!private_nh_.getParam("armlength", armlength_)) {
	  ROS_ERROR("armlength in nonlinear MPC controller is not loaded from ros parameter "
		  "server");
	  abort();
  }

 if (!private_nh_.getParam("antiwindup_ball", antiwindup_ball_)) {
    ROS_ERROR(
        "antiwindup_ball in nonlinear MPC controller is not loaded from ros parameter server");
    abort();
  }

  if (!private_nh_.getParam("position_error_integration_limit",
                            position_error_integration_limit_)) {
    ROS_ERROR(
        "position_error_integration_limit in nonlinear MPC is not loaded from ros parameter server");
    abort();
  }


  if (!private_nh_.getParam("sampling_time", sampling_time_)) {
    ROS_ERROR("sampling_time in nonlinear MPC is not loaded from ros parameter server");
    abort();
  }

  if (!private_nh_.getParam("prediction_sampling_time", prediction_sampling_time_)) {
    ROS_ERROR("prediction_sampling_time in nonlinear MPC is not loaded from ros parameter server");
    abort();
  }

  if (!private_nh_.getParam("record_to_csv", record_to_csv_)) {
      ROS_ERROR("record_to_csv in nonlinear MPC is not loaded from ros parameter server");
      abort();
    }

  if (!private_nh_.getParam("use_error_dynamics", use_error_dynamics_)) {
      ROS_ERROR("use_error_dynamics in nonlinear MPC is not loaded from ros parameter server");
      abort();
    }

  if (!private_nh_.getParam("force_constant", force_constant_)) {
    ROS_ERROR("force_constant in nonlinear MPC controller is not loaded from ros parameter "
              "server");
    abort();
  }

  for (int i = 0; i < ACADO_N + 1; i++) {
    acado_online_data_.block(i, 0, 1, ACADO_NOD) << mass_, inertia_(0), inertia_(1), inertia_(2), drag_coefficient_, armlength_, 0 ,0, 0, 0, 0; //zero for angVel_ref_, yaw_ref_ and 3 external forces
  }

  Eigen::Map<Eigen::Matrix<double, ACADO_NOD, ACADO_N + 1>>(const_cast<double*>(acadoVariables.od)) =
      acado_online_data_.transpose();

  if (verbose_) {
    std::cout << "acado online data: " << std::endl << acado_online_data_ << std::endl;
  }

  u_ref_ = mass_*kGravity / ACADO_NU;

  initialized_parameters_ = true;
  ROS_INFO("Nonlinear MPC: initialized correctly");
}

//applying dynamic config parameters
void NonlinearModelPredictiveControl::applyParameters()
{

  W_.block(0, 0, 3, 3) = q_position_.asDiagonal();
  W_.block(3, 3, 3, 3) = q_velocity_.asDiagonal();
  W_.block(6, 6, 3, 3) = q_attitude_.asDiagonal();
  W_.block(9, 9, 3, 3) = q_angVel_.asDiagonal();
  W_.block(12, 12, ACADO_NU, ACADO_NU) = r_command_.asDiagonal();

  WN_.block(0, 0, 3, 3) = q_position_.asDiagonal();
  WN_.block(3, 3, 3, 3) = q_velocity_.asDiagonal();

  Eigen::Map<Eigen::Matrix<double, ACADO_NY, ACADO_NY>>(const_cast<double*>(acadoVariables.W)) = W_.transpose();
  Eigen::Map<Eigen::Matrix<double, ACADO_NYN, ACADO_NYN>>(const_cast<double*>(acadoVariables.WN)) = WN_.transpose();

  for (size_t i = 0; i < ACADO_N*ACADO_NU; i++) {

	  acadoVariables.lbValues[i] = force_min_;
	  acadoVariables.ubValues[i] = force_max_;
  }

  for (size_t i = 0; i < ACADO_N+1; i++) {

	   acado_online_data_.block(i, ACADO_NOD - 5, 1, 1) << angVel_ref_;
  }

}

void NonlinearModelPredictiveControl::setOdometry(const mav_msgs::EigenOdometry& odometry)
{
  mpc_queue_.updateQueue();
  mpc_queue_.getQueue(position_ref_, velocity_ref_, acceleration_ref_, yaw_ref_, yaw_rate_ref_);

  //std::cout << "yaw_ref: \n" << yaw_ref_.front() << std::endl;
  //std::cout << "yaw_rate_ref: \n" << yaw_rate_ref_.front() << std::endl;

  for (int i = 0; i < ACADO_N + 1; i++) {
	 // acado_online_data_.block(i, ACADO_NOD-4, 1, 1) << yaw_ref_.front(); //change
	  acado_online_data_.block(i, ACADO_NOD-4, 1, 1) << yaw_ref_[i]; //change
  }

  static mav_msgs::EigenOdometry previous_odometry = odometry;

  Eigen::Quaterniond quat_desired(cos(yaw_ref_.front() / 2), 0, 0, sin(yaw_ref_.front() / 2));

  if (odometry.position_W.allFinite() == false) {
    odometry_.position_W = previous_odometry.position_W;
    ROS_WARN("Odometry.position has a non finite element");
  } else {
    odometry_.position_W = odometry.position_W;
    previous_odometry.position_W = odometry.position_W; //
  }

  if (odometry.velocity_B.allFinite() == false) {
    odometry_.velocity_B = previous_odometry.velocity_B;
    ROS_WARN("Odometry.velocity has a non finite element");
  } else {
    odometry_.velocity_B = odometry.velocity_B;
    previous_odometry.velocity_B = odometry.velocity_B;
  }

  if (odometry.angular_velocity_B.allFinite() == false) {
    odometry_.angular_velocity_B = previous_odometry.angular_velocity_B;
    ROS_WARN("Odometry.angular_velocity has a non finite element");
  } else {
    odometry_.angular_velocity_B = odometry.angular_velocity_B;
    previous_odometry.angular_velocity_B = odometry.angular_velocity_B;
  }

  odometry_.orientation_W_B = odometry.orientation_W_B;
  previous_odometry.orientation_W_B = odometry.orientation_W_B;

  position_error_ = odometry_.position_W - position_ref_.front();
  velocity_error_ = odometry_.getVelocityWorld() - velocity_ref_.front();

  if (use_error_dynamics_){
	  attitude_error_ = quat_desired.conjugate()*odometry_.orientation_W_B; //quaternion multiplication //change
  }
  else {
	  attitude_error_ = odometry_.orientation_W_B; //change
  }


  mrp_error_= attitude_error_.vec() / (1 + attitude_error_.w());

  mrp_error_magnitude_=  mrp_error_.transpose() *  mrp_error_;

 if  ( enable_shadow_mrp_ && mrp_error_magnitude_ > 1 )
 {
	 mrp_error_ = -(1/mrp_error_magnitude_) * mrp_error_;
	 ROS_INFO_STREAM("shadow mrp used");
 }

 // angVel_ref_vector_<< angVel_ref_, 0, yaw_rate_ref_; //change
    angVel_ref_vector_<< angVel_ref_, 0, 0; //change

  angVel_error_ = odometry_.angular_velocity_B - angVel_ref_vector_;

 /*if (verbose_) {

	ROS_INFO_STREAM("mrp error: " << mrp_error_.transpose());
	ROS_INFO_STREAM("orienW: " << odometry_.orientation_W_B.w() << " orien vec: " << odometry_.orientation_W_B.vec());
 }*/

  if (!received_first_odometry_) {

	  Eigen::VectorXd x0(ACADO_NX);

	  if (use_error_dynamics_){
		  x0 << position_error_, velocity_error_, mrp_error_, angVel_error_; // 0,0,0, 0,0,0, 0,0,0, 0,0,0; //change
	    }
	    else {
	  	  x0 << odometry_.position_W, odometry_.getVelocityWorld(), mrp_error_, odometry_.angular_velocity_B; //change
	    }





	    for (int i = 0; i < ACADO_N; i++) {
	    	input_.block(i, 0, 1, ACADO_NU) << Eigen::MatrixXd::Constant(1,ACADO_NU,u_ref_);
	    }

	  initializeAcadoSolver(x0);

	  received_first_odometry_ = true;
  }
}

void NonlinearModelPredictiveControl::setCommandTrajectoryPoint(
    const mav_msgs::EigenTrajectoryPoint& command_trajectory)
{
  mpc_queue_.insertReference(command_trajectory);
}

void NonlinearModelPredictiveControl::setCommandTrajectory(
    const mav_msgs::EigenTrajectoryPointDeque& command_trajectory)
{
  int array_size = command_trajectory.size();
  if (array_size < 1)
    return;

  mpc_queue_.insertReferenceTrajectory(command_trajectory);
}

void NonlinearModelPredictiveControl::initializeAcadoSolver(Eigen::VectorXd x0)
{

  for (int i = 0; i < ACADO_N+1; i++) {
    state_.block(i, 0, 1, ACADO_NX) << x0.transpose();
  }

  Eigen::Map<Eigen::Matrix<double, ACADO_NX, ACADO_N + 1>>(const_cast<double*>(acadoVariables.x)) =
      state_.transpose();
  Eigen::Map<Eigen::Matrix<double, ACADO_NU, ACADO_N>>(const_cast<double*>(acadoVariables.u)) =
      input_.transpose();
  Eigen::Map<Eigen::Matrix<double, ACADO_NY, ACADO_N>>(const_cast<double*>(acadoVariables.y)) =
      reference_.transpose();
  Eigen::Map<Eigen::Matrix<double, ACADO_NYN, 1>>(const_cast<double*>(acadoVariables.yN)) =
      referenceN_.transpose();
}

double NonlinearModelPredictiveControl::normalizeForce(double force) {
	return (2 * force - force_max_ - force_min_) / (force_max_ - force_min_);
}

void NonlinearModelPredictiveControl::calculateForcesCommand(
    Eigen::VectorXd* ref_normforces, Eigen::VectorXd* ref_forces)
{
//	ROS_INFO_STREAM("NMPC: calc Forces command");
 // assert(ref_forces != nullptr);
 // assert(ref_normforces != nullptr);
  assert(initialized_parameters_ == true);
  ros::WallTime starting_time = ros::WallTime::now();


  Eigen::Vector3d estimated_disturbances;
  Eigen::Matrix<double, ACADO_NX, 1> x_0;

  estimated_disturbances.setZero(kDisturbanceSize);


  if (enable_integrator_) {
    Eigen::Vector3d position_error = position_ref_.front() - odometry_.position_W;
    if (position_error.norm() < antiwindup_ball_) {
      position_error_integration_ += position_error * sampling_time_;
    } else {
      position_error_integration_.setZero();
    }

    position_error_integration_ = position_error_integration_.cwiseMax(
        Eigen::Vector3d(-position_error_integration_limit_, -position_error_integration_limit_,
                        -position_error_integration_limit_));

    position_error_integration_ = position_error_integration_.cwiseMin(
        Eigen::Vector3d(position_error_integration_limit_, position_error_integration_limit_,
                        position_error_integration_limit_));

    estimated_disturbances -= Eigen::Vector3d(Ki_xy_, Ki_xy_, Ki_altitude_).asDiagonal()
        * position_error_integration_;
  }


  for (size_t i = 0; i < ACADO_N; i++) {

	  if(use_error_dynamics_){
	 // reference_.block(i, 0, 1, ACADO_NY) << 0, 0, 0,   0, 0, 0,   0, 0, 0,   0, 0, 0,   Eigen::MatrixXd::Constant(1,ACADO_NU,u_ref_);//u_ref_, u_ref_, u_ref_, u_ref_;
	  reference_.block(i, 0, 1, ACADO_NY) << -position_ref_[0].transpose()+position_ref_[i].transpose(), -velocity_ref_[0].transpose()+velocity_ref_[i].transpose(),   0, 0, 0,  0, 0, 0,   Eigen::MatrixXd::Constant(1,ACADO_NU,0); //change
	  }
	  else {
		  reference_.block(i, 0, 1, ACADO_NY) << position_ref_[i].transpose(), velocity_ref_[i].transpose(),   0, 0, 0,  angVel_ref_vector_.transpose(),   Eigen::MatrixXd::Constant(1,ACADO_NU,0); //change
	  }

	  acado_online_data_.block(i, ACADO_NOD - 3, 1, 3) <<  acceleration_ref_[i].transpose();  // + estimated_disturbances.transpose();
  }

  acado_online_data_.block(ACADO_N, ACADO_NOD - 3, 1, 3) << acceleration_ref_[ACADO_N].transpose(); //  + estimated_disturbances.transpose();

  if(use_error_dynamics_){
 // referenceN_ << 0, 0, 0, 0, 0, 0;
  referenceN_ <<-position_ref_[0].transpose()+position_ref_[ACADO_N].transpose(), -velocity_ref_[0].transpose()+velocity_ref_[ACADO_N].transpose();//change
  }
  else{
	  referenceN_ << position_ref_[ACADO_N].transpose(), velocity_ref_[ACADO_N].transpose();//change
  }

  if(use_error_dynamics_){
	  x_0 << position_error_, velocity_error_, mrp_error_, angVel_error_;//change
  }
  else{
	    x_0 << odometry_.position_W, odometry_.getVelocityWorld(), mrp_error_, odometry_.angular_velocity_B;//change
  }

  Eigen::Map<Eigen::Matrix<double, ACADO_NX, 1>>(const_cast<double*>(acadoVariables.x0)) = x_0;
  Eigen::Map<Eigen::Matrix<double, ACADO_NY, ACADO_N>>(const_cast<double*>(acadoVariables.y)) =
      reference_.transpose();
  Eigen::Map<Eigen::Matrix<double, ACADO_NYN, 1>>(const_cast<double*>(acadoVariables.yN)) =
      referenceN_.transpose();
  Eigen::Map<Eigen::Matrix<double, ACADO_NOD, ACADO_N + 1>>(const_cast<double*>(acadoVariables.od)) =
      acado_online_data_.transpose();

  ros::WallTime time_before_solving = ros::WallTime::now();


  acado_preparationStep();

  int acado_status = acado_feedbackStep();


  solve_time_average_ += (ros::WallTime::now() - time_before_solving).toSec() * 1000.0;


  state_ = Eigen::Map<Eigen::Matrix<double, ACADO_N + 1, ACADO_NX, Eigen::RowMajor>>(
	      acadoVariables.x);

  input_ = Eigen::Map<Eigen::Matrix<double, ACADO_N, ACADO_NU, Eigen::RowMajor>>(
        acadoVariables.u);

  command_forces_ = input_.block(0, 0, 1, ACADO_NU).transpose();


  if (std::isnan(command_forces_(0)) || std::isnan(command_forces_(1)) || std::isnan(command_forces_(2)) || std::isnan(command_forces_(3))
      || acado_status != 0) {
    ROS_WARN_STREAM("Nonlinear MPC: Solver failed with status: " << acado_status);
    ROS_WARN("reinitializing...");

    if (record_to_csv_) {
		static std::ofstream failStatusFile("/home/barza/bagfiles/CSVfiles/fail_Status.csv",  std::ios_base::app);

		Eigen::VectorXd states(ACADO_NX);
		Eigen::VectorXd inputs(ACADO_NU);

		for (size_t i = 0; i < ACADO_N; i++) {

			  states = state_.block(i, 0, 1, ACADO_NX).transpose();
			  inputs=input_.block(i, 0, 1, ACADO_NU).transpose();

			  failStatusFile<< states(0) << "," << states(1) << "," << states(2) << "," << states(3) << "," << states(4) << "," << states(5) << "," << states(6) << "," << states(7) << "," << states(8) << "," << states(9) << "," << states(10) << "," << states(11) << "," << "," << inputs(0)<< "," << inputs(1)<< "," << inputs(2) << "," << inputs(3)<< std::endl;
			}

		states = state_.block(ACADO_N, 0, 1, ACADO_NX).transpose();
		failStatusFile<< states(0) << "," << states(1) << "," << states(2) << "," << states(3) << "," << states(4) << "," << states(5) << "," << states(6) << "," << states(7) << "," << states(8) << "," << states(9) << "," << states(10) << "," << states(11) << std::endl;

		failStatusFile << "end" << std::endl;
    }

    initializeAcadoSolver (x_0);

    double norm_uref = normalizeForce(u_ref_);

    *ref_forces<< Eigen::MatrixXd::Constant(1,ACADO_NU,u_ref_);
    *ref_normforces << Eigen::MatrixXd::Constant(1,ACADO_NU,norm_uref);

    return;
  }

  Eigen::VectorXd command_normforces(ACADO_NU);

  for (size_t i = 0; i < ACADO_NU; i++){
	  command_normforces(i)= normalizeForce(command_forces_(i));

  }
	  //ROS_INFO_STREAM("Reached this point: f1 " << command_forces_(0) << "\t" << "f2 : \t" << command_forces_(1) << "\t" << "f3 : \t" << command_forces_(2) << "\t" << "f4 \t" << command_forces_(3));
	  *ref_forces = command_forces_;
      *ref_normforces = command_normforces;


   //******************************************* Recording to csv **********************************

   if (record_to_csv_) {

   static std::ofstream qnormFile("/home/barza/bagfiles/CSVfiles/qnorm.csv",  std::ios_base::app);
   Eigen::Vector4d q;
   Eigen::Vector3d mrp;

   for (size_t i = 0; i < ACADO_N+1; i++) {
       mrp = state_.block(i, 6, 1, 3).transpose();

       q << ( 1-mrp.squaredNorm() )/(1 + mrp.squaredNorm()) , (2/(1 + mrp.squaredNorm())) * mrp;

       qnormFile << q.norm() << std::endl;
     }

   static double time_since_first_command_=0;
   static bool initialise_recording=true;
   static std::ofstream forcesFile("/home/barza/bagfiles/CSVfiles/forces.csv",  std::ios_base::app);
   static std::ofstream positionFile("/home/barza/bagfiles/CSVfiles/position.csv",  std::ios_base::app);
   static std::ofstream velocityFile("/home/barza/bagfiles/CSVfiles/velocity.csv",  std::ios_base::app);
   static std::ofstream angVelFile("/home/barza/bagfiles/CSVfiles/angVel.csv",  std::ios_base::app);
   static std::ofstream orientationFile("/home/barza/bagfiles/CSVfiles/orientation.csv",  std::ios_base::app);
   static std::ofstream timeFile("/home/barza/bagfiles/CSVfiles/time.csv",  std::ios_base::app);
   static std::ofstream predState("/home/barza/bagfiles/CSVfiles/pred_State.csv",  std::ios_base::app);
   static std::ofstream statesFile("/home/barza/bagfiles/CSVfiles/States.csv",  std::ios_base::app);
   static std::ofstream costFile("/home/barza/bagfiles/CSVfiles/cost.csv",  std::ios_base::app);

   Eigen::Vector3d worldvel= odometry_.getVelocityWorld();


if (initialise_recording)
{
	 if (ACADO_NU==4)
		   {forcesFile << 0 << "," << 0 << "," << 0  << "," << 0 << std::endl;}
	else
           {forcesFile << 0 << "," << 0 << "," << 0  << "," << 0 << ","<< 0 << "," << 0 << std::endl;}

	positionFile<< odometry_.position_W.x() << "," << odometry_.position_W.y() << "," << odometry_.position_W.z() << std::endl;
	velocityFile<< worldvel.x() << "," << worldvel.y() << "," << worldvel.z() << std::endl;
	angVelFile << odometry_.angular_velocity_B.x() << "," << odometry_.angular_velocity_B.y() << "," << odometry_.angular_velocity_B.z() << std::endl;
	orientationFile << odometry_.orientation_W_B.w() << "," << odometry_.orientation_W_B.x() << "," << odometry_.orientation_W_B.y()  << "," << odometry_.orientation_W_B.z()  << std::endl;
	timeFile << time_since_first_command_<< std::endl;

	initialise_recording=false;
}
else
{
	time_since_first_command_+=(ros::WallTime::now() - loop_start_time_).toSec();

	if (ACADO_NU==4)
		   {forcesFile << command_forces_(0) << "," << command_forces_(1) << "," << command_forces_(2)  << "," << command_forces_(3) << std::endl;}
    else
		  {forcesFile << command_forces_(0) << "," << command_forces_(1) << "," << command_forces_(2)  << "," << command_forces_(3) << "," << command_forces_(4) << "," << command_forces_(5) << std::endl;
	}

	positionFile<< odometry_.position_W.x() << "," << odometry_.position_W.y() << "," << odometry_.position_W.z() << std::endl;
	velocityFile<< worldvel.x() << "," << worldvel.y() << "," << worldvel.z() << std::endl;
	angVelFile << odometry_.angular_velocity_B.x() << "," << odometry_.angular_velocity_B.y() << "," << odometry_.angular_velocity_B.z() << std::endl;
	orientationFile << odometry_.orientation_W_B.w() << "," << odometry_.orientation_W_B.x() << "," << odometry_.orientation_W_B.y()  << "," << odometry_.orientation_W_B.z()  << std::endl;
	timeFile << time_since_first_command_<< std::endl;

}
loop_start_time_= ros::WallTime::now();

Eigen::VectorXd states(ACADO_NX);
Eigen::VectorXd inputs(ACADO_NU);
Eigen::Vector3d position_ref=  position_ref_.front();
Eigen::Vector3d velocity_ref=  velocity_ref_.front();
Eigen::Vector3d angVel_ref=  angVel_ref_vector_;

Eigen::VectorXd states_inputs(ACADO_NX+ACADO_NU);
Eigen::VectorXd statesN(6);
double cost=0;

double pred_time=time_since_first_command_;
for (size_t i = 0; i < ACADO_N; i++) {
	  states = state_.block(i, 0, 1, ACADO_NX).transpose();
      inputs=input_.block(i, 0, 1, ACADO_NU).transpose();

       if (i==0){
          statesFile<< pred_time << "," << states(0) + position_ref.x() << "," << states(1) + position_ref.y() << "," << states(2) + position_ref.z() << "," << states(3) + velocity_ref.x() << "," << states(4) + velocity_ref.y()<< "," << states(5) + velocity_ref.z()<< "," << states(6) << "," << states(7) << "," << states(8) << "," << states(9) + angVel_ref.x() << "," << states(10) + angVel_ref.y() << "," << states(11) + angVel_ref.z()<< std::endl;
       }
      predState<< pred_time << "," << states(0) + position_ref.x() << "," << states(1) + position_ref.y() << "," << states(2) + position_ref.z() << "," << states(3) + velocity_ref.x() << "," << states(4) + velocity_ref.y()<< "," << states(5) + velocity_ref.z()<< "," << states(6) << "," << states(7) << "," << states(8) << "," << states(9) + angVel_ref.x() << "," << states(10) + angVel_ref.y() << "," << states(11) + angVel_ref.z()<< std::endl;
      //predState<<  pred_time << "," << states(0) << "," << states(1) << "," << states(2) << "," << states(3) << "," << states(4) << "," << states(5) << "," << states(6) << "," << states(7) << "," << states(8) << "," << states(9) << "," << states(10) << "," << states(11) << std::endl;
      states_inputs<< states , inputs - Eigen::MatrixXd::Constant(ACADO_NU,1,u_ref_);
      cost+= states_inputs.transpose() * W_ * states_inputs;


      pred_time += prediction_sampling_time_;
    }

states = state_.block(ACADO_N, 0, 1, ACADO_NX).transpose();
predState<< pred_time << "," << states(0) + position_ref.x() << "," << states(1) + position_ref.y() << "," << states(2) + position_ref.z() << "," << states(3) + velocity_ref.x() << "," << states(4) + velocity_ref.y()<< "," << states(5) + velocity_ref.z()<< "," << states(6) << "," << states(7) << "," << states(8) << "," << states(9) + angVel_ref.x() << "," << states(10) + angVel_ref.y() << "," << states(11) + angVel_ref.z()<< std::endl;
//predState<<  pred_time << "," << states(0) << "," << states(1) << "," << states(2) << "," << states(3) << "," << states(4) << "," << states(5) << "," << states(6) << "," << states(7) << "," << states(8) << "," << states(9) << "," << states(10) << "," << states(11) << std::endl;
// predState << "end" << std::endl;
statesN = state_.block(ACADO_N, 0, 1, 6).transpose();
cost+=statesN.transpose() * WN_ * statesN;
costFile<<cost<<","<< angVel_ref_<<std::endl;

//******************************************* Recording solve time **********************************
   static std::ofstream avgSolveTimeFile("/home/barza/bagfiles/CSVfiles/avgSolveTime.csv",  std::ios_base::app);
   static std::ofstream solveTimeFile("/home/barza/bagfiles/CSVfiles/solveTime.csv",  std::ios_base::app);

   double diff_time = (ros::WallTime::now() - starting_time).toSec();

     static int counter = 0;
     if (counter > 100) {
       //ROS_INFO_STREAM("average solve time: " << solve_time_average_ / counter << " ms");
       avgSolveTimeFile << solve_time_average_ / counter << std::endl;
       //ROS_INFO_STREAM("Controller loop time : " << diff_time*1000.0 << " ms");

       solve_time_average_ = 0.0;
       counter = 0;
     }
     counter++;
     solveTimeFile<<diff_time*1000.0<< std::endl;

//******************************************* Recording to csv **********************************

   }




} 



bool NonlinearModelPredictiveControl::getCurrentReference(
    mav_msgs::EigenTrajectoryPoint* reference) const
{
  assert(reference != nullptr);

  (*reference).position_W = position_ref_.front();
  (*reference).velocity_W = velocity_ref_.front();
  (*reference).acceleration_W = acceleration_ref_.front();
  (*reference).setFromYaw(yaw_ref_.front());
     return true;
}

bool NonlinearModelPredictiveControl::getCurrentReference(
    mav_msgs::EigenTrajectoryPointDeque* reference) const
{
  assert(reference != nullptr);

  (*reference).clear();

  for (size_t i = 0; i < position_ref_.size(); i++) {
    mav_msgs::EigenTrajectoryPoint pnt;
    pnt.position_W = position_ref_.at(i);
    pnt.velocity_W = velocity_ref_.at(i);
    pnt.acceleration_W = acceleration_ref_.at(i);
    pnt.setFromYaw(yaw_ref_.at(i));
    (*reference).push_back(pnt);
  }

  return true;
}

bool NonlinearModelPredictiveControl::getPredictedState(
    mav_msgs::EigenTrajectoryPointDeque* predicted_state) const
{
  assert(predicted_state != nullptr);

  for (size_t i = 0; i < ACADO_N + 1; i++) {
    mav_msgs::EigenTrajectoryPoint pnt;
    pnt.position_W = state_.block(i, 0, 1, 3).transpose();
    (*predicted_state).push_back(pnt);
  }

  return true;
}

}



