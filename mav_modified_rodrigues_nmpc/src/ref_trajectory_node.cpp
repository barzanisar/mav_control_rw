/*
 Copyright (c) 2015, Barza Nisar, Mina Kamel, ASL, ETH Zurich, Switzerland

 You can contact the author at <nisarb@student.ethz.ch>

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

#include <thread>
#include <chrono>

#include <Eigen/Core>
#include <mav_msgs/conversions.h>
#include <mav_msgs/default_topics.h>
#include <ros/ros.h>
#include <std_srvs/Empty.h>
#include <trajectory_msgs/MultiDOFJointTrajectory.h>

int main(int argc, char** argv) {
  ros::init(argc, argv, "ref_trajectory_node");
  ros::NodeHandle nh;
  // Create a private node handle for accessing node parameters.
  ros::NodeHandle nh_private("~");
  ros::Publisher trajectory_pub =
      nh.advertise<trajectory_msgs::MultiDOFJointTrajectory>(
          mav_msgs::default_topics::COMMAND_TRAJECTORY, 1);
 
ROS_INFO("Start reference trajectory.");

double traj_time= 120; //5sec
double prediction_sampling_time=0.01;
int num_points=traj_time/prediction_sampling_time;
std::cout << " num_points: " <<  num_points << std::endl;
double radius=2; //meters
double omega=3.1416/4; //pi/4 rad/sec 
static const int64_t kNanoSecondsInSecond = 1000000000;
  

trajectory_msgs::MultiDOFJointTrajectoryPtr msg(new trajectory_msgs::MultiDOFJointTrajectory);
  msg->header.stamp = ros::Time::now();
  msg->points.resize(num_points);
  msg->joint_names.push_back("base_link");


  Eigen::Vector3d position;
  Eigen::Vector3d velocity;
  Eigen::Vector3d acceleration;
  double yaw=0;
  int64_t time_from_start_ns = 0;



  for (size_t i = 0; i < num_points; ++i) {

// 8 trajectory
   position.x()=radius*cos(omega*i*prediction_sampling_time);
   position.y()=radius*0.5*sin(2*omega*i*prediction_sampling_time);
   position.z()=1;

   velocity.x()= -radius*omega*sin(omega*i*prediction_sampling_time);
   velocity.y()= radius*omega*cos(2*omega*i*prediction_sampling_time);
   velocity.z()=0;

   acceleration.x()= -radius*omega*omega*cos(omega*i*prediction_sampling_time);
   acceleration.y()= -2*radius*omega*omega*sin(2*omega*i*prediction_sampling_time);
   acceleration.z()=0;

  /* //circle
   position.x()=radius*cos(omega*i*prediction_sampling_time) - radius;
   position.y()=radius*sin(omega*i*prediction_sampling_time);
   position.z()=1;

   velocity.x()= -radius*omega*sin(omega*i*prediction_sampling_time);
   velocity.y()= radius*omega*cos(omega*i*prediction_sampling_time);
   velocity.z()=0;

   acceleration.x()= -radius*omega*omega*cos(omega*i*prediction_sampling_time);
   acceleration.y()= -radius*omega*omega*sin(omega*i*prediction_sampling_time);
   acceleration.z()=0;
*/


    mav_msgs::EigenTrajectoryPoint trajectory_point;
    trajectory_point.position_W = position;
    trajectory_point.velocity_W = velocity;
	trajectory_point.acceleration_W = acceleration;
    trajectory_point.setFromYaw(yaw);
    trajectory_point.time_from_start_ns = time_from_start_ns;

    time_from_start_ns += static_cast<int64_t>(0.01 * kNanoSecondsInSecond);

    mav_msgs::msgMultiDofJointTrajectoryPointFromEigen(trajectory_point, &msg->points[i]);
  }

   ros::spinOnce();
   ros::Duration(2).sleep();

  trajectory_pub.publish(msg);
  	
  ros::spin();
  //ros::shutdown();

  return 0;
}
