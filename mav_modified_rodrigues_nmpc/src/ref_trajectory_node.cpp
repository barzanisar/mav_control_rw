/*
 * Copyright 2015 Fadri Furrer, ASL, ETH Zurich, Switzerland
 * Copyright 2015 Michael Burri, ASL, ETH Zurich, Switzerland
 * Copyright 2015 Mina Kamel, ASL, ETH Zurich, Switzerland
 * Copyright 2015 Janosch Nikolic, ASL, ETH Zurich, Switzerland
 * Copyright 2015 Markus Achtelik, ASL, ETH Zurich, Switzerland
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
  
 /* trajectory_msgs::MultiDOFJointTrajectory traj;
  traj.header.stamp = ros::Time::now();

  // Default desired position and yaw.

traj.points.resize(num_points);

for(int i=0; i<num_points;  i++) {

traj.points[i].transforms.translation.x()=radius*cos(omega*i*prediction_sampling_time) - radius; //to start from origin
traj.points[i].transforms.translation.y()=radius*sin(omega*i*prediction_sampling_time);
traj.points[i].transforms.translation.z()=5;


traj.points[i].transforms.rotation.x()=0;
traj.points[i].transforms.rotation.y()=0;
traj.points[i].transforms.rotation.z()=0;
traj.points[i].transforms.rotation.w()=1;

traj.points[i].velocities.linear.x()= -radius*omega*sin(omega*i*prediction_sampling_time);
traj.points[i].velocities.linear.y()= radius*omega*cos(omega*i*prediction_sampling_time);
traj.points[i].velocities.linear.z()=0;

traj.points[i].accelerations.linear.x()= -radius*omega*omega*cos(omega*i*prediction_sampling_time);
traj.points[i].accelerations.linear.y()= -radius*omega*omega*sin(omega*i*prediction_sampling_time);
traj.points[i].accelerations.linear.z()=0;

}*/

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
