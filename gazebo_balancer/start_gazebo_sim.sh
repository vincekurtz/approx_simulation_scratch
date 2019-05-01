#!/bin/bash

##
#
# A bash script to start gazebo from matlab, since matlab doesn't
# play well with ros launch files.
#
##

# Defines the launch file we want to run
ros_pkg="double_pendulum_gazebo"
launch_file="double_pendulum_world.launch"

# Matlab runs in it's own little world with a different LD_LIBRARY_PATH, so we need to set
# this to match the LD_LIBRARY_PATH of the real system.
export LD_LIBRARY_PATH="/usr/local/cuda-9.1/lib64:/home/vjkurtz/catkin_ws/devel/lib:/opt/ros/kinetic/lib:/opt/ros/kinetic/lib/x86_64-linux-gnu:/home/vjkurtz/stg/lib:/usr/lib/x86_64-linux-gnu/gazebo-7/plugins:/opt/gurobi801/linux64/lib:/opt/tomlab/shared"

# Run the simulation in the background.
# We'll pipe all output to /dev/null
roslaunch $ros_pkg $launch_file &> /dev/null &
PID=$!

echo $PID  # knowing the PID of the simulation will allow us to kill it from matlab
