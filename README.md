This is a scratch repository for code relating to approximate simulation for legged robots. 

## Example Usage:

For a simple matlab simulation of using the Linear Inverted Pendulum (LIP) model and a 
simulation interface to balance a double pendulum with a fixed ground contact, see `balancer/lip_dp_demo.m`.

For a more in-depth simulation using ROS/Gazebo to similarly balance a double pendulum mounted to a 
light foot/platform, see `gazebo_balancer/gazebo_balancer.m`.

## Dependencies

For the matlab simulation only:
- Matlab
- [spatial\_v2](http://royfeatherstone.org/spatial/v2/)

For the Gazebo simulation:
- Matlab 
- spatial\_v2
- ROS Kinetic
- Gazebo 7
- [This double pendulum model](https://github.com/vincekurtz/double_pendulum_gazebo)
