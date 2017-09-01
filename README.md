# CarND-Path-Planning
[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

This is the project repository for **Project No. 1 Path Planning**, part of Term 3 _Path Planning, Concentrations, and Systems_ of Udacity Self-Driving Car Nanodegree program, submitted by Albert Killer in August 2017. 

## Introduction
The goal of this project was to create a Path Planning algorithm which safely guides a self-driving vehicle through highway traffic within the [Udacity Term3 Simulator](https://github.com/udacity/self-driving-car-sim/releases). Driving on a highway also involves maneuvers like changing velocity or lanes in order to pass slower vehicles. At the same time safety is the number one priority for a self-driving car. To pass the test, the car has to avoid any collissions, is not allowed to exceed a maximum velocity of 50 mph and has to apply a comfortable driving style for potential human passengers. This means maximum acceleration of 10 m/s² and maximum jerk of 10 m/s³ are not exceeded. 

## Trajectory Generation
As an input for the simulator a trajectory is generated on basis of [frenet coordinates](https://en.wikipedia.org/wiki/Frenet%E2%80%93Serret_formulas), which are converted from the simulators output of the car’s position in cartesian coordinates. In order to **keep the car on the lane as smooth as possible** we  have to use the previous path points in order to assure more continuity between the old path and the new path. On this basis we then only create **anchor points** of the trajectory, using the simulated GPS coordinates. The space between those anchor points is filled up by interpolating the trajectory using a **Cubic Spline interpolation** ([spline.h](http://kluge.in-chemnitz.de/opensource/spline/)). The number of points filled up has to be defined on basis of the target velocity. The bigger the distance between the trajectory’s waypoints the faster the car will go as the simulator moves the car from point to point every 20 ms. Alternatively *Polynomial Trajectory Generation* (PTG) could by applied using a Quintic Polynomial Solver to get a *Jerk Minimized Trajectory* (JMT). 

## Behavior Planning:
Other than only following the lane our car should be able to execute certain maneuvers. If there is a slower car ahead on the same lane the self-driving car should change lanes or slow down. In order to do this in a safe way prediction have to be made using sensord fusion data from the simulator’s output stream. This involves checking for vehicles in car’s lane as well as potential target lanes. 

```c++
for(int i = 0; i < sensor_fusion.size(); i++){
   d = sensor_fusion[i][6];
   // Check if car is in ego vehicle's lane
   if(d < (2+4*lane+2) && d > (2+4*lane-2) ){
   ...
```
The decision to change to another lane is made by a transition function which uses position and velocity of sorounding cars relative to the ego car’s position in order to return the next best state (in this case lane) to transit to. The hereby implemented simple *Finite State Machine* (FSM) ignores useful states like “preparing lane change”. For the purpose of this simulation additional safety checks do the trick. 
To decide for the best lane **cost functions** evaluate the following basic situations:

1. Car ahead and far away
2. Car ahead but too close
3. Car behind and far away
4. Car behind, not too close but too fast
5. Car behind but too close 

To increase safety **additional and independent checks** were included, verifying that collisions are avoided under all circumstances by checking that:

* No vehicles on target lane are very close
* No vehicles on target lane are approaching from behind with a higher speed, before ...
* and during lange change
* If car keeps the lane, velocity is reduced to avoid collision with slower vehicles

![Screenshot of simulation result](Screenshot%20from%202017-09-01%2023-28-03.png?raw=true "Screenshot of simulation result")

## Discussion
As a result the implemented Path Planner is able to drive the car safely along the highway for several miles around the simulation’s track. 
In order to increase speed and safety while maneuvering through dense traffic, decision making of the Behavioural Planner can be improved. This is achieved by adding more detailed cost functions which for example involve the difference of several car’s velocities or smaller distance classifications.   
Safety while performing lane changes can be increased as well by switching to a more advanced FSM using additional states. 
By updating the relative velocity *rel_vel* within the trajectory generation while going through every point, the reaction time of speed adjustments could be improved.  
