# Unscented Kalman Filter Project 
Self-Driving Car Engineer Nanodegree Program


In this project implements an Unscented Kalman Filter to estimate the
state of a moving object of interest with noisy lidar and radar
measurements. 
In order to satisfy the project rubric requirements, 
UKF with RMSE below the pre-defined bound has been devised.

This project should be run in the Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

Cmake is used as a build system for the project. In order to compile the code:
1. mkdir build
2. cd build
3. cmake ../src/
4. make
5. ./UnscentedKF

## Important Dependencies

The project has been tested / run on Ubuntu.

* cmake >= 3.5
* make >= 4.1 (Linux, Mac),
* gcc/g++ >= 5.4

## RMSE

[image1]: ./img/ukf_pic.png

The Unscented Kalman filter is run against the test data set, obtaining the following error:

![EKF result pic][image1]


UKF bird's-eye-view:
* correct filter initialization;
* generating sigma points and predicting the state;
* updating filter state with radar and lidar measurements.
