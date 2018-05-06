Spectral Partially Observable Poisson Process
====

Spectral-POPP is a package representing the work described both in "A Poisson-spectral model for modelling temporal patterns in human data observed by a robot" (the Spectral part of the Poisson) and "Efficient Bayesian Methods for Counting Processes in Partially Observable Environments" (the partial observability part of the Poisson).

This package also contains the University of Birmingham dataset stored in folder "data".

Prerequisites
-------------

- pyyaml (>=3.12)
- matplotlib (>=1.3.1)
- scipy (>=0.19.1)


Getting started (general steps)
-------------------------------
1. git clone this repository

2. python-run the examples in folder "examples" such as activity_counter.py or detection_counter.py:
   ```
    $ python activity_counter.py -m 1 -l 1

    ```
   One can type ```-h``` to see other options.


The University of Birmingham Dataset
------------------------------------

The dataset is a collection of counts over time from three different person detectors attached to a mobile robot. The detectors are a leg detector (leg.yaml), an upper body detector (upper_body.yaml), and a change detector (scene.yaml). The dataset was gathered during a 69-day deployment from the lower-ground floor of a computer science building at the University of Birmingham. The mobile robot counts and observes the number of people passing by as it patrols around the perimeters of the lower-ground floor. Map is shown below.

![marker](https://github.com/ferdianjovan/spectral_popp/doc/map.png)

Each minute, each detector stores the timestamp of that minute where the detection was made (if there was a positive detection, i.e., there was an activity within that minute). The ground truth of each detection is stored in "present_activity.yaml" for the positive event (activity happened within a minute interval) and "absent_activity.yaml" for the negative event (no activity happened within a minute interval).  "observation_history.yaml" stores information where and when the robot was observing. 

As a mobile robot can not fully sense its environment, it can only perceive partial data at a particular time and place. Moreover, the robot’s patrol policy also affects where and when detections are perceived. Consequently, the detections are temporally and spatially scattered and they are not uniformly distributed across space. The detections are then organised according to
time/date and the spatial region where each detection was made. 
