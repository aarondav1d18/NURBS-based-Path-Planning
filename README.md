# Visualising NURBS Based Path Planning

## Background
This was made during my time in UGR and is an new approach to computing a path for a autonomous car to follow. This is a basic implementation that includes some data structures for cones etc.

The code will order the cones based on closest to the 'car' to furtherest away. From here it will produce a NURBS on both the inner and outer bounds of the track and then run a very simple averaging function that will genrate a middle line through the middle of the track. 

The idea of this was to produce a middle line through the centre of the track for the first lap that is very controllable and will also produce a nice line.


## Issues To Be Addressed
Currently the code does not perform as expected on the bigger_track.csv track. I have plans to fix my implementation to allow for this track to work as expected.

There currently is a `driveTrack.py` file that is not working very well but will drive certain tracks without editing. I intend to fix this to produce a much better visual, however this file is to show a basic idea of how it will perform in real world situations where the whole track is not visable. This file works well for the `big_track.csv` file and does work for `small_track.csv` however produces bad lines sometimes.