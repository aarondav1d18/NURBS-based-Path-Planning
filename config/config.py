# genearte config file that includes max distance from origin, max fake covariance uncertainty, origin, and the cone colours etc etc etc etc etc 
from core.Point import Point
import math
import numpy as np

MAX_DISTANCE_FROM_ORIGIN = 10
MAX_FAKE_COVARIANCE_UNCERTAINTY = 10
ORIGIN = Point()
ORIGIN.x, ORIGIN.y = -13.0, 10.3
HEADING = np.pi/2
MAX_DISTANCE_BETWEEN_CONES = 8
MAXIMUM_AVERAGE_CURVATURE_FOR_A_CONE_ARRAY = np.radians(35) / 1

