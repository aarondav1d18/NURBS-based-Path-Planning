import matplotlib.pyplot as plt
from geomdl import BSpline,convert,knotvector
import sys
import numpy as np
from typing import Any, List, Tuple
from core.Cone import Cone, ConeColour
from core.Point import Point
import config.config as config
from cone_ordering.cone_ordering import order_blue_and_yellow_cones


'''
This file will just produce a NURBS curve for the boundaries of the whole track and then average that
to produce a 'middle' line. However it does have some downfalls on ceratin corners etc etc which are being
worked on.
Bigger track does not work right now aswell and will be looked into
'''

def map_track(cones: List[Cone]) -> None:
    '''
    This function will take in a list of cones and plot it with its colour.
    This will also plot the last cones as red so you can see if the cone ordering
    is working as should. If not then just plot what the 'car' is taking in and mess
    around with how I reverse them etc etc

    Args:
        cones (List[Cone]): List of Cone to be plotted
    '''
    for cone in cones:
        rx,ry = cone.pos.x, cone.pos.y
        if cone == cones[-1]:
            plt.plot(rx,ry, 'o', c='red')
            continue
        if cone.colour == ConeColour.YELLOW:
            plt.plot(rx,ry, 'o', c='yellow')
        elif cone.colour == ConeColour.BLUE:
            plt.plot(rx,ry, 'o', c='blue')
        elif cone.colour == ConeColour.UNKNOWN:
            plt.plot(rx,ry, 'o', c='gray')
        elif cone.colour == ConeColour.ORANGE:
            plt.plot(rx,ry, 'o', c='orange')

def read_cones(filename: str) -> List[Cone]:
    '''
    Takes in a file from the tracks folder and reads the cones
    from it converting it to the type that is needed and appending
    to a list based on the colour of the cone.

    Args:
        filename (str): path to the file to read from

    Returns:
        list: list of blue and a list of yellow cones
    '''
    coords: List[Tuple[float, float, ConeColour]] = []

    with open(filename) as f:
        for l in f.readlines() :
            splat = l.split(":")
            colour, coord_list = splat[0], splat[1]
            
            coord_list = coord_list.replace("[", "").replace("]", "").replace("\s+", "")
            for coord_str in coord_list.split("),"):
                x, y = coord_str.split(",")
                pos = Point()
                pos.x, pos.y = float(x.replace("(", "")), float(y.replace(")", ""))
                coords.append((pos, ConeColour.YELLOW if colour == "yellow" else ConeColour.BLUE))

    yellow_cones = []
    blue_cones = []
    for coord in coords :
        pos = coord[0]
        dist = np.hypot(pos.x - config.ORIGIN.x, pos.y - config.ORIGIN.y)
        proportional_covariance = (dist / config.MAX_DISTANCE_FROM_ORIGIN) * config.MAX_FAKE_COVARIANCE_UNCERTAINTY
        if coord[1] == ConeColour.YELLOW:
            yellow_cones.append(
                Cone(coord[0], coord[1], [[proportional_covariance, 0.], [0., proportional_covariance]])
            )
        elif coord[1] == ConeColour.BLUE:
            blue_cones.append(
                Cone(coord[0], coord[1], [[proportional_covariance, 0.], [0., proportional_covariance]])
            )
    
    return yellow_cones, blue_cones

def generateNURBS(cones: List[Cone]) -> BSpline:
    '''
    This will take in a list of cones and calculate the weight of each cone
    based on the distance from the 'car' and the curvature that the surrounding
    cones make with the current cone. This will then generate a NURBS curve.
    can look into differetn ways to do this.

    Args:
        cones (List[Cone]): A list of cones to generate the NURBS curve from

    Returns:
        BSpline: The generated NURBS
    '''
    control_points = []
    crv = BSpline.Curve()
    for point in cones:
        # using a exponetial decay multiplier for weight
        dist = np.hypot(point.pos.x, point.pos.y)
        distWeight = 1 / (1 + dist / config.MAX_DISTANCE_FROM_ORIGIN)
        if point == cones[-1]:
            prev_point = cones[-2]
            next_point = point
        else:
            prev_point = cones[cones.index(point) - 1]
            next_point = cones[cones.index(point) + 1]

        # Calculate vectors
        vec1 = np.array([point.pos.x - prev_point.pos.x, point.pos.y - prev_point.pos.y])
        vec2 = np.array([next_point.pos.x - point.pos.x, next_point.pos.y - point.pos.y])

        # Calculate the angle between the vectors
        if np.linalg.norm(vec1) == 0 or np.linalg.norm(vec2) == 0:
            angle = 0
        else:
            angle = np.arccos(np.clip(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)), -1.0, 1.0))

        # Calculate curvature
        curvature = 2 * np.sin(angle / 2) / np.linalg.norm(vec1)
        # Calculate weight based on curvature
        curvateWeight = 1 / (1 + curvature)
        if curvature > 0.1:
            curvateWeight = 1
        if curvateWeight == np.nan:
            curvateWeight = 1
        weight = (distWeight + curvateWeight) / 2
        weighted_point = [point.pos.x * weight, point.pos.y * weight, weight]
        control_points.append(weighted_point)
    
    crv.degree = 4
    crv.ctrlpts = control_points
    knotvectors = knotvector.generate(crv.degree, len(control_points))
    crv.knotvector = knotvectors
    crv_rat = convert.bspline_to_nurbs(crv)
    return crv_rat

def averageNURBS(blueNURBS: List[List], yellowNURBS: List[List]) -> List[List]:
    '''
    A really simple function for averaging the two NURBS curves
    This may not be the best way to do this but it works and produces a curve
    as expected. Also the lenght of each NURBS will be the same so this does work
    for at least all the times i have tested it.

    Args:
        blueNURBS (List[List]): List of points for the blue NURBS
        yellowNURBS (List[List]): List of points for the yellow NURBS

    Returns:
        Bspline: The actual averaged NURBS in the BSpline format. used for plotting vel
        List[List]: A list of points that are the average of the two NURBS
    '''
    # will always be the same length. each will be 100 points no matter how many control points
    averagedNURBS = []
    for i,point in enumerate(blueNURBS):
        x = (point[0] + yellowNURBS[i][0]) / 2
        y = (point[1] + yellowNURBS[i][1]) / 2
        weight = (point[2] + yellowNURBS[i][2]) / 2
        averagedNURBS.append([x,y,weight])

    return averagedNURBS

def map_line(points: List[List], colour: str) -> None:
    '''
    This will take the points of a NURBS curve and plot it on the graph
    and take the weights back out of the points.
    '''
    x = []
    y = []
    for point in points:
        rx,ry = point[0],point[1]
        x.append(rx / point[2])
        y.append(ry / point[2])
    plt.plot(x, y, colour, '-')

def main(filename: str) -> None:
    '''
    This is the main entry point to the program.
    This will read the cone data and sort the cones based on colour and distance to the 'car'
    it will then plot the NURBS
    '''
    yellow_cones, blue_cones = read_cones(filename)
    orderedBlue, orderedYellow = order_blue_and_yellow_cones(blue_cones, yellow_cones)
    
    YellowNURBS = generateNURBS(orderedYellow)
    BlueNURBS = generateNURBS(orderedBlue)

    map_track(orderedBlue)
    map_track(orderedYellow)

    pointsYellow = YellowNURBS.evalpts
    pointsBlue = BlueNURBS.evalpts
    averaged = averageNURBS(pointsBlue, pointsYellow)

    map_line(pointsYellow, 'yellow')
    map_line(pointsBlue, 'blue')
    map_line(averaged, 'red')

    plt.show()

if __name__ == "__main__":
    if len(sys.argv) <= 1 :
        print("!!!!!!!!!!!!!!!!!!!! I NEED A FILE NAME !!!!!!!!!!!!!!!!!!!!")
        sys.exit(1)
    main(sys.argv[1])
    plt.show()
