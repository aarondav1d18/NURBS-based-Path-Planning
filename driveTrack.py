import matplotlib.pyplot as plt
from geomdl import BSpline,convert,knotvector
import sys
import numpy as np
from typing import Any, List, Tuple
from core.Cone import Cone, ConeColour
from core.Point import Point
import config.config as config
from core.PointVelocity import PointVelocity
from cone_ordering.cone_ordering import order_blue_and_yellow_cones
import time

PARTICLE_MASS = 300
MAX_FORCE = 1.1 * PARTICLE_MASS * 9.81
NUM_OF_CONES = 10
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
        rx,ry = roation_matrix(cone.pos.x,cone.pos.y, np.radians(90))
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

def roation_matrix(x: float, y: float, theta: float):
    '''
    Function to rotate a point around 0,0 by a certain angle. havnt done anything
    with this yet but could be nice for a local plot.
    '''
    theta = 0
    rx = x * np.cos(theta) - y * np.sin(theta)
    ry = x * np.sin(theta) + y * np.cos(theta)
    return rx,ry

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
            # print(splat)
            colour, coord_list = splat[0], splat[1]
            
            coord_list = coord_list.replace("[", "").replace("]", "").replace("\s+", "")
            for coord_str in coord_list.split("),"):
                x, y = coord_str.split(",")
                pos = Point()
                pos.x, pos.y = float(x.replace("(", "")), float(y.replace(")", ""))
                coords.append((pos, ConeColour.YELLOW if colour == "yellow" else ConeColour.BLUE))

    cones: List[Cone] = []
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
        # print(point.pos.x,point.pos.y)
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
        # print(weight)
        weighted_point = [point.pos.x * weight, point.pos.y * weight, weight]
        control_points.append(weighted_point)
    
    crv.degree = len(control_points) // 3
    # crv.degree = 4
    # if crv.degree < 3:
    #     crv.degree = 3
    # crv.degree = len(control_points) // 3
    crv.degree = 3
    print(crv.degree)
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
    print(len(blueNURBS))
    print(len(yellowNURBS))
    for i,point in enumerate(blueNURBS):
        x = (point[0] + yellowNURBS[i][0]) / 2
        y = (point[1] + yellowNURBS[i][1]) / 2
        weight = (point[2] + yellowNURBS[i][2]) / 2
        averagedNURBS.append([x,y,weight])
    
    crv = BSpline.Curve()
    crv.degree = len(averagedNURBS) // 3
    crv.ctrlpts = averagedNURBS
    knotvectors = knotvector.generate(crv.degree, len(averagedNURBS))
    crv.knotvector = knotvectors
    crv_rat = convert.bspline_to_nurbs(crv)

    return averagedNURBS, crv_rat

def convertNURBS(nurbs: BSpline ,points: List[List]) -> List[PointVelocity]:
    '''
    This is used to generate a PointVelocity for the averaged NURBS. This is only
    used if you want a plot of the velocity. But that plot is fucked atm so no real use.

    Args:
        nurbs (BSpline): The averaged NURBS in the BSpline format
        points (List[List]): The points of the averaged NURBS

    Returns:
        List[PointVelocity]: A list of the PointVelocity for the averaged NURBS
    '''
    ## This is raw velocity calculation would need to be inputted into a controller or something
    ## to generate a decent drivable vel but this is just showing the path so oh well
    Points = []

    nurbs.evaluate()
    u_min, u_max = nurbs.knotvector[nurbs.degree], nurbs.knotvector[-nurbs.degree-1]
    param_values = np.linspace(u_min, u_max, nurbs.sample_size)
    derivatives = [nurbs.derivatives(u, 1)[1] for u in param_values]
    velocities = [np.sqrt(d[0]**2 + d[1]**2 + d[2]**2) for d in derivatives]

    for i,point in enumerate(points):
        weight = point[-1]
        newPoint = PointVelocity()
        newPoint.x = point[0]/weight
        newPoint.y = point[1]/weight
        newPoint.vel = ((MAX_FORCE * (velocities[i] / weight))/ PARTICLE_MASS) ** 0.5
        Points.append(newPoint)

    return Points

def map_line(points: List[List], colour: str) -> None:
    '''
    This will take the points of a NURBS curve and plot it on the graph
    and take the weights back out of the points.
    '''
    x = []
    y = []
    for point in points:
        rx,ry = roation_matrix(point[0],point[1],np.radians(90))
        x.append(rx / point[2])
        y.append(ry / point[2])
    plt.plot(x, y, colour, '-')


def map_PointVel(points: List[PointVelocity]) -> None:
    ## This is me fucking around with generating colours based on velocity/acceleration
    x = []
    y = []
    colors = []
    for i, point in enumerate(points):
        rx, ry = roation_matrix(point.x, point.y, np.radians(90))
        x.append(rx)
        y.append(ry)
        if i > 0:
            dx = point.x - points[i-1].x
            dy = point.y - points[i-1].y
            dt = 1  # Assuming uniform time intervals for simplicity
            acceleration = (point.vel - points[i-1].vel) / dt
            if acceleration > 0:
                colors.append('green')  # Accelerating
            elif acceleration < 0:
                colors.append('red')  # Braking
            else:
                colors.append('yellow')  # Constant velocity
        else:
            colors.append('yellow')  # Starting point

    for i in range(len(x) - 1):
        plt.plot(x[i:i+2], y[i:i+2], color=colors[i])

def main(filename: str) -> None:
    '''
    This is the main entry point to the program.
    This pretty much deals with the logic of making the plot constantly update
    and calling function to order cones and generate the NURBS curves etc.
    
    TODO: implement a way to stop the sim once the car reaches the end of the track
    '''
    running = True
    i = 0
    j = 0
    yellow_cones, blue_cones = read_cones(filename)
    plt.ion()  # Turn on interactive mode

    # Will just order the cones to start then the code will cut the list of cones to take in only a ceratin amount
    orderedBlue, orderedYellow = order_blue_and_yellow_cones(blue_cones, yellow_cones)
    
    YellowNURBS = generateNURBS(orderedYellow)
    BlueNURBS = generateNURBS(orderedBlue)

    map_track(orderedBlue)
    map_track(orderedYellow)

    # Only so it will generate the points once
    pointsYellow = YellowNURBS.evalpts
    pointsBlue = BlueNURBS.evalpts
    averaged, averagedNURBS = averageNURBS(pointsBlue, pointsYellow)

    # print(np.array(averagedNURBS.evalpts))

    # map_line(pointsYellow, 'yellow')
    # map_line(pointsBlue, 'blue')
    # map_line(averaged, 'red')

    plt.show()

    
    totalLine = []
    # plt.show()

    while running:
        if i == 1:
            return
        if i > len(orderedYellow) - NUM_OF_CONES:
            i = 0
            orderedBlue, orderedYellow = order_blue_and_yellow_cones(blue_cones, yellow_cones)
            orderedYellow = orderedYellow[::-1]
            orderedBlue = orderedBlue[::-1]
            orderedYellow += [orderedYellow[0]]
            orderedBlue += [orderedBlue[0]]

        map_track(orderedBlue)
        map_track(orderedYellow)

        useableBlue = orderedBlue[i:i+NUM_OF_CONES]
        if i > 0:
            useableYellow = orderedYellow[(i -1) :i+(NUM_OF_CONES)]
        else:
            useableYellow = orderedYellow[i: i+NUM_OF_CONES]

        YellowNURBS = generateNURBS(useableYellow)
        BlueNURBS = generateNURBS(useableBlue)

        map_track(useableBlue)
        map_track(useableYellow)

        # Only so it will generate the points once
        pointsYellow = YellowNURBS.evalpts
        pointsBlue = BlueNURBS.evalpts
        averaged, averagedNURBS = averageNURBS(pointsBlue, pointsYellow)

        print(np.array(averagedNURBS.evalpts))

        map_line(pointsYellow, 'yellow')
        map_line(pointsBlue, 'blue')
        map_line(averaged, 'red')
        totalLine += averaged
        # map_line(totalLine, 'green')

        # test = convertNURBS(averagedNURBS, averaged)
        # map_PointVel(test)

        plt.draw()  # Update the plot
        plt.pause(0.1)  # Pause to allow for interaction
        time.sleep(2)
        
        # Save the plot to a file
        # add a variable to seperate the files
        # plt.savefig(f'visuals/small_track{j}.png')

        # need to look into this but in theory this will stop the sim once car reaches end of track
        # if np.abs(config.ORIGIN.x - averaged[-1][0]) < 5.5 and np.abs(config.ORIGIN.y - averaged[-1][1]) < 5.5:
        #     running = False

        config.ORIGIN.x = averaged[-1][0]
        config.ORIGIN.y = averaged[-1][1]

        dx = averaged[-2][0] - averaged[-1][0]
        dy = averaged[-2][1] - averaged[-1][1]
        config.HEADING += np.arctan2(dy, dx)

        plt.cla()
        i += NUM_OF_CONES

if __name__ == "__main__":
    if len(sys.argv) <= 1 :
        raise RuntimeError("I NEED A FUCKING FILE NAME!!!!!")

    main(sys.argv[1])
    plt.show()
