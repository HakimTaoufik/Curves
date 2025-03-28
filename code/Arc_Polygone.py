import json
import os
from util.processCurve import *
from util.saveFile import saveFile

############################################################################################################

def mainl(path : str):

    discretized_points = []

    with(open(path, 'r')) as f:
        data = json.load(f)

    N = data['N']

    if "3dr" in path:
        deg = 3
    elif "4dr" in path:
        deg = 4
    elif "5dr" in path:
        deg = 5
    elif "6dr" in path:
        deg = 6

    print(deg, path)
    a = 15

    for i in range(N):
        point1 = [data['x_ref'][i], data['y_ref'][i]]
        point2 = [data['x_ref'][i+1], data['y_ref'][i+1]]
        L = np.sqrt((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2)
        r = 1.1
        if i < N -1 and i > 0:
            # Discretizing the line
            pointM = [(point2[0] + point1[0])/2, (point2[1] + point1[1]) / 2]
            points1 = geometricalLine(pointM, point1, a, r)
            points2 = geometricalLine(pointM, point2, a, r)
            for p in points1:
                discretized_points.append(p.tolist())
            for p in points2:
                discretized_points.append(p.tolist())
        elif i < N - 1:
            points = geometricalLine(point1, point2, 2*a, r)
            for p in points:
                discretized_points.append(p.tolist())
        elif i > 0:
            points = geometricalLine(point2, point1, 2*a, r)
            for p in points:
                discretized_points.append(p.tolist())

    # Remove duplicate points
    seen = set()
    discretized_points = [point for point in discretized_points if tuple(point) not in seen and not seen.add(tuple(point))]
    
    return discretized_points

############################################################################################################

saveFile("Arc_Polygone", mainl)