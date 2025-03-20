import json
from util.processCurve import *
from util.saveFile import saveFile

############################################################################################################

def mainl(path : str):
    discretized_points = []

    with(open(path, 'r')) as f:
        data = json.load(f)


    point1 = [data['x_ref'][0], data['y_ref'][0]]
    point2 = [data['x_ref'][1], data['y_ref'][1]]
    point3 = [data['x_ref'][2], data['y_ref'][2]]

    points = geometricalLine(point1, point2, 55, r=1.05)
    for p in points:
        discretized_points.append(p.tolist())

    points = geometricalLine(point3, point2, 55, r=1.05)
    for p in points:
        discretized_points.append(p.tolist())

    # Remove duplicate points
    seen = set()
    discretized_points = [point for point in discretized_points if tuple(point) not in seen and not seen.add(tuple(point))]
    
    return discretized_points

############################################################################################################

saveFile("Rampe", mainl)