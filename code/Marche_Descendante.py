import json
from util.processCurve import *
from util.saveFile import saveFile

############################################################################################################

def mainl(path : str):

    discretized_points = []

    with(open(path, 'r')) as f:
        data = json.load(f)

    N = data['N']

    for i in range(N):
        point1 = [data['x_ref'][i], data['y_ref'][i]]
        point2 = [data['x_ref'][i+1], data['y_ref'][i+1]]
        L = np.sqrt((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2)
        r = 1 + (1.1 - 1) * np.exp(-L)
        print(r)
        if i < N -1 and i > 0:
            # Discretizing the line
            pointM = [(point2[0] + point1[0])/2, (point2[1] + point1[1]) / 2]
            points1 = geometricalLine(pointM, point1, 25, r)
            points2 = geometricalLine(pointM, point2, 25, r)
            for p in points1:
                discretized_points.append(p.tolist())
            for p in points2:
                discretized_points.append(p.tolist())
        elif i < N - 1:
            points = geometricalLine(point1, point2, 50, r)
            for p in points:
                discretized_points.append(p.tolist())
        elif i > 0:
            points = geometricalLine(point2, point1, 50, r)
            for p in points:
                discretized_points.append(p.tolist())

    # Remove duplicate points
    seen = set()
    discretized_points = [point for point in discretized_points if tuple(point) not in seen and not seen.add(tuple(point))]
    
    return discretized_points

############################################################################################################

saveFile("Marche_Descendante", mainl)