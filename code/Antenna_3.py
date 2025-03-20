from util.processCurve import *
from util.saveFile import saveFile

############################################################################################################

def mainRadome(path : str):
    # Number of points discretized
    t = 1000
    discretized_points = []


    with(open(path, 'r')) as f:
        data = json.load(f)

    N = data['N']
    curves = data['curves']
    types = [curves[i]['type'] for i in range(len(curves))]
    LT = 0

    # Extracting control points
    control_points = []
    point1 = [data['x_ref'][0], data['y_ref'][0]]
    point2 = [data['x_ref'][1], data['y_ref'][1]]

    control_points.append(point1)
    xref = curves[0]['x_ctr']
    yref = curves[0]['y_ctr']
    n = len(xref)
    for j in range(n):
        control_points.append([xref[j], yref[j]])
    control_points.append(point2)
    
    # Discretizing the Bezier curve
    points = discretize_bezier(t, control_points)
    #aux_points.append(points.tolist())
    flatAuxPoints = flatten_points_as_lists(points.tolist())
    npoints = np.array(geometrically_discretize_curve(points, 50, 1.05))
    discretized_points.append(npoints.tolist())

    # Discretizing the line
    point1 = [data['x_ref'][1], data['y_ref'][1]]
    point2 = [data['x_ref'][2], data['y_ref'][2]]
    pointM = [(point2[0] + point1[0])/2, (point2[1] + point1[1]) / 2]
    L = np.sqrt((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2)
    r = 1.02 - 0.01 * np.exp(-L**(1/10))
    points1 = geometricalLine(pointM, point1, 20, 1.1)
    points2 = geometricalLine(pointM, point2, 20, 1.1)
    for p in points1:
        discretized_points.append(p.tolist())
    for p in points2[::-1]:
        discretized_points.append(p.tolist())
    # Uncomment the following block if you want to handle additional line discretization cases
    # point1 = [data['x_ref'][i], data['y_ref'][i]]
    # point2 = [data['x_ref'][i+1], data['y_ref'][i+1]]
    # L = np.sqrt((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2)
    # r = (1.05 - 1) * np.exp(-L)
    # points = discretize_line(point1, point2, t)
    # aux.append(points.tolist())
    # flatAuxPoints = flatten_points_as_lists(aux)
    # newPoints = harmonize(flatAuxPoints, r)
    # discretized_points.append(newPoints)
    # aux = []

    # Extracting control points
    control_points = []
    point1 = [data['x_ref'][2], data['y_ref'][2]]
    point2 = [data['x_ref'][3], data['y_ref'][3]]

    control_points.append(point1)
    xref = curves[2]['x_ctr']
    yref = curves[2]['y_ctr']
    n = len(xref)
    for j in range(n):
        control_points.append([xref[j], yref[j]])
    control_points.append(point2)
    
    # Discretizing the Bezier curve
    points = discretize_bezier(t, control_points)
    #aux_points.append(points.tolist())
    flatAuxPoints = flatten_points_as_lists(points.tolist())
    npoints = np.array(geometrically_discretize_curve(points[::-1], 50, 1.05))
    discretized_points.append(npoints.tolist()[::-1])
    discretized_points = flatten_points_as_lists(discretized_points)    # Removing duplicates
    seen = set()
    discretized_points = [point for point in discretized_points if tuple(point) not in seen and not seen.add(tuple(point))]

    return discretized_points

############################################################################################################

saveFile("Antenna_3", mainRadome)