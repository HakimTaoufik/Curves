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


    # Bezier curves
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
    discretized_points.append(points.tolist())

    # Bezier curves
    control_points = []
    point1 = [data['x_ref'][1], data['y_ref'][1]]
    point2 = [data['x_ref'][2], data['y_ref'][2]]

    control_points.append(point1)
    xref = curves[1]['x_ctr']
    yref = curves[1]['y_ctr']
    n = len(xref)
    for j in range(n):
        control_points.append([xref[j], yref[j]])
    control_points.append(point2)

    # Discretizing the Bezier curve
    points = discretize_bezier(t, control_points)
    discretized_points.append(points.tolist())

    flatten_points = flatten_points_as_lists(discretized_points)

    seen = set()
    flatten_points = [point for point in flatten_points if tuple(point) not in seen and not seen.add(tuple(point))]

    x,y = np.array(flatten_points)[:,0], np.array(flatten_points)[:,1]
    x_new, y_new = resample_by_curvature(x, y, num_points=100, max_dist=0.04, min_dist=0.02)
    discretized_points = np.array([x_new, y_new]).T.tolist()
    
    point1 = [data['x_ref'][2], data['y_ref'][2]]
    point2 = [data['x_ref'][3], data['y_ref'][3]]

    L = np.sqrt((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2)
    r = 1 + (1.1 - 1) * np.exp(-L)
    print(math.floor(math.sqrt(L*10)))
    pointM = [(point2[0] + point1[0])/2, (point2[1] + point1[1]) / 2]
    points1 = geometricalLine(pointM, point1, math.floor(math.sqrt(L*250)), 1.2)
    points2 = geometricalLine(pointM, point2, math.floor(math.sqrt(L*250)), 1.2)
    
    for p in points1:
        discretized_points.append(p.tolist())
    for p in points2:
        discretized_points.append(p.tolist())



    # Bezier curves
    control_points = []
    point1 = [data['x_ref'][3], data['y_ref'][3]]
    point2 = [data['x_ref'][4], data['y_ref'][4]]

    control_points.append(point1)
    xref = curves[3]['x_ctr']
    yref = curves[3]['y_ctr']
    n = len(xref)
    for j in range(n):
        control_points.append([xref[j], yref[j]])
    control_points.append(point2)

    # Discretizing the Bezier curve
    points = discretize_bezier(t, control_points)
    x = points
    x = x.tolist()

    a,y = np.array(x)[:,0], np.array(x)[:,1]
    x_new, y_new = resample_by_curvature(a, y, num_points=70, max_dist=0.04, min_dist=0.015)
    discretized_points.append(np.array([x_new, y_new]).T.tolist())

    # Line
    point1 = [data['x_ref'][4], data['y_ref'][4]]
    point2 = [data['x_ref'][5], data['y_ref'][5]]

    L = np.sqrt((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2)
    r = 1 + (1.2 - 1) * np.exp(-L)
    print(L, r)

    points = geometricalLine(point2, point1, math.floor(math.sqrt(L)*30), 1.1)
    for p in points:
        discretized_points.append(p.tolist())

    flatten_points = flatten_points_as_lists(discretized_points)
    # Remove duplicate points
    seen = set()
    discretized_points = [point for point in flatten_points if tuple(point) not in seen and not seen.add(tuple(point))]

    return discretized_points





############################################################################################################

saveFile("Radome_2", mainRadome)