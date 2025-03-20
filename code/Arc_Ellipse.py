from util.processCurve import *
from util.saveFile import saveFile

############################################################################################################

def mainl(path : str):
    t = 1000
    discretized_points = []

    with(open(path, 'r')) as f:
        data = json.load(f)

    curves = data['curves']

    x_c = curves[0]['x_c']
    y_c = curves[0]['y_c']
    a = curves[0]['a']
    b = curves[0]['b']
    theta = curves[0]['theta']
    t1 = curves[0]['t1']
    t2 = curves[0]['t2']
    x, y = generate_ellipse_points(x_c, y_c, a, b, theta, t1, t2, num_points=1000)
    xn, yn = resample_by_curvature(x, y, num_points=80, max_dist=None, min_dist=0.02)
    discretized_points.append(np.array([xn, yn]).T.tolist())

    discretized_points = discretized_points[0]
    
    # Remove duplicate points
    seen = set()
    discretized_points = [point for point in discretized_points if tuple(point) not in seen and not seen.add(tuple(point))]
    
    return discretized_points

############################################################################################################

saveFile("Arc_Ellipse", mainl)