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

    for i in range(N):
        if types[i] == 'bezier':
            # Extracting control points
            control_points = []
            point1 = [data['x_ref'][i], data['y_ref'][i]]
            point2 = [data['x_ref'][i+1], data['y_ref'][i+1]]

            control_points.append(point1)
            xref = curves[i]['x_ctr']
            yref = curves[i]['y_ctr']
            n = len(xref)
            for j in range(n):
                control_points.append([xref[j], yref[j]])
            control_points.append(point2)
            
            # Discretizing the Bezier curve
            points = discretize_bezier(t, control_points)
            discretized_points.append(points.tolist())

        elif types[i] == 'droite':
            # Discretizing the line
            point1 = [data['x_ref'][i], data['y_ref'][i]]
            point2 = [data['x_ref'][i+1], data['y_ref'][i+1]]
            points = discretize_line(point1, point2, t)
            discretized_points.append(points.tolist())
        
        elif types[i] == 'cercle':
            # Discretizing the circle
            z1 = (data['x_ref'][i], data['y_ref'][i])
            z2 = (data['x_ref'][i+1], data['y_ref'][i+1])
            trigo = curves[i]['trigo']
            c = (curves[i]['c'][0], curves[i]['c'][1])
            points = generate_half_circle_points(z1, z2, c, trigo)
            discretized_points.append(points)

    flat = flatten_points_as_lists(discretized_points)
    x, y = np.array(flat)[:,0], np.array(flat)[:,1]
    xn, yn = resample_by_curvature(x, y, num_points=1000, max_dist=0.05, min_dist=0.02)
    xn = xn.tolist()
    yn = yn.tolist()
    discretized_points = np.array([xn, yn]).T.tolist()
    return discretized_points

############################################################################################################

saveFile("EOS_1", mainRadome)