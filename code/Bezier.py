from util.processCurve import *
from util.saveFile import saveFile

############################################################################################################

def mainl(path : str):
    # Number of points discretized
    t = 1000
    discretized_points = []


    with(open(path, 'r')) as f:
        data = json.load(f)

    N = data['N']
    curves = data['curves']

    for i in range(N):
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

    flattened_points = flatten_points_as_lists(discretized_points)
    
    # Remove duplicate points
    seen = set()
    flattened_points = [point for point in flattened_points if tuple(point) not in seen and not seen.add(tuple(point))]

    x,y = np.array(flattened_points)[:,0], np.array(flattened_points)[:,1]
    x_new, y_new = resample_by_curvature(x, y, num_points=100, max_dist=0.05, min_dist=0.01)
    discretized_points = np.array([x_new, y_new]).T.tolist()

    return discretized_points

############################################################################################################

saveFile("Bezier", mainl)