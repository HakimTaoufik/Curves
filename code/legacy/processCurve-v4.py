import numpy as np
import os
import yaml
import json
import matplotlib.pyplot as plt
from scipy.special import comb
import math
from scipy.interpolate import interp1d

def generate_ellipse_points(x_c, y_c, a, b, theta, t1, t2, num_points=100):
    """
    Generate points along an elliptical arc.
    
    Parameters:
        x_c, y_c : float
            The center coordinates of the ellipse.
        a, b : float
            The semi-axis lengths. (Note: Negative values will reverse the parameter direction.)
        theta : float
            The rotation angle of the ellipse in radians.
        t1, t2 : float
            The starting and ending parameter values.
        num_points : int, optional
            The number of points to generate along the arc (default is 100).
    
    Returns:
        x_points, y_points : ndarray
            Arrays of x and y coordinates of the ellipse arc.
    """
    # Generate parameter values between t1 and t2
    t_values = np.linspace(t1, t2, num_points)
    
    # Precompute cosine and sine of the rotation angle for efficiency
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    
    # Compute the x and y coordinates using the parametric equations of the ellipse
    x_points = x_c + a * np.cos(t_values) * cos_theta - b * np.sin(t_values) * sin_theta
    y_points = y_c + a * np.cos(t_values) * sin_theta + b * np.sin(t_values) * cos_theta
    
    return x_points, y_points

def discretize_bezier(N, control_points):
    # Compute the Bezier curve for a given set of control points.
    control_points = np.array(control_points)
    n = len(control_points) - 1

    def bernstein_poly(i, n, t):
        return comb(n, i) * (t ** i) * ((1 - t) ** (n - i))

    t = np.linspace(0, 1, N)
    curve = np.zeros((N, control_points.shape[1]))

    for i in range(n + 1):
        curve += np.outer(bernstein_poly(i, n, t), control_points[i])

    return curve

def discretize_line(point1, point2, N):
    # Discretize a line given 2 points and N number of points.
    point1 = np.array(point1)
    point2 = np.array(point2)
    return np.linspace(point1, point2, N)


def relative_line_size(point1, point2):
    """
    Calculate the relative size of the line between two points,
    scaled between 1.05 and 1.2.
    """
    x1, y1 = point1
    x2, y2 = point2

    # Calculate Euclidean distance
    distance = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    # Normalize the distance to fit within the range [1.05, 1.2]
    min_dist, max_dist = 0.05, 0.9  # Define a reasonable distance range
    min_scale, max_scale = 1.05, 1.2
    
    # Scale the distance proportionally within [1.05, 1.2]
    scaled_size = min_scale + (max_scale - min_scale) * ((distance - min_dist) / (max_dist - min_dist))
    
    # Ensure it stays within bounds
    return max(min_scale, min(scaled_size, max_scale))


def geometricalLine(A, B, N, r):
    """
    Discretizes the line segment from A to B into N points with a geometric distribution clustered near B.
    
    Parameters:
    A (array-like): First point coordinates.
    B (array-like): Last point coordinates.
    N (int): Number of points.
    r (float): Ratio of geometric progression (>1).
    
    Returns:
    list: List of points along the line.
    """
    A = np.array(A)
    B = np.array(B)
    delta = B - A
    
    # Compute sum of geometric series
    if r == 1:
        S = N - 1
    else:
        S = (r**(N-1) - 1) / (r - 1)
    d = 1.0 / S
    
    # Generate t values
    t_values = [1.0]
    current_t = 1.0
    current_d = d
    for _ in range(N-2):
        current_t -= current_d
        t_values.append(current_t)
        current_d *= r
    # The last step should reach 0, but due to precision, we enforce it
    t_values.append(0.0)
    t_values = t_values[::-1]  # Reverse to go from 0 to 1
    
    # Compute points
    points = [A + t * delta for t in t_values]
    return points

def flatten_points_as_lists(nested_list):
    # Flattens a deeply nested list of lists into a single list of points, 
    # keeping points as lists.
    flat_points = []
    for item in nested_list:
        if isinstance(item, list):
            # If the item is a list, recursively process it
            if len(item) == 2 and all(isinstance(coord, (int, float)) for coord in item):
                # Append directly if it's a point-like structure
                flat_points.append(item)
            else:
                # Otherwise, recursively flatten the sublist
                flat_points.extend(flatten_points_as_lists(item))
    return flat_points

def interp_s(x, y, s_old, s_target):
    '''Interpolate x and y at s_target based on s_old'''
    i = 0
    while i < len(s_old) - 1 and not (s_old[i] <= s_target <= s_old[i+1]):
        i += 1
    if i == len(s_old) - 1:
        return x[-1], y[-1]
    seg_length = s_old[i+1] - s_old[i]
    if seg_length < 1e-14:
        return x[i], y[i]
    t = (s_target - s_old[i]) / seg_length
    xx = x[i] + t * (x[i+1] - x[i])
    yy = y[i] + t * (y[i+1] - y[i])
    return xx, yy


def curvature(points):
    ds_list = []
    x, y = np.array(points)[:, 0], np.array(points)[:, 1]
    for i in range(len(x)-1):
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
        ds_list.append(math.sqrt(dx*dx + dy*dy))

    s = [0.0]
    for dist in ds_list:
        s.append(s[-1] + dist)

    L = s[-1]

    return s, L

def findN(L, alpha, beta):
    # Finding the number of points N for the geometrical law
    a = 1 - L/alpha * (1 - beta)
    if a <= 0 or beta == 1 :
        raise ValueError("Invalid values for alpha, beta and L")
    return round(1 + math.log(a) / math.log(beta))

def harmonize_points(points, delta, epsilon=0.0005, max_points=1000):
    '''Uniformly sample points along the curve'''
    s = curvature(points)[0]
    N = max_points


    s_new = [0 for i in range(N)]
    for j in range(1,N):
        s_new[j] = s_new[j-1] + delta

    x, y = np.array(points)[:, 0], np.array(points)[:, 1]
    npoints = []

    for val in s_new:
        xx, yy = interp_s(x, y, s, val)
        npoints.append([xx, yy])

    return npoints

def curvatureBased(x, y, num_points):
    x = np.asarray(x)
    y = np.asarray(y)
    
    # Compute arc length
    dx = np.diff(x)
    dy = np.diff(y)
    ds = np.sqrt(dx**2 + dy**2)
    s = np.zeros(len(x))
    s[1:] = np.cumsum(ds)
    
    # Calculate curvature at each point
    dx_dt = np.gradient(x)
    dy_dt = np.gradient(y)
    d2x_dt2 = np.gradient(dx_dt)
    d2y_dt2 = np.gradient(dy_dt)
    curvature = np.abs(dx_dt * d2y_dt2 - dy_dt * d2x_dt2) / (dx_dt**2 + dy_dt**2)**1.5
    curvature = np.nan_to_num(curvature, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Compute weights for each segment (average curvature * segment length)
    avg_curvature = (curvature[:-1] + curvature[1:]) / 2
    weights = avg_curvature * ds
    total_weight = np.sum(weights)
    
    if total_weight <= 0:
        # Fallback to uniform sampling if all weights are zero
        s_targets = np.linspace(0, s[-1], num_points)
    else:
        # Compute cumulative weights
        cumulative_weights = np.zeros(len(s))
        cumulative_weights[1:] = np.cumsum(weights)
        
        # Generate target weights
        target_weights = np.linspace(0, total_weight, num_points)
        
        # Find corresponding segment indices and fractions
        indices = np.searchsorted(cumulative_weights, target_weights) - 1
        indices = np.clip(indices, 0, len(weights) - 1)
        
        # Calculate fraction within each segment
        valid = (cumulative_weights[indices + 1] > cumulative_weights[indices])
        remaining = target_weights - cumulative_weights[indices]
        segment_weights = cumulative_weights[indices + 1] - cumulative_weights[indices]
        fraction = np.divide(remaining, segment_weights, where=valid, out=np.zeros_like(remaining))
        fraction = np.clip(fraction, 0, 1)
        
        # Compute target arc lengths
        s_targets = s[indices] + fraction * (s[indices + 1] - s[indices])
    
    # Interpolate x and y at the target arc lengths
    f_x = interp1d(s, x, kind='linear', fill_value='extrapolate')
    f_y = interp1d(s, y, kind='linear', fill_value='extrapolate')
    x_resampled = f_x(s_targets)
    y_resampled = f_y(s_targets)
    
    return x_resampled, y_resampled


def resample_by_curvature(x, y, num_points, max_dist=None, min_dist=None):
    x = np.asarray(x)
    y = np.asarray(y)
    
    # Calculate arc length
    dx = np.diff(x)
    dy = np.diff(y)
    ds = np.sqrt(dx**2 + dy**2)
    s = np.zeros(len(x))
    s[1:] = np.cumsum(ds)
    
    # Compute curvature
    dx_dt = np.gradient(x)
    dy_dt = np.gradient(y)
    d2x_dt2 = np.gradient(dx_dt)
    d2y_dt2 = np.gradient(dy_dt)
    curvature = np.abs(dx_dt * d2y_dt2 - dy_dt * d2x_dt2) / (dx_dt**2 + dy_dt**2)**1.5
    curvature = np.nan_to_num(curvature, 0) #replace NaN with 0
    
    # Calculate segment weights
    avg_curvature = (curvature[:-1] + curvature[1:]) / 2
    weights = avg_curvature * ds
    total_weight = np.sum(weights)
    
    # Generate target arc lengths
    if total_weight <= 0:
        s_targets = np.linspace(0, s[-1], num_points)
    else:
        cumulative_weights = np.insert(np.cumsum(weights), 0, 0)
        target_weights = np.linspace(0, total_weight, num_points)
        indices = np.clip(np.searchsorted(cumulative_weights, target_weights) - 1, 0, len(weights)-1)
        t_frac = (target_weights - cumulative_weights[indices]) / (cumulative_weights[indices+1] - cumulative_weights[indices])
        s_targets = s[indices] + t_frac * (s[indices+1] - s[indices])
    
    # Apply maximum and minimum distance constraints
    if (max_dist is not None and max_dist > 0) or (min_dist is not None and min_dist > 0):
        s_processed = [s_targets[0]]
        for current_s in s_targets[1:]:
            prev_s = s_processed[-1]
            delta_s = current_s - prev_s
            
            # Enforce maximum distance
            if max_dist is not None and delta_s > max_dist:
                n_segments = int(np.ceil(delta_s / max_dist))
                new_s = np.linspace(prev_s, current_s, n_segments + 1)[1:]
                s_processed.extend(new_s)
            
            # Enforce minimum distance
            elif min_dist is not None and delta_s < min_dist:
                # Skip this point if it's too close to the previous one
                continue
            
            else:
                s_processed.append(current_s)
        s_targets = np.array(s_processed)
    
    
    # Ensure s_targets are within the valid range of the original curve
    s_targets = np.clip(s_targets, 0, s[-1])
    
    # Interpolate new points
    f_x = interp1d(s, x, kind='linear', fill_value='extrapolate')
    f_y = interp1d(s, y, kind='linear', fill_value='extrapolate')
    return f_x(s_targets), f_y(s_targets)

# Geometrical law
def geometrical(points, alpha, beta):
    s = curvature(points)[0]
    L = curvature(points)[1]
    N = findN(L, alpha, beta)

    s_new = [0 for i in range(N)]
    for j in range(1,N):
        s_new[j] = s_new[j-1] + alpha*(beta**(N-j-1))

    x, y = np.array(points)[:, 0], np.array(points)[:, 1]
    npoints = []

    for val in s_new:
        xx, yy = interp_s(x, y, s, val)
        npoints.append([xx, yy])

    return npoints

def harmonize(points, delta, epsilon=0.0005, max_points=1000):
    # Uniformly sample points along the curve
    s = curvature(points)[0]
    N = max_points


    s_new = [0 for i in range(N)]
    for j in range(1,N):
        s_new[j] = s_new[j-1] + delta

    x, y = np.array(points)[:, 0], np.array(points)[:, 1]
    npoints = []

    for val in s_new:
        xx, yy = interp_s(x, y, s, val)
        npoints.append([xx, yy])

    return npoints
import numpy as np

def geometrically_discretize_curve(points, num_points, ratio):
    # Compute cumulative arc lengths
    arc_lengths = [0.0]
    for i in range(1, len(points)):
        dx = points[i][0] - points[i-1][0]
        dy = points[i][1] - points[i-1][1]
        arc_lengths.append(arc_lengths[-1] + np.hypot(dx, dy))
    
    total_length = arc_lengths[-1]
    if total_length == 0:
        return points  # All points are the same
    
    # Generate geometrically spaced t values
    M = num_points
    r = ratio
    if r <= 1:
        raise ValueError("Ratio must be greater than 1.")
    
    delta_1 = (r - 1) / (r ** (M - 1) - 1)
    t_values = [0.0]
    current_t = 0.0
    for i in range(M - 2, -1, -1):
        current_t += delta_1 * (r ** i)
        t_values.append(current_t)
    t_values[-1] = 1.0  # Ensure last point is exactly 1 due to floating point issues
    
    # Convert t_values to arc lengths
    s_values = [t * total_length for t in t_values]
    
    # Interpolate points
    new_points = []
    for s in s_values:
        # Find the segment containing s
        idx = 0
        while idx < len(arc_lengths) - 1 and arc_lengths[idx + 1] <= s:
            idx += 1
        if idx == len(arc_lengths) - 1:
            new_points.append(points[-1])
            continue
        # Interpolate between points[idx] and points[idx + 1]
        seg_start = arc_lengths[idx]
        seg_end = arc_lengths[idx + 1]
        seg_length = seg_end - seg_start
        if seg_length == 0:
            frac = 0.0
        else:
            frac = (s - seg_start) / seg_length
        x = points[idx][0] + frac * (points[idx + 1][0] - points[idx][0])
        y = points[idx][1] + frac * (points[idx + 1][1] - points[idx][1])
        new_points.append((x, y))
    
    return new_points


def main(path : str):
    # Number of points discretized
    t = 1000
    discretized_points = []

    with(open(path, 'r')) as f:
        data = json.load(f)

    name = file_path.replace("json/", "").replace(".json", "")
    N = data['N']
    curves = data['curves']
    types = [curves[i]['type'] for i in range(len(curves))]
    
    if "Antenna_1" not in name:
        for i in range(N):
            aux = []
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
                #aux_points.append(points.tolist())
                flatAuxPoints = flatten_points_as_lists(points.tolist())
                x,y = np.array(flatAuxPoints)[:,0], np.array(flatAuxPoints)[:,1]
                xn, yn = resample_by_curvature(x, y, num_points=100, max_dist=0.05, min_dist=0.02)
                discretized_points.append(np.array([xn, yn]).T.tolist())

            elif types[i] == 'droite':
                if i < N -1 and i > 0 and types[i-1] == 'droite' and types[i+1] == 'droite':
                    # Discretizing the line
                    point1 = [data['x_ref'][i], data['y_ref'][i]]
                    point2 = [data['x_ref'][i+1], data['y_ref'][i+1]]

                    pointM = [(point2[0] + point1[0])/2, (point2[1] + point1[1]) / 2]
                    points1 = geometricalLine(pointM, point1, 25, r = 1.15)
                    points2 = geometricalLine(pointM, point2, 25, r = 1.15)
                    for p in points1:
                        discretized_points.append(p.tolist())
                    for p in points2:
                        discretized_points.append(p.tolist())
                    aux = []
                elif i < N - 1 and types[i+1] == 'droite':
                    # Discretizing the line
                    point1 = [data['x_ref'][i], data['y_ref'][i]]
                    point2 = [data['x_ref'][i+1], data['y_ref'][i+1]]

                    points = geometricalLine(point1, point2, 50, r=1.15)
                    for p in points:
                        discretized_points.append(p.tolist())
                    aux = []
                elif i > 0 and types[i-1] == 'droite':
                    point1 = [data['x_ref'][i], data['y_ref'][i]]
                    point2 = [data['x_ref'][i+1], data['y_ref'][i+1]]

                    points = geometricalLine(point2, point1, 50, r=1.15)
                    for p in points:
                        discretized_points.append(p.tolist())
                    aux = []
                else:
                    # Discretizing the line
                    point1 = [data['x_ref'][i], data['y_ref'][i]]
                    point2 = [data['x_ref'][i+1], data['y_ref'][i+1]]

                    points = discretize_line(point1, point2, t)
                    aux.append(points.tolist())
                    flatAuxPoints = flatten_points_as_lists(aux)
                    newPoints = harmonize(flatAuxPoints, 0.02)
                    discretized_points.append(newPoints)
                    aux = []
            elif types[i] == "ellipse":
                x_c = curves[i]['x_c']
                y_c = curves[i]['y_c']
                a = curves[i]['a']
                b = curves[i]['b']
                theta = curves[i]['theta']
                t1 = curves[i]['t1']
                t2 = curves[i]['t2']
                x, y = generate_ellipse_points(x_c, y_c, a, b, theta, t1, t2, num_points=1000)
                xn, yn = resample_by_curvature(x, y, num_points=70, max_dist=None, min_dist=None)
                discretized_points.append(np.array([xn, yn]).T.tolist())
    else:
        for i in range(N):
            aux = []
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
                #aux_points.append(points.tolist())
                flatAuxPoints = flatten_points_as_lists(points.tolist())
                mid_index = len(points) // 2
                left_side = points[:mid_index]
                left_side = left_side[::-1]
                right_side = points[mid_index:]
                npointsL = np.array(geometrically_discretize_curve(left_side, 25, 1.1))
                npointsR = np.array(geometrically_discretize_curve(right_side, 25, 1.1))
                discretized_points.append(npointsL.tolist()[::-1])
                discretized_points.append(npointsR.tolist())
            elif types[i] == 'droite':
                if i != 3:
                    # Discretizing the line
                    point1 = [data['x_ref'][i], data['y_ref'][i]]
                    point2 = [data['x_ref'][i+1], data['y_ref'][i+1]]
                    pointM = [(point2[0] + point1[0])/2, (point2[1] + point1[1]) / 2]
                    points1 = geometricalLine(pointM, point1, 5, r = 1.3)
                    points2 = geometricalLine(pointM, point2, 5, r = 1.3)
                    for p in points1:
                        discretized_points.append(p.tolist())
                    for p in points2:
                        discretized_points.append(p.tolist())
                else:
                    # Discretizing the line
                    point1 = [data['x_ref'][i], data['y_ref'][i]]
                    point2 = [data['x_ref'][i+1], data['y_ref'][i+1]]
                    points = discretize_line(point1, point2, t)
                    aux.append(points.tolist())
                    flatAuxPoints = flatten_points_as_lists(aux)
                    newPoints = harmonize(flatAuxPoints, 0.02)
                    discretized_points.append(newPoints)
                    aux = []

    return discretized_points

# Main Code

def mainRadone(path : str):
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

    flat = flatten_points_as_lists(discretized_points)
    x, y = np.array(flat)[:,0], np.array(flat)[:,1]
    xn, yn = resample_by_curvature(x, y, num_points=1000, max_dist=0.05, min_dist=0.02)
    xn = xn.tolist()
    yn = yn.tolist()
    discretized_points = np.array([xn, yn]).T.tolist()
    return discretized_points

############################################################################################################

def get_file_paths(folder):
    file_paths = []
    for file in os.listdir(folder):
        full_path = os.path.join(folder, file)
        if os.path.isfile(full_path):  # Ensure it's a file, not a folder
            file_paths.append(full_path)
    return file_paths

# Example usage
folder = "json"
all_file_paths = get_file_paths(folder)
for file_path in all_file_paths:
    try:
        if "Radome_2" in file_path:
            discretized_points = mainRadone(file_path)
        else:
            discretized_points = main(file_path)
        name = file_path.replace("json/", "").replace(".json", "")
        flattened = flatten_points_as_lists(discretized_points)
        if flattened[0] != [0,0]:
            flattened.insert(0, [0,0])
        x, y = np.array(flattened)[:,0], np.array(flattened)[:,1]

        plt.figure(figsize=(8, 6))  # Set the figure size
        plt.plot(x, y, marker="o", linestyle="-", color="b", label="y = f(x)")  # Plot the curve
        plt.title("Plot of y vs x")  # Add a title
        plt.xlabel("x")  # Label the x-axis
        plt.ylabel("y")  # Label the y-axis
        plt.grid(True)  # Add a grid
        plt.legend()  # Show the legend

        plt.savefig("generated/imgs/plot_" + name + ".png")

        points = [{"x":x, "y":y} for x, y in zip(x.tolist(), y.tolist())]

        with open("generated/data/points_" + name + ".yaml", "w") as file:
            yaml.dump(points, file, default_flow_style=False)
    except:
        continue