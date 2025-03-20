import numpy as np
import os
import yaml
import json
import matplotlib.pyplot as plt
from scipy.special import comb
import math
from scipy.interpolate import interp1d

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

# Main Code
def exctractFromPath(path : str):
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

    return discretized_points

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
    # Interpolate x and y at s_target based on s_old
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

def harmonize_points(points, delta, epsilon=0.0005, max_points=1000):
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

def resample_by_curvature(x, y, num_points):
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




def sample_reg(X_ctr, Y_ctr ,N=1000):
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

    points = np.array([X_ctr, Y_ctr]).T
    s = curvature(points)[0]
    L = curvature(points)[1]

    s_new = [0 for i in range(N)]
    for j in range(1,N):
        s_new[j] = s_new[j-1] + L/N

    x, y = np.array(points)[:, 0], np.array(points)[:, 1]
    npoints = []

    for val in s_new:
        xx, yy = interp_s(x, y, s, val)
        npoints.append([xx, yy])

    return npoints



############################################################################################################





def get_file_paths(folder):
    file_paths = []
    for file in os.listdir(folder):
        full_path = os.path.join(folder, file)
        if os.path.isfile(full_path):  # Ensure it's a file, not a folder
            file_paths.append(full_path)
    return file_paths

# Example usage
folder = "tests/radone1"
all_file_paths = get_file_paths(folder)
i = 0
for file_path in all_file_paths:
    discretized_points = exctractFromPath(file_path)
    flattened_points = flatten_points_as_lists(discretized_points)
    x,y = np.array(flattened_points)[:,0], np.array(flattened_points)[:,1]
    x_new, y_new = resample_by_curvature(x, y, num_points=1000,max_dist=0.05, min_dist=0.02)

    plt.figure(figsize=(8, 6))  # Set the figure size
    plt.plot(x_new, y_new, marker="o", linestyle="-", color="b", label="y = f(x)")  # Plot the curve
    plt.title("Plot of y vs x")  # Add a title
    plt.xlabel("x")  # Label the x-axis
    plt.ylabel("y")  # Label the y-axis
    plt.grid(True)  # Add a grid
    plt.legend()  # Show the legend

    plt.savefig("generated/imgs/radone1/plot_" + str(i) + ".png")

    points = [{"x":x, "y":y} for x, y in zip(x_new.tolist(), y_new.tolist())]

    with open("generated/data/radone1/points_" + str(i) +".yaml", "w") as file:
        yaml.dump(points, file, default_flow_style=False)

    i += 1