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
                # Calculate the length of the curve represented by left_side
                Ll = 0
                for i in range(len(left_side) - 1):
                    p1 = left_side[i]
                    p2 = left_side[i+1]
                    Ll += np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)
                Lr = 0
                for i in range(len(right_side) - 1):
                    p1 = right_side[i]
                    p2 = right_side[i+1]
                    Lr += np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)
                LT += Ll + Lr
                rl = 1.05 - 0.1 * np.exp(-Ll**(1/10))
                rr = 1.05 - 0.1 * np.exp(-Lr**(1/10))
                npointsL = np.array(geometrically_discretize_curve(left_side, math.floor(math.sqrt(Ll*800)), rl))
                npointsR = np.array(geometrically_discretize_curve(right_side, math.floor(math.sqrt(Lr*800)), rr))
                discretized_points.append(npointsL.tolist()[::-1])
                discretized_points.append(npointsR.tolist())
            elif types[i] == 'droite':
                # Discretizing the line
                point1 = [data['x_ref'][i], data['y_ref'][i]]
                point2 = [data['x_ref'][i+1], data['y_ref'][i+1]]
                pointM = [(point2[0] + point1[0])/2, (point2[1] + point1[1]) / 2]
                L = np.sqrt((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2)
                r = 1.02 - 0.01 * np.exp(-L**(1/10))
                points1 = geometricalLine(pointM, point1, math.floor(np.exp(-L) * 20), r)
                points2 = geometricalLine(pointM, point2, math.floor(np.exp(-L) * 20), r)
                for p in points1:
                    discretized_points.append(p.tolist())
                for p in points2:
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

    return discretized_points

############################################################################################################

saveFile("Antenna_4", mainRadome)