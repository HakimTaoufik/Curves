import os
import yaml
import numpy as np
import matplotlib.pyplot as plt
from util.processCurve import *

def saveFile(family, main):
    def get_file_paths(folder):
        file_paths = []
        for file in os.listdir(folder):
            full_path = os.path.join(folder, file)
            if os.path.isfile(full_path):  # Ensure it's a file, not a folder
                file_paths.append(full_path)
        return file_paths

    # Main
    folder = "json/" + family
    all_file_paths = get_file_paths(folder)
    for file_path in all_file_paths:
        try:
            discretized_points = main(file_path)
            name = file_path.replace(f"json/{family}/", "").replace(".json", "")
            flattened = flatten_points_as_lists(discretized_points)
            if flattened[0] != [0,0]:
                flattened.insert(0, [0,0])
            x, y = np.array(flattened)[:,0], np.array(flattened)[:,1]

            plt.figure(figsize=(8, 8))  # Set square figure for consistent size
            plt.axis('equal')  # Make plot orthonormal (equal scaling on both axes)
            plt.plot(x, y, marker="o", linestyle="-", color="b", label="y = f(x)")  # Plot the curve
            plt.title("Plot of y vs x")  # Add a title
            plt.xlabel("x")  # Label the x-axis
            plt.ylabel("y")  # Label the y-axis
            plt.grid(True)  # Add a grid
            plt.legend()  # Show the legend

            plt.savefig(f"generated/imgs/{family}/plot_" + name + ".png")
            plt.close()

            points = [{"x":x, "y":y} for x, y in zip(x.tolist(), y.tolist())]

            with open(f"generated/data/{family}/points_" + name + ".yaml", "w") as file:
                yaml.dump(points, file, default_flow_style=False)
        except:
            continue