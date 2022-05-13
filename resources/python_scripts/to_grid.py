import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Courtesy of github.com/iaraota
def to_grid():
    columns=["x", "y", "z", "data"]
    data = np.genfromtxt("/home/lucas/GRLensing/Release/Superimposed Kerr binary in Kerr-Schild coordinates_lapse_0.ascii")
    df = pd.DataFrame(data, columns=columns)

    # Pick a z value to plot
    z_filtered_df = df[df["z"] == 0.0]

    plt.tricontourf(z_filtered_df["x"], z_filtered_df["y"], z_filtered_df["data"], 100, vmin=0, vmax=1)
    plt.colorbar()
    plt.show()

to_grid()
