import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

file_name = "trajectory.ascii"
M = 1

font_size = 20

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'Latin Modern Roman'
plt.rcParams['figure.figsize'] = [9, 9]
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size

def lapse(x, y, z):
    return 1 - (2 * M) / (M + 2 * np.sqrt(x**2 + y**2 + z**2))

def global_energy(local_energy, x, y, z):
    return lapse(x, y, z) * local_energy

vars = [
    "time",
    "V1",
    "V2",
    "V3",
    "X1",
    "X2",
    "X3",
    "EN"
]

data = pd.read_csv(file_name, delim_whitespace=True, names=vars)

data["diff"] = data.apply(lambda row : row["EN"] - global_energy(row["EN"], row["X1"], row["X2"], row["X3"]), axis = 1)

plt.close("all")

plt.plot(data["time"], data["diff"])

plt.xlabel("$t$", fontsize = font_size)
plt.ylabel("Deviation to conserved", fontsize = font_size)

plt.show()