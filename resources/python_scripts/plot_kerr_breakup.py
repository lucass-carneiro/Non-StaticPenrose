import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import yaml
import sys

# MPL settings.
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'Latin Modern Roman'
plt.rcParams['figure.figsize'] = [9, 9]

color_1 = "red"
color_2 = "black"
color_3 = "blue"

font_size = 18

line_width = 2.5

# Datafile vars, for pandas
vars = [
    "time",
    "V1",
    "V2",
    "V3",
    "X1",
    "X2",
    "X3",
    "El",
    "Eg",
    "residual_Eg"
]

# Config and trajectory files
config_file_name = sys.argv[1]
trajectory_1 = sys.argv[2]
trajectory_2 = sys.argv[3]
trajectory_3 = sys.argv[4]

# Config data
with open(config_file_name, "r") as file:
    config_file = yaml.safe_load(file)

plot_radius = float(config_file["background_radius"])

# Trajectory data
data_1 = pd.read_csv(trajectory_1, delim_whitespace=True, names=vars)
data_2 = pd.read_csv(trajectory_2, delim_whitespace=True, names=vars)
data_3 = pd.read_csv(trajectory_3, delim_whitespace=True, names=vars)

# Trajectory plot
plt.plot(data_1["X1"], data_1["X2"], color=color_1, linewidth=line_width)
plt.plot(data_2["X1"], data_2["X2"], color=color_2, linewidth=line_width)
plt.plot(data_3["X1"], data_3["X2"], color=color_3, linewidth=line_width)

# Event Horizon data
M = float(config_file["KerrSchild_Kerr_Settings"]["M"])
a = float(config_file["KerrSchild_Kerr_Settings"]["a"])

if a*a > M*M:
    warn("Naked singularity detected. Drawing a unit radius horizon.")
    bh_radius = 1
else:
    bh_radius = M + np.sqrt(M**2 - a**2)

region_calc_radius = 2.0 * bh_radius

x = np.arange(-region_calc_radius, region_calc_radius,
              2.0 * region_calc_radius / 100.0)
y = x
X, Y = np.meshgrid(x, y)

rho2ma2 = X**2 + Y**2 - a**2

R = np.sqrt((rho2ma2 + np.abs(rho2ma2))/2)
gtt = 2.0 * M/R - 1.0

# Draw ergosphere
plt.contour(X, Y, gtt, [0.0], colors="black", linewidths=line_width)
plt.contourf(X, Y, gtt, [0.0, 1.0], colors=[
             (209.0/255.0, 231.0/255.0, 207.0/255.0)])

# Draw Event horizon
plt.contour(X, Y, R, [bh_radius], colors="black", linewidths=line_width)
plt.contourf(X, Y, R, [0.0, bh_radius], colors="white")

# Axes
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)

plt.xlabel("$x$", fontsize=font_size)
plt.ylabel("$y$", fontsize=font_size)

# Plot range
plt.xlim([-plot_radius, plot_radius])
plt.ylim([-plot_radius, plot_radius])

# Display/save
plt.tight_layout()
plt.show()
