#!/bin/python

"""GRLensing data plotter.

Usage:
  data_plotter.py trajectory <trajectory_config_file> <trajectory_output_file> [--font_size=<size>] [--delta=<value>] [--plot_radius=<radius>] [--color=<color>]
  data_plotter.py energy (local|global|residual) <trajectory_output_file> [--font_size=<size>]
  data_plotter.py penrose <penrose_config_file> <trajectory_1> <trajectory_2> <trajectory_3> [--font_size=<size>] [--color_1=<color1>] [--color_2=<color2>] [--color_3=<color3>] [--plot_radius=<radius>]
  data_plotter.py gridfunction <grid_function_data_file>
  data_plotter.py (-h | --help)
  data_plotter.py --version

Options:
  -h --help               Show this screen.
  --version               Show version.
  --font_size=<size>      The size of the font in the plots [default: 20].
  --delta=<value>         The grid spacing used in internal grid generation [default: 0.1].
  --plot_radius=<radius>  The radius of the background.
  --color=<color>         The color of the trajectory. [default: black]
  --color_1=<color1>      The color of the first breakup trajectory. [default: red]
  --color_2=<color2>      The color of the first breakup trajectory. [default: black]
  --color_3=<color3>      The color of the first breakup trajectory. [default: blue]

"""

from docopt import docopt
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import yaml

# MPL settings.
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'Latin Modern Roman'
plt.rcParams['figure.figsize'] = [9, 9]

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

def plot_horizons(plt, ax, arguments, config_file):
  background_radius = float(config_file["background_radius"])
  background_metric = config_file["background_metric"]

  if arguments["--plot_radius"] == None:
    plot_radius = background_radius
  else:
    plot_radius = float(arguments["--plot_radius"])

  delta = float(arguments["--delta"])
  x = np.arange(-plot_radius, plot_radius, delta)
  y = x
  X, Y = np.meshgrid(x, y)

  if background_metric == "Isotropic Schwarzschild":
    M = float(config_file["Isotropic_Schwarzschild_Settings"]["M"])
    bh_radius = 2 * M
    R = np.sqrt(X**2 + Y**2)
  elif background_metric == "Kerr-Schild Kerr":
    M = float(config_file["KerrSchild_Kerr_Settings"]["M"])
    a = float(config_file["KerrSchild_Kerr_Settings"]["a"])

    if a*a > M*M:
      warn("Naked singularity detected. Drawing a unit radius horizon.")
      bh_radius = 1
    else:
      bh_radius = M + np.sqrt(M**2 - a**2)

    part1 = X**2 + Y**2 - a**2
    part2 = np.abs(part1)
    
    R = np.sqrt((part1 + part2)/2)
    ergo = M * np.sqrt(2*(part1 + part2))/part2 - 1
    
    # Draw ergosphere
    ax.contour(X, Y, ergo, [0.0], colors="black", linestyles="--")
  else:
    raise Exception("Cannot plot data due to unrecognized metric: " + background_metric)
    
  # Draw Event horizon
  ax.contour(X, Y, R, [bh_radius], colors="black")

  # Draw Background
  background = plt.Circle((0, 0), background_radius, color="red", fill=False)
  ax.add_patch(background)

  plt.xlim([-plot_radius, plot_radius])
  plt.ylim([-plot_radius, plot_radius])


def plot_single_trajectory(plt, ax, clr, output_file_name):  
  data = pd.read_csv(output_file_name, delim_whitespace=True, names=vars)
  ax.plot(data["X1"], data["X2"], color=clr)

  plt.xlabel("$x$", fontsize = font_size);
  plt.ylabel("$y$", fontsize = font_size);

def plot_trajectory(arguments):
  config_file_name = arguments["<trajectory_config_file>"]
  output_file_name = arguments["<trajectory_output_file>"]

  with open(config_file_name, "r") as file:
    config_file = yaml.safe_load(file)
  
  plt.close("all")
  fig, ax = plt.subplots()

  plot_single_trajectory(plt, ax, arguments["--color"], output_file_name)
  plot_horizons(plt, ax, arguments, config_file)

  plt.show()

def plot_penrose(arguments):
  config_file_name = arguments["<penrose_config_file>"]
  trajectory_1 = arguments["<trajectory_1>"]
  trajectory_2 = arguments["<trajectory_2>"]
  trajectory_3 = arguments["<trajectory_3>"]

  with open(config_file_name, "r") as file:
    config_file = yaml.safe_load(file)
  
  color_1 = arguments["--color_1"]
  color_2 = arguments["--color_2"]
  color_3 = arguments["--color_3"]

  plt.close("all")
  fig, ax = plt.subplots()

  plot_single_trajectory(plt, ax, color_1, trajectory_1)
  plot_single_trajectory(plt, ax, color_2, trajectory_2)
  plot_single_trajectory(plt, ax, color_3, trajectory_3)
  plot_horizons(plt, ax, arguments, config_file)

  plt.show()

def plot_energy(arguments):
  output_file_name = arguments["<trajectory_output_file>"]
  font_size = int(arguments["--font_size"])

  data = pd.read_csv(output_file_name, delim_whitespace=True, names=vars)

  plt.close("all")

  plt.xlabel("$t$", fontsize = font_size)

  if arguments["local"]:
    plt.plot(data["time"], data["El"])
    plt.ylabel("$E_l$", fontsize = font_size)
  elif arguments["global"]:
    plt.plot(data["time"], data["Eg"])
    plt.ylabel("$E_g$", fontsize = font_size)
  elif arguments["residual"]:
    plt.loglog(data["time"], data["residual_Eg"])
    plt.ylabel("$|E_g(t) - E_g(0)|$", fontsize = font_size)

  plt.show()

def plot_gf(arguments):
    columns=["x", "y", "z", "data"]
    data = np.genfromtxt(arguments["<grid_function_data_file>"])
    df = pd.DataFrame(data, columns=columns)

    # Pick a z value to plot
    z_filtered_df = df[df["z"] == 0.0]

    plt.close("all")

    plt.xlabel("$x$", fontsize = font_size)
    plt.ylabel("$y$", fontsize = font_size)

    plt.tricontourf(z_filtered_df["x"], z_filtered_df["y"], z_filtered_df["data"], levels=np.linspace(0, 4.0, 100))
    plt.colorbar()
    plt.show()

# Main
if __name__ == '__main__':
  arguments = docopt(__doc__, version="GRLensing data plotter 1.0")
  
  # Global mpl settings
  font_size = int(arguments["--font_size"])
  mpl.rcParams['xtick.labelsize'] = font_size
  mpl.rcParams['ytick.labelsize'] = font_size

  if arguments["trajectory"]:
    plot_trajectory(arguments)
  elif arguments["energy"]:
    plot_energy(arguments)
  elif arguments["penrose"]:
    plot_penrose(arguments)
  elif arguments["gridfunction"]:
    plot_gf(arguments)
