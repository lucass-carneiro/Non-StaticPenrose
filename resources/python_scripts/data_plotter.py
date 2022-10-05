#!/bin/python

"""GRLensing data plotter.

Usage:
  data_plotter.py trajectory <trajectory_config_file> <trajectory_output_file> [--font_size=<size>] [--delta=<value>] [--plot_radius=<radius>] [--color=<color>]
  data_plotter.py animated-trajectory <trajectory_config_file> <trajectory_output_file> <animation_output_file> [--workers=<workers>] [--delay=<delay>] [--loop] [--keep_frames] [--intermediate_format=<format>] [--font_size=<size>] [--delta=<value>] [--plot_radius=<radius>] [--color=<color>]
  data_plotter.py energy (local|global|residual) <trajectory_output_file> [--font_size=<size>]
  data_plotter.py penrose <penrose_config_file> <trajectory_1> <trajectory_2> <trajectory_3> [--font_size=<size>] [--color_1=<color1>] [--color_2=<color2>] [--color_3=<color3>] [--plot_radius=<radius>]
  data_plotter.py penrose-energy (local|global) <trajectory_1> <trajectory_3> [--font_size=<size>] [--color_1=<color1>] [--color_3=<color3>] [--plot_radius=<radius>]
  data_plotter.py penrose-energy-diff <trajectory_1> <trajectory_3>
  data_plotter.py gridfunction <grid_function_data_file>
  data_plotter.py (-h | --help)
  data_plotter.py --version

Options:
  -h --help                         Show this screen.
  --version                         Show version.
  --workers=<workers>               Number of parallel workers to use while producing frames. [default: 4]
  --delay=<delay>                   The delay between animation frames (in ms). [default: 0].
  --loop                            Whether or not to create a looping gif.
  --keep_frames                     Whether or not to keep the intermediate animation frames.
  --intermediate_format=<format>    The intermediate format to use while creating animation frames. [default: pdf]
  --font_size=<size>                The size of the font in the plots [default: 20].
  --delta=<value>                   The grid spacing used in internal grid generation. [default: 0.1].
  --plot_radius=<radius>            The radius of the background.
  --color=<color>                   The color of the trajectory. [default: black]
  --color_1=<color1>                The color of the first breakup trajectory. [default: red]
  --color_2=<color2>                The color of the first breakup trajectory. [default: black]
  --color_3=<color3>                The color of the first breakup trajectory. [default: blue]

"""

from docopt import docopt
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import yaml
import os
from multiprocess import Pool

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

def plot_horizons(plt, ax, arguments, config_file, t):
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
    
    # Draw Event horizon
    ax.contour(X, Y, R, [bh_radius], colors="black")

  elif background_metric == "Kerr-Schild Kerr":
    M = float(config_file["KerrSchild_Kerr_Settings"]["M"])
    a = float(config_file["KerrSchild_Kerr_Settings"]["a"])

    if a*a > M*M:
      warn("Naked singularity detected. Drawing a unit radius horizon.")
      bh_radius = 1
    else:
      bh_radius = M + np.sqrt(M**2 - a**2)

    rho2ma2 = X**2 + Y**2 - a**2
    
    R = np.sqrt((rho2ma2 + np.abs(rho2ma2))/2)
    gtt = 2*M/R - 1
    
    # Draw ergosphere
    ax.contour(X, Y, gtt, [0.0], colors="black", linestyles="--")

    # Draw Event horizon
    ax.contour(X, Y, R, [bh_radius], colors="black")
    
  elif background_metric == "Superimposed Kerr binary in Kerr-Schild coordinates":
    M1 = float(config_file["SKS_Settings"]["M1"])
    M2 = float(config_file["SKS_Settings"]["M2"])
    
    a1 = float(config_file["SKS_Settings"]["a1"])
    a2 = float(config_file["SKS_Settings"]["a2"])
    
    b = float(config_file["SKS_Settings"]["b"])

    Omega = np.sqrt((M1 + M2)/b**3)
    
    partA_1 = X**2 + Y**2 - a1**2
    partA_2 = X**2 + Y**2 - a1**2
    
    partB_1 = np.abs(partA_1)
    partB_2 = np.abs(partA_2)
    
    R_1 = np.sqrt((partA_1 + partB_1)/2)
    R_2 = np.sqrt((partA_2 + partB_2)/2)

    bh_radius_1 = M1 + np.sqrt(M1**2 - a1**2)
    bh_radius_2 = M2 + np.sqrt(M2**2 - a2**2)
    
    # Draw Event horizons
    ax.contour(X - b/2 * np.cos(Omega * t), Y - b/2 * np.sin(Omega * t), R_1, [bh_radius_1], colors="black")
    ax.contour(X + b/2 * np.cos(Omega * t), Y + b/2 * np.sin(Omega * t), R_2, [bh_radius_2], colors="black")
  else:
    raise Exception("Cannot plot data due to unrecognized metric: " + background_metric)

  # Draw Background
  background = plt.Circle((0, 0), background_radius, color="red", fill=False)
  ax.add_patch(background)

  plt.xlim([-plot_radius, plot_radius])
  plt.ylim([-plot_radius, plot_radius])


def plot_single_trajectory(plt, ax, clr, output_file_name, config_file, draw_horizons=True):
  data = pd.read_csv(output_file_name, delim_whitespace=True, names=vars)
  ax.plot(data["X1"], data["X2"], color=clr)

  plt.xlabel("$x$", fontsize = font_size);
  plt.ylabel("$y$", fontsize = font_size);

  if draw_horizons:
    plot_horizons(plt, ax, arguments, config_file, data["time"].iloc[-1])

def plot_trajectory(arguments):
  config_file_name = arguments["<trajectory_config_file>"]
  output_file_name = arguments["<trajectory_output_file>"]

  with open(config_file_name, "r") as file:
    config_file = yaml.safe_load(file)
  
  plt.close("all")
  fig, ax = plt.subplots()

  plot_single_trajectory(plt, ax, arguments["--color"], output_file_name, config_file)

  plt.tight_layout()
  plt.show()

def plot_instant(plt, ax, clr, data, index, config_file, draw_horizons=True):
  ax.plot(data["X1"].iloc[index], data["X2"].iloc[index], marker="o", markersize=6, markeredgecolor=clr, markerfacecolor=clr)

  plt.xlabel("$x$", fontsize = font_size);
  plt.ylabel("$y$", fontsize = font_size);

  if draw_horizons:
    plot_horizons(plt, ax, arguments, config_file, data["time"].iloc[index])

def plot_animated_trajectory(arguments):
  config_file_name = arguments["<trajectory_config_file>"]
  output_file_name = arguments["<trajectory_output_file>"]

  with open(config_file_name, "r") as file:
    config_file = yaml.safe_load(file)
  
  data = pd.read_csv(output_file_name, delim_whitespace=True, names=vars)
  
  def frame_output_kernerl(index):
    plt.close("all")
    fig, ax = plt.subplots()
    plot_instant(plt, ax, arguments["--color"], data, index, config_file)
    plt.savefig("frame_" + (f"{index}").rjust(4, "0") + ".pdf")
  
  process_pool = Pool(int(arguments["--workers"]))
  process_pool.map_async(frame_output_kernerl, range(0, data["time"].size))
  process_pool.close()
  process_pool.join()
  
  # Call imagemagick to convert stills to gif
  os.system("convert -delay " + arguments["--delay"] + " -loop " + ("0" if arguments["--loop"] else "1") + " *." + arguments["--intermediate_format"] + " " + arguments["<animation_output_file>"] + ".gif")
  
  if not arguments["--keep_frames"]:
    os.system("rm *." + arguments["--intermediate_format"])

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

  plot_single_trajectory(plt, ax, color_1, trajectory_1, config_file, False)
  plot_single_trajectory(plt, ax, color_2, trajectory_2, config_file, False)
  plot_single_trajectory(plt, ax, color_3, trajectory_3, config_file)

  plt.tight_layout()
  plt.show()

def compute_penrose_energy_diff(arguments):
  #data_plotter.py penrose-renegy-diff (local|global) <trajectory_1> <trajectory_3>
  trajectory_1 = arguments["<trajectory_1>"]
  trajectory_3 = arguments["<trajectory_3>"]

  data_1 = pd.read_csv(trajectory_1, delim_whitespace=True, names=vars)
  data_3 = pd.read_csv(trajectory_3, delim_whitespace=True, names=vars)

  local_energy_1 = data_1["El"].iloc[-1]
  local_energy_3 = data_3["El"].iloc[-1]

  global_energy_1 = data_1["Eg"].iloc[-1]
  global_energy_3 = data_3["Eg"].iloc[-1]

  print("Local energy difference (3 - 1): ", local_energy_3 - local_energy_1)
  print("Global energy difference (3 - 1): ", global_energy_3 - global_energy_1)
  print("Global/localenergy absolute difference (1): ", np.abs(global_energy_1 - local_energy_1))
  print("Global/localenergy absolute difference (3): ", np.abs(global_energy_3 - local_energy_3))

def plot_penrose_energy(arguments):
  trajectory_1 = arguments["<trajectory_1>"]
  trajectory_3 = arguments["<trajectory_3>"]
  
  color_1 = arguments["--color_1"]
  color_3 = arguments["--color_3"]

  data_1 = pd.read_csv(trajectory_1, delim_whitespace=True, names=vars)
  data_3 = pd.read_csv(trajectory_3, delim_whitespace=True, names=vars)

  plt.close("all")

  plt.xlabel("$t$", fontsize = font_size)

  if arguments["local"]:
    plt.plot(data_1["time"], data_1["El"], color=color_1)
    plt.plot(data_3["time"], data_3["El"], color=color_3)
    plt.ylabel("$E_l$", fontsize = font_size)
  elif arguments["global"]:
    plt.plot(data_1["time"], data_1["Eg"], color=color_1)
    plt.plot(data_3["time"], data_3["Eg"], color=color_3)
    plt.ylabel("$E_g$", fontsize = font_size)

  plt.tight_layout()
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

  plt.tight_layout()
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

    plt.tricontourf(z_filtered_df["y"], z_filtered_df["x"], z_filtered_df["data"], levels=np.linspace(-7, 1, 100))
    plt.colorbar()
    plt.tight_layout()
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
  elif arguments["animated-trajectory"]:
    plot_animated_trajectory(arguments)
  elif arguments["energy"]:
    plot_energy(arguments)
  elif arguments["penrose"]:
    plot_penrose(arguments)
  elif arguments["penrose-energy"]:
    plot_penrose_energy(arguments)
  elif arguments["penrose-energy-diff"]:
    compute_penrose_energy_diff(arguments)
  elif arguments["gridfunction"]:
    plot_gf(arguments)
