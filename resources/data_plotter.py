"""GRLensing data plotter.

Usage:
  data_plotter.py trajectory <trajectory_config_file> <trajectory_output_file> [--font_size=<size>]
  data_plotter.py energy (local|global) <trajectory_output_file>
  data_plotter.py (-h | --help)
  data_plotter.py --version

Options:
  -h --help           Show this screen.
  --version           Show version.
  --font_size=<size>  The size of the font in the plots [default: 20].

"""

from docopt import docopt
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
  "Eg"
]

def plot_trajectory(arguments):
  config_file_name = arguments["<trajectory_config_file>"]
  output_file_name = arguments["<trajectory_output_file>"]
  font_size = int(arguments["--font_size"])
  
  # Open the trajectory config file and extract relevant data
  with open(config_file_name, "r") as file:
    config_file = yaml.safe_load(file)
  
  background_radius = float(config_file["background_radius"])
  background_metric = config_file["background_metric"]

  if background_metric == "Isotropic Schwarzschild":
    M = float(config_file["Isotropic_Schwarzschild_Settings"]["M"])
    bh_radius = 2 * M
  else:
    raise Exception("Cannot plot data due to unrecognized metric: " + background_metric)

  # mpl options
  mpl.rcParams['xtick.labelsize'] = font_size
  mpl.rcParams['ytick.labelsize'] = font_size

  data = pd.read_csv(output_file_name, delim_whitespace=True, names=vars)

  plt.close("all")

  horizon = plt.Circle((0, 0), bh_radius, color="black", fill=False)
  background = plt.Circle((0, 0), background_radius, color="red", fill=False)

  fig, ax = plt.subplots()
  ax.add_patch(horizon)
  ax.add_patch(background)


  ax.plot(
      data["X1"],
      data["X2"]
  )

  plt.xlabel("$x$", fontsize = font_size);
  plt.ylabel("$y$", fontsize = font_size);

  plt.xlim([-background_radius, background_radius])
  plt.ylim([-background_radius, background_radius])

  plt.show()

def plot_energy(arguments):
  output_file_name = arguments["<trajectory_output_file>"]
  font_size = int(arguments["--font_size"])

  mpl.rcParams['xtick.labelsize'] = font_size
  mpl.rcParams['ytick.labelsize'] = font_size

  data = pd.read_csv(output_file_name, delim_whitespace=True, names=vars)

  plt.close("all")

  plt.xlabel("$t$", fontsize = font_size)

  if arguments["local"]:
    plt.plot(
      data["time"],
      data["El"]
    )

    plt.ylabel("$E_l$", fontsize = font_size)

  else:
    plt.plot(
      data["time"],
      data["Eg"]
    )

    plt.ylabel("$E_g$", fontsize = font_size)

  plt.show()

# Main
if __name__ == '__main__':
  arguments = docopt(__doc__, version="GRLensing data plotter 1.0")
  
  if arguments["trajectory"]:
    plot_trajectory(arguments)
  elif arguments["energy"]:
    plot_energy(arguments)
