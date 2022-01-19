import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Courtesy of github.com/iaraota
def to_grid():
	columns=["x", "y", "z", "data"]
	data = np.genfromtxt("Isotropic Schwarzschild_lapse_0.ascii")
	df = pd.DataFrame(
			data,
			columns=columns,
			).drop_duplicates().reset_index(drop=True)	

	x_unique = df["x"].drop_duplicates().sort_values().reset_index(drop=True)	
	y_unique = df["y"].drop_duplicates().sort_values().reset_index(drop=True)	

	df = df.sort_values(by=["y"])

	x_grid, y_grid = np.meshgrid(x_unique, y_unique)
	z_grid = np.array([df["data"][df["x"] == x_unique[i]].values for i in range(len(x_unique))])

	plt.contourf(x_grid, y_grid, z_grid)
	plt.colorbar()

	# plt.contour(x_grid, y_grid, z_grid, cmap = "Greys_r")

	plt.show()


to_grid()