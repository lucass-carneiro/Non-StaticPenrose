import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

file_name = "trajectory.ascii"
font_size = 20
range = 30

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'Latin Modern Roman'
plt.rcParams['figure.figsize'] = [9, 9]
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size

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

data = pd.read_csv(file_name, delim_whitespace=True, names=vars);

plt.close("all")

horizon = plt.Circle((0, 0), 10.0, color="black", fill=False)

fig, ax = plt.subplots()
ax.add_patch(horizon)


ax.plot(
    data["X1"],
    data["X2"]
)

plt.xlabel("$x$", fontsize = font_size);
plt.ylabel("$y$", fontsize = font_size);

plt.xlim([-range, range])
plt.ylim([-range, range])

plt.show()