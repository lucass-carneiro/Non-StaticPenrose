import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

file_name = "trajectory.ascii"
font_size = 20

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

plt.plot(
    data["time"],
    data["EN"]
)

plt.xlabel("$t$", fontsize = font_size);
plt.ylabel("$EN$", fontsize = font_size);

plt.show()