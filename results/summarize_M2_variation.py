import os
import pandas as pd
import yaml

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

def main():
    data_dirs = make_file_list(0, 59)

    energy_gains = {}

    for data_dir in data_dirs:
        compute_extracted_energy(data_dir, energy_gains)

    save_results(energy_gains)

def make_file_list(start, end):
    list = []
    for i in range(start, end + 1):
        list.append("sks_breakup_M2_" +  str(i).zfill(3))
    return list

def compute_extracted_energy(data_dir, energy_gains):
    os.chdir(data_dir)

    with open("../../configs/sks/" + data_dir + ".yaml", "r") as file:
        config_file = yaml.safe_load(file)
    
    M2 = float(config_file["SKS_Settings"]["M2"])

    data_1 = pd.read_csv("trajectory_1.ascii", delim_whitespace=True, names=vars)
    data_3 = pd.read_csv("trajectory_3.ascii", delim_whitespace=True, names=vars)

    energy_gains[data_dir] = {
        "M2": M2,
        "El": data_3["El"].iloc[-1] - data_1["El"].iloc[-1],
        "Eg": data_3["Eg"].iloc[-1] - data_1["Eg"].iloc[-1],
        "Eff_l": (data_3["El"].iloc[-1] - data_1["El"].iloc[-1])/data_1["El"].iloc[-1],
        "Eff_g": (data_3["Eg"].iloc[-1] - data_1["Eg"].iloc[-1])/data_1["Eg"].iloc[-1]
    }

    os.chdir("../")

def save_results(energy_gains):
    output_file = open("M2_variation_sumary.txt", "w")

    output_file.write("# 1:M2 2:Local energy gain 3:Global energy gain 4:Eff. with global energy, Eff. with local energy\n")
    
    for run in sorted(energy_gains.keys()):
        output_file.write("%.16f" % energy_gains[run]["M2"])
        output_file.write("    ")
        output_file.write("%.16f" % energy_gains[run]["El"])
        output_file.write("    ")
        output_file.write("%.16f" % energy_gains[run]["Eg"])
        output_file.write("    ")
        output_file.write("%.16f" % energy_gains[run]["Eff_l"])
        output_file.write("    ")
        output_file.write("%.16f" % energy_gains[run]["Eff_g"])
        output_file.write("\n")

    output_file.close()

main()
