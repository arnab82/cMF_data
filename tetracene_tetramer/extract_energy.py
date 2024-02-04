import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors


# Function to extract energy values from an output file
def extract_energy_values(file_path):
    energy_values = {}
    found_keywords = ['BFGS', 'CG', 'GD','DIIS_GD', 'DIIS_WOP', 'DIIS_WP', 'NEWTON_WP','NEWTON_WOP']

    for keyword in found_keywords:
        energy_values[keyword] = []

    with open(file_path, 'r') as file:
        for line in file:
            for keyword in found_keywords:
                if keyword in line:
                    while True:
                        match1 = re.search(r'E:\s*(-?\d+\.\d+)', line)
                        print(match1)
                        match2 = re.search(r'Total=\s*(-?\d+\.\d+)', line)
                        print(match2)
                        if match1:
                            energy_n = float(match1.group(1))
                            energy_values[keyword].append(energy_n)
                        elif match2:
                            energy_n = float(match2.group(1))
                            energy_values[keyword].append(energy_n)
                        line = next(file, None)
                        if line is None or any(keyword in line for keyword in found_keywords):
                            break
                else:
                    continue

    return energy_values
# Replace 'your_output_file.txt' with the path to your output file
output_file_path = '/Users/arnab/arnab/workspace/plot/cmf.out'

# Extract energy values from the output file
energy_values = extract_energy_values(output_file_path)
print(energy_values)

# Convert energy_values to a 1D array
energy_bfgs = np.array(energy_values['BFGS']).flatten()
energy_cg = np.array(energy_values['CG']).flatten()
energy_gd = np.array(energy_values['GD']).flatten()
energy_diis_gd = np.array(energy_values['DIIS_GD']).flatten()
energy_diis_wop = np.array(energy_values['DIIS_WOP']).flatten()
energy_diis_wp = np.array(energy_values['DIIS_WP']).flatten()
energy_newton_wp = np.array(energy_values['NEWTON_WP']).flatten()
energy_newton_wop = np.array(energy_values['NEWTON_WOP']).flatten()

# Stack the arrays vertically
energy_lists = [energy_bfgs, energy_cg, energy_gd, energy_diis_gd, energy_diis_wop, energy_diis_wp, energy_newton_wp, energy_newton_wop]

with open('energy_values.txt', 'w') as file:
    for i, energy_list in enumerate(energy_lists):
        list_name = ['BFGS', 'CG', 'GD', 'DIIS_GD', 'DIIS_WOP', 'DIIS_WP', 'NEWTON_WP', 'NEWTON_WOP'][i]
        file.write(f"{list_name}:\n")
        for energy in energy_list:
            file.write(f"{energy:.4f},")
        file.write('\n')

