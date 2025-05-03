import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append('../code')  # To import from code directory

from msd_analysis_protein import simulate_overdamped, compute_msd, compute_log_slope

# === Load dataset ===
data = pd.read_csv("../examples/histone.csv")

# Assumed structure: time, x, y
time_data = data['time'].values
x_data = data['x'].values
y_data = data['y'].values

dt = time_data[1] - time_data[0]
T = time_data[-1] - time_data[0]
x0, y0 = x_data[0], y_data[0]

# === Simulate using extracted parameters ===
inputs = {
    'mass': 1e-20,  # dummy mass (not used in overdamped)
    'radius': 1e-7,  # assumed particle size
    'x0': x0,
    'y0': y0,
    'dt': dt,
    'T': T,
    'is_underdamped': False
}

# From code/msd_analysis_protein.py
from msd_analysis_protein import initialize_parameters

gamma, D, N, time, kB, T_kelvin = initialize_parameters(inputs)
x_sim, y_sim = simulate_overdamped(inputs, D, N)

# === MSD from real data ===
msd_x_data = compute_msd(x_data)
msd_y_data = compute_msd(y_data)
msd_data_total = msd_x_data + msd_y_data

# === MSD from simulation ===
msd_x_sim = compute_msd(x_sim)
msd_y_sim = compute_msd(y_sim)
msd_sim_total = msd_x_sim + msd_y_sim

# === Slope comparisons ===
slope_data, _ = compute_log_slope(time[1:], msd_data_total[1:], fit_range=(time[1], time[-1]/10))
slope_sim, _ = compute_log_slope(time[1:], msd_sim_total[1:], fit_range=(time[1], time[-1]/10))

print(f"Slope from real data   : {slope_data:.4f}")
print(f"Slope from simulation  : {slope_sim:.4f}")
print(f"Difference              : {abs(slope_sim - slope_data):.4e}")

# === Test pass condition ===
assert np.isclose(slope_data, slope_sim, rtol=1e-2), "Slopes do not match within tolerance!"

print("✅ Test passed: MSD slopes match closely between real data and simulation.")
