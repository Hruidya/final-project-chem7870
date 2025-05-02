"""
msd_analysis_protein.py

This script computes and plots the mean squared displacement (MSD) from single-molecule tracking (SMT) data
for a protein undergoing underdamped Brownian motion in water. It uses numerical integration of the velocity
autocorrelation function to calculate MSD in both x and y directions.

The SMT dataset must contain three columns: 't' (time in seconds), 'x', and 'y' (positions in meters).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def load_data(filepath):
    """Load the CSV file and extract time, x, and y arrays."""
    df = pd.read_csv(filepath)
    return df['t'].values, df['x'].values, df['y'].values

def get_physical_constants():
    """Prompt the user for protein mass and radius, and return all physical constants."""
    mass = float(input("Enter the mass of the protein (in kilograms): "))
    radius = float(input("Enter the radius of the protein (in meters): "))
    kB = 1.380649e-23  # J/K
    T = 298.15  # K
    eta = 1e-3  # Pa·s
    gamma = 6 * np.pi * eta * radius
    return mass, radius, kB, T, eta, gamma

def compute_velocity(x, t):
    """Compute velocity using finite differences."""
    return np.gradient(x, t)

def compute_vacf(v):
    """Compute the velocity autocorrelation function."""
    N = len(v)
    return np.array([np.mean(v[:N-lag] * v[lag:]) for lag in range(N)])

def compute_msd_from_vacf(vacf, dt):
    """Compute MSD from the velocity autocorrelation function using double integration."""
    return 2 * np.cumsum([np.sum(vacf[:i]) for i in range(1, len(vacf)+1)]) * dt**2

def plot_results(time, msd_x, msd_y):
    """Plot log(MSD) versus time."""

    plt.figure(figsize=(8, 5))
    plt.plot(time, np.log10(msd_x), label='log10(MSD_x)', linewidth=2)
    plt.plot(time, np.log10(msd_y), label='log10(MSD_y)', linewidth=2)
    plt.xlabel("Time (s)")
    plt.ylabel("log10(MSD) (m²)")
    plt.title("Log-MSD vs Time")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

def main():
    filepath = 'filename.csv'  # Update as necessary
    t, x, y = load_data(filepath)
    mass, radius, kB, T, eta, gamma = get_physical_constants()
    
    v_x = compute_velocity(x, t)
    v_y = compute_velocity(y, t)
    dt = np.mean(np.diff(t))

    vacf_x = compute_vacf(v_x)
    vacf_y = compute_vacf(v_y)

    msd_x = compute_msd_from_vacf(vacf_x, dt)
    msd_y = compute_msd_from_vacf(vacf_y, dt)

    plot_results(t, msd_x, msd_y)

if __name__ == "__main__":
    main()
