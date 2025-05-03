"""
msd_analysis_protein.py

This Python script simulates Brownian motion of a protein-sized particle in a fluid medium, using either 
the underdamped or overdamped Langevin equation depending on user input. It models stochastic motion due to 
thermal fluctuations and viscous drag, calculates the mean squared displacement (MSD), and plots the log-log plot of MSD versus time.

import numpy as np
import matplotlib.pyplot as plt

def get_user_inputs():
    """
    Prompt the user for simulation parameters.

    Returns:
        dict: Dictionary of input parameters.
    """
    inputs = {
        'mass': float(input("Enter the mass of the particle (kg): ")),
        'radius': float(input("Enter the radius of the particle (m): ")),
        'x0': float(input("Enter the initial x position (m): ")),
        'y0': float(input("Enter the initial y position (m): ")),
        'dt': float(input("Enter the timestep (s): ")),
        'T': float(input("Enter total simulation time (s): ")),
        'is_underdamped': input("Use underdamped Langevin? (yes/no): ").lower().startswith('y')
    }
    return inputs

def initialize_parameters(inputs):
    """
    Initialize derived physical and simulation parameters.

    Args:
        inputs (dict): User input parameters.

    Returns:
        tuple: Derived parameters (gamma, D, N, time array).
    """
    kB = 1.380649e-23  # Boltzmann constant (J/K)
    T_kelvin = 298.15  # Room temperature (K)
    eta = 1e-3         # Water viscosity (Pa·s)
    
    gamma = 6 * np.pi * eta * inputs['radius']
    D = kB * T_kelvin / gamma
    N = int(inputs['T'] / inputs['dt'])
    time = np.linspace(0, inputs['T'], N)

    return gamma, D, N, time, kB, T_kelvin

def simulate_underdamped(inputs, gamma, N, kB, T_kelvin):
    """
    Simulate underdamped Langevin dynamics.

    Returns:
        tuple: Position arrays x, y
    """
    np.random.seed(42)
    dt = inputs['dt']
    mass = inputs['mass']

    x, y = np.zeros(N), np.zeros(N)
    vx, vy = np.zeros(N), np.zeros(N)
    x[0], y[0] = inputs['x0'], inputs['y0']

    for i in range(1, N):
        fx = -gamma * vx[i-1] + np.sqrt(2 * gamma * kB * T_kelvin / dt) * np.random.randn()
        fy = -gamma * vy[i-1] + np.sqrt(2 * gamma * kB * T_kelvin / dt) * np.random.randn()
        ax, ay = fx / mass, fy / mass
        vx[i] = vx[i-1] + ax * dt
        vy[i] = vy[i-1] + ay * dt
        x[i] = x[i-1] + vx[i] * dt
        y[i] = y[i-1] + vy[i] * dt

    return x, y

def simulate_overdamped(inputs, D, N):
    """
    Simulate overdamped Langevin dynamics.

    Returns:
        tuple: Position arrays x, y
    """
    np.random.seed(42)
    dt = inputs['dt']

    x, y = np.zeros(N), np.zeros(N)
    x[0], y[0] = inputs['x0'], inputs['y0']

    for i in range(1, N):
        x[i] = x[i-1] + np.sqrt(2 * D * dt) * np.random.randn()
        y[i] = y[i-1] + np.sqrt(2 * D * dt) * np.random.randn()

    return x, y

def compute_msd(r):
    """
    Compute mean squared displacement (MSD) of a trajectory.

    Args:
        r (np.ndarray): Position array.

    Returns:
        np.ndarray: MSD as a function of time lag.
    """
    N = len(r)
    msd = np.zeros(N)
    for lag in range(1, N):
        diffs = r[lag:] - r[:-lag]
        msd[lag] = np.mean(diffs**2)
    return msd

def compute_log_slope(time, msd_total, fit_range=(0.1, 1.0)):
    """
    Fit a line to log-log MSD data in a specified time range to compute the slope.

    Args:
        time (np.ndarray): Time array.
        msd_total (np.ndarray): Total MSD array.
        fit_range (tuple): (t_min, t_max) for fitting in seconds.

    Returns:
        float: Slope of log10(MSD) vs log10(time) in the specified range.
    """
    # Filter based on time range
    mask = (time > fit_range[0]) & (time < fit_range[1])
    log_time = np.log10(time[mask])
    log_msd = np.log10(msd_total[mask])

    # Linear fit
    slope, intercept = np.polyfit(log_time, log_msd, 1)
    return slope, (log_time, log_msd)

def plot_log_msd(time, msd_total, slope, fit_data):
    """
    Plot log10(MSD) versus log10(time) and annotate slope.

    Args:
        time (np.ndarray): Time array.
        msd_total (np.ndarray): Total MSD array.
        slope (float): Fitted slope of log-log MSD.
        fit_data (tuple): (log_time, log_msd) used in fitting.
    """
    log_time_full = np.log10(time[1:])
    log_msd_full = np.log10(msd_total[1:])
    
    plt.figure(figsize=(8, 5))
    plt.plot(log_time_full, log_msd_full, label='log10(MSD_total)', linewidth=2)
    plt.plot(fit_data[0], np.poly1d(np.polyfit(fit_data[0], fit_data[1], 1))(fit_data[0]), 'r--', label=f'Fit Slope = {slope:.2f}')
    plt.xlabel("log10(Time) [s]")
    plt.ylabel("log10(MSD) [m²]")
    plt.title("Simulated Brownian Motion: log(MSD) vs log(Time)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def main():
    inputs = get_user_inputs()
    gamma, D, N, time, kB, T_kelvin = initialize_parameters(inputs)

    if inputs['is_underdamped']:
        x, y = simulate_underdamped(inputs, gamma, N, kB, T_kelvin)
    else:
        x, y = simulate_overdamped(inputs, D, N)

    msd_x = compute_msd(x)
    msd_y = compute_msd(y)
    msd_total = msd_x + msd_y

    slope, fit_data = compute_log_slope(time[1:], msd_total[1:], fit_range=(time[1], time[-1] / 10))
    print(f"Estimated slope of log-log MSD: {slope:.3f}")

    plot_log_msd(time, msd_total, slope, fit_data)
