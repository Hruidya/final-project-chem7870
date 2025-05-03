This is a GitHub repository created by Mr. Aakash Gupta, Ms. Hruidya Babu and Ms. Jyotirsika Dalal for the purpose of their final project on the course CHEM7870 offered by C&amp;CB, Cornell University. The code package aims to simulate the underdamped Brownian motion of proteins in water with the help of experimental Single-Molecule Tracking data.

# Brownian Motion Simulator: Overdamped & Underdamped Langevin Dynamics

This Python package simulates the Brownian motion of a protein-like particle suspended in a fluid, governed by either **overdamped** or **underdamped Langevin equations**. It computes the **mean squared displacement (MSD)** and visualizes the result on a log-log plot to study diffusive behavior.

---

## Features

- Simulates both **overdamped** (diffusive) and **underdamped** (inertial) Langevin motion.
- Computes **mean squared displacement (MSD)**.
- Plots **log(MSD) vs log(time)**.
- Computes and overlays the **slope of the log-log MSD** plot to characterize the regime (slope ≈ 1 for diffusive regime and ≈ 2 for ballistic regime).

---
## Input Simulation Parameters:

    1.Particle mass (kg)
    2.Particle radius (m)
    3.Initial position (x and y, in meters)
    4.Timestep dt (seconds)
    5.Total simulation time T (seconds)
    6.Damping regime (yes for underdamped, no for overdamped)

## Requirements

- Python 3.x
- `numpy`
- `matplotlib`

## Runnning

Run the script using Python:

```bash
python msd_analysis_protein.py
