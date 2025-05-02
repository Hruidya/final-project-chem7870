This is a GitHub repository created by Mr. Aakash Gupta, Ms. Hruidya Babu and Ms. Jyotirsika Dalal for the purpose of their final project on the course CHEM7870 offered by C&amp;CB, Cornell University. The code package aims to simulate the underdamped Brownian motion of proteins in water with the help of experimental Single-Molecule Tracking data.

# MSD Analysis for Single-Molecule Protein Tracking

This project analyzes the **mean squared displacement (MSD)** of a protein undergoing **underdamped Brownian motion** based on time-resolved single-molecule tracking data. It uses numerical integration of the velocity autocorrelation function (VACF) to compute MSD along both x and y directions.

---

## Input

- A CSV file named `filename.csv` containing:
  - `t`: time values (in seconds)
  - `x`, `y`: positions of the tracked molecule (in meters)

---

## How to Run

1. Make sure you have Python 3 and the following packages:
   - `numpy`
   - `pandas`
   - `matplotlib`

2. Prepare your dataset in the correct format.

3. Run the script:

```bash
python msd_analysis_protein.py
```

4. You will be prompted to enter:
   - Mass of the protein (in kg)
   - Radius of the protein (in meters)

---

## Output

A plot will be generated:

- **log₁₀(MSD) vs Time**

---

## 📌 Notes

- The analysis assumes zero initial position and velocity: \( x(0) = 0, v(0) = 0 \)
- Temperature is set to 298.15 K (room temperature)
- Solvent is assumed to be water (viscosity (eta) = 1*10^{-3} Pa·s)
