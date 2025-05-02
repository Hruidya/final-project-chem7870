import numpy as np
import pandas as pd
import pytest
import sys
from pathlib import Path

# Add the code directory to sys.path
sys.path.append(str(Path(__file__).resolve().parents[1] / "code"))

from msd_analysis_protein import (
    load_data, compute_velocity, compute_vacf, compute_msd_from_vacf
)

@pytest.fixture(scope="module")
def sample_data():
    # Load dataset from example directory
    data_path = Path(__file__).resolve().parents[1] / "example" / "filename.csv"
    time, x, y = load_data(str(data_path))
    dt = np.mean(np.diff(time))
    return time, x, y, dt

def test_data_loaded_correctly(sample_data):
    time, x, y, _ = sample_data
    assert len(time) == len(x) == len(y)
    assert np.all(np.diff(time) > 0), "Time values must be strictly increasing."

def test_velocity_calculation(sample_data):
    time, x, _, _ = sample_data
    v_x = compute_velocity(x, time)
    assert len(v_x) == len(x)
    assert not np.isnan(v_x).any(), "Velocity contains NaNs."

def test_vacf_output(sample_data):
    time, _, y, _ = sample_data
    v_y = compute_velocity(y, time)
    vacf = compute_vacf(v_y)
    assert len(vacf) == len(v_y)
    assert np.all(np.isfinite(vacf)), "VACF contains non-finite values."

def test_msd_monotonicity(sample_data):
    time, x, _, dt = sample_data
    v_x = compute_velocity(x, time)
    vacf_x = compute_vacf(v_x)
    msd_x = compute_msd_from_vacf(vacf_x, dt)
    assert np.all(msd_x >= 0), "MSD must be non-negative."
    assert np.all(np.diff(msd_x) >= 0), "MSD must be non-decreasing."
