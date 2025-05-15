from utils import is_cuda_available
import pytest
import os

# Skip all tests in this file if no GPU is available
if not is_cuda_available():
    pytest.skip("Skipping tests because no GPU is available.", allow_module_level=True)

import pyUAMMD
import numpy as np

atol = 1e-5
rtol = 1e-5


def bond3_base_simulation():
    temperature = 0
    pos0 = [-1.1, -5.5, 2]
    pos1 = [4.5, 2.2, -4.3]
    pos2 = [3.5, 6.1, 3.3]
    L = 20
    k = 200
    theta0 = 1.23
    box = [L, L, L]

    base_simulation = {
        "system": {
            "info": {
                "type": ["Simulation", "Information"],
                "parameters": {"name": "testHessian"},
            }
        },
        "global": {
            "types": {
                "type": ["Types", "Basic"],
                "labels": ["name", "radius", "mass", "charge"],
                "data": [["A", 1, 0, 0]],
            },
            "ensemble": {
                "type": ["Ensemble", "NVT"],
                "labels": ["box", "temperature"],
                "data": [[box, temperature]],
            },
        },
        "integrator": {
            "BBK": {
                "type": ["Langevin", "BBK"],
                "parameters": {"timeStep": 0.001, "frictionConstant": 1},
            },
            "schedule": {
                "type": ["Schedule", "Integrator"],
                "labels": ["order", "integrator", "steps"],
                "data": [[1, "BBK", 1]],
            },
        },
        "state": {
            "labels": ["id", "position"],
            "data": [[0, pos0], [1, pos1], [2, pos2]],
        },
        "topology": {
            "forceField": {
                "Bond": {
                    "type": ["Bond3", "HarmonicAngular"],
                    "labels": ["id_i", "id_j", "id_k", "K", "theta0"],
                    "data": [[0, 1, 2, k, theta0]],
                    "parameters": {},
                }
            },
            "structure": {
                "labels": ["id", "type"],
                "data": [[0, "A"], [1, "A"], [2, "A"]],
            },
        },
    }
    return base_simulation


@pytest.mark.parametrize("hessian_mode", ["Analytical", "Numerical"])
def test_bond3_harmonic_symmetric(tmp_path, hessian_mode):
    """Tests that the HessianMeasure is symmetric."""

    hessian_path = os.path.join(tmp_path, "Hessian.out")
    base_simulation = bond3_base_simulation()
    base_simulation["simulationStep"] = {
        "Hessian": {
            "type": ["MechanicalMeasure", "HessianMeasure"],
            "parameters": {
                "mode": hessian_mode,
                "outputFilePath": hessian_path,
                "intervalStep": 1,
            },
        }
    }

    simulation = pyUAMMD.simulation(base_simulation)
    simulation.run()
    assert os.path.exists(hessian_path)

    hessian = np.loadtxt(hessian_path)
    # Hessian file has shape (npairs, 12), first two columns are the pair indices
    # Transform to a (N, N, 3, 3) array
    n = int(np.sqrt(hessian.shape[0]))
    hessian = hessian[:, 2:].reshape((n, n, 3, 3))
    # Check symmetry, hessian[i,j] == hessian[j,i].T
    hessian_t = np.transpose(hessian, (1, 0, 2, 3))
    max_diff = np.max(np.abs(hessian - hessian_t) / np.mean(hessian))
    np.allclose(
        hessian, hessian_t, rtol=rtol, atol=atol
    ), f"The HessianMeasure is not symmetric. Max diff: {max_diff}"


def test_bond3_harmonic_self_consistent(tmp_path):
    """Tests that the HessianMeasure is self-consistent between the analytical and numerical methods."""

    analytical_hessian_path = os.path.join(tmp_path, "HessianAnalytical.out")
    numerical_hessian_path = os.path.join(tmp_path, "HessianNumerical.out")
    base_simulation = bond3_base_simulation()
    base_simulation["simulationStep"] = {
        "Hessian_analytical": {
            "type": ["MechanicalMeasure", "HessianMeasure"],
            "parameters": {
                "mode": "Analytical",
                "outputFilePath": analytical_hessian_path,
                "intervalStep": 1,
            },
        },
        "Hessian_numerical": {
            "type": ["MechanicalMeasure", "HessianMeasure"],
            "parameters": {
                "mode": "Numerical",
                "outputFilePath": numerical_hessian_path,
                "intervalStep": 1,
            },
        },
    }

    simulation = pyUAMMD.simulation(base_simulation)
    simulation.run()
    assert os.path.exists(analytical_hessian_path)
    assert os.path.exists(numerical_hessian_path)

    numerical = np.loadtxt(analytical_hessian_path)
    analytical = np.loadtxt(numerical_hessian_path)

    assert numerical.shape == analytical.shape
    max_diff = np.max(np.abs(numerical - analytical) / np.mean(analytical))
    assert np.allclose(
        numerical, analytical, rtol=rtol, atol=atol
    ), f"The HessianMeasure results are not the same between the analytical and numerical methods. Max diff: {max_diff}"
