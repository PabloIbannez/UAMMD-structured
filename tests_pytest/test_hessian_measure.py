from utils import is_cuda_available
import pytest
import os

# Skip all tests in this file if no GPU is available
if not is_cuda_available():
    pytest.skip("Skipping tests because no GPU is available.", allow_module_level=True)

import pyUAMMD
import numpy as np

atol = 1e-5
rtol = 1e-4


def create_state(pos):
    state = {"labels": ["id", "position"], "data": list(enumerate(pos))}

    return state


def create_state_and_topology_bond1_fixedharmonic():
    pos0 = [-1.1, 2.2, -3.5]
    posL = [-0.6, 1.6, -2]
    k = 123
    r0 = 0.7

    state = create_state([pos0])
    topology = {
        "forceField": {
            "Bond": {
                "type": ["Bond1", "FixedHarmonic"],
                "labels": ["id_i", "K", "r0", "position"],
                "data": [[0, k, r0, posL]],
                "parameters": {},
            }
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_bond2_harmonic():
    pos0 = [-1.1, -5.5, 2]
    pos1 = [4.5, 2.2, -4.3]
    k = 123
    r0 = 2.13

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "Bond": {
                "type": ["Bond2", "Harmonic"],
                "labels": ["id_i", "id_j", "K", "r0"],
                "data": [[0, 1, k, r0]],
                "parameters": {},
            }
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_bond2_harmonicanisotropic():
    pos0 = [-1.1, 2.2, -3.5]
    pos1 = [1, 4, 1.3]
    k = [123, 5, 321]
    r0 = [0.3, 1.5, 2.2]

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "Bond": {
                "type": ["Bond2", "HarmonicAnisotropic"],
                "labels": ["id_i", "id_j", "K", "r0"],
                "data": [[0, 1, k, r0]],
                "parameters": {},
            }
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_bond2_lennardjones():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    epsilon = 3.1
    sigma = 2.2

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "Bond": {
                "type": ["Bond2", "LennardJonesType1"],
                "labels": ["id_i", "id_j", "epsilon", "sigma"],
                "data": [[0, 1, epsilon, sigma]],
                "parameters": {},
            }
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_bond2_steric():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    epsilon = 3.1
    sigma = 2.2

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "Bond": {
                "type": ["Bond2", "Steric6"],
                "labels": ["id_i", "id_j", "epsilon", "sigma"],
                "data": [[0, 1, epsilon, sigma]],
                "parameters": {},
            }
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_bond2_fene():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    r0 = 1.2
    R0 = 2.3
    k = 102

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "Bond": {
                "type": ["Bond2", "Fene"],
                "labels": ["id_i", "id_j", "K", "r0", "R0"],
                "data": [[0, 1, k, r0, R0]],
                "parameters": {},
            }
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_clashed():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1.8, 2.6, -2.1]
    gamma = 1.2
    lambd = 2.3

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "Clashed"],
                "parameters": {"lambda": lambd, "gamma": gamma, "condition": "all"},
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_debyehuckel():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1.8, 2.6, -2.1]
    epsilon = 0.94
    lambd = 3.42

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "DebyeHuckel"],
                "parameters": {
                    "dielectricConstant": epsilon,
                    "debyeLength": lambd,
                    "condition": "all",
                    "cutOffFactor": 2,
                },
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_dlvo():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    eps = 2.94
    debLength = 7.33
    dielConst = 3.12
    sig = 4.75

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "DLVOType1"],
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", eps, sig]],
                "parameters": {
                    "dielectricConstant": dielConst,
                    "debyeLength": debLength,
                    "condition": "all",
                    "cutOffNPFactor": 2 ** (1.0 / 6.0),
                    "cutOffDHFactor": 2 ** (1.0 / 6.0),
                },
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_lennardjonestype1():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1.8, 2.6, -2.1]
    epsilon = -2.33
    sigma = 4.75

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "GeneralLennardJonesType1"],
                "parameters": {"condition": "all", "cutOffFactor": 2 ** (1.0 / 6)},
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_lennardjonestype2():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1.8, 2.6, -2.1]
    epsilon = -2.33
    sigma = 4.75

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "GeneralLennardJonesType2"],
                "parameters": {"condition": "all", "cutOffFactor": 2 ** (1.0 / 6)},
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_lennardjonestype3():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1.8, 2.6, -2.1]
    epsilon = -2.33
    sigma = 4.75

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "GeneralLennardJonesType3"],
                "parameters": {"condition": "all", "cutOffFactor": 2 ** (1.0 / 6)},
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_kimhummer():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [0.5, 4, 1.4]
    epsilon = 2.5
    sigma = 3
    dielConst = 3.94
    debyeLen = 3.42

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "KimHummer"],
                "parameters": {
                    "condition": "all",
                    "cutOffNPFactor": 2,
                    "cutOffDHFactor": 2,
                    "dielectricConstant": dielConst,
                    "debyeLength": debyeLen,
                },
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_pairlennardjonestype1():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    epsilon = 3
    sigma = 4.75

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "LennardJonesType1"],
                "parameters": {"condition": "all", "cutOffFactor": 2 ** (1.0 / 6)},
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_pairlennardjonestype2():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    epsilon = 3
    sigma = 4.75

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "GeneralLennardJonesType2"],
                "parameters": {"condition": "all", "cutOffFactor": 2 ** (1.0 / 6)},
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_pairlennardjonestype3():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    epsilon = 3
    sigma = 4.75

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "LennardJonesType3"],
                "parameters": {"condition": "all", "cutOffFactor": 2 ** (1.0 / 6)},
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_pairsteric12():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    epsilon = 3
    sigma = 4.75

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "Steric12"],
                "parameters": {"condition": "all", "cutOffFactor": 2 ** (1.0 / 6)},
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_pairsteric6():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    epsilon = 3
    sigma = 4.75

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "Steric6"],
                "parameters": {"condition": "all", "cutOffFactor": 2 ** (1.0 / 6)},
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_splitlennardjones():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    epsilon = 3
    sigma = 4.75
    eps_r = 3.1
    eps_a = 1.8

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "SplitLennardJones"],
                "parameters": {
                    "condition": "all",
                    "cutOffFactor": 2 ** (1.0 / 6),
                    "epsilon_r": eps_r,
                    "epsilon_a": eps_a,
                },
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_wcatype1():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    epsilon = 3
    sigma = 4.75

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "WCAType1"],
                "parameters": {"condition": "all", "cutOffFactor": 2 ** (1.0 / 6)},
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_wcatype2():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    epsilon = 3
    sigma = 4.75

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "WCAType2"],
                "parameters": {"condition": "all", "cutOffFactor": 2 ** (1.0 / 6)},
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_nonbonded_wcatype3():
    pos0 = [1.5, 2.2, -3.3]
    pos1 = [1, 4, 1]
    epsilon = 3
    sigma = 4.75

    state = create_state([pos0, pos1])
    topology = {
        "forceField": {
            "nl": {
                "type": ["VerletConditionalListSet", "all"],
                "parameters": {"cutOffVerletFactor": 1.2},
            },
            "Bond": {
                "type": ["NonBonded", "WCAType3"],
                "parameters": {"condition": "all", "cutOffFactor": 2 ** (1.0 / 6)},
                "labels": ["name_i", "name_j", "epsilon", "sigma"],
                "data": [["A", "A", epsilon, sigma]],
            },
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"]],
        },
    }

    return state, topology


def create_state_and_topology_bond3_harmonic():
    pos0 = [-1.1, -5.5, 2]
    pos1 = [4.5, 2.2, -4.3]
    pos2 = [3.5, 6.1, 3.3]

    k = 123
    theta0 = 1.23
    state = create_state([pos0, pos1, pos2])

    topology = {
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
    }

    return state, topology


def create_state_and_topology_bond4_dihedral():

    pos0 = [1.5, 2.1, 4.342]
    pos1 = [0.1231, 2.12, -1.65]
    pos2 = [0.65, 2, 2]
    pos3 = [1, 5, 2.1]

    k = 123
    phi0 = 1.321
    n = 3

    state = create_state([pos0, pos1, pos2, pos3])

    topology = {
        "forceField": {
            "Bond": {
                "type": ["Bond4", "Dihedral"],
                "labels": ["id_i", "id_j", "id_k", "id_l", "n", "K", "phi0"],
                "data": [[0, 1, 2, 3, k, n, phi0]],
                "parameters": {},
            }
        },
        "structure": {
            "labels": ["id", "type"],
            "data": [[0, "A"], [1, "A"], [2, "A"], [3, "A"]],
        },
    }

    return state, topology


INTERACTOR_BUILDERS = {
    "Bond1_FixedHarmonic": create_state_and_topology_bond1_fixedharmonic,
    "Bond2_Harmonic": create_state_and_topology_bond2_harmonic,
    "Bond2_HarmonicAnisotropic": create_state_and_topology_bond2_harmonicanisotropic,
    "Bond2_LennardJones": create_state_and_topology_bond2_lennardjones,
    "Bond2_Steric": create_state_and_topology_bond2_steric,
    "Bond2_FENE": create_state_and_topology_bond2_fene,
    "NonBonded_Clashed": create_state_and_topology_nonbonded_clashed,
    "NonBonded_DebyeHuckel": create_state_and_topology_nonbonded_debyehuckel,
    "NonBonded_DLVO": create_state_and_topology_nonbonded_dlvo,
    "NonBonded_LennardJonesType1": create_state_and_topology_nonbonded_lennardjonestype1,
    "NonBonded_LennardJonesType2": create_state_and_topology_nonbonded_lennardjonestype2,
    "NonBonded_LennardJonesType3": create_state_and_topology_nonbonded_lennardjonestype3,
    "NonBonded_KimHummer": create_state_and_topology_nonbonded_kimhummer,
    "NonBonded_PairLennardJonesType1": create_state_and_topology_nonbonded_pairlennardjonestype1,
    "NonBonded_PairLennardJonesType2": create_state_and_topology_nonbonded_pairlennardjonestype2,
    "NonBonded_PairLennardJonesType3": create_state_and_topology_nonbonded_pairlennardjonestype3,
    "NonBonded_PairSteric12": create_state_and_topology_nonbonded_pairsteric12,
    "NonBonded_PairSteric6": create_state_and_topology_nonbonded_pairsteric6,
    "NonBonded_SplitLennardJones": create_state_and_topology_nonbonded_splitlennardjones,
    "NonBonded_WCAType1": create_state_and_topology_nonbonded_wcatype1,
    "NonBonded_WCAType2": create_state_and_topology_nonbonded_wcatype2,
    "NonBonded_WCAType3": create_state_and_topology_nonbonded_wcatype3,
    "Bond3_Harmonic": create_state_and_topology_bond3_harmonic,
    "Bond4_Dihedral": create_state_and_topology_bond4_dihedral,
}


def create_base_simulation(interactor):
    temperature = 0
    L = 20
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
                "data": [["A", 1, 0, 1]],
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
    }

    if interactor in INTERACTOR_BUILDERS:
        state, topology = INTERACTOR_BUILDERS[interactor]()
    else:
        raise KeyError(f"Builder function not defined for interactor: {interactor}")

    base_simulation["state"] = state
    base_simulation["topology"] = topology
    return base_simulation


def read_hessian_file(file_path):
    hessian_f = np.atleast_2d(np.loadtxt(file_path))
    # Hessian file has shape (npairs, 11), first two columns are the pair indices
    assert (
        hessian_f.shape[1] == 11
    ), f"Hessian file has unexpected shape {hessian_f.shape}"
    # Transform to a (n, n, 3, 3) array
    n = int(np.sqrt(hessian_f.shape[0]))
    assert (
        n * n == hessian_f.shape[0]
    ), f"Unexpected number of pairs {hessian_f.shape[0]}"
    i = hessian_f[:, 0].astype(int)
    j = hessian_f[:, 1].astype(int)
    matrices = hessian_f[:, 2:].reshape(-1, 3, 3)
    hessian = np.empty((n, n, 3, 3))
    hessian[i, j] = matrices
    return hessian


@pytest.mark.parametrize("interactor", INTERACTOR_BUILDERS.keys())
@pytest.mark.parametrize("hessian_mode", ["Analytical", "Numerical"])
def test_is_hessian_symmetric(tmp_path, hessian_mode, interactor):
    """Tests that the HessianMeasure is symmetric."""

    hessian_path = os.path.join(tmp_path, "Hessian.out")
    base_simulation = create_base_simulation(interactor)

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
    hessian = read_hessian_file(hessian_path)
    # Check symmetry, hessian[i,j] == hessian[j,i].T
    assert np.any(hessian != 0), "Hessian is entirely zero"

    hessian_t = np.transpose(hessian, (1, 0, 3, 2))
    max_diff = np.max(np.abs(hessian - hessian_t) / np.abs(np.mean(hessian)))
    assert np.allclose(
        hessian, hessian_t, rtol=rtol, atol=atol
    ), f"The HessianMeasure is not symmetric. Max diff: {max_diff}"


@pytest.mark.parametrize("interactor", INTERACTOR_BUILDERS.keys())
@pytest.mark.parametrize("hessian_mode", ["Analytical", "Numerical"])
def test_hessian_row_sums_are_zero(tmp_path, hessian_mode, interactor):
    """Tests that the sum of all elements in each row of the Hessian is zero.

    This verifies that sum_j^a H_ij^{ab} = 0 for each particle i and coordinate α,
    which is a consequence of Newton's First Law (translational invariance)."""

    if "Bond1" in interactor:
        pytest.skip("Bond1 does not satisfy Newton's First Law")

    hessian_path = os.path.join(tmp_path, "Hessian.out")
    base_simulation = create_base_simulation(interactor)

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
    hessian = read_hessian_file(hessian_path)
    assert np.any(hessian != 0), "Hessian is entirely zero"

    tol = 1e-8 if hessian_mode == "Analytical" else 1e-4

    nRows = hessian.shape[0]  # Three times the number of particles
    for i in range(nRows):
        # Sum over all j (axis=1) and over block columns (axis=3)
        net_sum = np.sum(hessian[i, :, :, :], axis=(0, 2))  # shape (3,)
        assert np.allclose(
            net_sum, 0.0, atol=tol
        ), f"Particle {i} violates translational invariance: sum = {net_derivative}"


@pytest.mark.parametrize("interactor", INTERACTOR_BUILDERS.keys())
def test_hessian_self_consistent(tmp_path, interactor):
    """Tests that the HessianMeasure is self-consistent between the analytical and numerical methods."""

    analytical_hessian_path = os.path.join(tmp_path, "HessianAnalytical.out")
    numerical_hessian_path = os.path.join(tmp_path, "HessianNumerical.out")
    base_simulation = create_base_simulation(interactor)
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
