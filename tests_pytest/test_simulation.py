import pyUAMMD


def test_can_construct():
    sim = pyUAMMD.simulation()


def test_can_construct_with_dict():

    dict_data = {
        "system": {
            "info": {
                "type": ["Simulation", "Information"],
                "parameters": {"name": "testConstruct"},
            }
        }
    }

    sim = pyUAMMD.simulation(dict_data)
    assert sim is not None
