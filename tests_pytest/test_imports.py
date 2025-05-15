def test_import_pyuammd():
    import pyUAMMD


def test_has_simulation():
    # Check pyUAMMD has a simulation module
    import pyUAMMD

    assert hasattr(pyUAMMD, "simulation")
