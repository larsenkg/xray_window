import pytest
from xraywindow.xray_data import *

def test_import_xray_data_csv():
    mat = "test_xray_data"
    energies, transmissions, thickness, density = import_xray_data_csv(mat)
    assert energies[0] == 10 # First energy should be 10 in 'test_xray_data'

def test_import_xray_data_failed():
    mat = "does_not_exist"
    with pytest.raises(FileNotFoundError):
        _ = import_xray_data_csv(mat)

def test_xray_data_transmission():
    mat = "test_xray_data"
    xd  = XRayData(mat)
    # Transmission at 10 eV should be 0.1 for 100 nm
    assert xd.transmission(10, 100e-9) == 0.1

def test_xray_data_transmission_interpolation():
    mat = "test_xray_data"
    xd  = XRayData(mat)
    # Transmission at 10.5 eV should be 0.105 for 100 nm
    assert xd.transmission(10.5, 100e-9) == pytest.approx(0.105)