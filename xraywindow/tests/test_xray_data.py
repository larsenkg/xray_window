from xraywindow.xray_data import *

def test_import_xray_data_csv():
    mat = "test_xray_data"
    energies, transmissions, thickness, density = import_xray_data_csv(mat)
    assert energies[0] == 10 # First energy should be 10 in 'test_xray_data'