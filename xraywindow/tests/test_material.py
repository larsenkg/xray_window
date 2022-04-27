from xraywindow.material import Material

def test_xray_data_import():
    """Test import of xray data. xray_data should not be None."""
    mat = Material("silicon", 200e9, 7e9, 0.2, 20e-9)
    mat.get_xray_data()
    assert mat.xray_data is not None