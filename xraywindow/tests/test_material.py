from xraywindow.material import Material, import_materials

def test_xray_data_import():
    """Test import of xray data. xray_data should not be None."""
    mat = Material("silicon", 200e9, 7e9, 0.2, 20e-9)
    mat.get_xray_data()
    assert mat.xray_data is not None

def test_import_materials():
    """Import materials.yml and confirm TestMat matches expections."""
    materials = import_materials()
    assert materials['TestMat'].modulus == 123.456e9