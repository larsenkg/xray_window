from xraywindow.mechanical import MechanicalWindow, BeamLayer, RectangularMembraneLayer
from xraywindow.material   import Material

def test_membrane_transmission():
    silicon    = Material("silicon", 150e9, 7000e6, 0.17)
    membrane   = RectangularMembraneLayer("Membrane", silicon, 100e-6, 100e-9)
    xray_layer = membrane.to_xray_window_layer()

    # Transmission should not be the same at vastly different energies
    assert xray_layer.transmission(50) != xray_layer.transmission(5000)