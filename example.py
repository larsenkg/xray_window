import matplotlib.pyplot as plt
from scipy.optimize        import differential_evolution, Bounds
from xraywindow.plotting   import plot_transmission
from xraywindow.xray_data  import XRayData
from xraywindow.mechanical import MechanicalWindow, BeamLayer, RectangularMembraneLayer
from xraywindow.material   import import_materials


materials = import_materials()
polymer   = materials["polymer"]
silicon   = materials["silicon"]
aluminum  = materials["aluminum"]
boron     = materials["boron"]
SiN       = materials["SiNx"]

# Energies for which transmission will be optimized
opt_energies = [54.3, 108.5, 183.3, 277, 392.4, 524.9, 676.8, 1041, 1740]

def light_block(window):
    """Visible light blocking and charge dissipation layer."""
    window.add_layer(RectangularMembraneLayer("Light Block", material=aluminum, thickness=30e-9))

    return window

def gas_barrier(window):
    """Gas barrier layer."""
    window.add_layer(RectangularMembraneLayer("Gas Barrier", material=boron, thickness=20e-9))

    return window

def common_layers(window):
    """Both light blocking and gas barrier layers."""
    window = light_block(window)
    window = gas_barrier(window)
    
    return window

def get_ap3_df():
    """Create AP3 model and return Pandas dataframe."""
    prim_spacing = 190e-6
    prim_width   = 60e-6
    prim_length  = 10.2e-3
    prim_thick   = 380e-6

    tert_width   = prim_spacing
    tert_thick   = 300e-9

    prim        = BeamLayer("Primary", silicon, prim_spacing, prim_width, prim_length, prim_thick)
    tert        = RectangularMembraneLayer("Tertiary", polymer, tert_width, tert_thick)
    light_block = RectangularMembraneLayer("Light Block", material=aluminum, thickness=30e-9)
    gas_barrier = RectangularMembraneLayer("Gas Barrier", material=boron,    thickness=20e-9)
    
    mech_win = MechanicalWindow()
    mech_win.add_layer(prim)
    mech_win.add_layer(tert)
    mech_win.add_layer(light_block)
    mech_win.add_layer(gas_barrier)
    
    return mech_win.to_xray_window().spectrum().df()

def make_three_layer_window(
    p,
    prim_width  = 60e-6,
    prim_length = 10.2e-3,
    prim_thick  = 380e-6,
    sec_thick   = 45e-6,
    prim_mat    = silicon,
    sec_mat     = silicon,
    tert_mat    = polymer,
    use_light   = True,
    use_gas     = True
):
    """Return the MechanicalWindow consisting of the three specified layers."""
    prim_spacing, sec_spacing, sec_width = p

    sec_length = prim_spacing
    tert_width = sec_spacing

    prim = BeamLayer("Primary", prim_mat, prim_spacing, prim_width, prim_length, prim_thick)
    sec  = BeamLayer("Secondary", sec_mat, sec_spacing,  sec_width,  sec_length, sec_thick)
    tert = RectangularMembraneLayer("Membrane", tert_mat, width=tert_width)
    
    tert.thickness = tert.calc_min_thickness()*1.01  # Increase by one percent to not straddle the minimum
    tert.calc_stress()
    
    if prim.failure() or sec.failure() or tert.failure():
        return 0

    mech_win = MechanicalWindow()
    mech_win.add_layer(prim)
    mech_win.add_layer(sec)
    mech_win.add_layer(tert)
    
    # Add light blocking and gas barrier layers
    if use_light:
        mech_win = light_block(mech_win)
        
    if use_gas:
        mech_win = gas_barrier(mech_win)
    
    return mech_win

def make_two_layer_window(
    p,
    prim_width  = 60e-6,
    prim_length = 10.2e-3,
    prim_thick  = 380e-6,
    prim_mat    = silicon,
    tert_mat    = polymer,
    use_light   = True,
    use_gas     = True
):
    """Return the MechanicalWindow consisting of the two specified layers."""
    (prim_spacing,) = p

    tert_width = prim_spacing

    prim = BeamLayer("Primary", prim_mat, prim_spacing, prim_width, prim_length, prim_thick)
    tert = RectangularMembraneLayer("Membrane", tert_mat, width=tert_width)
    
    tert.thickness = tert.calc_min_thickness()*1.01  # Increase by one percent to not straddle the minimum
    tert.calc_stress()
    
    if prim.failure() or tert.failure():
        return 0

    mech_win = MechanicalWindow()
    mech_win.add_layer(prim)
    mech_win.add_layer(tert)
    
    # Add light blocking and gas barrier layers
    if use_light:
        mech_win = light_block(mech_win)
        
    if use_gas:
        mech_win = gas_barrier(mech_win)
    
    return mech_win



# These `trans_*` functions are the evaluation function passed to `differential_evolution`.

def trans_two_layer_polymer(p):
    """Calculate the transmission through a two-layer window. Returns 0 if parameters yield a 
    window that mechanically fails (cannot withstand pressure)."""
    # Returns 0 (not MechanicalWindow) if stress in any layer exceeds material strength
    win = make_two_layer_window(p)
    
    if not isinstance(win, MechanicalWindow):
        return 0
    
    return -win.to_xray_window().spectrum().integrate(energies=opt_energies)

def trans_two_layer_polymer_54(p):
    """Calculate the transmission at 54 eV through a two-layer window. Returns 0 if parameters 
    yield a window that mechanically fails (cannot withstand pressure)."""
    win = make_two_layer_window(p)
    
    if not isinstance(win, MechanicalWindow):
        return 0
    
    return -win.to_xray_window().spectrum().integrate(energies=[54.3])

def trans_three_layer_polymer(p):
    """Calculate the transmission through a three-layer window. Returns 0 if parameters yield a 
    window that mechanically fails (cannot withstand pressure)."""
    win = make_three_layer_window(p)
    
    if not isinstance(win, MechanicalWindow):
        return 0
    
    return -win.to_xray_window().spectrum().integrate(energies=opt_energies)


def trans_two_layer_si3n4(p):
    """Calculate the transmission through a two-layer SiN window. Returns 0 if parameters yield a 
    window that mechanically fails (cannot withstand pressure)."""
    win = make_two_layer_window(p, tert_mat=SiN, use_gas=False)
    
    if not isinstance(win, MechanicalWindow):
        return 0
    
    return -win.to_xray_window().spectrum().integrate(energies=opt_energies)

def trans_three_layer_si3n4(p):
    """Calculate the transmission through a three-layer SiN window. Returns 0 if parameters yield a 
    window that mechanically fails (cannot withstand pressure)."""
    win = make_three_layer_window(p, tert_mat=SiN, use_gas=False)
    
    if not isinstance(win, MechanicalWindow):
        return 0
    
    return -win.to_xray_window().spectrum().integrate(energies=opt_energies)


# Polymer-based window like AP3
lb   = [100e-6]
ub   = [2000e-6]
bnds = Bounds(lb, ub)

two_layer_polymer = differential_evolution(trans_two_layer_polymer, bnds)

# Polymer-based window like AP3, but with extra support layer
lb   = [100e-6, 1e-6, 5e-6]
ub   = [2000e-6, 200e-6, 30e-6]
bnds = Bounds(lb, ub)

three_layer_polymer = differential_evolution(trans_three_layer_polymer, bnds)

# Two layer Si3N4-based window
lb   = [100e-6]
ub   = [2000e-6]
bnds = Bounds(lb, ub)

two_layer_si3n4 = differential_evolution(trans_two_layer_si3n4, bnds)

# Three layer Si3N4-based window
lb   = [100e-6, 1e-6, 5e-6]
ub   = [2000e-6, 200e-6, 30e-6]
bnds = Bounds(lb, ub)

three_layer_si3n4 = differential_evolution(trans_three_layer_si3n4, bnds)

print(two_layer_polymer)
print(three_layer_polymer)
print(two_layer_polymer)
print(three_layer_polymer)


light_gas_only_xray_win       = common_layers(MechanicalWindow()).to_xray_window("light/gas barrier")
two_layer_polymer_xray_win    = make_two_layer_window(two_layer_polymer.x).to_xray_window("2-layer polymer")
#two_layer_polymer_54_xray_win = make_two_layer_window(two_layer_polymer_54.x).to_xray_window("2-layer polymer (54 eV)")
three_layer_polymer_xray_win  = make_three_layer_window(three_layer_polymer.x).to_xray_window("3-layer polymer")
two_layer_si3n4_xray_win      = make_two_layer_window(two_layer_si3n4.x, tert_mat=SiN, use_gas=False).to_xray_window("2-layer Si3N4")
three_layer_si3n4_xray_win    = make_three_layer_window(three_layer_si3n4.x, tert_mat=SiN, use_gas=False).to_xray_window("3-layer Si3N4")

ap3_sin = make_two_layer_window([200e-6], tert_mat=SiN, use_gas=False).to_xray_window("AP3 SiN")

print(light_gas_only_xray_win)
print(two_layer_polymer_xray_win)
#print(two_layer_polymer_54_xray_win)
print(three_layer_polymer_xray_win)
print(two_layer_si3n4_xray_win)
print(three_layer_si3n4_xray_win)
print(ap3_sin)

ap3_df = get_ap3_df()

plt.close("all")

fig, ax = plt.subplots(figsize=(8,3))
ax = plot_transmission(
    data = [
        ap3_df, 
        two_layer_polymer_xray_win.spectrum().df(),
#        two_layer_polymer_54_xray_win.spectrum().df(),
        three_layer_polymer_xray_win.spectrum().df(), 
        two_layer_si3n4_xray_win.spectrum().df(), 
        three_layer_si3n4_xray_win.spectrum().df(),
        ap3_sin.spectrum().df()
    ], 
    labels = [
        "AP3", 
        "2-layer polymer", 
#        "2-layer polymer (54 eV)", 
        "3-layer polymer", 
        "2-layer Si$_3$N$_4$", 
        "3-layer Si$_3$N$_4$",
        "AP3 SiN"
    ],  
    use_log    = True,
    add_lines  = True,
    min_energy = 10,
    max_energy = 10000,
    ax         = ax
)

plt.show()

print("Polymer x-ray data ", polymer.xray_data)