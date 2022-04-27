import numpy as np
from xraywindow.transmission import XRayWindow, XRayWindowLayer
from math import sqrt

ATM_PRESSURE  = 101.3e3         # Pa
TEST_PRESSURE = 2*ATM_PRESSURE

class MechanicalWindowLayer:
    '''Defines mechanical support layer of x-ray detector window.'''
    def __init__(self, material, name):
        self.max_stress  = np.nan
        self.fail_stress = material.stress
        self.modulus     = material.modulus
        self.poisson     = material.poisson
        self.material    = material
        self.name        = name
        
    def xray_thickness(self):
        '''Get the thickness that x-rays must pass through. This method allows child classes to use 
        semantically named properties (e.g. thicknes or height).'''
        return 0
    
    def calc_stress(self):
        '''Should set and return self.max_stress'''
        return 0
    
    def open_area(self):
        '''Open area is used for x-ray transmission. Should be a fraction from 0 to 1 
        that indicates what percent of collimated light would not hit a feature in the layer.'''
        return 1
    
    def failure(self):
        '''If max_stress is greater than fail_stress, then the material has failed.'''
        if np.isnan(self.max_stress):
            raise ValueError("self.max_stress not set. Be sure to include all necessary parameters.")
        return self.max_stress > self.fail_stress
    
    def to_xray_window_layer(self):
        return XRayWindowLayer(
            self.name, 
            self.material.get_xray_data(), 
            self.xray_thickness(), 
            self.open_area(),
            self,
            self.material.name,
        )
            
    
class BeamLayer(MechanicalWindowLayer):
    '''A support structure layer consiting of an Euler-Bernoulli fixed-fixed beam.'''
    def __init__(self, name, material, spacing=np.nan, width=np.nan, length=np.nan, height=np.nan, pressure = TEST_PRESSURE):
        MechanicalWindowLayer.__init__(self, material=material, name=name)
        
        self.spacing  = spacing
        self.width    = width
        self.length   = length
        self.height   = height
        self.pressure = pressure
        
        self.slenderness_ratio()
        
        self.calc_stress()
        
    def xray_thickness(self):
        return self.height
    
    def calc_stress(self):
        '''
        $\sigma_max = Mc/I$
        $M = pL^2/12$
        $c = h/2$
        $I = wh^3/12$
        $\sigma_max = pL^2/(2wh^2)$
        $p = linear pressure [N/m]$
        $P = pressure [N/m^3]$
        $s = spacing between ribs$
        $A = (w + s) * L$
        $P = F/A$
        $p = F/L = P*A/L = P*(w + s)*L/L = P*(w + s)$
        '''
        dist_load       = (self.spacing + self.width) * self.pressure
        self.max_stress = dist_load * (self.length**2)/(2.0*self.width*self.height**2)
        return self.max_stress
    
    def calc_max_deflection(self):
        dist_load           = (self.spacing + self.width) * self.pressure
        self.max_deflection = dist_load * self.length**4/(32 * self.modulus * self.width * self.height**3)
        return self.max_deflection
    
    def open_area(self):
        return self.spacing / (self.spacing + self.width)
    
    def slenderness_ratio(self):
        '''This was a test and is not correct/useful.'''
        I = self.width*self.height**3/12
        A = self.width*self.height
        r = sqrt(I/A)   # radius of gyration: TO DO-> Does width not matter? Why?
        self.slenderness_ratio = self.length/r
        
        return self.slenderness_ratio
    
    def calc_max_spacing(self):
        '''Calculate the maximum spacing.'''
        s = self.fail_stress
        w = self.width
        h = self.height
        p = self.pressure
        L = self.length
        
        spacing = 2*s*w*h**2/(p*L**2)
        
        return spacing
    
    def __repr__(self):
        msg  = f"{self.name}: BeamLayer | OA {self.open_area()*100:4.1f}%\n"
        msg += f"  Spacing:  \t{self.spacing * 1e6:7.1f} µm\n"
        msg += f"  Width:    \t{self.width   * 1e6:7.1f} µm\n"
        msg += f"  Length:   \t{self.length  * 1e6:7.1f} µm\n"
        msg += f"  Height:   \t{self.height  * 1e6:7.1f} µm\n"
        msg += f"  Material: \t{self.material.name}\n"
        
        return msg
        
class RectangularMembraneLayer(MechanicalWindowLayer):
    '''A support structure layer consiting of a rectangular membrane. Assumes length is more than
    six times the width (i.e. an infitely long membrane).'''
    def __init__(self, name, material, width=np.nan, thickness=np.nan, pressure = TEST_PRESSURE):
        MechanicalWindowLayer.__init__(self, material=material, name=name)
        
        self.width     = width
        self.thickness = thickness
        self.pressure  = pressure
        
        self.calc_stress()
        
    def xray_thickness(self):
        return self.thickness
        
    def calc_stress(self):
        E = self.modulus
        v = self.poisson
        t = self.thickness
        p = self.pressure
        a = self.width / 2.0
        
        self.max_stress = (E * p**2 * a**2 /(6.0 * t**2 * (1 - v**2)))**(1/3.0)
        return self.max_stress
    
    def open_area(self):
        return 0
    
    def calc_max_width(self):
        '''Calculate the maximum width a long rectangular membrane can be and still
        withstand the specified pressure.'''
        E = self.modulus     #* (1 - porosity)
        s = self.fail_stress #* (1 - porosity)
        v = self.poisson
        t = self.thickness
        p = self.pressure

        w = 2 * 2.4495 * t / p * sqrt(((1-v**2) * s**3) / E)
        
        if np.isnan(w):
            raise ValueError("Missing parameter, possibly self.thickness.")
            
        return w
    
    def calc_min_thickness(self):
        '''Calculate the minimum required thickness for a rectangular membrane to 
        withstand the given pressure.'''
        E = self.modulus     #* (1 - porosity)
        s = self.fail_stress #* (1 - porosity)
        v = self.poisson
        p = self.pressure
        a = self.width / 2.0

        t = p * a / 2.4495 * sqrt(E / ((1-v**2) * s**3))
        
        if np.isnan(t):
            raise ValueError("Missing parameter, possibly self.width.")
            
        if t < self.material.min_thickness:
            t = self.material.min_thickness
            
        return t
    
    def __repr__(self):
        msg  = f"{self.name}: RectangularMembraneLayer | OA {self.open_area()*100:4.1f}%\n"
        msg += f"  Width:     \t{self.width     * 1e6:7.1f} µm\n"
        msg += f"  Thickness: \t{self.thickness * 1e6:7.1f} µm\n"
        msg += f"  Material:  \t{self.material.name}\n"
        
        return msg

class MechanicalWindow:
    '''Defines entire mechanical support structure of x-ray detector window.'''
    def __init__(self):
        self.layers = []
    
    def add_layer(self, layer:MechanicalWindowLayer):
        self.layers.append(layer)
        
    def to_xray_window(self, name=""):
        window = XRayWindow()
        window.name = name
        for layer in self.layers:
            window.add_layer(layer.to_xray_window_layer())
        return window
    
    def __repr__(self):
        #msg = f"Name: {self.name}\n"
        msg = ""
        for L in self.layers:
            msg += str(L) + "\n"
            
        return msg