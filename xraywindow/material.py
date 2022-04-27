from xraywindow.xray_data import XRayData

class Material:
    '''Contains material properties for a given material.'''
    def __init__(self, name, modulus, stress, poisson, min_thickness=0):
        self.name          = name
        self.modulus       = modulus
        self.stress        = stress
        self.poisson       = poisson
        self.min_thickness = min_thickness
        
        # X-ray data
        self.xray_data = None
        
    def biaxial_modulus(self):
        return self.modulus / (1 - self.poisson)
    
    def get_xray_data(self):
        if self.xray_data is None: 
            self.xray_data = XRayData(self.name)
        return self.xray_data