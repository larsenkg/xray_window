import os
import yaml

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

def import_materials(material_data_dir = None, material_filename = "materials.yml"):
    '''Load data from the materials.yml file and create new Material 
    objects for each loaded material.'''
    if material_data_dir is None:
        material_data_dir = os.path.join("data")

    filename  = os.path.join(material_data_dir, material_filename)
    materials = {}

    with open(filename, encoding="utf-8") as f:
        yml = yaml.load(f, yaml.SafeLoader)
        for mat in yml:
            materials[mat['name']] = Material(mat['name'], mat['youngs_modulus'], mat['ultimate_stress'], mat['poisson_ratio'], mat['min_thickness'])
    
    return materials
