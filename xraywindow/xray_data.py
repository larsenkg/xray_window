import os
from scipy.interpolate import interp1d
import pandas as pd

class XRayData:
    '''This object holds the x-ray transmission data for a given material. It 
    can calculate the transmission at a specified energy and material thickness.'''
    def __init__(self, material_name):
        self.energies, self.transmissions, self.thickness, self.density = import_xray_data_csv(material_name)
        self.material_name = material_name
        self.interp_trans = interp1d(self.energies, self.transmissions)

    # def transmission(self, energy, thickness):
    #     '''Calculate the transmission of the material with `thickness` and `energy`.'''
    #     #Interpolation
    #     # Find closest value that is smaller than specified energy
    #     #closest = min(self.energies, key=lambda x:abs(x-energy))
    #     #index = self.energies.index(closest)
    #     index = 0
    #     for ind, E in enumerate(self.energies):
    #         if E > energy:
    #             index = ind - 1
    #             break

    #     T1 = self.transmissions[index]
    #     T2 = self.transmissions[index + 1]
    #     E1 = self.energies[index]
    #     E2 = self.energies[index + 1]

    #     T = T1 + (energy - E1) * (T2 - T1) / (E2 - E1)

    #     return T ** (thickness / self.thickness)

    def transmission(self, energy, thickness):
        '''Use the interpolated transmission values to estimate the transmission at a certain energy
        through a given material thickness.
        TODO: Use this interpolated method instead of what I currently do in `transmission()`'''
        return self.interp_trans(energy) ** (thickness / self.thickness)

def import_xray_data_csv(material_name, xray_data_dir = None):
    '''Load the x-ray data from the appropriate .csv file.
    Expected columns: Energy, Transmission, Density, Thickness
    Energy should be in eV, Transmission should be a number less than 1,
    Thickness should be in meters, and density should be in kg/m^3.'''
    if xray_data_dir is None:
        xray_data_dir = os.path.join("data", "xray")

    filename = os.path.join(xray_data_dir, material_name + ".csv")
    df       = pd.read_csv(filename)

    return df["Energy"].values, df["Transmission"].values, df["Thickness"].values[0], df["Density"].values[0]

    