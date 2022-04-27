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

    def transmission(self, energy, thickness):
        '''Calculate the transmission of the material with `thickness` and `energy`.'''
        #Interpolation
        # Find closest value that is smaller than specified energy
        #closest = min(self.energies, key=lambda x:abs(x-energy))
        #index = self.energies.index(closest)
        index = 0
        for ind, E in enumerate(self.energies):
            if E > energy:
                index = ind - 1
                break

        T1 = self.transmissions[index]
        T2 = self.transmissions[index + 1]
        E1 = self.energies[index]
        E2 = self.energies[index + 1]

        T = T1 + (energy - E1) * (T2 - T1) / (E2 - E1)

        return T ** (thickness / self.thickness)

    def transmission_all(self, energy, thickness):
        return self.interp_trans(energy) ** (thickness / self.thickness)

def import_xray_data(material_name, xray_data_dir = None):
    '''Load the x-ray data from the appropriate .dat file. Expects xray data to be under data/xray.
    TO DO: Describe x-ray data format.'''
    if xray_data_dir is None:
        xray_data_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), os.pardir, "data", "xray")

    filename = os.path.join(xray_data_dir, material_name + ".dat")
    print("Trying to load " + filename)
    # Get our arrays ready
    energies, transmissions = [], []
    try:
        with open (filename) as datfile:
            first = datfile.readline()

            try:
                thickness = float(first.split(" ")[3].split("=")[1])
            except:
                thickness = 0.005 # Set default to 0.005 micron

            try:
                density = float(first.split(" ")[2].split("=")[1])
            except:
                density = 1 # Set default to 1 g/cm^3

            for line in datfile:
                e = line[0:13].strip()
                t = line[13:].strip()
                try:
                    e = float(e)
                    t = float(t)
                except ValueError:
                    continue
                energies.append(e)
                transmissions.append(t)
            return energies, transmissions, thickness, density
    except EnvironmentError:
        print("Could not read file.")

def import_xray_data_csv(material_name):
    '''Load the x-ray data from the appropriate .csv file.'''
    filename = os.path.join("data", "xray", material_name + ".csv")
    df       = pd.read_csv(filename)

    return df["Energy"].values, df["Transmission"].values, df["Thickness"].values, df["Density"].values

    