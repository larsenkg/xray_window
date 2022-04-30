import numpy as np
import pandas as pd
import scipy.integrate as integrate
from scipy.interpolate import interp1d

#from xraywindow.mechanical import BeamLayer

class XRaySpectrum:
    '''X-ray spectrum information.'''
    def __init__(self, energy=None, transmission=None, spectrum=None):
        if spectrum is not None:
            energy       = spectrum[:, 0]
            transmission = spectrum[:, 1]
            
        self.energy       = np.array(energy)
        self.transmission = np.array(transmission)
        
        self.interp_trans = interp1d(self.energy, self.transmission)
        
        self.last_integration = None
        self.min_energy       = -np.inf
        self.max_energy       = np.inf
    
    def spectrum(self):
        return np.stack([self.energy, self.transmission]).T
    
    def df(self):
        df = pd.DataFrame(columns=['Energy', 'Transmission'], data=self.spectrum())
        return df
    
    def integrate(self, min_energy=-np.inf, max_energy=np.inf, energies=None):
        # Check if we need to integrate again
        #if self.last_integration is not None and self.min_energy == min_energy and self.max_energy == max_energy:
        #    return self.last_integration
        
        if energies is None:
            self.min_energy, self.max_energy = min_energy, max_energy
            e = self.energy[(self.energy >= min_energy) & (self.energy <= max_energy)]
            t = self.transmission[(self.energy >= min_energy) & (self.energy <= max_energy)]
            
            self.last_integration = integrate.simps(t, e)
        else: # Assume energies is a list
            e = energies
            t = self.interp_trans(e)
            
            self.last_integration = t.sum()   
        
        return self.last_integration

class XRayWindowLayer:
    '''Class that represents a layer in an x-ray window.'''
    def __init__(self, layer_name, xray_data, thickness, open_area, mech_layer, mat_name=""):
        self.layer_name = layer_name
        self.xray_data  = xray_data
        # TODO: Old method used microns for distance. We're using meters. Fix to be consistent (use meters).
        self.thickness  = thickness#*1e6
        self.open_area  = open_area
        self.mat_name   = mat_name
        self.mech_layer = mech_layer
        
    # def transmission(self, energy):
    #     '''Calculates sum of transmission through material and open area.'''
    #     return (1 - self.open_area) * self.xray_data.transmission(energy, self.thickness) + self.open_area

    def transmission(self, energy):
        '''Calculates sum of transmission through material and open area.'''
        # TODO: Vectorize 'energy'!
        return (1 - self.open_area) * self.xray_data.transmission(energy, self.thickness) + self.open_area
    
    def thick_str(self, t = None):
        if t is None:
            t   = self.thickness*1e6
        msg = ""
        
        if t < 1:
            msg = f"{t*1000:6.2f} nm"
        elif t > 1000:
            msg = f"{t/1000:6.2f} mm"
        else:
            msg = f"{t:6.2f} µm"
        return msg
    
    def __repr__(self):
        return f"{self.layer_name} ({self.mat_name}):\t{self.thick_str()} | OA {self.open_area*100:4.1f}%"

class XRayWindow:
    '''Class that represents a collection of layers that together make up an 
    x-ray window.'''
    def __init__(self):
        self.layers = []
        self.name = ""
        self.plot = True  # Set to false if this window shouldn't be plotted
        self.line_color = None # Line color for plot
        return
    
    def add_layer(self, layer):
        '''Add a layer to the window.'''
        self.layers.append(layer)
        
#     def add_layer(self, layer_name="", xray_data=None, thickness=0, open_area=0, layer=None, mat_name=""):
#         '''Add a layer to the window.'''
#         # If xray_data is just a string, as opposed to a XRayData object,
#         # then try and fetch the xray data
#         if xray_data is not None and isinstance(xray_data, str):
#             # TO DO: Add importer or get rid of importing
#             xray_data = importer.get_data(xray_data)
            
#         if not isinstance(layer, XRayWindowLayer):
#             layer = XRayWindowLayer(layer_name, xray_data, thickness, open_area, mat_name)
            
#         self.layers.append(layer)

    # def transmission(self, energy):
    #     '''Calculates the transmission of the x-ray window stack as a fraction.'''
    #     total = 1
    #     for layer in self.layers:
    #         total *= layer.transmission(energy)
    #     return total
    
    def transmission(self, energy=range(10, 10000)):
        '''Calculates the transmission of the x-ray window stack as a fraction.'''
        if isinstance(energy, int):
            total = 1
        else:
            total = np.ones(len(energy))
        
        for layer in self.layers:
            total *= layer.transmission(energy)
        return total

    def spectrum(self, energy=range(10,10000)):
        '''Return the transmission spectrum of the window. Energy should be in eV.'''
        #return energy,list(map(self.transmission, energy))
        #return list(zip(list(energy),list(map(self.transmission, energy))))
        #trans    = list(map(self.transmission, energy))
        trans = self.transmission(energy)
        return XRaySpectrum(list(energy), trans)
        #return np.stack([np.array(energy), trans]).T
        
    def __repr__(self):
        msg = f"Name: {self.name}\n"
        for L in self.layers:
            msg += str(L) + "\n"
            
        return msg
    
    # def draw_window(
    #     self,
    #     frame_width  = 2e-3,
    #     ax           = None
    # ):
    #     # Assumes that layer 0 is primary support structure and layer 1 is secondary support structure
    #     if isinstance(self.layers[0].mech_layer, BeamLayer):
    #         prim_layer = self.layers[0].mech_layer
    #     else:
    #         raise ValueError("First layer is not a BeamLayer.")
        
    #     if isinstance(self.layers[1].mech_layer, BeamLayer):
    #         sec_layer = self.layers[1].mech_layer
    #         draw_sec  = True
    #     else:
    #         draw_sec  = False
                

    #     if ax is None:
    #         ax  = plt.gca()

    #     prim_width   = prim_layer.width
    #     prim_spacing = prim_layer.spacing
    #     prim_length  = prim_layer.length
        
    #     if draw_sec:
    #         sec_width    = sec_layer.width
    #         sec_spacing  = sec_layer.spacing
            
    #     #prim_x = np.arange(prim_spacing, prim_length-prim_spacing, prim_spacing)
    #     prim_x = np.arange(prim_length/2, prim_length, prim_spacing)
    #     prim_x = np.hstack((np.flip(prim_length-prim_x[1:]), prim_x))

    #     if draw_sec:
    #         sec_y = np.arange(prim_length/2, prim_length, sec_spacing)
    #         sec_y = np.hstack((np.flip(prim_length-sec_y[1:]), sec_y))

    #     prim_support = [plt.Rectangle((x, 0), prim_width, prim_length, fc='black') for x in prim_x]
        
    #     if draw_sec:
    #         sec_support  = [plt.Rectangle((0, y), prim_length, sec_width, fc='black') for y in sec_y]
        
    #     # Add frame
    #     ax.add_patch(plt.Rectangle((0, 0), -frame_width, prim_length + frame_width, fc='black'))                       # Left
    #     ax.add_patch(plt.Rectangle((0, prim_length), prim_length + frame_width, frame_width, fc='black'))              # Top
    #     ax.add_patch(plt.Rectangle((prim_length, prim_length), frame_width, -(prim_length + frame_width), fc='black')) # Right
    #     ax.add_patch(plt.Rectangle((prim_length, 0), -(prim_length + frame_width), -frame_width, fc='black'))          # Bottom

    #     for p in prim_support:
    #         ax.add_patch(p)

    #     if draw_sec:
    #         for p in sec_support:
    #             ax.add_patch(p)

    #     ax.axis('scaled')
    #     ax.axis('off')

    #     return ax
    
    # def draw_window_circ(
    #     self,
    #     frame_width = 2e-3,
    #     ax          = None
    # ):
    #     # Assumes that layer 0 is primary support structure and layer 1 is secondary support structure
    #     if isinstance(self.layers[0].mech_layer, BeamLayer):
    #         prim_layer = self.layers[0].mech_layer
    #     else:
    #         raise ValueError("First layer is not a BeamLayer.")
        
    #     if isinstance(self.layers[1].mech_layer, BeamLayer):
    #         sec_layer = self.layers[1].mech_layer
    #         draw_sec  = True
    #     else:
    #         draw_sec  = False
                
    #     if ax is None:
    #         ax  = plt.gca()

    #     prim_width   = prim_layer.width
    #     prim_spacing = prim_layer.spacing
    #     prim_length  = prim_layer.length
        
    #     if draw_sec:
    #         sec_width    = sec_layer.width
    #         sec_spacing  = sec_layer.spacing

    #     #prim_x = np.arange(prim_spacing, prim_length-prim_spacing, prim_spacing)
    #     prim_x = np.arange(prim_spacing/2, prim_length/2, prim_spacing)
    #     prim_x = np.hstack((np.flip(-prim_x), prim_x))

    #     if draw_sec:
    #         sec_y = np.arange(sec_spacing/2, prim_length/2, sec_spacing)
    #         sec_y = np.hstack((np.flip(-sec_y), sec_y))

    #     prim_support = []
    #     sec_support  = []

    #     for x in prim_x:
    #         r     = prim_length/2 + frame_width/2
    #         chord = np.sqrt(r**2 - x**2)
    #         xy    = (x - prim_width/2, -chord)
    #         w     = prim_width
    #         h     = 2*chord
    #         rect  = plt.Rectangle(xy, w, h, fc='black')
    #         prim_support.append(rect)

    #     if draw_sec:
    #         for y in sec_y:
    #             r     = prim_length/2 + frame_width/2
    #             chord = np.sqrt(r**2 - y**2)
    #             xy    = (-chord, y - sec_width/2)
    #             w     = 2*chord
    #             h     = sec_width
    #             rect  = plt.Rectangle(xy, w, h, fc='black')
    #             sec_support.append(rect)    

    #     # Add frame
    #     center = (0, 0)
    #     ax.add_patch(plt.Circle(center, prim_length/2 + frame_width, fc='black'))
    #     ax.add_patch(plt.Circle(center, prim_length/2, fc='white'))

    #     for p in prim_support:
    #         ax.add_patch(p)

    #     if draw_sec:
    #         for p in sec_support:
    #             ax.add_patch(p)

    #     ax.axis('scaled')
    #     ax.axis('off')

    #     return ax