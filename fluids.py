""" Model classes for liquids and gases. """


import numpy as np
import logging as log


class Silicone():
    
    def __init__(self, weight, viscosity):
        self.weight = weight
        self.viscosity = viscosity
    
    def volume(self, t):
        #  Volume formula form Wacker Silicone Fluids AK
        return self.weight * (1/0.97) * (1 + 0.00092*(t-25) + 0.00000045*(t-25)**2)

    @classmethod
    def solubility_coeficient(cls, gas, t):
        if gas == 'co2':
            temps = [-25.,  0.,  25.,  50.,  75., 100., 125., 150., 175., 200.]
            solub = [2.7, 2.3, 1.9, 1.55, 1.22, 0.95,  0.8,  0.7,  0.6, 0.55]
        else:
            raise NotImplementedError(f'Gas typ not implemented: {gas}')
        t_array = np.array(t, ndmin=1)
        if any(t_array < temps[0]) or any(t_array > temps[-1]):
            raise LookupError(f'Temperature out of range: {temps[0]} to {temps[-1]}')
        return np.interp(t, temps, solub)
    
    def __repr__(self):
        return f'<Silicone Oil, {self.weight}g, {self.viscosity}vsc.>'


class CO2():
    __typ__ = 'co2'
    name = 'co2'
    r_j_kg = 188.92
    r_j_g = r_j_kg#*0.001
    density = 1.98
    density_g_ml = density/1000
    
    def __repr__(self):
        return f'<Gas {self.name}>'
