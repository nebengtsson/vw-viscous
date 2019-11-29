import numpy as np
import logging as log


KELVIN = 273.15


class Viscous():
    """ Model of dynamics in a Viscous coupling in a VW T3 syncro. """
    
    atm_pressure = 1.025
    
    # VC dynamic volume, pressure dependent
    # Coefficent estimated from extensive messurements with dail indicator on VC under pressure 0-9 bar
    lid_play = 1.21  # cm^3 att positiv pressure
    c_ring_play_exp = [0.49, -1]  # [A, B], V = A*e^B*p
    lid_bending_coef = 0.136  # A, V = A*p
    # TODO, determine bottom
    bottom_bending_coef = 0.04  # A, V = A*p
    
    def __init__(self, volume):
        self.fix_volume = volume
        
        self.temp = 25  # Celcius
        self.pressure = self.atm_pressure  # bar
        self.slip_speed = 0  # rpm
        self.filled = False
        
    def fill(self, oil, gas, temp, pressure, solved_gas):
        self.oil = oil
        self.gas = gas
        
        k = temp + KELVIN
        gas_solved = (oil.solubility_coeficient(gas.name, temp) * gas.density_g_ml * 
                           oil.weight * solved_gas * pressure)
        gas_gas = (pressure * (self.volume(pressure) - self.oil.volume(temp)) / 
                            (gas.r_j_g * k))
        
        self.gas_weight = gas_solved + gas_gas
        self.gas_solved = gas_solved
        self.filled = True
        log.debug(f'Init gas, solved {gas_solved:.4}g, gas {gas_gas:.4}g.')
        print(f'Viscous filled with {oil} and {gas} {self.gas_weight:.4}g')
        print(f'Filling: {100*self.oil.volume(25)/self.fix_volume:0.1f}% of fix volume @25c')
        print(f'Filling: {100*self.oil.volume(25)/self.volume(10):0.1f}% of fix volume @10bar')
        
    def volume(self, pressure):
        vol = self.fix_volume
        rel_p = pressure - self.atm_pressure
        if rel_p >= 0:
            vol += self.lid_play
            vol += self.c_ring_play_exp[0]*(1 - np.exp(self.c_ring_play_exp[1]*rel_p))
            vol += self.lid_bending_coef * rel_p
            vol += self.bottom_bending_coef * rel_p
        return vol
        
    def pressure_gas_desolved(self, temp):
        ALMOST_VACUME = 0.01
        oil_volume = self.oil.volume(temp)
        vc_volume_no_expassion = self.volume(self.atm_pressure)
        oil_extra_volume = oil_volume - vc_volume_no_expassion
        if oil_extra_volume > 0:
            #return self.atm_pressure + oil_extra_volume/self.expantion_coef
            return self._pressure(oil_volume)
        return ALMOST_VACUME
        
    def pressure_equalibrium(self, temp):
        gas_solvability = self.oil.solubility_coeficient(self.gas.name, temp) * self.gas.density_g_ml
        oil_volume = self.oil.volume(temp)
        # Iterate solution
        
        MAX_ITERATION = 20
        LIMIT_ERROR = 0.001
        
        pressure_gas_desolved = self.pressure_gas_desolved(temp)
        pressure = None
        pressure_last_iter = pressure_gas_desolved        
        for i in range(MAX_ITERATION):
            gas_volume = self.volume(pressure_last_iter) - oil_volume
            pressure_equalibrium = (
                self.gas_weight / 
                (gas_solvability * self.oil.weight + (gas_volume/(self.gas.r_j_g * (temp + KELVIN))))
            )
            pressure = max(pressure_equalibrium, pressure_gas_desolved)
            log.debug(f'Pressure {pressure}, pressure_last_iter {pressure_last_iter}, pressure_equalibrium: {pressure_equalibrium}')
            if abs(pressure - pressure_last_iter) < LIMIT_ERROR:
                break
            pressure_last_iter = pressure
            
        else:
            raise ValueError('Did not converge')
        return pressure
    
    def saturated_solved_gas(self):
        gas_solvability = self.oil.solubility_coeficient(self.gas.name, self.temp) * self.gas.density_g_ml
        return gas_solvability * self.oil.weight * self.pressure
    
    #def 
    
    def delta_gas_solve(self, seconds):
        pass
        # defussion * vicosity / temperature(k) = konstant
        
        
    def _pressure(self, volume):
        ALMOST_VACUME = 0.01
        if volume <= self.fix_volume:
            return ALMOST_VACUME
        elif volume <= self.volume(self.atm_pressure):
            return self.atm_pressure
        
        MAX_ITERATION = 20
        LIMIT_ERROR = 0.001
        delta_p = 0.01
        init_step_size = 1.0
        p = self.atm_pressure + init_step_size
        
        for i in range(MAX_ITERATION):
            v_ = self.volume(p)
            error = volume - v_
            if error < LIMIT_ERROR:
                return p
            gradient = (self.volume(p + delta_p) - v_)/delta_p
            p = p + error/gradient
            p = max(p, self.atm_pressure)  # dont use p where volume(p) is discontinious
        else:
            raise ValueError('Pressure by volume did not converge')
        
