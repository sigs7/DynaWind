# DC-link class
# Contains DC-link voltage, capacitance, and chopper

# from tops.cosim_models.machine_side_converter import MachineSideConverter
# from tops.cosim_models.grid_side_converter import GridSideConverter
# from tops.cosim_models.pmsm import PMSM
import numpy as np

class DClink():
    def __init__(self, params : dict):

        # self.msc = msc
        # self.gsc = gsc
        # self.pmsm = pmsm

        self.vdc = params.get("vdc_0", 1.0)
        self.vdc_ref = params.get("vdc_ref", 1.0)
        self.Cdc = params.get("Cdc", 0.2)
        self.K_p_dc = params.get("K_p_dc", 1.0)
        self.T_i_dc = params.get("T_i_dc", 0.2)
        # self.duty = 0.0
        self.chopper_threshold = params.get("chopper_threshold", 1.2)
        self.chopper_resistance = params.get("chopper_resistance", 10)

        # Controller signals
        self.x_pref_adj = 0.0


    def connect_dclink(self, msc, pmsm):
        self.msc = msc
        self.pmsm = pmsm

    def derivatives(self, p_gsc):
        dX = {}
        dX["vdc"] = self.vdc*(self.msc.i_dc() - self.i_chopper()) - p_gsc / (self.Cdc * self.vdc)
        dX["x_pref_adj"] = (self.vdc - self.vdc_ref) * (self.K_p_dc / self.T_i_dc)

        return dX       

    def i_chopper(self):
        if self.vdc > self.chopper_threshold:
            return (self.vdc * self.duty()) / self.chopper_resistance
        else:
            return 0.0
        
    def gsc_p_ref(self):
        return self.vdc*self.msc.i_dc() + self.p_adjust()
        
        
    def duty(self):
        if self.vdc > 1*self.chopper_threshold:
            return 1.0
        else:
            return 0.0

    def step_dclink(self, time, step_size, p_gsc):
        # Update all the states from other models

        dX = self.derivatives(p_gsc)
        self.vdc += dX["vdc"] * step_size
        self.x_pref_adj += dX["x_pref_adj"] * step_size

    def p_adjust(self):
        return self.x_pref_adj + (self.vdc - self.vdc_ref) * self.K_p_dc





