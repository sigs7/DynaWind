# DC-link class
# Contains DC-link voltage, capacitance, and chopper

# from tops.cosim_models.machine_side_converter import MachineSideConverter
# from tops.cosim_models.grid_side_converter import GridSideConverter
# from tops.cosim_models.pmsm import PMSM
import numpy as np

class DClink():
    def __init__(self, params : dict):

        """          
        dclink_params = {
                "vdc_0": 1.0,
                "vdc_ref": 1.0,
                "chopper_threshold": 1.1,
                "Cdc": 0.3,
                "K_p_dc": 4.0,
                "T_i_dc": 0.1,
                "chopper_resistance": 1,
            }
        """

        # Initiate states
        self.vdc = params.get("vdc_0", 1.0)
        self.x_pref_adj = 0.0
        self.vdc_ref = params.get("vdc_ref", 1.0)
        self.chopper_on = False
        self.chopper_timer = 0

        # DC-link parameters
        self.params = params


    def connect_dclink(self, msc, pmsm):
        self.msc = msc
        self.pmsm = pmsm

    def derivatives(self, p_gsc):
        dX = {}
        p = self.params
        dX["vdc"] = self.msc.p_e_dq() - self.vdc*self.i_chopper() - p_gsc / (p["Cdc"] * self.vdc)
        dX["x_pref_adj"] = (self.vdc - self.vdc_ref) * (p["K_p_dc"] / p["T_i_dc"])

        return dX       

    def i_chopper(self):

        if self.vdc >= self.params["chopper_threshold"]:
            self.chopper_on = True
            self.chopper_timer = 50  # Number of timesteps to remain on         50 steps * 5ms = 0.25s

        if self.chopper_on:
            self.chopper_timer -= 1
            if self.chopper_timer <= 0:
                self.chopper_on = False

        if self.vdc < 1:
            self.chopper_on = False

        if self.chopper_on:
            return (self.vdc * self.duty()) / self.params["chopper_resistance"]
        else:
            return 0.0
        
    # def i_chopper(self):
    #     if self.vdc >= self.params["chopper_threshold"]:
    #         return (self.vdc * self.duty()) / self.params["chopper_resistance"]
    #     else:
    #         return 0.0
        
    def gsc_p_ref(self):
        return max(0.0, self.msc.p_e_dq() - (self.vdc * self.i_chopper()) + self.p_adjust())
        
        
    def duty(self):
        if self.vdc > 1.2*self.params["chopper_threshold"]:
            return 1
        if self.vdc > 1.05*self.params["chopper_threshold"]:
            return 0.75
        if self.vdc > 1.0*self.params["chopper_threshold"]:
            return 0.5
        
        else:
            return 0.0

    def step_dclink(self, time, step_size, p_gsc):

        # Update all the states from other models
        dX = self.derivatives(p_gsc)
        self.vdc += dX["vdc"] * step_size

        # Anti-windup for the integral term
        self.x_pref_adj += dX["x_pref_adj"] * step_size
        self.x_pref_adj = np.clip(self.x_pref_adj, -1, 1)

    def p_adjust(self):
            return np.clip(self.x_pref_adj + (self.vdc - self.vdc_ref) * self.params["K_p_dc"], -1, 1)





