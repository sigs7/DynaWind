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

        dX["vdc"] = (self.msc.p_e_dq() - self.vdc*self.i_chopper() - p_gsc) / (p["Cdc"] * self.vdc)
        dX["x_pref_adj"] = (self.vdc - self.vdc_ref) * (p["K_p_dc"] / p["T_i_dc"])

        return dX

    def i_chopper(self):
        if self.chopper_on:
            return (self.vdc * self.duty()) / self.params["chopper_resistance"]
        else:
            return 0.0
        
        
    def gsc_p_ref(self):
        return max(0.0, self.msc.p_e_dq() - (self.vdc * self.i_chopper()) + self.p_adjust())
        
        
    def duty(self):
        if self.vdc > 1.2*self.params["chopper_threshold"]:
            return 1
        elif self.vdc > 1.05*self.params["chopper_threshold"]:
            return 0.75
        elif self.vdc > 1.02*self.params["chopper_threshold"]:
            return 0.5
        elif self.vdc > 1.0*self.params["chopper_threshold"]:
            return 0.25
        else:
            return 0.0

        
    # def update_chopper_status(self, time, step_size):
    #     chopper_activation_duration = 5  # seconds
    #     chopper_threshold_on = self.params["chopper_threshold"]
    #     chopper_threshold_off = 1.001  # Off threshold set to 1.001

    #     # Debug prints
    #     if 3.99 < time < 4.5:
    #         print(f"vdc: {self.vdc}, chopper_on: {self.chopper_on}, chopper_timer: {self.chopper_timer}")

    #     # Turn on the chopper if the voltage exceeds the threshold
    #     if self.vdc >= chopper_threshold_on and self.chopper_on == False:
    #         self.chopper_on = True
    #         self.chopper_timer = chopper_activation_duration / step_size  # Number of timesteps to remain on

    #     # Turn off the chopper if the timer expires or the voltage drops below the off threshold
    #     if self.chopper_on:
    #         self.chopper_timer -= 1
    #         if self.chopper_timer <= 0 or self.vdc <= chopper_threshold_off:
    #             self.chopper_on = False

    #     # Optional safety check to turn off the chopper if the voltage is too low
    #     if self.vdc < 1:
    #         self.chopper_on = False

    #     # # Debug prints after update
    #     # if 3.99 < time < 4.5:
    #     #     print(f"Updated chopper_on: {self.chopper_on}, chopper_timer: {self.chopper_timer}")

    def update_chopper_status(self, time, step_size):

        min_on_time = 0.5  # Ensure chopper stays on for at least 50 ms

        if self.vdc >= self.params["chopper_threshold"] and not self.chopper_on:
            self.chopper_on = True
            self.chopper_timer = max(self.chopper_timer, min_on_time / step_size)  

        if self.chopper_on:
            self.chopper_timer -= step_size
            if self.chopper_timer <= 0 or self.vdc <= 1.02:
                self.chopper_on = False

        # # Debug prints
        # if 3.99 < time < 4.5:
        #     print(f"vdc: {self.vdc}, chopper_on: {self.chopper_on}, chopper_timer: {self.chopper_timer}")


    def step_dclink(self, time : float, step_size : float, p_gsc : float):

        # Update chopper status
        self.update_chopper_status(time, step_size)

        # Update all the states from other models
        dX = self.derivatives(p_gsc)
        
        self.vdc += dX["vdc"] * step_size

        # Anti-windup for the integral term
        self.x_pref_adj += dX["x_pref_adj"] * step_size
        self.x_pref_adj = np.clip(self.x_pref_adj, -1, 1)

    def p_adjust(self):
            # return np.clip(self.x_pref_adj + (self.vdc - self.vdc_ref) * self.params["K_p_dc"], -1, 1)
            return self.x_pref_adj + (self.vdc - self.vdc_ref) * self.params["K_p_dc"]





