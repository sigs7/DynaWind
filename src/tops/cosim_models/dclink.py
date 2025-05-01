# DC-link class

# Contains DC-link voltage, controll, capacitance, and chopper
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
        self.vdc_ref = params.get("vdc_ref", 1.0)

        self.x_pref_adj = 0.0
        
        self.chopper_on = False
        self.chopper_timer = 0

        # DC-link parameters
        self.params = params

    def connect_dclink(self, pmsm):
        self.pmsm = pmsm

    def derivatives(self, p_gsc):
        dX = {}
        p = self.params

        dX["vdc"] = (self.pmsm.p_e() - self.vdc*self.i_chopper() - p_gsc) / (p["Cdc"] * self.vdc)
        dX["x_pref_adj"] = (self.vdc - self.vdc_ref) * (p["K_p_dc"] / p["T_i_dc"])

        return dX

    def i_chopper(self):
        if self.chopper_on:
            return (self.vdc * self.duty()) / self.params["chopper_resistance"]
        else:
            return 0.0

        
    def gsc_p_ref(self):
        # if self.pmsm.state == "fault" or self.pmsm.state == "hold":     # During fault, the power balance is handled by the PMSM
        #     return max(0.0, self.pmsm.p_e())
        # else:
        #     return max(0.0, self.pmsm.p_e() + self.p_adjust()) # - (self.vdc * self.i_chopper())
        return max(0.0, self.pmsm.p_e() + self.p_adjust()) # - (self.vdc * self.i_chopper())

        
    def duty(self):
        vdc = self.vdc
        vth = self.params["chopper_threshold"]  # Should be 1.05
        max_duty = 1                            # To limit current spikes
        k_slope = 15.0                          # Duty increases 0 → max over 0.1 pu (e.g., 1.05 to 1.15)

        if vdc > vth:
            return np.clip(k_slope * (vdc - vth), 0.0, max_duty)
        else:
            return 0.0

    
    def update_chopper_status(self, time, step_size):

        min_on_time = 5e-5  # Ensure chopper stays on for at least 50 mikroseconds

        if self.vdc >= self.params["chopper_threshold"] and not self.chopper_on:
            self.chopper_on = True
            self.chopper_timer = max(self.chopper_timer, min_on_time / step_size)

        if self.chopper_on:
            self.chopper_timer -= step_size
            if self.chopper_timer <= 0 or self.vdc <= self.params["chopper_threshold"]:
                self.chopper_on = False

    #     # # Debug prints
    #     # if 3.99 < time < 4.5:
    #     #     print(f"vdc: {self.vdc}, chopper_on: {self.chopper_on}, chopper_timer: {self.chopper_timer}")

    def step_dclink(self, time : float, step_size : float, p_gsc : float):

        # Update chopper status
        self.update_chopper_status(time, step_size)

        # Update all the states from other models
        dX = self.derivatives(p_gsc)
        
        self.vdc += dX["vdc"] * step_size

        # Anti-windup for the integral term
        self.x_pref_adj += dX["x_pref_adj"] * step_size
        # self.x_pref_adj = np.clip(self.x_pref_adj, -0.2, 0.2)


    def p_adjust(self):
            # return np.clip(self.x_pref_adj + (self.vdc - self.vdc_ref) * self.params["K_p_dc"], -0.1, 0.1)
            if self.chopper_on:
                self.x_pref_adj = 0.0
            return self.x_pref_adj + (self.vdc - self.vdc_ref) * self.params["K_p_dc"]

