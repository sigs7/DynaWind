
# from tops.cosim_models.dclink import DClink
# from tops.cosim_models.pmsm import PMSM
import numpy as np

# region Converter class

class MachineSideConverter:

    def __init__(self, params):

        self.params = params

        # Initialize voltages and time constant
        self.v_d = self.params.get("v_d_0", 0.0)
        self.v_q = self.params.get("v_q_0", 0.0)
        self.T = self.params.get("T_conv", 1.0)

        self.v_d_ref = self.v_d
        self.v_q_ref = self.v_q



    def connect_msc(self, dclink, pmsm):
        self.dclink = dclink
        self.pmsm = pmsm


    def derivatives(self):
        dX = {}

        # Get rotor electrical angle
        theta_r = self.pmsm.theta_elec

        # Scale control voltages by DC-link voltage
        if self.dclink.vdc_ref == 0:
            raise ValueError("DC-link reference voltage cannot be zero.")
        
        v_d_scaled = (self.pmsm.v_d_ctrl * self.dclink.vdc) / self.dclink.vdc_ref
        v_q_scaled = (self.pmsm.v_q_ctrl * self.dclink.vdc) / self.dclink.vdc_ref

        # V = Vdc * k0 * m ???

        # Converter boys forslag
        # v_d_scaled = self.pmsm.v_d_ctrl * (self.dclink.vdc_ref - self.dclink.vdc) * K_DC
        # v_q_scaled = self.pmsm.v_q_ctrl * (self.dclink.vdc_ref - self.dclink.vdc) * K_DC

        # # Convert control voltages from dq to αβ (stator reference frame)
        # v_alpha = v_d_scaled * np.cos(theta_r) - v_q_scaled * np.sin(theta_r)
        # v_beta = v_d_scaled * np.sin(theta_r) + v_q_scaled * np.cos(theta_r)

        # # Transform back to dq-frame using current rotor position
        # v_d_transformed = v_alpha * np.cos(theta_r) + v_beta * np.sin(theta_r)
        # v_q_transformed = -v_alpha * np.sin(theta_r) + v_beta * np.cos(theta_r)

        # First-order dynamics for voltage updates
        dX["v_d"] = (v_d_scaled - self.v_d) * (1 / self.T)
        dX["v_q"] = (v_q_scaled - self.v_q) * (1 / self.T)

        return dX


    def step_msc(self, time, step_size):

        # self.update_voltages(self.pmsm.v_d_ctrl, self.pmsm.v_q_ctrl, step_size)

        dX = self.derivatives()
        # Update the voltages
        self.v_d += dX["v_d"] * step_size
        self.v_q += dX["v_q"] * step_size


    def i_dc(self):
        # return (-3/2 * (self.pmsm.i_d*self.v_d + self.pmsm.i_q*self.v_q)) / self.dclink.vdc
        return self.p_e_dq() / self.dclink.vdc
        # return 1
    

    def clamp_value(self, value, min_value=None, max_value=None):
        """
        Clamp the value within the specified minimum and maximum limits.

        Parameters:
        value (float): The value to be clamped.
        min_value (float, optional): The minimum limit. Defaults to None.
        max_value (float, optional): The maximum limit. Defaults to None.

        Returns:
        float: The clamped value.
        """
        if min_value is not None and max_value is not None:
            # Clamp value between min_value and max_value
            return max(min_value, min(value, max_value))
        elif min_value is not None:
            # Clamp value to be at least min_value
            return max(min_value, value)
        elif max_value is not None:
            # Clamp value to be at most max_value
            return min(value, max_value)
        else:
            # No clamping needed
            return value

    
    def safe_arctan(self, x, y, epsilon=1e-10):
        if abs(y) < epsilon:
            return 0.0
        return np.arctan(x / y)
    
    def p_e_dc(self):
        # return self.i_dc() * self.dclink.vdc
        return 1
    
    def p_e_dq(self):
        return -3/2 * (self.pmsm.i_d*self.v_d + self.pmsm.i_q*self.v_q)

    
