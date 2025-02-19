# region Converter class

class MachineSideConverter:

    def __init__(self, params):

        self.params = params

        self.v_d = self.params["v_d_0"]
        self.v_q = self.params["v_q_0"]
        self.T = self.params["T_conv"]

        self.v_d_ref = self.v_d
        self.v_q_ref = self.v_q


        
    def clamp_value(self, value, limit):
        if value > limit:
            value = limit
        if value < -limit:
            value = -limit
        return value

    def set_reference_voltages(self, v_d_ref, v_q_ref):
        self.v_d_ref = v_d_ref
        self.v_q_ref = v_q_ref

    def update_voltages(self, v_d_ctrl, v_q_ctrl, dt):
        # Apply first-order filter to v_d and v_q
        self.v_d += (1/self.T) * (v_d_ctrl - self.v_d) * dt        # Tsw = 2/3*fsw 
        self.v_q += (1/self.T) * (v_q_ctrl - self.v_q) * dt

        # Limit the voltages
        self.v_d = self.clamp_value(self.v_d, 2)
        self.v_q = self.clamp_value(self.v_q, 2)

    def get_voltages(self):
        self.update_voltages()
        return self.v_d, self.v_q
    
# endregion