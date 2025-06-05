#Ideal Generator model

import numpy as np

class IdealGenerator:
    """
    Ideal Generator simplifies the electric transients in the turbine. A simpler model enables faster computation.
    """

    def __init__(self, ideal_generator_parameters : dict):

        self.parameters = ideal_generator_parameters

        self.torque_ref = 0
        self.torque = 0
        
        self.speed = 0
        self.SPEED = 0

    def update_ideal_generator(self, fast):
        self.torque_ref = -1 * fast.fmu.getReal([fast.vrs['GenTq']])[0] / self.parameters["T_r"]
        self.speed = fast.fmu.getReal([fast.vrs['GenSpeed']])[0] / self.parameters["rpm_n"]
        self.SPEED = self.speed * self.parameters["rpm_n"]


    def derivatives(self):
        dX = {}
        p = self.parameters

        dX["torque"] = (self.torque_ref - self.torque) / p["Tg"]

        return dX
    
    def step_ideal_generator(self, fast, step_size):

        self.update_torque_ref(fast)

        dX = self.derivatives()

        self.torque += dX["torque"] * step_size


    def p_e(self):
        return self.torque * self.speed
    
