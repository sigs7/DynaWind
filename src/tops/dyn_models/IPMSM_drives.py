# Setting up a IPMSM as done in electrical drives.
# The goal is to include the dynamics of electrical torque and include the inverter control

import numpy as np
from tops.dyn_models.blocks import *
from tops.dyn_models.utils import DAEModel, output
import numpy as np
from tops.dyn_models.utils import DAEModel
import tops.utility_functions as dps_uf
import time
import matplotlib.pyplot as plt

class PIController_LAH:
    def __init__(self, kp, ki):
        self.kp = kp
        self.ki = ki
        self.integral = 0.0

    def compute(self, error, dt):
        self.integral += error * dt
        return self.kp * error + self.ki * self.integral


class MSIC:
    """
    Machine Side Inverter Controller

    Inputs: 
        P_ref
        V_ref
    
    Outputs:
        Vd_ref
        Vq_ref

    """
class PrimeMover:
    def __init__(self, torque=1, speed=1, alpha=0.5):
        self.torque = torque
        self.speed = speed
        self.torque_ref = torque
        self.speed_ref = speed
        self.alpha = alpha

    def set_reference_values(self, torque_ref, speed_ref):
        self.torque_ref = torque_ref
        self.speed_ref = speed_ref

    def update_values(self):
        # Apply first-order filter to torque and speed
        self.torque += self.alpha * (self.torque_ref - self.torque)
        # self.speed += self.alpha * (self.speed_ref - self.speed)

    def get_values(self):
        self.update_values()
        return self.torque, self.speed

class Converter:
    def __init__(self, vd=0.0, vq=0.0, alpha=0.95):
        self.vd = vd
        self.vq = vq
        self.vd_ref = vd
        self.vq_ref = vq
        self.alpha = alpha

    def set_reference_voltages(self, vd_ref, vq_ref):
        self.vd_ref = vd_ref
        self.vq_ref = vq_ref
        # converter.update_voltages()

    def update_voltages(self):
        # Apply first-order filter to vd and vq
        self.vd += self.alpha * (self.vd_ref - self.vd)
        self.vq += self.alpha * (self.vq_ref - self.vq)

    def get_voltages(self):
        self.update_voltages()
        return self.vd, self.vq

# Tror jeg kan unngå å bruke DAEModel ettersom at MSC er dekoblet fra det dynamiske nettet
class IPMSM(Converter, PrimeMover, MSIC):
    """
    Internally Permanent Magnet Synchrnous Machine
    
    """
    
    def __init__(self, params : dict, converter : Converter, prime_mover : PrimeMover, i_d0=0.0, i_q0=1.0):
        """
        params = {
            "x_q": 0.533,         q-axis reactance [pu]
            "x_d": 1.07,         d-axis reactancce
            "speed": 1.0,       nominal speed
            "rs": 0.01,         stator resistance
            "Psi_m": 0.8,       Permanent magnet flux linkage
            "T_m" : 3.7         Mechanical time constant -> Tm =(J*omega**2)/Sn
            
            # Other parameters...
            }
        
        """
        # Assigning the basic parameters of the IPMSM
        self.params = params

        # Assigning the converter and prime mover class. They should be initated before initiating the IPMSM
        self.converter = converter
        self.primemover = prime_mover

        # Initiating some basic initial values of operation, should be updated later based of load flow solution
        self.i_d = i_d0
        self.i_q = i_q0
        self.speed = self.primemover.speed

        # Initialize PI controllers for i_d and i_q
        self.pi_controller_id = PIController_LAH(kp=1.0, ki=0.1)
        self.pi_controller_iq = PIController_LAH(kp=20.0, ki=0.2)
        self.pi_controller_speed = PIController_LAH(kp=5.0, ki=0.1)

        # Desired currents and speed
        self.i_d_ref = 0.0          # Må kanskje endres senere
        # self.i_q_ref = self.i_q
        # self.speed_ref = self.speed

    def derivatives(self, dx, x, v):
        """
        From block diagram electric drives IPMSM with currents as state variables, dynamic analysis state space model forelesning
        
        """
        dX = {}
        p = self.params
        
        Tq = p["x_q"]/(p["wn"]*p["rs"])
        Td = p["x_d"]/(p["wn"]*p["rs"])
        
        Te = p["Psi_m"]*self.i_q - (p["x_q"]-p["x_d"])*self.i_d*self.i_q

        vd = self.converter.vd
        vq = self.converter.vq
        
        # dX["i_d"] = -1/Td * self.i_d  + (self.speed*p["x_q"]/p["x_d"])*self.i_q + (self.speed/p["x_d"])*vd
        # dX["i_q"] = -1/Tq * self.i_q - (self.speed*p["x_d"]/p["x_q"])*self.i_d - (self.speed/p["x_q"])*vq - self.speed*p["Psi_m"]/p["x_q"] + (self.speed/p["x_q"])*vq

        dX["speed"] = 1/p["Tm"] * (self.primemover.torque - Te)
        dX["i_d"] = -1/Td * self.i_d  + (self.speed/p["x_d"])*vd
        dX["i_q"] = -1/Tq * self.i_q  + (self.speed/p["x_q"])*vq

        return dX

    # These should be set from the machine side controllers    
    def set_current_references(self, i_d_ref, i_q_ref):
        self.i_d_ref = i_d_ref
        self.i_q_ref = i_q_ref
    
    # This should come from the prime mover speed controller
    def set_speed_reference(self, speed_ref):
        self.speed_ref = speed_ref

    def update_reference_voltages(self, dt):
        # Compute errors
        error_id = self.i_d_ref - self.i_d
        error_iq = self.i_q_ref - self.i_q

        # Compute control signals
        v_d = self.pi_controller_id.compute(error_id, dt)
        v_q = self.pi_controller_iq.compute(error_iq, dt)

        # Set converter reference voltages
        self.set_converter_voltages(v_d, v_q)
    
    def update_reference_torque(self, dt):
        # Compute speed error
        error_speed = self.primemover.speed_ref - self.speed

        # Compute control signal for torque
        torque_ref = self.pi_controller_speed.compute(error_speed, dt)

        # # Set prime mover reference torque
        # self.set_prime_mover_reference(torque_ref, self.speed_ref)

        # Calculate the required i_q_ref to produce the needed torque
        p = self.params
        self.i_q_ref = torque_ref / p["Psi_m"]

    def update_states(self, dt):

        # Update reference torque using PI controller
        self.update_reference_torque(dt)

        # Update reference voltages using PI controllers
        self.update_reference_voltages(dt)



        # Create a dictionary to hold the derivatives
        dx = {"speed": 0.0, "i_d": 0.0, "i_q": 0.0}
        x = {"speed": self.speed, "i_d": self.i_d, "i_q": self.i_q}
        v = {}  # Add any additional inputs if needed

        # Calculate the derivatives
        dX = self.derivatives(dx, x, v)

        # Update the states using Euler integration
        self.speed += dX["speed"] * dt
        self.i_d += dX["i_d"] * dt
        self.i_q += dX["i_q"] * dt

    def set_converter_voltages(self, vd_ref, vq_ref):
        self.converter.set_reference_voltages(vd_ref, vq_ref)
        self.converter.update_voltages()

    def set_prime_mover_reference(self, torque_ref, speed_ref):
        self.primemover.set_reference_values(torque_ref, speed_ref)
        self.primemover.update_values()

    # region Utility functions

    def get_Te(self):
        p = self.params
        return p["Psi_m"]*self.i_q - (p["x_q"]-p["x_d"])*self.i_d*self.i_q

    def get_vd(self):
        return self.converter.vd
    
    def get_vq(self):
        return self.converter.vq

    def get_speed(self):
        return self.speed
    
    def get_Pe(self):
        p = self.params
        return self.Vq*self.i_q + self.Vd*self.i_d + (self.i_d**2 + self.i_q**2)*p["rs"]
    
    def get_Pm(self):
        return self.speed * self.primemover.torque

    # endregion




######################   


    # def state_list(self):
    #     """
    #     Using a model with currents as state variables

    #     vd = rs*id + xd/wn * di_d/dt - n*xq*iq
    #     vq = rs*iq + xd/wn * di_q/dt - n*xd*id + n*Psi_m

    #     Psi_d = xd * id + Psi_m
    #     Psi_q = xq*iq

    #     te = Psi_m*iq - (xq-xd)*id*iq

    #     Adjusting the equations to fit a current controller

    #     Tq = xq/(wn*rs)
    #     Td = xd/(wn*rs)

    #     di_d/dt = -1/Td * id + wn/xd * ud
    #     di_q/dt = -1/Tq * iq + wn/xq * uq

    #     """
    #     return ["i_d", "i_q", "speed", "angle"]
    
        # def input_list(self):
    #     return ['v_d', 'v_q']


    # def states(self):
    #     """
    #     Using a model with currents as state variables

    #     vd = rs*id + xd/wn * di_d/dt - n*xq*iq
    #     vq = rs*iq + xd/wn * di_q/dt - n*xd*id + n*Psi_m

    #     Psi_d = xd * id + Psi_m
    #     Psi_q = xq*iq

    #     te = Psi_m*iq - (xq-xd)*id*iq

    #     Adjusting the equations to fit a current controller

    #     Tq = xq/(wn*rs)
    #     Td = xd/(wn*rs)

    #     di_d/dt = -1/Td * id + wn/xd * ud
    #     di_q/dt = -1/Tq * iq + wn/xq * uq

    #     """

    #     X = {"i_d" : self.i_d,
    #          "i_q" : self.i_q,
    #          "speed" : self.speed
    #         }
    #     return 

        
    # def derivatives(self, dx, x, v):
    #     """
    #     From block diagram electric drives IPMSM with currents as state variables
        
    #     """
    #     dX = {}
    #     p = self.params

    #     if all(key in p for key in ["x_q", "rs", "x_d", "Psi_m", "Tm"]):
    #         Tq = p["x_q"]/(self.speed*p["rs"])
    #         Td = p["x_d"]/(self.speed*p["rs"])

    #         # Te = p["Psi_m"]*self.i_q - (p["xq"]-p["xd"])*self.i_d*self.i_q
    #         Psi_d = self.i_d*p["x_d"] + p["Psi_m"]
    #         Psi_q = self.i_q*p["x_q"]
  
    #         Te = Psi_d*self.i_q - Psi_q*self.i_d

    #         vd = self.converter.vd
    #         vq = self.converter.vq
            
    #         dX["speed"] = 1/p["Tm"] * (self.primemover.torque - Te)
    #         dX["i_d"] = -1/Td * self.i_d + (self.speed/p["rs"])*vd
    #         dX["i_q"] = -1/Tq * self.i_q + (self.speed/p["rs"])*vq
    #     else:
    #         raise KeyError("One or more required parameters are missing in the params dictionary")
    #     return dX