# Setting up a IPMSM as done in electrical drives.
# The goal is to include the dynamics of electrical torque and include the inverter control

import numpy as np
import matplotlib.pyplot as plt

# region Controller class
class PIController_LAH:
    def __init__(self, kp, ti):
        self.kp = kp
        self.ti = ti
        self.integral = 0.0

    def compute(self, error, dt):
        self.integral += error * dt

        # Anti-windup
        if self.integral > 1.5:
            self.integral = 1.5
        if self.integral < -1.5:
            self.integral = -1.5
        
        return self.kp * error + (self.kp/self.ti) * self.integral

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

# region PrimeMover class

class PrimeMover:
    def __init__(self, torque=1, speed=1, alpha=0.5):
        self.torque = torque
        self.speed = speed
        self.torque_ref = torque
        self.speed_ref = speed
        self.alpha = alpha
        self.T = 1e-1           # Time constant for the first-order filter

    def set_reference_values(self, torque_ref, speed_ref):
        self.torque_ref = torque_ref
        self.speed_ref = speed_ref

    def update_values(self, dt):
        # Apply first-order filter to torque and speed
        self.torque += (1/self.T)*(self.torque_ref - self.torque)*dt

    def get_values(self):
        self.update_values()
        return self.torque, self.speed

# endregion

# region Converter class

class Converter:
    def __init__(self, vd=0.0, vq=0.0, alpha=0.95):
        self.vd = vd
        self.vq = vq
        self.vd_ref = vd
        self.vq_ref = vq
        self.alpha = alpha
        self.T = 1e-4           # Time constant for the first-order filter

    def clamp_value(self, value, limit):
        if value > limit:
            value = limit
        if value < -limit:
            value = -limit
        return value

    def set_reference_voltages(self, vd_ref, vq_ref):
        self.vd_ref = vd_ref
        self.vq_ref = vq_ref
        # converter.update_voltages()

    def update_voltages(self, vd_ctrl, vq_ctrl, dt):
        # Apply first-order filter to vd and vq
        self.vd += (1/self.T) * (vd_ctrl - self.vd) * dt
        self.vq += (1/self.T) * (vq_ctrl - self.vq) * dt

        # Using instantaneous values and limiting the voltages
        # self.vd = vd_ctrl
        # self.vq = vq_ctrl

        # Limit the voltages
        self.vd = self.clamp_value(self.vd, 5)
        self.vq = self.clamp_value(self.vq, 5)


    def get_voltages(self):
        self.update_voltages()
        return self.vd, self.vq
    
# endregion

# region IMPSM class
# Tror jeg kan unngå å bruke DAEModel ettersom at MSC er dekoblet fra det dynamiske nettet
class IPMSM(Converter, PrimeMover, PIController_LAH):
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

        self.speed_measured = self.speed

        # Initialize PI controllers for i_d and i_q
        self.pi_controller_id = PIController_LAH(kp=1.27, ti=0.14)
        self.pi_controller_iq = PIController_LAH(kp=1.27, ti=0.14)
        self.pi_controller_speed = PIController_LAH(kp=56, ti=0.04)

        # Desired currents and speed
        self.i_d_ref = 0.0          # Må kanskje endres senere
        self.torque_ref = 0.0       # Will be asigned from speed controller
        self.i_q_ref = 0.0
        # self.speed_ref = self.speed

        # refernce initial values
        self.target_speed_ref = self.primemover.speed_ref
        self.target_torque_ref = self.primemover.torque_ref
        self.ramp_duration = 0
        self.ramp_start_time = 0

    # region derivatives

    def derivatives(self):
        """
        V3
        From block diagram electric drives IPMSM SIMULINK, currents as state variables

        parmas = {% Electrical
                    r_s   = 0.03; [pu] stator resistance
                    x_s   = 0.4;  [pu] stator inductance
                    x_d   = 0.4;  [pu] stator inductance
                    x_q   = 0.4;  [pu] stator inductance
                    psi_m = 0.66; [pu] magnetic flux 

                    f_n = 50;       [Hz] nominal frequency
                    w_n = 2*pi*f_n; [rad/s] nominal frequency

                    % Mechanical
                    T_m = 0.8; [s] mechanical time constant
        
        """
        dX = {}
        p = self.params
        
        psi_q = self.i_q*p["x_q"]
        psi_d = self.i_d*p["x_d"] + p["Psi_m"]

        dX["i_d"] = (self.converter.vd - p["rs"]*self.i_d + psi_q*self.speed) * (p["w_n"]/p["x_d"])
        dX["i_q"] = (self.converter.vq - p["rs"]*self.i_q - psi_d*self.speed) * (p["w_n"]/p["x_d"])

        Te = psi_d*self.i_q - psi_q*self.i_d

        dX["speed"] = 1/p["Tm"] * (self.primemover.torque - Te)

        return dX    
       
    # endregion
    
    # region Speed controll functions

    def update_speed_control(self, dt):
        """
        Method to update the reference torque based on the speed error
        Input: dt - time step
        
        Output: self.i_q_ref - updated reference q-axis current
        """

        # Compute speed error
        error_speed = -(self.primemover.speed_ref - self.primemover.speed)

        # Compute control signal for torque
        self.torque_ref = self.pi_controller_speed.compute(error_speed, dt)

        # Limit torque reference
        self.torque_ref = self.clamp_value(self.torque_ref, 1.5)

        # Calculate the required i_q_ref to produce the needed torque (algebraic)
        p = self.params
        self.i_q_ref = self.torque_ref / p["Psi_m"]

        self.i_q_ref = self.clamp_value(self.i_q_ref, 1.2)

    # endregion


    # region Current controll functions
    def update_current_control(self, dt):
        p = self.params

        ### Decoupling and current control ###

        # Compute errors
        error_id = self.i_d_ref - self.i_d
        error_iq = self.i_q_ref - self.i_q

        # Compute decoupled voltage control signals
        v_dII = self.i_q * p["x_q"] * self.speed
        v_qII = self.i_d * p["x_d"] + p["Psi_m"] * self.speed

        # I_d current control
        Kd = error_id*self.pi_controller_id.kp
        Ti_d = Kd/self.pi_controller_id.ti
        Ti_d = self.clamp_value(Ti_d, 5)
        v_d_ctrl = v_dII + Kd + Ti_d
        v_d_ctrl = self.clamp_value(v_d_ctrl, 3)

        # I_q current control
        Kq = error_iq*self.pi_controller_iq.kp
        Ti_q = Kq/self.pi_controller_iq.ti
        Ti_q = self.clamp_value(Ti_q, 5)
        v_q_ctrl = v_qII + Kq + Ti_q
        v_q_ctrl = self.clamp_value(v_q_ctrl, 3)

        # Input voltage control signal to converter and update the voltages
        self.set_converter_voltages(v_d_ctrl, v_q_ctrl, dt)
    
    # endregion

    def update_states(self, t, dt):

        # Update the reference values of the prime mover, this is because of any changes in the reference values (ramping)
        self.set_primemover_ramp_reference_values(current_time=t)

        # Applying torque ref change in the prime mover with the first order filter
        self.primemover.update_values(dt)

        # Update reference torque using PI controller
        self.update_speed_control(dt)

        # Update reference voltages using PI controllers
        self.update_current_control(dt)

        # Calculate the derivatives
        dX = self.derivatives()

        # Update the states using Euler integration
        self.speed += dX["speed"] * dt
        self.primemover.speed += dX["speed"] * dt
        self.i_d += dX["i_d"] * dt
        self.i_q += dX["i_q"] * dt

    def set_converter_voltages(self, v_d_ctrl, v_q_ctrl, dt):
        self.converter.set_reference_voltages(v_d_ctrl, v_q_ctrl)
        self.converter.update_voltages(v_d_ctrl, v_q_ctrl, dt)


    def set_prime_mover_reference(self, speed_ref, torque_ref, ramp_time, current_time, dt):
        self.target_speed_ref = speed_ref
        self.target_torque_ref = torque_ref
        self.ramp_duration = ramp_time / dt
        self.ramp_start_time = current_time


    def set_primemover_ramp_reference_values(self, current_time):
        if self.ramp_duration > 0:
            elapsed_time = current_time - self.ramp_start_time
            if elapsed_time < self.ramp_duration:
                ramp_factor = elapsed_time / self.ramp_duration
                self.primemover.speed_ref = (1 - ramp_factor) * self.primemover.speed_ref + ramp_factor * self.target_speed_ref
                self.primemover.torque_ref = (1 - ramp_factor) * self.primemover.torque_ref + ramp_factor * self.target_torque_ref
            else:
                self.primemover.speed_ref = self.target_speed_ref
                self.primemover.torque_ref = self.target_torque_ref
                self.ramp_duration = 0  # Ramp completed


    # region Utility functions
    def clamp_value(self, value, limit):
        if value > limit:
            value = limit
        if value < -limit:
            value = -limit
        return value

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
# 
# 
 # """
        # V2
        # From block diagram electric drives IPMSM with currents as state variables, dynamic analysis state space model forelesning
        
        # """
        # dX = {}
        # p = self.params
        
        # Tq = p["x_q"]/(p["wn"]*p["rs"])
        # Td = p["x_d"]/(p["wn"]*p["rs"])
        
        # Te = (p["Psi_m"]*self.i_q - (p["x_q"]-p["x_d"])*self.i_d*self.i_q)

        # vd = self.converter.vd
        # vq = self.converter.vq
        
        # # dX["i_d"] = -1/Td * self.i_d  + (self.speed*p["x_q"]/p["x_d"])*self.i_q + (self.speed/p["x_d"])*vd
        # # dX["i_q"] = -1/Tq * self.i_q - (self.speed*p["x_d"]/p["x_q"])*self.i_d - (self.speed/p["x_q"])*vq - self.speed*p["Psi_m"]/p["x_q"] + (self.speed/p["x_q"])*vq

        # dX["speed"] = 1/p["Tm"] * (self.primemover.torque - Te)
        # dX["i_d"] = (-1/Td * self.i_d)  + (self.speed/p["x_d"])*vd
        # dX["i_q"] = (-1/Tq * self.i_q)  + (self.speed/p["x_q"])*vq
        # return dX    

        # """
        # V1
        # From block diagram electric drives IPMSM with currents as state variables
        
        # """
        # dX = {}
        # p = self.params

        # if all(key in p for key in ["x_q", "rs", "x_d", "Psi_m", "Tm"]):
        #     Tq = p["x_q"]/(self.speed*p["rs"])
        #     Td = p["x_d"]/(self.speed*p["rs"])

        #     # Te = p["Psi_m"]*self.i_q - (p["xq"]-p["xd"])*self.i_d*self.i_q
        #     Psi_d = self.i_d*p["x_d"] + p["Psi_m"]
        #     Psi_q = self.i_q*p["x_q"]
  
        #     Te = Psi_d*self.i_q - Psi_q*self.i_d

        #     vd = self.converter.vd
        #     vq = self.converter.vq
            
        #     dX["speed"] = 1/p["Tm"] * (self.primemover.torque - Te)
        #     dX["i_d"] = -1/Td * self.i_d + (self.speed/p["rs"])*vd
        #     dX["i_q"] = -1/Tq * self.i_q + (self.speed/p["rs"])*vq
        # else:
        #     raise KeyError("One or more required parameters are missing in the params dictionary")
        # return dX

   


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