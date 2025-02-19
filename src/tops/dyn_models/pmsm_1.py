# Setting up a PMSM as done in electrical drives.
# The goal is to include the dynamics of electrical torque and include the inverter control

import numpy as np
import matplotlib.pyplot as plt
from tops.dyn_models.utils import DAEModel

def rkf45(f, y0, t0, dt):
    c = [0, 1/4, 3/8, 12/13, 1, 1/2]
    a = [
        [],
        [1/4],
        [3/32, 9/32],
        [1932/2197, -7200/2197, 7296/2197],
        [439/216, -8, 3680/513, -845/4104],
        [-8/27, 2, -3544/2565, 1859/4104, -11/40]
    ]
    b = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55]
    b_star = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0]

    k = [f(t0, y0)]
    for i in range(1, 6):
        yi = y0 + dt * sum(a[i][j] * k[j] for j in range(i))
        k.append(f(t0 + c[i] * dt, yi))

    y_new = y0 + dt * sum(b[i] * k[i] for i in range(6))
    y_star = y0 + dt * sum(b_star[i] * k[i] for i in range(6))
    error = np.linalg.norm(y_new - y_star)
    return y_new, error

# Adaptive timestep control
def adaptive_rkf45(f, y0, t0, dt, tol=1e-6):
    safety_factor = 0.9
    max_factor = 5.0
    min_factor = 0.2

    while True:
        y_new, error = rkf45(f, y0, t0, dt)
        if error < tol:
            break
        dt = max(min_factor, min(max_factor, safety_factor * (tol / error)**0.25)) * dt

    return y_new, dt

# region PI Controller class
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


# region PrimeMover class

class PrimeMover:
    def __init__(self, **params):

        self.params = params

        if "T_pm" in params:
            self.Tpm = params["T_pm"]
        else:
            self.Tpm = 1e-1

        self.torque = params["torque_0"] if "torque_0" in params else 0.6
        self.speed = params["speed_0"] if "speed_0" in params else 0.6
        self.torque_ref = self.torque
        self.speed_ref = self.speed
        

    def set_reference_values(self, torque_ref, speed_ref):
        self.torque_ref = torque_ref
        self.speed_ref = speed_ref

    def update_values(self, dt):
        # Apply first-order filter to torque and speed
        self.torque += (1 / self.Tpm) * (self.torque_ref - self.torque) * dt


    def get_values(self):
        self.update_values()
        return self.torque, self.speed

# endregion

# region Converter class

class MachineSideConverter:

    def __init__(self, **params):
        self.params = params
        if "T_conv" in params:
            self.T = params["T_conv"]
        else:
            self.T = 1e-4

        self.vd = params["vd_0"] if "vd_0" in params else 0.0       # Same code as above only in one line :D
        self.vq = params["vq_0"] if "vq_0" in params else 0.0

        self.vd_ref = self.vd
        self.vq_ref = self.vq
        

    def clamp_value(self, value, limit):
        if value > limit:
            value = limit
        if value < -limit:
            value = -limit
        return value

    def set_reference_voltages(self, vd_ref, vq_ref):
        self.vd_ref = vd_ref
        self.vq_ref = vq_ref

    def update_voltages(self, vd_ctrl, vq_ctrl, dt):
        # Apply first-order filter to vd and vq
        self.vd += (1/self.T) * (vd_ctrl - self.vd) * dt        # Tsw = 2/3*fsw 
        self.vq += (1/self.T) * (vq_ctrl - self.vq) * dt

        # Limit the voltages
        self.vd = self.clamp_value(self.vd, 2)
        self.vq = self.clamp_value(self.vq, 2)


    def get_voltages(self):
        self.update_voltages()
        return self.vd, self.vq
    
# endregion

# region PMSM class

# Tror jeg kan unngå å bruke DAEModel ettersom at MSC er dekoblet fra det dynamiske nettet
class PMSM(MachineSideConverter, PrimeMover, PIController_LAH):
    """
    Internally Permanent Magnet Synchrnous Machine

    """
    
    def __init__(self, pmsm_params : dict, MSC_params : dict, prime_mover_params : dict):
        
        """
        pmsm_params = {
            "s_n" : 10,     # MVA
            "U_n" : 730,    # V
            "rs": 0.03,     # Stator resistance
            "x_d": 0.4,     # Stator d-axis inductance
            "x_q": 0.4,     # Stator q-axis inductance
            "Psi_m": 0.9,   # Magnetic flux
            "w_n" : 2*np.pi*50  # nominal rad
    }
        
        """
        # Assigning the basic parameters of the PMSM
        self.params = pmsm_params

        # Assigning the converter and prime mover class. They should be initated before initiating the PMSM
        self.converter = MachineSideConverter(**MSC_params)
        self.primemover = PrimeMover(**prime_mover_params)

        # Initiating some basic initial values of operation, should be updated later based of load flow solution
        self.i_d = 0.0
        self.i_q = 0.0
        self.speed = self.primemover.speed

        # Initialize PI controllers for i_d, i_q, and speed with parameters adjusted for a larger timestep
        self.pi_controller_id = PIController_LAH(kp=0.5, ti=0.2)
        self.pi_controller_iq = PIController_LAH(kp=0.5, ti=0.2)
        self.pi_controller_speed = PIController_LAH(kp=20, ti=0.1)

        # Initiate reference values
        self.i_d_ref = 0.0          # Må kanskje endres senere
        self.torque_ref = 0.0       # Will be asigned from speed controller
        self.i_q_ref = 0.0
        # self.speed_ref = self.speed

        # reference initial values
        self.target_speed_ref = self.primemover.speed_ref
        self.target_torque_ref = self.primemover.torque_ref
        self.ramp_duration = 0
        self.ramp_start_time = 0

    # region Derivatives

    def derivatives(self):
        """
        V3
        From block diagram electric drives PMSM SIMULINK, currents as state variables

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
        Te = psi_d*self.i_q - psi_q*self.i_d

        # Motor convention
        dX["i_d"] = (self.converter.vd - p["rs"]*self.i_d + psi_q*self.speed) * (p["w_n"]/p["x_d"])
        dX["i_q"] = (self.converter.vq - p["rs"]*self.i_q - psi_d*self.speed) * (p["w_n"]/p["x_q"])

        # change in speed T_t - T_e
        dX["speed"] = 1/p["Tm"] * (self.primemover.torque + Te)

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
        error_speed = (self.primemover.speed_ref - self.primemover.speed)

        # Compute control signal for torque
        self.torque_ref = self.pi_controller_speed.compute(error_speed, dt)

        # Limit torque reference
        self.torque_ref = self.clamp_value(self.torque_ref, max_value=0, min_value=-2)        #, min_value=0.0

        # Calculate the required i_q_ref to produce the needed torque (algebraic)
        p = self.params
        self.i_q_ref = self.torque_ref / p["Psi_m"]

        # self.i_q_ref = self.clamp_value(self.i_q_ref, 1.2)

    # endregion


    # region Current controll functions
    def update_current_control(self, dt):
        p = self.params

        # Compute errors
        error_id = self.i_d_ref - self.i_d
        error_iq = self.i_q_ref - self.i_q

        ### Decoupling and current control ###
        # Compute decoupled voltage control signals
        v_dII = -self.i_q * p["x_q"] * self.speed
        v_qII = self.i_d * p["x_d"] + p["Psi_m"] * self.speed

        # I_d current control
        v_d_ctrl = v_dII + self.pi_controller_id.compute(error_id, dt)
        v_d_ctrl = self.clamp_value(v_d_ctrl, max_value=2, min_value=-2)
        
        # I_q current control
        v_q_ctrl = v_qII + self.pi_controller_iq.compute(error_iq, dt)
        v_q_ctrl = self.clamp_value(v_q_ctrl, max_value=2, min_value=-2)

        # Input voltage control signal to converter and update the voltages
        self.set_converter_voltages(v_d_ctrl, v_q_ctrl, dt)
    
    # endregion

    def update_states(self, t, dt, tol=1e-6):

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

        # # Update the states using Euler integration
        self.primemover.speed += dX["speed"] * dt
        self.speed += dX["speed"] * dt
        self.i_d += dX["i_d"] * dt
        self.i_q += dX["i_q"] * dt

        # region RKF45
        # Calculate the derivatives using adaptive RKF45 method
        # Define the Runge-Kutta-Fehlberg method (RKF45)

        # region Traapezoidal rule
        # # Update the states using the trapezoidal rule
        # dX1 = self.derivatives()
        
        # # Predict the next state
        # temp_speed = self.primemover.speed + dX1["speed"] * dt
        # temp_i_d = self.i_d + dX1["i_d"] * dt
        # temp_i_q = self.i_q + dX1["i_q"] * dt
        
        # # Calculate the derivatives at the predicted next state
        # self.primemover.speed = temp_speed
        # self.i_d = temp_i_d
        # self.i_q = temp_i_q
        # dX2 = self.derivatives()
        
        # # Correct the state using the average of the derivatives
        # self.primemover.speed += 0.5 * (dX1["speed"] + dX2["speed"]) * dt
        # self.speed = self.primemover.speed
        # self.i_d += 0.5 * (dX1["i_d"] + dX2["i_d"]) * dt
        # self.i_q += 0.5 * (dX1["i_q"] + dX2["i_q"]) * dt
        # endregion

        # region Runge-Kutta
        # # Update the states using Runge-Kutta 4th order method
        # Runge-Kutta 4th order method
        # k1 = self.derivatives()
        # self.i_d += k1["i_d"] * dt / 2
        # self.i_q += k1["i_q"] * dt / 2
        # self.speed += k1["speed"] * dt / 2

        # k2 = self.derivatives()
        # self.i_d += (k2["i_d"] - k1["i_d"]) * dt / 2
        # self.i_q += (k2["i_q"] - k1["i_q"]) * dt / 2
        # self.speed += (k2["speed"] - k1["speed"]) * dt / 2

        # k3 = self.derivatives()
        # self.i_d += (k3["i_d"] - k2["i_d"]) * dt
        # self.i_q += (k3["i_q"] - k2["i_q"]) * dt
        # self.speed += (k3["speed"] - k2["speed"]) * dt

        # k4 = self.derivatives()
        # self.i_d += (k4["i_d"] - k3["i_d"]) * dt / 6
        # self.i_q += (k4["i_q"] - k3["i_q"]) * dt / 6
        # self.speed += (k4["speed"] - k3["speed"]) * dt / 6
        # endregion

        # region Modified Euler
        # Update the states using Modified Euler method
        # k1 = self.derivatives()
        # temp_i_d = self.i_d + k1["i_d"] * dt
        # temp_i_q = self.i_q + k1["i_q"] * dt
        # temp_speed = self.speed + k1["speed"] * dt

        # # Temporarily update the states
        # self.i_d = temp_i_d
        # self.i_q = temp_i_q
        # self.speed = temp_speed

        # k2 = self.derivatives()

        # # Final update of the states
        # self.i_d += 0.5 * (k1["i_d"] + k2["i_d"]) * dt
        # self.i_q += 0.5 * (k1["i_q"] + k2["i_q"]) * dt
        # self.speed += 0.5 * (k1["speed"] + k2["speed"]) * dt
        # endregion

    # Machine side converter voltage updates
    # -> The references are first updated by setting the reference voltages from the current controll
    # -> The voltages are then updated by the converter via a first order filter
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
    # def clamp_value(self, value, limit, min_value=None):
    #     if min_value is not None:
    #         if value > limit:
    #             value = limit
    #         if value < min_value:
    #             value = min_value
    #     else:
    #         if value > limit:
    #             value = limit
    #         if value < -limit:
    #             value = -limit
        
    #     return value

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


    def get_Te(self):
        p = self.params
        psi_q = self.i_q*p["x_q"]
        psi_d = self.i_d*p["x_d"] + p["Psi_m"]

        return psi_d*self.i_q - psi_q*self.i_d

    def get_vd(self):
        return self.converter.vd
    
    def get_vq(self):
        return self.converter.vq

    def get_speed(self):
        return self.speed
    
    def get_Pe(self):
        return (-3/2)*(self.converter.vd*self.i_d + self.converter.vq*self.i_q)      # -3/2 fac to adjust to three-phase power and change direction as injection
        # p = self.params
        # return self.Vq*self.i_q + self.Vd*self.i_d + (self.i_d**2 + self.i_q**2)*p["rs"]

    def get_Qe(self):
        return (-3/2)*(self.converter.vq*self.i_d - self.converter.vd*self.i_q)      # 3/2 fac?
        # p = self.params
        # return self.Vq*self.i_q + self.Vd*self.i_d + (self.i_d**2 + self.i_q**2)*p["rs"]
    
    def get_Pm(self):
        return self.speed * self.primemover.torque

    # endregion
    def get_pmsm_results(self):
        return {
            'speed': self.speed,
            'i_d': self.i_d,
            'i_q': self.i_q,
            'vd': self.converter.vd,
            'vq': self.converter.vq,
            'Pe': self.get_Pe()
        }

    def store_pmsm_states(self, res, t):
        pmsm_results = self.get_pmsm_results()
        for key, value in pmsm_results.items():
            res[key].append(value)
        res['t'].append(t)




#### TESTING ####

# # region GridSideConverter class
# class GridSideConverter(DAEModel):
#     """
#     VSC PQ model - ref to VSC1 in tops/dyn_models/vsc1.py - cite Sjur Føyen

#     Model of the DC link and Grid Side Converter.

#     Its purpose is to control P and Q to the grid and to control the DC link voltage. 

#     'vsc': {
#         'VSC_PQ': [
#             ['name',   'bus', 'S_n',   "pref",   "qref",   'Cdc',   'k_p',    'k_q',  'T_p',  'T_q',     'k_pll',   'T_pll',    'T_i',    'i_max'],
#             ['VSC1',   'B1',    50,      1,         0,      0.1,       1,       1,     0.1,     0.1,        5,        1,         0.01,      1.2],
#         ],
#     }

#     Parameters:
#     w_n: Nominal frequency
#     Cdc: DC link capacitance

#     References:
#     Vdc_ref
#     Pref
#     Qref
#     id_ref
#     iq_ref

#     States:
#     i_d
#     i_q
#     vdc
#     angle
#     xpll    integral i pll

#     """

#     def __init__(self, params : dict, pmsm : PMSM, vsc_params : dict, *args, **kwargs):
#         super().__init__(*args, **kwargs)

#         self.bus_idx = np.array(np.zeros(self.n_units), dtype=[(key, int) for key in self.bus_ref_spec().keys()])
#         self.bus_idx_red = np.array(np.zeros(self.n_units), dtype=[(key, int) for key in self.bus_ref_spec().keys()])

#         # self.pref = pmsm.get_Pe()              # Active power reference will come from the PMSM
#         # self.qref = 0                           # Reactive power reference can/should come from grid side

#         self.vdc = 1.0
#         self.vdc_ref = self.vdc

#         self.vdc_controller = PIController_LAH(kp=0.3, ti=0.003)

#         self.par = params

#     # region Definitions

#     def init_from_load_flow(self, x_0, v_0, S):
#         X = self.local_view(x_0)

#         self._input_values['p_ref'] = self.par['p_ref']       # NESCESERRY
#         self._input_values['q_ref'] = self.par['q_ref']

#         vg = v_0[self.bus_idx_red['terminal']]

#         X['angle'] = np.angle(vg)
#         X['i_d'] = self.par['p_ref']/abs(vg)
#         X['i_q'] = self.par['q_ref']/abs(vg)
#         X['x_pll'] = 0
#         X["vdc"] = 1.0
        
#     def state_derivatives(self, dx, x, v):
#         dX = self.local_view(dx)
#         X = self.local_view(x)
#         par = self.par

#         self.pref = PMSM.get_Pe() - self.vdc_controller.compute((self.vdc_ref - X["vdc"]))*X["vdc"]
#         self.qref = 0.0     # Should be updated from grid side command

#         i_d_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (self.pref * self.v_d(x,v) - self.qref * self.v_q(x,v))
#         i_q_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (self.pref * self.v_q(x,v) + self.qref * self.v_d(x,v))

#         # Limiters
#         i_ref = (i_d_ref+1j*i_q_ref)
#         i_ref = i_ref*par['i_max']/np.maximum(par['i_max'],abs(i_ref))

#         dX['i_d'][:] = 1 / (par['T_i']) * (i_ref.real - X['i_d'])
#         dX['i_q'][:] = 1 / (par['T_i']) * (i_ref.imag - X['i_q'])
#         dX['x_pll'][:] = par['k_pll'] / (par['T_pll']) * (self.v_q(x,v))
#         dX['angle'][:] = X['x_pll']+par['k_pll']*self.v_q(x,v)
#         dX["vdc"][:] = (self.pref - self.p_e(x,v)) / (par["Cdc"] * X["vdc"])
#         #dX['angle'][:] = 0
#         return    

#     def load_flow_pq(self):
#         return self.bus_idx['terminal'], -self.par['p_ref']*self.par['S_n'], -self.par['q_ref']*self.par['S_n']

#     def int_par_list(self):
#         return ['f']

#     def bus_ref_spec(self):
#         return {'terminal': self.par['bus']}

#     # endregion

#     def state_list(self):
#         """
#         All states in pu.

#         i_d: d-axis current, first-order approximation of EM dynamics
#         i_q: q-axis current -------------""--------------
#         x_p: p-control integral
#         x_q: q-control integral
#         x_pll: pll q-axis integral
#         angle: pll angle
#         """
#         return ['i_d', 'i_q', 'x_p', 'x_q', 'x_pll', 'angle']

#     def input_list(self):
#         """
#         All values in pu.
#         p_ref: outer loop active power setpoint
#         q_ref: outer loop reactive power setpoint
#         """
#         return ['p_ref', 'q_ref']

#     def current_injections(self, x, v):
#         i_n_r = self.par['S_n'] / self.sys_par['s_n']
#         return self.bus_idx_red['terminal'], self.i_inj(x, v) * i_n_r

#     # region Utility methods
#     def i_inj(self, x, v):
#         X = self.local_view(x)
#         #v_t = self.v_t(x,v)
#         return (X['i_d'] + 1j * X['i_q']) * np.exp(1j*X['angle'])

#     def v_t(self, x, v):
#         return v[self.bus_idx_red['terminal']]

#     def s_e(self, x, v):
#         # Apparent power in p.u. (generator base units)
#         return self.v_t(x, v)*np.conj(self.i_inj(x, v))

#     def v_q(self,x,v):
#         return (self.v_t(x,v)*np.exp(-1j*self.local_view(x)['angle'])).imag
    
#     def v_d(self,x,v):
#         return (self.v_t(x,v)*np.exp(-1j*self.local_view(x)['angle'])).real     # Added
    
#     def p_e(self, x, v):
#         return self.s_e(x,v).real

#     def q_e(self, x, v):
#         return self.s_e(x,v).imag

#     # endregion


# class WindPowerSystem(GridSideConverter):         # Nødvendig?




    














# endregion


######################
# 
# 
 # """
        # V2
        # From block diagram electric drives PMSM with currents as state variables, dynamic analysis state space model forelesning
        
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
        # From block diagram electric drives PMSM with currents as state variables
        
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
    #     From block diagram electric drives PMSM with currents as state variables
        
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