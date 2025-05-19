# region PMSM class
from tops.cosim_models.pi_controller import Controller
import numpy as np

class PMSM:

    """
    Permanent Magnet Synchrnous Machine

    """
    
    def __init__(self, pmsm_params : dict, fast):

        """
        # pmsm_params = {
        #     "s_n" : 15e6,     # 15 MW
        #     "w_r" : 0.792*(60/(2*np.pi)),  # Rated speed in rpm
        #     "U_n" : 4770.34,    # Rated phase voltage
        #     "I_n" : 1084.55,    # Nominal phase current
        #     "T_r" : 21.03e3,   # kNm torque        
        #     "r_s": 0.03,     # Stator resistance
        #     "x_d": 0.4,     # Stator d-axis inductance
        #     "x_q": 0.4,     # Stator q-axis inductance
        #     "Psi_m": 0.9,   # Magnetic flux
        #     "p": 200,         # Number of poles
        # }

        # https://www.nrel.gov/docs/fy20osti/75698.pdf

        """
    
        ### Might need to adress an eventuall gearbox ratio
        # Not all parameters are used in the current implementation, but they are included for future reference and expansion of the model.
        
        # Assigning the basic parameters of the PMSM
        self.params = pmsm_params
        p = self.params

        # Initiating some basic initial values of operation
        self.speed = fast.fmu.getReal([fast.vrs['GenSpeed']])[0] / self.params["rpm_n"]
        self.SPEED = self.speed * self.params["rpm_n"]        # Mulig drit kodepraksis menmen     (rpm)
        self.torque_ref_raw = fast.fmu.getReal([fast.vrs['GenTq']])[0] / self.params["T_r"]       # Torque reference from FAST FMU in pu
        self.torque_ref = self.torque_ref_raw

        self.torque_ref_applied = 0.0  # This will hold the delayed reference


        # Initiate starting state of ELectric Drive
        self.state = "normal"

        # Initiate reference values
        self.i_d_ref = 0.0          
        self.i_q_ref = self.torque_ref / self.params["Psi_m"]

        self.i_d = self.i_d_ref
        self.i_q = self.i_q_ref

        self.v_d = 0.0
        self.v_q = 0.0

        self.v_dII = self.i_q * p["x_q"] * self.speed               
        self.v_qII = ((self.i_d*p["x_d"]) + p["Psi_m"]) * self.speed

        # Possible expansion of future model, not used in current configuration
        self.theta_mech = 0.0
        self.theta_elec = 0.0

        # Initiate controllers
        self.pi_controller_id = Controller(kp=3, ki=3/0.02, td=0.0, kaw=15.0, output_limit=10)
        self.pi_controller_iq = Controller(kp=3, ki=3/0.02, td=0.0, kaw=15.0, output_limit=10)
        self.pi_controller_vdc = Controller(kp=0.5, ki=5, td=0.0, kaw=15.0, output_limit=10)

        # Add a filtered electric power to be fed into turbine controller??

    # Connect the PMSM to the DC-link
    def connect_pmsm(self, dclink):
        self.dclink = dclink

    # region Step

    def step_pmsm(self,dclink, time : float, step_size_elec : float):

            # Update reference voltages using PI controllers
            self.update_current_control(dclink=dclink, dt=step_size_elec, time=time)

            # Calculate the derivatives
            dX = self.derivatives()

            # # Update the states using Euler nunmerical integration
            self.i_d += dX["i_d"] * step_size_elec
            self.i_q += dX["i_q"] * step_size_elec

    # endregion

    # region Derivatives

    def derivatives(self):
        """
        V3
        From block diagram electric drives PMSM SIMULINK, currents as state variables

        """

        dX = {}
        p = self.params

        # Old version
        # psi_q = self.i_q*p["x_q"]
        # psi_d = self.i_d*p["x_d"] + p["Psi_m"]

        # # Motor convention
        # dX["i_d"] = (self.v_d - p["r_s"]*self.i_d + psi_q*self.speed) * (p["w_n"]/p["x_d"])
        # dX["i_q"] = (self.v_q - p["r_s"]*self.i_q - psi_d*self.speed) * (p["w_n"]/p["x_q"])

        w_e = self.w_e()        # Electrical per unit rotational speed is essentially the same as pu mechanical speed
        L_d = p["L_d"]
        L_q = p["L_q"]
        R_s = p["r_s"]

        dX["i_d"] = (self.v_d - R_s*self.i_d + L_q*self.i_q*w_e) / L_d
        dX["i_q"] = (self.v_q - R_s*self.i_q - (L_d*self.i_d + p["Psi_m"])*w_e) / L_q

        # Mechanical speed is here exluded from the derivatives as all mechanical states are solved in the OpenFAST FMU

        return dX
       
    # endregion
    
    # region Current control

    def update_current_control(self, dclink, dt, time):
        p = self.params

        # Run the torque reference through the Electric Drive control and externally clamp the value
        self.torque_ref = self.drive_control(time, dclink, self.torque_ref_raw, dt)
        self.torque_ref = self.clamp_value(self.torque_ref, max_value=0, min_value=-1.5)

        # Remove drive control for now, as it is not implemented yet
        # self.torque_ref = self.torque_ref_raw

        # Calculate the required i_q_ref to produce the needed torque (algebraic)
        self.i_q_ref = self.torque_ref / self.params["Psi_m"]

        # Number of poles in torque equation confuses me..
        # self.i_q_ref = (2*self.torque_ref) / (3*self.params["Psi_m"]*(self.params["Poles"]/2))

        # Clamping current references, this should be done based on I = sqrt(id^2 + iq^2) and not on id and iq individually, for now it is a quick fix
        self.i_q_ref = self.clamp_value(self.i_q_ref, max_value=1.2, min_value=-1.2)
        self.i_d_ref = self.clamp_value(self.i_d_ref, max_value=1.2, min_value=-1.2)

        # Compute errors
        error_id = (self.i_d_ref - self.i_d)
        error_iq = (self.i_q_ref - self.i_q)

        ### Decoupling and current control ###
        self.v_dII = p["L_q"]*self.i_q*self.w_e()
        self.v_qII = (p["L_d"]*self.i_d + p["Psi_m"])*self.w_e()

        # I_d current control
        self.v_d_ctrl = self.pi_controller_id.compute(error_id, dt, time) - self.v_dII
        # self.v_d_ctrl = self.clamp_value(self.v_d_ctrl, max_value=1.5, min_value=0)
        
        # I_q current control
        self.v_q_ctrl = self.pi_controller_iq.compute(error_iq, dt, time) + self.v_qII
        # self.v_q_ctrl = self.clamp_value(self.v_q_ctrl, max_value=1.5, min_value=0)

        # Skipping the MSC voltage delay, applying voltages imidiately
        self.v_d = self.v_d_ctrl
        self.v_q = self.v_q_ctrl

    # endregion

    # region Drive control
    def drive_control(self, time, dclink, torque_ref_raw, dt):
        """
        Method to control the drive based on current state of the system

        Input: dclink - DC-link voltage
               Frequency - grid frequency? 
        
        Output: torque_ref - updated torque reference

        Goal:
          - Generator controlls the DC-link voltage during fault conditions via regulating the torque
          - Add ramp rate limiter for the torque reference
          - Grid supporting capabilities? Inertia emulation - would need grid measurements?


        First step: Find the state of the system to detect faults and determine drive behavior
        Second step: Alter torque reference based on the state of the system

        """

        # Configurables
        vdc = dclink.vdc                                # DC-link voltage [pu]
        vdc_ref = dclink.params["vdc_ref"]              # Reference voltage for DC-link [pu]
        vdc_th = dclink.params["chopper_threshold"]     # Threshold voltage for fault detection [pu]
        hold_duration = 3                               # Duration to hold the reduced torque after fault is cleared [s]
        ramp_rate_fault = 100                           # Maximum ramp rate for torque reference during fault [pu/s]
        ramp_rate_hold = 1                              # Maximum ramp rate for torque reference during hold [pu/s]
        ramp_rate_normal = 2                            # Maximum ramp rate for torque reference during normal operation [pu/s]


        # Start by adressing the state of the system
        # --- Fault detection ---
        if self.state == "normal":
            if vdc > vdc_th:
                self.state = "fault"
                self.torque_at_fault = torque_ref_raw
                self.fault_timestamp = time
                print(f"[{time:.3f}s] FAULT DETECTED")


        # --- Fault cleared, enter hold ---    
        elif self.state == "fault":
            if vdc <= vdc_th:
                self.state = "hold"
                self.hold_start_time = time
                print(f"[{time:.3f}s] FAULT CLEARED → HOLD")


        # --- Hold back to normal operation torque ---
        elif self.state == "hold":
            if time - self.hold_start_time > hold_duration:
                self.state = "normal"
                print(f"[{time:.3f}s] HOLD FINISHED → NORMAL")
                

        # Now describe the drive behavior based on state:
        if self.state == "normal":
            return self.ramp_to(self.torque_ref, torque_ref_raw, ramp_rate_normal, dt)
        
        elif self.state == "fault":
            # Initiate VDC regulation
            # error_vdc = -(vdc_ref - vdc)
            # self.torque_ref_pivdc = self.pi_controller_vdc.compute(error_vdc, dt, time) + self.torque_ref_raw
            # return self.ramp_to(current_value=self.torque_ref, target_value=self.torque_ref_pivdc, max_rate=ramp_rate_fault, dt=dt)
            return self.ramp_to(self.torque_ref, 0.66*self.torque_at_fault, ramp_rate_fault, dt)

        elif self.state == "hold":
            # Initiate VDC regulation
            # error_vdc = -(vdc_ref - vdc)
            # self.torque_ref_pivdc = self.pi_controller_vdc.compute(error_vdc, dt, time) + self.torque_ref_raw
            # return self.ramp_to(current_value=self.torque_ref, target_value=self.torque_ref_pivdc, max_rate=ramp_rate_hold, dt=dt)
            return self.ramp_to(self.torque_ref, 0.85*self.torque_at_fault, ramp_rate_hold, dt)

    # endregion



    # region Utility functions

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

    def ramp_to(self, current_value: float, target_value: float, max_rate: float, dt: float) -> float:
        """
        Linearly ramp a value toward a target at a limited rate.

        Parameters:
        - current_value: the current value (e.g., torque_ref)
        - target_value: the desired value to reach
        - max_rate: maximum change per second
        - dt: simulation timestep in seconds

        Returns:
        - new_value: updated value after applying the ramp
        """
        delta = target_value - current_value
        max_delta = max_rate * dt

        if abs(delta) <= max_delta:
            return target_value
        else:
            return current_value + max_delta * (1 if delta > 0 else -1)

    # def W_e(self):
    #     return self.params["Poles"] * self.params["w_r"] / 2
    
    def w_e(self):
        return self.speed       # equal to speed in pu?

    # def L_d(self):
    #     return self.params["x_d"] / self.W_e()

    # def L_q(self):
    #     return self.params["x_q"] / self.W_e() 

    def t_e(self):
        return ((self.params["Psi_m"] + (self.w_e()*self.params["L_d"] - self.w_e()*self.params["L_q"])*self.i_d))*self.i_q
    
    def T_e(self):
        return self.t_e()*self.params["T_r"]
    
    def p_e(self):
        return -1*(self.v_d*self.i_d + self.v_q*self.i_q)      # -3/2 fac?? to adjust to three-phase power and change direction as injection

    def P_e(self):
        return self.p_e()*self.params["s_n"]            # kW
    
    def q_e(self):
        return -1*(self.v_q*self.i_d - self.v_d*self.i_q)      # 3/2 fac? Not when entire system is in pu
    
    def Q_e(self):
        return self.q_e()*self.params["s_n"]

    def P_mech(self):
        return -self.t_e() * self.params["T_r"] * self.speed * self.params["w_r"]       # Mechanical power in kW

















    # endregion

    # region Old code
    # def step_pmsm(self, fast, time : float, step_size_elec : float):

    #     from tops.cosim_models.fast import FAST     # Såkalt "lazy import" for å unngå sirkulær import

    #     # Update states from fast FMU
    #     self.speed = fast.fmu.getReal([fast.vrs['GenSpeed']])[0] / self.params["rpm_n"]
    #     self.SPEED = self.speed * self.params["rpm_n"]        # Mulig drit kodepraksis menmen     (rpm)

    #     # Torque control
    #     self.update_torque_control(fast)

    #     # Update reference voltages using PI controllers
    #     self.update_current_control(dt=step_size_elec, time=time)

    #     # Calculate the derivatives
    #     dX = self.derivatives()
    #     p = self.params

    #     # # Update the states using Euler integration
    #     self.i_d += dX["i_d"] * step_size_elec
    #     self.i_q += dX["i_q"] * step_size_elec
    #     self.theta_mech += dX["theta_mech"] * step_size_elec


    #     # Normalize the mechanical angle to be within 0 and 2*pi
    #     self.theta_mech = self.theta_mech % (2 * np.pi)

    #     if self.theta_mech < 0:
    #         self.theta_mech += 2 * np.pi

    #     # Calculate the electrical angle
    #     self.theta_elec = (self.theta_mech * p["Poles"] / 2 ) % (2 * np.pi)

    # def set_converter_voltages(self, v_d_ctrl, v_q_ctrl, dt):
    # # def set_converter_voltages(self, v_d_ctrl, v_q_ctrl, dt):
    #     # self.converter.set_reference_voltages(v_d_ctrl, v_q_ctrl)
    #     # self.msc.v_d_ref = v_d_ctrl
    #     # self.msc.v_q_ref = v_q_ctrl
    #     # self.msc.apply_voltages(v_d_ctrl, v_q_ctrl)
    #     self.msc.update_voltages(v_d_ctrl, v_q_ctrl, dt)