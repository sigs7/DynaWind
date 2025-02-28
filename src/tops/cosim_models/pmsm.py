# region PMSM class
from tops.cosim_models.pi_controller import PIController_LAH
from tops.cosim_models.machine_side_converter import MachineSideConverter



# Tror jeg kan unngå å bruke DAEModel ettersom at MSC er dekoblet fra det dynamiske nettet, stemmer
class PMSM:

    """
    Internally Permanent Magnet Synchrnous Machine

    """
    
    def __init__(self, pmsm_params : dict):

        """
        # pmsm_params = {
        #     "s_n" : 15e6,     # 15 MW
        #     "w_r" : 0.792*(60/(2*np.pi)),  # Rated speed in rpm
        #     "U_n" : 4770.34,    # Rated phase voltage
        #     "I_n" : 1084.55,    # Nominal phase current
        #     "T_r" : 21.03e3,   # kNm torque        
        #     "rs": 0.03,     # Stator resistance
        #     "x_d": 0.4,     # Stator d-axis inductance
        #     "x_q": 0.4,     # Stator q-axis inductance
        #     "Psi_m": 0.9,   # Magnetic flux
        # }

        # https://www.nrel.gov/docs/fy20osti/75698.pdf

        """
        

        ### Might need to adress an eventuall gearbox ratio
        
        # Assigning the basic parameters of the PMSM
        self.params = pmsm_params

        # Initiating some basic initial values of operation, should be updated later based of load flow solution
        self.speed = 0.0
        self.torque_ref = 0.0

        # Initialize PI controllers for i_d, i_q, and speed with parameters adjusted for a larger timestep
        self.pi_controller_id = PIController_LAH(kp=1, ti=0.1)
        self.pi_controller_iq = PIController_LAH(kp=1, ti=0.1)

        # Initiate reference values
        self.i_d_ref = 0.0          # Må kanskje endres senere
        self.i_q_ref = self.torque_ref / self.params["Psi_m"]

        self.i_d = 0.0
        self.i_q = self.i_q_ref
        self.theta = 0.0

    def connect_pmsm(self, msc):
        self.msc = msc

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

        # Motor convention
        dX["i_d"] = (self.msc.v_d - p["rs"]*self.i_d + psi_q*self.speed) * (p["w_n"]/p["x_d"])
        dX["i_q"] = (self.msc.v_q - p["rs"]*self.i_q - psi_d*self.speed) * (p["w_n"]/p["x_q"])
        dX["theta"] = self.speed

        # Speed is here exluded from the derivatives, as it is updated from the FAST FMU

        return dX
       
    # endregion
    
    # region Controll functions

    def update_torque_control(self, fast):
        """
        Method to update the reference torque based on the speed error
        Input: dt - time step
        
        Output: self.i_q_ref - updated reference q-axis current
        """
        # fast.fmu: FMU2Slave

        self.torque_ref = -1 * fast.fmu.getReal([fast.vrs['GenTq']])[0]/self.params["T_r"]       # Torque reference from FAST FMU in pu

        # Limit torque reference
        self.torque_ref = self.clamp_value(self.torque_ref, max_value=0, min_value=-1.5)        

        # Calculate the required i_q_ref to produce the needed torque (algebraic)
        self.i_q_ref = self.torque_ref / self.params["Psi_m"]


    def update_current_control(self, dt):
        p = self.params

        # Compute errors
        error_id = (self.i_d_ref - self.i_d)
        error_iq = (self.i_q_ref - self.i_q)

        ### Decoupling and current control ###
        # Compute decoupled voltage control signals
        self.v_dII = self.i_q * p["x_q"] * self.speed               
        self.v_qII = ((self.i_d*p["x_d"]) + p["Psi_m"]) * self.speed

        # I_d current control
        self.v_d_ctrl = self.pi_controller_id.compute(error_id, dt) - self.v_dII
        self.v_d_ctrl = self.clamp_value(self.v_d_ctrl, max_value=2, min_value=-2)
        
        # I_q current control
        self.v_q_ctrl = self.pi_controller_iq.compute(error_iq, dt) + self.v_qII
        self.v_q_ctrl = self.clamp_value(self.v_q_ctrl, max_value=2, min_value=-2)
    
    # endregion

    def step_pmsm(self, fast, time : float, step_size : float):

        from tops.cosim_models.fast import FAST     # Såkalt "lazy import" for å unngå sirkulær import

        # Update states from fast FMU
        self.speed = fast.fmu.getReal([fast.vrs['RotSpeed']])[0] / self.params["w_n"]
        self.SPEED = self.speed * self.params["w_n"]        # Mulig drit kodepraksis menmen     (rpm)

        # Torque control
        self.update_torque_control(fast)

        # Update reference voltages using PI controllers
        self.update_current_control(step_size)

        # Calculate the derivatives
        dX = self.derivatives()

        # # Update the states using Euler integration
        self.i_d += dX["i_d"] * step_size
        self.i_q += dX["i_q"] * step_size
        self.theta += dX["theta"] * step_size



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


    def get_t_e(self):
        p = self.params
        psi_q = self.i_q*p["x_q"]
        psi_d = self.i_d*p["x_d"] + p["Psi_m"]

        return psi_d*self.i_q - psi_q*self.i_d
    
    def get_T_e(self):
        return self.get_t_e()*self.params["T_r"]

    def get_v_d(self):
        return self.msc.v_d
    
    def get_v_q(self):
        return self.msc.v_q

    def get_speed(self):
        return self.speed
    
    def get_p_e(self):
        return (-3/2)*(self.msc.v_d*self.i_d + self.msc.v_q*self.i_q)      # -3/2 fac to adjust to three-phase power and change direction as injection

    def get_P_e(self):
        return self.get_p_e()*self.params["s_n"]

    def get_q_e(self):
        return (-3/2)*(self.msc.v_q*self.i_d - self.msc.v_d*self.i_q)      # 3/2 fac?
    
    def get_Q_e(self):
        return self.get_q_e()*self.params["s_n"]

    def get_Pm(self):
        return self.speed * self.primemover.torque

    # endregion


    def set_converter_voltages(self, v_d_ctrl, v_q_ctrl, dt):
    # def set_converter_voltages(self, v_d_ctrl, v_q_ctrl, dt):
        # self.converter.set_reference_voltages(v_d_ctrl, v_q_ctrl)
        # self.msc.v_d_ref = v_d_ctrl
        # self.msc.v_q_ref = v_q_ctrl
        # self.msc.apply_voltages(v_d_ctrl, v_q_ctrl)
        self.msc.update_voltages(v_d_ctrl, v_q_ctrl, dt)