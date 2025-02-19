# WindTurbine class to combine all the classes and initate insteances as structured as possible

from tops.cosim_models.fast import FAST
from tops.cosim_models.pmsm import PMSM
import numpy as np

class WindTurbine(FAST, PMSM):
    def __init__(self, name : str, index : int, pmsm_params : dict = None, MSC_params : dict = None, FAST_params : dict = None):

        if FAST_params is None:
            FAST_params = {
                'fmu_filename' : 'fast.fmu',
                'start_time' : 0.0,
                'mode' : 3
            }

        if pmsm_params is None:
            pmsm_params = {
            "s_n" : 15e3,     # 15 000 kW
            "w_n" : 0.792*(60/(2*np.pi)),  # Rated speed rpm       # 0.792 rad/s
            "U_n" : 4770.34,    # Rated line voltage
            "I_n" : 1084.55,    # Nominal phase current
            "T_r" : 21.03e3,   # kNm torque         #### NOTE #### Might need to change to Nm      30e3, 
            "rs": 0.03,     # Stator resistance
            "x_d": 0.4,     # Stator d-axis inductance
            "x_q": 0.4,     # Stator q-axis inductance
            "Psi_m": 0.9,   # Magnetic flux

            # Converter parameters
            "T_conv" : 1/(300),
              "v_q_0" : 0.0,
              "v_d_0" : 0.0
            }

        self.fast = FAST(FAST_params)
        self.pmsm = PMSM(pmsm_params)

        # Assigning the name and index of the wind turbine to be used in co-simulation storage and plotting of results
        self.name = name
        self.index = index
        
    def step_windturbine(self, time : float, step_size : float):
        self.fast.step_fmu(self.pmsm, time, step_size)
        self.pmsm.step_pmsm(self.fast, time, step_size)
