# WindTurbine class to combine all the classes and initate insteances as structured as possible

from tops.cosim_models.fast import FAST
from tops.cosim_models.pmsm import PMSM
from tops.cosim_models.dclink import DClink
from tops.cosim_models.machine_side_converter import MachineSideConverter
from tops.cosim_models.grid_side_converter import GridSideConverter

import numpy as np

class WindTurbine(FAST, PMSM, DClink, MachineSideConverter, GridSideConverter):
    def __init__(self, name : str, index : int, pmsm_params : dict = None, msc_params : dict = None, FAST_params : dict = None, dclink_params : dict = None):

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
            }
        
        if msc_params is None:
            msc_params = {
            # Converter parameters
            "T_conv" : 1/(300),
            "v_q_0" : 0.0,
            "v_d_0" : 0.0
            }

        if dclink_params is None:
            dclink_params = {
                "vdc_0": 1.0,
                "vdc_ref": 1.0,
                "chopper_threshold": 1.3,
                "Cdc": 2,
                "K_p_dc": 1.0,
                "T_i_dc": 0.2,
                "chopper_resistance": 10,
            }

        self.fast = FAST(FAST_params)
        self.msc = MachineSideConverter(msc_params)
        self.dclink = DClink(params=dclink_params)
        self.pmsm = PMSM(pmsm_params)

        # Connect the models to eachother
        self.connect_models()
        

        # Assigning the name and index of the wind turbine to be used in co-simulation storage and plotting of results
        self.name = name
        self.index = index

    def connect_models(self):
        self.msc.connect_msc(self.dclink, self.pmsm)
        self.dclink.connect_dclink(self.msc, self.pmsm)
        self.pmsm.connect_pmsm(self.msc)

        
    def step_windturbine(self, time : float, step_size : float, p_gsc : float):
        self.fast.step_fmu(self.pmsm, time, step_size)
        self.pmsm.step_pmsm(self.fast, time, step_size)
        self.msc.step_msc(time, step_size)
        self.dclink.step_dclink(time, step_size, p_gsc)

    def set_gsc_power_reference(self):
        return self.dclink.gsc_p_ref()
