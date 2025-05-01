# WindTurbine class to combine all the classes and initate insteances as structured as possible

from tops.cosim_models.fast import FAST
from tops.cosim_models.pmsm import PMSM
from tops.cosim_models.dclink import DClink


import numpy as np

class WindTurbine(FAST, PMSM, DClink):
    def __init__(self, name : str, index : int, gsc_control : str, pmsm_params : dict = None, ideal_generator_params = None, FAST_params : dict = None, dclink_params : dict = None):

        if FAST_params is None:
            FAST_params = {
                'fmu_filename' : 'fast.fmu',
                'start_time' : 0.0,
                'mode' : 3
            }

        if pmsm_params is None:
            pmsm_params = {
            "s_n"   : 15e3*1.11,     # 15 000 kW 15e3, # 
            "rpm_n" : 0.792*(60/(2*np.pi)),  # Rated speed rpm       # 0.792 rad/s
            "w_r"   : 0.792,  # Rated speed rad/s
            "f_n"   : 12.6,    # Rated frequency
            "U_n"   : 4770.34,    # Rated phase voltage
            "I_n"   : 1084.55,    # Nominal phase current
            "T_r"   : 21.03e3,  #30e3,    #21.03e3,   # kNm torque         #### NOTE #### Might need to change to Nm      30e3, 
            "r_s"   : 0.03,       # Stator resistance 0.03
            "x_d"   : 0.3,        # 0.4,#0.608,     # Stator d-axis inductance
            "x_q"   : 0.3,        #1,#0.608,     # Stator q-axis inductance
            "Psi_m" : 1, #1.5,#0.9255,     # Magnetic flux
            "Poles" : 200,        # Number of poles
            }
        

        if dclink_params is None:
            dclink_params = {
                "vdc_0"              : 1.5,
                "vdc_ref"            : 1.5,
                "chopper_threshold"  : 1.55,
                "Cdc"                : 5e-3,  #0.6, # 1.2, #0.45, #1.2,   #0.1,      0.45,             # Zbase = (sqrt(3)*4770.34)**2 / 15e6 = 4.55      Cbase = 1 / 2pi*f_n*Zbase = 0.00275      Cdc = 0.01 F / 0.0082 = 1.2 
                "K_p_dc"             : 2,
                "T_i_dc"             : 0.2,
                "chopper_resistance" : 2, #2,
            }

        # if ideal_generator_params is None:
        #     ideal_generator_params = {
        #         "T_r" : 21.03e3,   # kNm torque
        #         "T_g" : 0.1e-3,
        #         "rpm_n" : 0.792*(60/(2*np.pi)),  # Rated speed in rpm
        #         "s_n" : 15e6,
        #     }

        self.fast = FAST(FAST_params)
        self.dclink = DClink(dclink_params)
        self.pmsm = PMSM(pmsm_params, self.fast)

        # Connect the models to eachother
        self.connect_models()
        
        # Assigning the name and index of the wind turbine to be used in co-simulation storage and plotting of results
        self.name = name
        self.index = index

        self.gsc_control = gsc_control
        self.gsc_sn = None

    def connect_models(self):
        self.dclink.connect_dclink(self.pmsm)
        self.pmsm.connect_pmsm(self.dclink)
        self.fast.connect_fmu(self)


    def step_windturbine_multirate(self, ps, x, v, time : float, step_size_mech : float, step_size_elec : float, results):
        p_gsc = self.calculate_p_gsc(ps, x, v)
        self.p_gsc = p_gsc      # Debugging printing

        elec_steps_per_mech_step = int(step_size_mech / step_size_elec)


        for _ in range(elec_steps_per_mech_step):

            # results.store_time_elec(time)
            # results.store_pmsm_results(self)
            # results.store_dclink_results(self, ps, x, v)
            # results.store_fmu_results(self)

            self.pmsm.step_pmsm(dclink=self.dclink, time=time, step_size_elec=step_size_elec)
            self.dclink.step_dclink(time, step_size_elec, self.calculate_p_gsc(ps, x, v))

            # results.store_pmsm_results(self)
            # results.store_dclink_results(self, ps, x, v)
            # results.store_time_elec(time)
            # results.store_fmu_results(self)
            
            time += step_size_elec
        
        self.fast.step_fmu(self.pmsm, ps, v, x, time, step_size_mech)



        if self.gsc_control == "PQ":
            ps.vsc["GridSideConverter_PQ"].set_pref_grid(x, v, index=self.index, pref=self.dclink.gsc_p_ref(), power_rating=self.pmsm.params["s_n"])
        if self.gsc_control == "PV":
            ps.vsc["GridSideConverter_PV"].set_pref_grid(x, v, index=self.index, pref=self.dclink.gsc_p_ref(), power_rating=self.pmsm.params["s_n"])

        self.pref_dclink_pu = self.dclink.gsc_p_ref()

        if self.gsc_sn == None:
            self.get_converter_rating(ps)

    # ps info
    def get_converter_rating(self, ps):
        if self.gsc_control == "PQ":
            self.gsc_sn = ps.vsc["GridSideConverter_PQ"].par["S_n"][self.index]
        if self.gsc_control == "PV":
            self.gsc_sn = ps.vsc["GridSideConverter_PV"].par["S_n"][self.index]

    def calculate_p_gsc(self, ps, x, v):
        if self.gsc_control == "PQ":
            return ps.vsc["GridSideConverter_PQ"].p_e(x, v)[self.index] * ps.vsc["GridSideConverter_PQ"].par["S_n"][self.index] / (self.pmsm.params["s_n"] * 1e-3)
        if self.gsc_control == "PV":
            return ps.vsc["GridSideConverter_PV"].p_e(x, v)[self.index] * ps.vsc["GridSideConverter_PV"].par["S_n"][self.index] / (self.pmsm.params["s_n"] * 1e-3)






    # def step_windturbine(self, ps, time : float, step_size : float, x, v):
    #     p_gsc = self.calculate_p_gsc(ps, x, v)
    #     self.p_gsc = p_gsc      # Debugging printing
    #     self.fast.step_fmu(self.pmsm, time, step_size)
    #     self.pmsm.step_pmsm(fast=self.fast, time=time, step_size_elec=step_size)
    #     self.dclink.step_dclink(time, step_size, p_gsc)

    #     if self.gsc_control == "PQ":
    #         ps.vsc["GridSideConverter_PQ"].set_pref_grid(x, v, index=self.index, pref=self.dclink.gsc_p_ref(), power_rating=self.pmsm.params["s_n"])
    #     if self.gsc_control == "PV":
    #         ps.vsc["GridSideConverter_PV"].set_pref_grid(x, v, index=self.index, pref=self.dclink.gsc_p_ref(), power_rating=self.pmsm.params["s_n"])

    #     self.pref_dclink_pu = self.dclink.gsc_p_ref()

    #     if self.gsc_sn == None:
    #         self.get_converter_rating(ps)


