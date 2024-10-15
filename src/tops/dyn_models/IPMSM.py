# Setting up a IPMSM as done in electrical drives.
# The goal is to include the dynamics of electrical torque and include the inverter control

import numpy as np
from tops.dyn_models.blocks import *
from tops.dyn_models.utils import DAEModel, output
import numpy as np
from tops.dyn_models.utils import DAEModel
import tops.utility_functions as dps_uf


# Må tror jeg kan unngå å bruke DAEModel ettersom at MSC er dekoblet fra det dynamiske nettet

class IMPMSM():
    
    def __init__(self, params):
        self.p = params
    

    def state_list(self):
        """
        Using a model with currents as state variables

        vd = rs*id + xd/wn * di_d/dt - n*xq*iq
        vq = rs*iq + xd/wn * di_q/dt - n*xd*id + n*Psi_m

        Psi_d = xd * id + Psi_m
        Psi_q = xq*iq

        te = Psi_m*iq - (xq-xd)*id*iq

        Adjusting the equations to fit a current controller

        Tq = xq/(wn*rs)
        Td = xd/(wn*rs)

        di_d/dt = -1/Td * id + wn/xd * ud
        di_q/dt = -1/Tq * iq + wn/xq * uq

        

        """
        return ["i_d", "i_q", "speed", "angle"]
    
    def input_list(self):
        return ['v_d', 'v_q']


    def derivatives(self, dx, x, v):

        dX = {}
        X = self.local_view(x)
        par = self.par

        Tq = par["x_q"]/(par["wn"]*par["rs"])
        Td = par["x_d"]/(par["wn"]*par["rs"])
        
        # Need a MSI model to acsess vd and vq
        # vd = par["rs"]*X["i_d"] + (par["xd"]/par["rs"]) * dX["i_d"] - n*xq*iq
        # vq = par["rs"]*X["i_q"] + xd/wn * di_q/dt - n*xd*id + n*Psi_m
        vd = 0
        vq = 1

        dX["i_d"][:] = -1/Td * X["i_d"] + (par["wn"]/par["rs"])*vd
        dX["i_q"][:] = -1/Td * X["i_q"] + (par["wn"]/par["rs"])*vq

        return


        def v_d(self, x, v):
            return 





