from tops.dyn_models.utils import DAEModel
import numpy as np

# region GridSideConverter class
class GridSideConverter(DAEModel):
    """
    VSC PQ model - ref to VSC1 in tops/dyn_models/vsc1.py - cite Sjur Føyen

    Model of the DC link and Grid Side Converter.

    Its purpose is to control P and Q to the grid and to control the DC link voltage. 

    'vsc': {
        'VSC_PQ': [
            ['name',   'bus', 'S_n',   "p_ref",   "q_ref",   'Cdc',   'k_p',    'k_q',  'T_p',  'T_q',     'k_pll',   'T_pll',    'T_i',    'i_max'],
            ['VSC1',   'B1',    50,      1,         0,      0.1,       1,       1,     0.1,     0.1,        5,        1,         0.01,      1.2],
        ],
    }

    Parameters:
    w_n: Nominal frequency
    Cdc: DC link capacitance

    References:
    Vdc_ref
    p_ref
    q_ref
    id_ref
    iq_ref

    States:
    i_d
    i_q
    angle   PLL angle
    xpll    integral i pll
    vdc     DC link voltage

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.bus_idx = np.array(np.zeros(self.n_units), dtype=[(key, int) for key in self.bus_ref_spec().keys()])
        self.bus_idx_red = np.array(np.zeros(self.n_units), dtype=[(key, int) for key in self.bus_ref_spec().keys()])

    def connect_gsc(self, dclink):
        self.dclink = dclink

    # region Definitions

    def init_from_load_flow(self, x_0, v_0, S):
        X = self.local_view(x_0)

        self._input_values['p_ref_grid'] = self.par['p_ref_grid']       # unescassary?
        self._input_values['q_ref_grid'] = self.par['q_ref_grid']

        vg = v_0[self.bus_idx_red['terminal']]

        X['angle'] = np.angle(vg)
        X['x_pll'] = 0
        X['i_d'] = self.par['p_ref_grid']/abs(vg)
        X['i_q'] = self.par['q_ref_grid']/abs(vg)
        # X["vdc"] = 1.0
        # X['x_pref_adj'] = 0.0

        
    def state_derivatives(self, dx, x, v):
        dX = self.local_view(dx)
        X = self.local_view(x)
        par = self.par

        # self.p_ref = IPMSM.get_Pe() - self.vdc_controller.compute((self.vdc_ref - X["vdc"]))*X["vdc"]
        # self.q_ref = 0.0     # Should be updated from grid side command

        i_d_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (self.par["p_ref_grid"] * self.v_d(x,v) - self.par["q_ref_grid"] * self.v_q(x,v))
        i_q_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (self.par["p_ref_grid"] * self.v_q(x,v) + self.par["q_ref_grid"] * self.v_d(x,v))

        # i_d_ref = max(0, i_d_ref)
        # i_q_ref = max(0, i_q_ref)

        # i_d_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (X["p_ref"] * self.v_d(x,v) - X["q_ref"] * self.v_q(x,v))
        # i_q_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (X["p_ref"] * self.v_q(x,v) + X["q_ref"] * self.v_d(x,v))

        # Limiters
        i_ref = (i_d_ref+1j*i_q_ref)
        i_ref = i_ref*par['i_max']/np.maximum(par['i_max'],abs(i_ref))

        dX['i_d'][:] = 1 / (par['T_i']) * (i_ref.real - X['i_d'])
        dX['i_q'][:] = 1 / (par['T_i']) * (i_ref.imag - X['i_q'])

        dX['x_pll'][:] = par['k_pll'] / (par['T_pll']) * (self.v_q(x,v))
        dX['angle'][:] = X['x_pll']+par['k_pll']*self.v_q(x,v)
        
        # dX["vdc"][:] = (par["p_ref_gen"] - self.p_e(x,v)) / (par["Cdc"] * X["vdc"])
        # dX["x_pref_adj"][:] = (X["vdc"] - par["vdc_ref"]) * (par["K_p_dc"] / par["T_i_dc"])


        return    

    def load_flow_pq(self):
        return self.bus_idx['terminal'], -self.par['p_ref_grid']*self.par['S_n'], -self.par['q_ref_grid']*self.par['S_n']

    def int_par_list(self):
        return ['f']

    def bus_ref_spec(self):
        return {'terminal': self.par['bus']}

    # endregion

    def state_list(self):
        """
        All states in pu.

        i_d: d-axis current, first-order approximation of EM dynamics
        i_q: q-axis current -------------""--------------
        x_p: p-control integral
        x_q: q-control integral
        x_pll: pll q-axis integral
        angle: pll angle
        x_pref_adj: DC link voltage control integral
        pref_adj: DC link voltage control output
        vdc: DC link voltage
        """
        # return ['i_d', 'i_q', 'x_p', 'x_q', 'x_pll', 'angle', 'vdc', 'vdc_ref', 'x_pref_adj', 'pref_adj' ,"p_ref", "q_ref"]
        return ['i_d', 'i_q', 'x_pll', 'angle']#, 'vdc', 'x_pref_adj']

    def input_list(self):
        """
        All values in pu.
        p_ref: outer loop active power setpoint
        q_ref: outer loop reactive power setpoint
        """
        return ['p_ref_grid', 'q_ref_grid']

    def current_injections(self, x, v):
        i_n_r = self.par['S_n'] / self.sys_par['s_n']
        return self.bus_idx_red['terminal'], self.i_inj(x, v) * i_n_r

    # region Utility methods
    def i_inj(self, x, v):
        X = self.local_view(x)
        #v_t = self.v_t(x,v)
        return (X['i_d'] + 1j * X['i_q']) * np.exp(1j*X['angle'])

    def v_t(self, x, v):
        return v[self.bus_idx_red['terminal']]

    def s_e(self, x, v):
        # Apparent power in p.u. (generator base units)
        return self.v_t(x, v)*np.conj(self.i_inj(x, v))

    def v_q(self,x,v):
        return (self.v_t(x,v)*np.exp(-1j*self.local_view(x)['angle'])).imag
    
    def v_d(self,x,v):
        return (self.v_t(x,v)*np.exp(-1j*self.local_view(x)['angle'])).real     # Added
    
    def p_e(self, x, v):
        return self.s_e(x,v).real

    def q_e(self, x, v):
        return self.s_e(x,v).imag

    def vdc(self, x, v):
        return self.local_view(x)["vdc"]
    
    def i_dc(self, x, v):
        # X = self.local_view(x)
        return abs(self.i_inj(x,v)) ## dette er en rask løsning, kommer tilbake til det
    
    # def set_pref_grid(self, x, v, index):
    #     X = self.local_view(x)
    #     self.par["p_ref_grid"][index] = self.par["p_ref_gen"][index] + self.p_ref_adj(x,v)[index]

    def set_pref_grid(self, x, v, index, pref):
        X = self.local_view(x)
        self.par["p_ref_grid"][index] = pref
        # self.par["p_ref_gen"][index] = pref
        # self.par['p_ref_grid'][index] = self.par["p_ref_gen"][index] + self.p_ref_adj(x,v)[index]

    def p_ref_adj(self, x, v):
        X = self.local_view(x)
        return X["x_pref_adj"] + (X["vdc"] - self.par["vdc_ref"]) * self.par["K_p_dc"] #* X["vdc"]

    def set_qref_grid(self, Qref, index):
        self.par['q_ref_grid'][index] = Qref
    


    
    # endregion