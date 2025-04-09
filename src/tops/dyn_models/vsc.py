# from tops.dyn_models.blocks import *
# from .pll import PLL1
from tops.dyn_models.utils import DAEModel
import numpy as np

# from IPMSM_drives import IPMSM

# region GridSideConverter class
class GridSideConverter_PQ(DAEModel):
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

        # pe = 3/2 /vq*iq + vd*id)
        i_d_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (self.par["p_ref_grid"] * self.v_d(x,v) - self.par["q_ref_grid"] * self.v_q(x,v))  # ---      # Byttet om + og - ?
        i_q_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (self.par["p_ref_grid"] * self.v_q(x,v) + self.par["q_ref_grid"] * self.v_d(x,v))  # +++

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
        return abs(self.i_inj(x,v)) ## dette er en rask løsning, kommer tilbake til det

    def set_pref_grid(self, x, v, index, pref, power_rating):
        X = self.local_view(x)
        self.par["p_ref_grid"][index] = pref * power_rating*1e-3 / self.par["S_n"]          # MVA / MVA

    def p_ref_adj(self, x, v):
        X = self.local_view(x)
        return X["x_pref_adj"] + (X["vdc"] - self.par["vdc_ref"]) * self.par["K_p_dc"] #* X["vdc"]

    def set_qref_grid(self, Qref, index):
        self.par['q_ref_grid'][index] = Qref
    
    
    # endregion

class GridSideConverter_PV(DAEModel):
    """
    Instantiate:
    'vsc': {
            'VSC_PV': [
                ['name', 'bus', 'S_n', 'p_ref', 'V', 'k_p', 'k_v', 'T_p', 'T_v', 'k_pll', 'T_pll', 'T_i', 'i_max'],
                ['VSC1', 'B1',    50,   0.8,    0.93,    1,      1,   0.1,   0.1,     5,      1,      0.01,    1.2],
            ],
        }
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.bus_idx = np.array(np.zeros(self.n_units), dtype=[(key, int) for key in self.bus_ref_spec().keys()])
        self.bus_idx_red = np.array(np.zeros(self.n_units), dtype=[(key, int) for key in self.bus_ref_spec().keys()])

    # region Definitions

    def load_flow_pv(self):
        return self.bus_idx['terminal'], -self.par['p_ref_grid']*self.par['S_n'], self.par['v_ref_grid']

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
        x_v: v-control integral
        x_pll: pll q-axis integral
        angle: pll angle
        """
        return ['i_d', 'i_q', 'x_p', 'x_v', 'x_pll', 'angle']

    def input_list(self):
        """
        All values in pu.

        p_ref: outer loop active power setpoint
        v_ref: outer loop voltage setpoint
        """
        return ['p_ref_grid', 'v_ref_grid']

    def state_derivatives(self, dx, x, v):
        dX = self.local_view(dx)
        X = self.local_view(x)
        par = self.par

        # Active and reactive power errors
        dp = self.par["p_ref_grid"] - self.p_e(x, v)
        dv = self.par["v_ref_grid"] - abs(self.v_t(x, v))

        # # Gain Scheduling: Adjust k_p dynamically based on dp magnitude
        # k_p_eff = par["k_p"] * (1 + 0.5 * np.abs(dp))  

        # # Improved reference current calculation
        # i_d_ref = np.clip(k_p_eff * dp + X['x_p'], 0, par['i_max'])
        # i_q_ref = -dv * par['k_v'] - X['x_v']

        self.q_ref_grid = - self.par["k_v"] * dv - X['x_v']

        # pe = vq*iq + vd*id)
        i_d_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (self.par["p_ref_grid"] * self.v_d(x,v) - self.q_ref_grid * self.v_q(x,v))  # ---      # Byttet om + og - ?
        i_q_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (self.par["p_ref_grid"] * self.v_q(x,v) + self.q_ref_grid * self.v_d(x,v))  # +++




        # Apply current limiters
        i_ref = (i_d_ref + 1j * i_q_ref)
        i_ref = i_ref * par['i_max'] / np.maximum(par['i_max'], abs(i_ref))

        # Limit integrator values
        X['x_p'] = np.clip(X['x_p'], -par['i_max'], par['i_max'])
        X['x_v'] = np.clip(X['x_v'], -par['i_max'], par['i_max'])

        # Update state derivatives
        dX['i_d'][:] = 1 / (par['T_i']) * (i_ref.real - X['i_d'])
        dX['i_q'][:] = 1 / (par['T_i']) * (i_ref.imag - X['i_q'])
        dX['x_p'][:] = par['k_p'] / par['T_p'] * dp
        dX['x_v'][:] = par['k_v'] / par['T_v'] * dv
        dX['x_pll'][:] = par['k_pll'] / par['T_pll'] * self.v_q(x, v)
        dX['angle'][:] = X['x_pll'] + par['k_pll'] * self.v_q(x, v)

        return

        # dp = self.par["p_ref_grid"] - self.p_e(x, v)
        # dv = self.par["v_ref_grid"] - abs(self.v_t(x, v))

        # # i_d_ref = np.clip(dp * par['k_p'] + X['x_p'], 0, par['i_max'])      # Do not allow negative d-axis current
        # # i_d_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (self.par["p_ref_grid"] * self.v_d(x,v) - self.par["q_ref_grid"] * self.v_q(x,v))  # ---      # Byttet om + og - ?

        # i_d_ref = np.clip(dp * par['k_p'] + X['x_p'], 0, par['i_max'])
        # i_q_ref = -dv * par['k_v'] - X['x_v']

        # # Limiters
        # i_ref = (i_d_ref+1j*i_q_ref)
        # i_ref = i_ref*par['i_max']/np.maximum(par['i_max'],abs(i_ref))

        # # X['x_p'] = np.maximum(np.minimum(X['x_p'],par['i_max']),-par['i_max'])
        # # X['x_v'] = np.maximum(np.minimum(X['x_v'],par['i_max']),-par['i_max'])

        # X['x_p'] = np.clip(X['x_p'], -par['i_max'], par['i_max'])
        # X['x_v'] = np.clip(X['x_v'], -par['i_max'], par['i_max'])

        # dX['i_d'][:] = 1 / (par['T_i']) * (i_ref.real - X['i_d'])
        # dX['i_q'][:] = 1 / (par['T_i']) * (i_ref.imag - X['i_q'])
        # dX['x_p'][:] = par['k_p']/par['T_p'] * dp
        # dX['x_v'][:] = par['k_v']/par['T_v'] * dv
        # dX['x_pll'][:] = par['k_pll']/par['T_pll'] * (self.v_q(x,v))
        # dX['angle'][:] = X['x_pll']+par['k_pll']*self.v_q(x,v)
        # return

    def init_from_load_flow(self, x_0, v_0, S):
        X = self.local_view(x_0)

        self._input_values['p_ref_grid'] = self.par['p_ref_grid']
        self._input_values['v_ref_grid'] = self.par['v_ref_grid']

        v0 = v_0[self.bus_idx_red['terminal']]
        # s0 = S/ self.par["S_n"]
        # /self.par['S_n']
        # s00 = s0 / (self.par["S_n"])
        X['x_p'] = 0
        X['x_v'] = 0
        X['i_d'] = X['x_p']
        X['i_q'] = -X['x_v']
        X['x_pll'] = 0
        X['angle'] = np.angle(v0)

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
    
    def set_pref_grid(self, x, v, index, pref, power_rating):
        X = self.local_view(x)
        self.par["p_ref_grid"][index] = pref * power_rating*1e-3 / self.par["S_n"]          # MVA / MVA

    # endregion


# class VSC(DAEModel):
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)

#         self.bus_idx = np.array(np.zeros(self.n_units), dtype=[(key, int) for key in self.bus_ref_spec().keys()])
#         self.bus_idx_red = np.array(np.zeros(self.n_units), dtype=[(key, int) for key in self.bus_ref_spec().keys()])

#     def bus_ref_spec(self):
#         return {'terminal': self.par['bus']}
#         return {'terminal': self.par['bus']}

#     def add_blocks(self):

#         p = self.par
#         self.pll = PLL1(T_filter=self.par['T_pll'], bus=p['bus'])

#         self.pi_p = PIRegulator(K_p=p['P_K_p'], K_i=p['P_K_i'])
#         self.pi_p.input = lambda x, v: self.P_setp(x, v) - self.P(x, v)

#         self.pi_q = PIRegulator(K_p=p['Q_K_p'], K_i=p['Q_K_i'])
#         self.pi_q.input = lambda x, v: self.Q_setp(x, v) - self.Q(x, v)

#         self.lag_p = TimeConstant(T=p['T_i'])
#         self.lag_p.input = self.pi_p.output
#         self.lag_q = TimeConstant(T=p['T_i'])
#         self.lag_q.input = self.pi_q.output

#         self.I_d = self.lag_p.output
#         self.I_q = self.lag_q.output

#     def I_inj(self, x, v):
#         return (self.I_d(x, v) - 1j*self.I_q(x, v))*np.exp(1j*self.pll.output(x, v))

#     def input_list(self):
#         return ['P_setp', 'Q_setp']

#     def P(self, x, v):
#         v_n = self.sys_par['bus_v_n'][self.bus_idx_red['terminal']]
#         V = abs(v[self.bus_idx_red['terminal']])*v_n
#         return np.sqrt(3)*V*self.I_d(x, v)

#     def Q(self, x, v):
#         v_n = self.sys_par['bus_v_n'][self.bus_idx_red['terminal']]
#         V = abs(v[self.bus_idx_red['terminal']])*v_n
#         return np.sqrt(3)*V*self.I_q(x, v)

#     def load_flow_pq(self):
#         return self.bus_idx['terminal'], -self.par['P_setp'], -self.par['Q_setp']

#     def init_from_load_flow(self, x_0, v_0, S):
#         self._input_values['P_setp'] = self.par['P_setp']
#         self._input_values['Q_setp'] = self.par['Q_setp']

#         v_n = self.sys_par['bus_v_n'][self.bus_idx_red['terminal']]

#         V_0 = v_0[self.bus_idx_red['terminal']]*v_n

#         I_d_0 = self.par['P_setp']/(abs(V_0)*np.sqrt(3))
#         I_q_0 = self.par['Q_setp']/(abs(V_0)*np.sqrt(3))

#         self.pi_p.initialize(
#             x_0, v_0, self.lag_p.initialize(x_0, v_0, I_d_0)
#         )

#         self.pi_q.initialize(
#             x_0, v_0, self.lag_q.initialize(x_0, v_0, I_q_0)
#         )

#     def current_injections(self, x, v):
#         i_n = self.sys_par['s_n'] / (np.sqrt(3) * self.sys_par['bus_v_n'])
#         # self.P(x, v)
#         return self.bus_idx_red['terminal'], self.I_inj(x, v)/i_n[self.bus_idx_red['terminal']]
    

# class VSC_PQ(DAEModel):
#     """
#     Instantiate:
#     'vsc': {
#             'VSC_PQ': [
#                 ['name', 'bus', 'S_n', 'p_ref', 'q_ref',  'k_p', 'k_q', 'T_p', 'T_q', 'k_pll','T_pll', 'T_i', 'i_max'],
#                 ['VSC1', 'B1',    50,     1,       0,       1,      1,    0.1,   0.1,     5,      1,      0.01,    1.2],
#             ],
#         }
#     """
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)

#         self.bus_idx = np.array(np.zeros(self.n_units), dtype=[(key, int) for key in self.bus_ref_spec().keys()])
#         self.bus_idx_red = np.array(np.zeros(self.n_units), dtype=[(key, int) for key in self.bus_ref_spec().keys()])

#     # region Definitions

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

#     def state_derivatives(self, dx, x, v):
#         dX = self.local_view(dx)
#         X = self.local_view(x)
#         par = self.par

#         dp = self.p_ref(x,v) - self.p_e(x, v)
#         dq = -self.q_ref(x,v) + self.q_e(x, v)
#         i_d_ref = dp * par['k_p'] + X['x_p']
#         i_q_ref = dq * par['k_q'] + X['x_q']

#         # Limiters
#         i_ref = (i_d_ref+1j*i_q_ref)
#         i_ref = i_ref*par['i_max']/np.maximum(par['i_max'],abs(i_ref))

#         X['x_p'] = np.maximum(np.minimum(X['x_p'],par['i_max']),-par['i_max'])
#         X['x_q'] = np.maximum(np.minimum(X['x_q'],par['i_max']),-par['i_max'])

#         dX['i_d'][:] = 1 / (par['T_i']) * (i_ref.real - X['i_d'])
#         dX['i_q'][:] = 1 / (par['T_i']) * (i_ref.imag - X['i_q'])
#         dX['x_p'][:] = par['k_p'] / (par['T_p']) * dp
#         dX['x_q'][:] = par['k_q'] / (par['T_q']) * dq
#         dX['x_pll'][:] = par['k_pll'] / (par['T_pll']) * (self.v_q(x,v))
#         dX['angle'][:] = X['x_pll']+par['k_pll']*self.v_q(x,v)
#         #dX['angle'][:] = 0
#         return

#     def init_from_load_flow(self, x_0, v_0, S):
#         X = self.local_view(x_0)

#         self._input_values['p_ref'] = self.par['p_ref']
#         self._input_values['q_ref'] = self.par['q_ref']

#         v0 = v_0[self.bus_idx_red['terminal']]

#         X['i_d'] = self.par['p_ref']/abs(v0)
#         X['i_q'] = self.par['q_ref']/abs(v0)
#         X['x_p'] = X['i_d']
#         X['x_q'] = X['i_q']
#         X['x_pll'] = 0
#         X['angle'] = np.angle(v0)

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
    
#     def p_e(self, x, v):
#         return self.s_e(x,v).real

#     def q_e(self, x, v):
#         return self.s_e(x,v).imag

#     # endregion






# #################### OLD VERSION #############################

# # # region GridSideConverter class
# # class GridSideConverter(DAEModel):
# #     """
# #     VSC PQ model - ref to VSC1 in tops/dyn_models/vsc1.py - cite Sjur Føyen

# #     Model of the DC link and Grid Side Converter.

# #     Its purpose is to control P and Q to the grid and to control the DC link voltage. 

# #     'vsc': {
# #         'VSC_PQ': [
# #             ['name',   'bus', 'S_n',   "p_ref",   "q_ref",   'Cdc',   'k_p',    'k_q',  'T_p',  'T_q',     'k_pll',   'T_pll',    'T_i',    'i_max'],
# #             ['VSC1',   'B1',    50,      1,         0,      0.1,       1,       1,     0.1,     0.1,        5,        1,         0.01,      1.2],
# #         ],
# #     }

# #     Parameters:
# #     w_n: Nominal frequency
# #     Cdc: DC link capacitance

# #     References:
# #     Vdc_ref
# #     p_ref
# #     q_ref
# #     id_ref
# #     iq_ref

# #     States:
# #     i_d
# #     i_q
# #     angle   PLL angle
# #     xpll    integral i pll
# #     vdc     DC link voltage

# #     """

# #     def __init__(self, *args, **kwargs):
# #         super().__init__(*args, **kwargs)

# #         self.bus_idx = np.array(np.zeros(self.n_units), dtype=[(key, int) for key in self.bus_ref_spec().keys()])
# #         self.bus_idx_red = np.array(np.zeros(self.n_units), dtype=[(key, int) for key in self.bus_ref_spec().keys()])

#         # self.p_ref = IPMSM.get_Pe()              # Active power reference will come from the IPMSM
#         # self.q_ref = 0                           # Reactive power reference can/should come from grid side

#         # self.vdc = 1.0
#         # self.vdc_ref = self.vdc

#         # self.vdc_controller = PIController_LAH(kp=0.3, ti=0.003)

#     # region Definitions

#     def init_from_load_flow(self, x_0, v_0, S):
#         X = self.local_view(x_0)

#         self._input_values['p_ref_grid'] = self.par['p_ref_grid']       # unescassary?
#         self._input_values['q_ref_grid'] = self.par['q_ref_grid']

#         vg = v_0[self.bus_idx_red['terminal']]

#         X['angle'] = np.angle(vg)
#         X['x_pll'] = 0
#         X['i_d'] = self.par['p_ref_grid']/abs(vg)
#         X['i_q'] = self.par['q_ref_grid']/abs(vg)
#         X["vdc"] = 1.0
#         X['x_pref_adj'] = 0.0

#         # X["pref_adj"] = 0.0
#         # X['vdc_ref'] = self.par['vdc_ref']
#         # X['p_ref'] = self.par['p_ref']
#         # X['q_ref'] = self.par['q_ref']
        
#     def state_derivatives(self, dx, x, v):
#         dX = self.local_view(dx)
#         X = self.local_view(x)
#         par = self.par

#         # self.p_ref = IPMSM.get_Pe() - self.vdc_controller.compute((self.vdc_ref - X["vdc"]))*X["vdc"]
#         # self.q_ref = 0.0     # Should be updated from grid side command

#         i_d_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (self.par["p_ref_grid"] * self.v_d(x,v) - self.par["q_ref_grid"] * self.v_q(x,v))
#         i_q_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (self.par["p_ref_grid"] * self.v_q(x,v) + self.par["q_ref_grid"] * self.v_d(x,v))

#         # i_d_ref = max(0, i_d_ref)
#         # i_q_ref = max(0, i_q_ref)

#         # i_d_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (X["p_ref"] * self.v_d(x,v) - X["q_ref"] * self.v_q(x,v))
#         # i_q_ref = 1/(self.v_d(x,v)**2 + self.v_q(x,v)**2) * (X["p_ref"] * self.v_q(x,v) + X["q_ref"] * self.v_d(x,v))

#         # Limiters
#         i_ref = (i_d_ref+1j*i_q_ref)
#         i_ref = i_ref*par['i_max']/np.maximum(par['i_max'],abs(i_ref))

#         dX['i_d'][:] = 1 / (par['T_i']) * (i_ref.real - X['i_d'])
#         dX['i_q'][:] = 1 / (par['T_i']) * (i_ref.imag - X['i_q'])

#         dX['x_pll'][:] = par['k_pll'] / (par['T_pll']) * (self.v_q(x,v))
#         dX['angle'][:] = X['x_pll']+par['k_pll']*self.v_q(x,v)
        
#         dX["vdc"][:] = (par["p_ref_gen"] - self.p_e(x,v)) / (par["Cdc"] * X["vdc"])
#         dX["x_pref_adj"][:] = (X["vdc"] - par["vdc_ref"]) * (par["K_p_dc"] / par["T_i_dc"])


#         return    

#     def load_flow_pq(self):
#         return self.bus_idx['terminal'], -self.par['p_ref_grid']*self.par['S_n'], -self.par['q_ref_grid']*self.par['S_n']

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
#         x_pref_adj: DC link voltage control integral
#         pref_adj: DC link voltage control output
#         vdc: DC link voltage
#         """
#         # return ['i_d', 'i_q', 'x_p', 'x_q', 'x_pll', 'angle', 'vdc', 'vdc_ref', 'x_pref_adj', 'pref_adj' ,"p_ref", "q_ref"]
#         return ['i_d', 'i_q', 'x_pll', 'angle', 'vdc', 'x_pref_adj']

#     def input_list(self):
#         """
#         All values in pu.
#         p_ref: outer loop active power setpoint
#         q_ref: outer loop reactive power setpoint
#         """
#         return ['p_ref_grid', 'q_ref_grid']

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

#     def vdc(self, x, v):
#         return self.local_view(x)["vdc"]
    
#     # def set_pref_grid(self, x, v, index):
#     #     X = self.local_view(x)
#     #     self.par["p_ref_grid"][index] = self.par["p_ref_gen"][index] + self.p_ref_adj(x,v)[index]

#     def set_pref_grid(self, x, v, index, pref):
#         X = self.local_view(x)
#         self.par["p_ref_gen"][index] = pref
#         self.par['p_ref_grid'][index] = self.par["p_ref_gen"][index] + self.p_ref_adj(x,v)[index]

#     def set_pref_gen(self, ref, index):
#         self.par['p_ref_gen'][index] = ref

#     def p_ref_adj(self, x, v):
#         X = self.local_view(x)
#         return (X["x_pref_adj"] + (X["vdc"] - self.par["vdc_ref"]) * self.par["K_p_dc"]) #* X["vdc"]

#     def set_qref_grid(self, Qref, index):
#         self.par['q_ref_grid'][index] = Qref
    


    
#     # endregion