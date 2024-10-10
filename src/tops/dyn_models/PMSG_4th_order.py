"""
Machine Side Converter, controls the tourque/power from the generator, 

Inputs:
Vd_ref
Vq_ref

Outputs/States:
Vq
Vd


PMSG - Permanent Magnet Synchronous Generator - fourth-order model

States:
Eq
Ed
Iq
Id
Tm
Pe

Derivetives:
M


"""

class PIController_LAH:
    def __init__(self, Kp, Ki, setpoint=0):
        self.Kp = Kp
        self.Ki = Ki
        self.setpoint = setpoint
        self.integral = 0

    def update(self, measurement, dt):
        error = self.setpoint - measurement
        self.integral += error * dt
        output = self.Kp * error + self.Ki * self.integral
        return output


import numpy as np
from tops.dyn_models.utils import DAEModel
import tops.utility_functions as dps_uf

class PMSGFourthOrder:
    """
Instantiate:
'WT-T4': {
        'PMSGFourthOrder': [
            ['name',   'S_n',  'V_n',   'P',    'H',    'D',    'X_d',    'X_q',    'X_d_t',    'X_q_t',    'X_d_st',    'X_q_st',    'T_d0_t',    'T_q0_t',    'T_d0_st',    'T_q0_st'],
            ['PMSG',    10,     4.0,    10,      3.0,      0,      0.91,     0.66,       1,        0.01,        1.2,         1.2,          0.1,         0.1,        0.01,         0.01],
        ],
    }
"""

    def __init__(self, params):
        self.params = params
        self.state = {
            'angle': 0.0,
            'speed': 0.0,
            'e_q_t': 0.0,
            'e_d_t': 0.0
        }
        self.input_values = {
            'P_m': 0.0,
            'E_f': 0.0
        }

        self.speed_controller = PIController_LAH(Kp=1.0, Ki=0.1)
        self.iq_controller = PIController_LAH(Kp=1.0, Ki=0.1)

    def initialize(self, X_0, v_g_dq, I_d, I_q, s_pu):
        p = self.params

        angle = np.angle(v_g_dq)
        speed = np.zeros_like(v_g_dq)

        v_d = v_g_dq.real
        v_q = v_g_dq.imag

        e_q_t = v_q + p['X_d_t'] * I_d
        e_d_t = v_d - p['X_q_t'] * I_q
        e_t = e_q_t * np.cos(angle) + e_d_t * np.sin(angle)

        e_q_st = v_q + p['X_d_st'] * I_d
        e_d_st = v_d - p['X_q_st'] * I_q
        e_st = e_q_st * np.cos(angle) + e_d_st * np.sin(angle)

        e_q = e_q_t + I_d * (p['X_d'] - p['X_d_t'])
        e = e_q * np.exp(1j * angle)
        e_q_0 = e_q

        PF_n = 1        #p['PF_n'] if 'PF_n' in p.dtype.names else 1
        self.input_values['P_m'] = s_pu.real / PF_n
        self.input_values['E_f'] = e_q_0

        X_0['speed'][:] = speed
        X_0['angle'][:] = angle
        X_0['e_d_t'][:] = e_d_t
        X_0['e_q_t'][:] = e_q_t

    def state_derivatives(self, dx, x, v):
        p = self.params
        dX = self.local_view(dx)
        X = self.state

        T_m = self.P_m(x, v)/(1 + X['speed'])
        P_e = self.p_e(x, v)

        PF_n = p['PF_n'] if 'PF_n' in p.dtype.names else 1
        H = p['H']/PF_n

        dX['speed'][:] = 1 / (2 * H) * (T_m - P_e/PF_n - p['D'] * X['speed'])
        dX['angle'][:] = X['speed'] * 2 * np.pi * self.sys_par['f_n']
        dX['e_q_t'][:] = 1 / (p['T_d0_t']) * (self.E_f(x, v)  - X['e_q_t'] - self.i_d(x, v) * (p['X_d'] - p['X_d_t']))
        dX['e_d_t'][:] = 1 / (p['T_q0_t']) * (-X['e_d_t'] - self.i_q(x, v) * (p['X_q'] - p['X_q_t']))


    def output(self, x, v):
        X = self.state
        p = self.params

        # Terminal voltage calculation
        V_d = X["e_d_t"] - p["R"]*v["I_d"] - p["X_q_t"]*v["I_q"]
        V_q = X["e_q_t"] + p["X_d_t"]*v["I_d"] - p["R"]*v["I_q"]

        # Electrical power calculation
        P_e = X["e_q_t"]*v["I_q"] + X["e_d_t"]*v["I_d"] + (p["X_d_t"] - p["X_q_t"])**v["I_d"]*v["I_q"]

        # Torque calculation
        T_e = P_e / p['omega_base']

        # Generator voltage (magnitude)
        V_g = np.sqrt(V_d**2 + V_q**2)

        # Generator current (magnitude)
        I_g = np.sqrt(v['I_d']**2 + v['I_q']**2)

        return {
            'V_g': V_g,
            'I_g': I_g,
            'I_d': v['I_d'],
            'I_q': v['I_q'],
            'V_d': V_d,
            'V_q': V_q,
            'P_e': P_e,
            'T_e': T_e,
            'angle': X['angle'],
            'speed': X['speed']
        }
    
    def control(self, desired_speed, desired_iq, dt):
        # Update speed controller
        self.speed_controller.setpoint = desired_speed
        iq_ref = self.speed_controller.update(self.state['speed'], dt)

        # Update Iq controller
        self.iq_controller.setpoint = iq_ref
        vq_ref = self.iq_controller.update(self.state['psi_q'], dt)

        return {
            'V_q_ref': vq_ref,
            'I_q_ref': iq_ref
        }




# Define realistic parameters for the PMSG model
params = {
    'S_n': 10,  # Rated apparent power (MVA)
    'V_n': 5,   # Rated voltage (kV)
    'P': 4,     # Number of poles
    'H': 3.5,   # Inertia constant (s)
    'D': 0.01,  # Damping factor
    'X_d': 1.8, # Direct axis synchronous reactance (pu)
    'X_q': 1.7, # Quadrature axis synchronous reactance (pu)
    'X_d_t': 0.3, # Direct axis transient reactance (pu)
    'X_q_t': 0.55, # Quadrature axis transient reactance (pu)
    'X_d_st': 0.2, # Direct axis sub-transient reactance (pu)
    'X_q_st': 0.25, # Quadrature axis sub-transient reactance (pu)
    'T_d0_t': 6.0, # Direct axis transient open circuit time constant (s)
    'T_q0_t': 0.4, # Quadrature axis transient open circuit time constant (s)
    'T_d0_st': 0.03, # Direct axis sub-transient open circuit time constant (s)
    'T_q0_st': 0.05, # Quadrature axis sub-transient open circuit time constant (s)
    'R_s': 0.01, # Stator resistance (pu)
    'omega_base': 1.0, # Base angular speed (pu)
    'J': 0.1, # Moment of inertia (kg*m^2)
    'PF_n': 0.9 # Power factor
}

# Initial conditions
X_0 = {
    'angle': np.array([0.0]),
    'speed': np.array([1.0]),
    'psi_d': np.array([0.0]),
    'psi_q': np.array([0.0])
}

# Initial voltages and currents
v_g_dq = 1.0 + 0.0j
I_d = 0.0
I_q = 0.0
s_pu = 1.0 + 0.0j

# Create PMSG model instance
pmsg = PMSGFourthOrder(params)
pmsg.initialize(X_0, v_g_dq, I_d, I_q, s_pu)

# Simulation loop
dt = 0.01  # Time step
time = np.arange(0, 10, dt)
outputs_list = []

for t in time:
    control_signals = pmsg.control(desired_speed=1.0, desired_iq=0.5, dt=dt)
    # Update system inputs with control signals
    v = {
        'V_d': 0.0,  # Assuming V_d is controlled externally
        'V_q': control_signals['V_q_ref'],
        'I_d': I_d,
        'I_q': control_signals['I_q_ref']
    }
    # Update state derivatives and outputs
    dx = {key: 0.0 for key in X_0}
    pmsg.state_derivatives(dx, X_0, v)
    outputs = pmsg.output(X_0, v)
    outputs_list.append(outputs)
    # Update state variables
    for key in X_0:
        X_0[key] += dx[key] * dt

# Convert outputs to numpy arrays for easier plotting
outputs_array = {key: np.array([output[key] for output in outputs_list]) for key in outputs_list[0]}

# Plot results
import matplotlib.pyplot as plt

plt.figure(figsize=(12, 8))

plt.subplot(3, 2, 1)
plt.plot(time, outputs_array['speed'])
plt.title('Rotor Speed')
plt.xlabel('Time (s)')
plt.ylabel('Speed (pu)')

plt.subplot(3, 2, 2)
plt.plot(time, outputs_array['P_e'])
plt.title('Electrical Power')
plt.xlabel('Time (s)')
plt.ylabel('Power (pu)')

plt.subplot(3, 2, 3)
plt.plot(time, outputs_array['V_g'])
plt.title('Generator Voltage')
plt.xlabel('Time (s)')
plt.ylabel('Voltage (pu)')

plt.subplot(3, 2, 4)
plt.plot(time, outputs_array['I_g'])
plt.title('Generator Current')
plt.xlabel('Time (s)')
plt.ylabel('Current (pu)')

plt.subplot(3, 2, 5)
plt.plot(time, outputs_array['V_d'], label='V_d')
plt.plot(time, outputs_array['V_q'], label='V_q')
plt.title('Direct and Quadrature Axis Voltages')
plt.xlabel('Time (s)')
plt.ylabel('Voltage (pu)')
plt.legend()

plt.subplot(3, 2, 6)
plt.plot(time, outputs_array['I_d'], label='I_d')
plt.plot(time, outputs_array['I_q'], label='I_q')
plt.title('Direct and Quadrature Axis Currents')
plt.xlabel('Time (s)')
plt.ylabel('Current (pu)')
plt.legend()

plt.tight_layout()
plt.show()