# Lets try to test the entire WT system

import sys
from collections import defaultdict
import matplotlib.pyplot as plt
import time
import tops.dynamic as dps
import tops.solvers as dps_sol
import importlib
importlib.reload(dps)
import importlib

from tops.dyn_models.IPMSM_drives import *
from tops.dyn_models.vsc1 import VSC_PV

# Define parameters
ipmsm_params = {
    "rs": 0.03,
    "x_d": 0.4,
    "x_q": 0.4,
    "Psi_m": 0.9,
    "Tm": 4,
    "w_n" : 2*np.pi*50  # nominal rad
}

MSC_params = {"T_conv" : 1e-4,
              "vq_0" : 0.5,
              "vd_0" : 0.0
}

prime_mover_params = {"T_pm" : 1e-1,
                      "speed_0" : 0.5,
                      "torque_0" : 0.7
}

# Create IPMSM instance
ipmsm = IPMSM(ipmsm_params, MSC_params = MSC_params, prime_mover_params = prime_mover_params)         # Initiate the IPMSM

# Create GridSideConverter instance
# WindTurbine = GridSideConverter(DAEModel, ipmsm = ipmsm, vsc_params = vsc_params)               # Initiate the grid side converter / Wind turbine to be connected to TOPS


if __name__ == '__main__':

    # Load model
    import tops.ps_models.user_ps_models.k2a_vsc as model_data
    model = model_data.load()

    model["vsc"] = {"GridSideConverter": [
        ['name',   'bus', 'S_n',   "p_ref",   "q_ref",   'Cdc',   'k_p',    'k_q',  'T_p',  'T_q',     'k_pll',   'T_pll',    'T_i',    'i_max'],
        ['WT1',    'B1',    50,      1,         0,      0.1,       1,       1,     0.1,     0.1,        5,        1,         0.01,      1.2],
    ]}

    # Power system model
    ps = dps.PowerSystemModel(model=model)


    ps.init_dyn_sim()
    print(max(abs(ps.state_derivatives(0, ps.x_0, ps.v_0))))

    t_end = 3

    x_0 = ps.x_0.copy()

    dt = 7e-4

    # Solver
    sol = dps_sol.ModifiedEulerDAE(ps.state_derivatives, ps.solve_algebraic, 0, x_0, t_end, max_step=dt)

    # Initialize simulation
    t = 0
    res = defaultdict(list)
    t_0 = time.time()

    event_flag1 = True
    event_flag2 = True

    # Run simulation
    while t < t_end:
        sys.stdout.write("\r%d%%" % (t/(t_end)*100))

        # Update the states
        if t > 1 and event_flag1:
            ipmsm.set_prime_mover_reference(speed_ref=0.5, torque_ref=0.4, ramp_time=1, dt=dt, current_time=t)
            event_flag1 = False

        if t > 20 and event_flag2:
            ipmsm.set_prime_mover_reference(speed_ref=0.7, torque_ref=0.8, ramp_time=10, dt=dt, current_time=t)
            event_flag2 = False        

        ipmsm.update_states(t=t, dt=sol.dt)
        Pref = ipmsm.get_Pe()

        ps.vsc['GridSideConverter'].set_input('p_ref', Pref)


        # Simulate next step
        result = sol.step()
        x = sol.y
        v = sol.v
        t = sol.t

        dx = ps.ode_fun(0, ps.x_0)

        # Store result
        res['t'].append(t)

        res['gen_speed'].append(ps.gen['GEN'].speed(x, v).copy())

        res['ipmsm_speed'].append(ipmsm.speed)
        res['ipmsm_torque'].append(ipmsm.get_Te())
        res['speed_ref'].append(ipmsm.primemover.speed_ref)
        res['speed_error'].append(-(ipmsm.primemover.speed_ref - ipmsm.primemover.speed))
        res['ipmsm_i_d'].append(ipmsm.i_d)
        res['ipmsm_i_q'].append(ipmsm.i_q)

        res['primemover_speed'].append(ipmsm.primemover.speed)
        res['primemover_torque'].append(ipmsm.primemover.torque)

        res['electrical_power'].append(ipmsm.get_Pe())
        res['mechanical_power'].append(ipmsm.get_Pm())

        res['VSC_p_e'].append(ps.vsc["GridSideConverter"].p_e(x, v).copy())
        res['VSC_p_ref'].append(ps.vsc["GridSideConverter"].p_ref(x, v).copy())


        res['VSC_q'].append(ps.vsc["GridSideConverter"].q_e(x, v).copy())
        res['VSC_q_ref'].append(ps.vsc["GridSideConverter"].q_ref(x, v).copy())

        res["WT i_inj"].append(abs(ps.vsc["GridSideConverter"].i_inj(x, v).copy()))


    print('\n Simulation completed in {:.2f} seconds.'.format(time.time() - t_0))

    # Plot the results
    fig, axs = plt.subplots(5, 1, figsize=(15, 20), sharex=True)

    # Mechanical results
    axs[0].plot(res['t'], res['ipmsm_speed'], label='IPMSM Speed')
    axs[0].plot(res['t'], res['primemover_speed'], label='Prime Mover Speed')
    axs[0].plot(res['t'], res['speed_ref'], label='Speed Reference', linestyle='--')
    axs[0].set_ylabel('Speed (pu)')
    axs[0].legend()
    axs[0].grid(True)

    axs[1].plot(res['t'], res['ipmsm_torque'], label='IPMSM Torque')
    axs[1].plot(res['t'], res['primemover_torque'], label='Prime Mover Torque')
    axs[1].plot(res['t'], res['speed_error'], label='Speed Error', linestyle='--')
    axs[1].set_ylabel('Torque (pu)')
    axs[1].legend()
    axs[1].grid(True)

    # Electrical results
    axs[2].plot(res['t'], res['ipmsm_i_d'], label='i_d')
    axs[2].plot(res['t'], res['ipmsm_i_q'], label='i_q')
    axs[2].set_ylabel('Current (pu)')
    axs[2].legend()
    axs[2].grid(True)

    axs[3].plot(res['t'], res['electrical_power'], label='Electrical Power')
    axs[3].plot(res['t'], res['mechanical_power'], label='Mechanical Power')
    axs[3].set_ylabel('Power (pu)')
    axs[3].legend()
    axs[3].grid(True)

    # VSC results
    axs[4].plot(res['t'], res['VSC_p_e'], label='VSC P_e')
    axs[4].plot(res['t'], res['VSC_p_ref'], label='VSC P_ref')
    axs[4].plot(res['t'], res['VSC_q'], label='VSC Q')
    axs[4].plot(res['t'], res['VSC_q_ref'], label='VSC Q_ref')
    axs[4].plot(res['t'], res["WT i_inj"], label='WT i_inj')
    axs[4].set_ylabel('VSC (pu)')
    axs[4].legend()
    axs[4].grid(True)

    plt.xlabel('Time (s)')
    plt.tight_layout()
    plt.savefig("WT_system_test.png")