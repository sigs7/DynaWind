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

from tops.dyn_models.pmsm import *
from tops.dyn_models.vsc1 import VSC_PV

# Define parameters
pmsm_params = {
    "rs": 0.03,
    "x_d": 0.4,
    "x_q": 0.4,
    "Psi_m": 0.9,
    "Tm": 4,
    "w_n" : 2*np.pi*50  # nominal rad
}

MSC_params = {"T_conv" : 1e-2,
              "vq_0" : 0.5,
              "vd_0" : 0.0
}

prime_mover_params = {"T_pm" : 1e-1,
                      "speed_0" : 0.5,
                      "torque_0" : 0.7
}

# Create PMSM instance
pmsm = PMSM(pmsm_params, MSC_params = MSC_params, prime_mover_params = prime_mover_params)         # Initiate the PMSM


if __name__ == '__main__':

    # Load model
    import tops.ps_models.user_ps_models.k2a_vsc as model_data
    model = model_data.load()

    model["vsc"] = {"GridSideConverter": [
        ['name',   'bus', 'S_n',   "p_ref",   "q_ref",   'Cdc',   'k_p',    'k_q',  'T_p',  'T_q',     'k_pll',   'T_pll',    'T_i',    'i_max'],
        ['WT1',    'B1',    50,      0.0,         0,      0.1,       5,       1,     0.1,     0.1,        5,        1,         0.01,      1.2],
    ]}

    # Power system model
    ps = dps.PowerSystemModel(model=model)
    ps.init_dyn_sim()
    print(max(abs(ps.state_derivatives(0, ps.x_0, ps.v_0))))

    x_0 = ps.x_0.copy()
    t_end = 20
    dt = 1e-3

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
        if t > 5 and event_flag1:
            pmsm.set_prime_mover_reference(speed_ref=0.8, torque_ref=0.3, ramp_time=3, dt=dt, current_time=t)
            event_flag1 = False

        if t > 15 and event_flag2:
            pmsm.set_prime_mover_reference(speed_ref=0.4, torque_ref=0.7, ramp_time=3, dt=dt, current_time=t)
            event_flag2 = False        

        pmsm.update_states(t=t, dt=sol.dt)
        Pref = pmsm.get_Pe()

        ps.vsc['GridSideConverter'].set_pref(Pref)

        # Simulate next step
        result = sol.step()
        x = sol.y
        v = sol.v
        t = sol.t

        dx = ps.ode_fun(0, ps.x_0)

        # Store result
        res['t'].append(t)


        ### Generator results ###
        # res['gen_speed'].append(ps.gen['GEN'].speed(x, v).copy())
        res["G2 speed"].append(ps.gen['GEN'].speed(x, v)[1].copy())
        res["Generator speeds"].append(ps.gen['GEN'].speed(x, v).copy())
        res["Generator current injections"].append(abs(ps.gen['GEN'].i(x, v).copy()))
        res["G2_i_inj"].append(abs(ps.gen['GEN'].i(x, v)[1].copy()))


        ### PMSM results ###
        res['pmsm_speed'].append(pmsm.speed)
        res['pmsm_torque'].append(pmsm.get_Te())
        res['speed_ref'].append(pmsm.primemover.speed_ref)
        res['speed_error'].append(-(pmsm.primemover.speed_ref - pmsm.primemover.speed))
        res['pmsm_i_d'].append(pmsm.i_d)
        res['pmsm_i_q'].append(pmsm.i_q)

        res['primemover_speed'].append(pmsm.primemover.speed)
        res['primemover_torque'].append(pmsm.primemover.torque)

        res['electrical_power'].append(pmsm.get_Pe())
        res['mechanical_power'].append(pmsm.get_Pm())


        ### VSC results ###
        res['VSC_p_e'].append(ps.vsc["GridSideConverter"].p_e(x, v).copy())
        # res['VSC_p_ref'].append(ps.vsc["GridSideConverter"].p_ref(x, v).copy())

        res['VSC_q'].append(ps.vsc["GridSideConverter"].q_e(x, v).copy())
        # res['VSC_q_ref'].append(ps.vsc["GridSideConverter"].q_ref(x, v).copy())

        res["WT i_inj"].append(abs(ps.vsc["GridSideConverter"].i_inj(x, v).copy()))


    print('\n Simulation completed in {:.2f} seconds.'.format(time.time() - t_0))

    # Plot the results
    fig, axs = plt.subplots(5, 1, figsize=(15, 20), sharex=True)

    # Mechanical results
    axs[0].plot(res['t'], res['pmsm_speed'], label='pmsm Speed')
    axs[0].plot(res['t'], res['primemover_speed'], label='Prime Mover Speed')
    axs[0].plot(res['t'], res['speed_ref'], label='Speed Reference', linestyle='--')
    axs[0].set_ylabel('Speed (pu)')
    axs[0].legend()
    axs[0].grid(True)

    axs[1].plot(res['t'], res['pmsm_torque'], label='pmsm Torque')
    axs[1].plot(res['t'], res['primemover_torque'], label='Prime Mover Torque')
    axs[1].plot(res['t'], res['speed_error'], label='Speed Error', linestyle='--')
    axs[1].set_ylabel('Torque (pu)')
    axs[1].legend()
    axs[1].grid(True)

    # Electrical results
    axs[2].plot(res['t'], res['pmsm_i_d'], label='i_d')
    axs[2].plot(res['t'], res['pmsm_i_q'], label='i_q')
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
    # axs[4].plot(res['t'], res['VSC_p_ref'], label='VSC P_ref')
    axs[4].plot(res['t'], res['VSC_q'], label='VSC Q')
    # axs[4].plot(res['t'], res['VSC_q_ref'], label='VSC Q_ref')
    axs[4].plot(res['t'], res["WT i_inj"], label='WT i_inj')
    axs[4].set_ylabel('VSC (pu)')
    axs[4].legend()
    axs[4].grid(True)

    plt.xlabel('Time (s)')
    plt.tight_layout()
    plt.savefig("Figures/TOPS_GridSideConverter_response.png")

    fig, axs = plt.subplots(3, 1, figsize=(15, 10), sharex=True)

    # Generator Speeds
    axs[0].plot(res['t'], res["Generator speeds"], label='Generator speeds')
    axs[0].set_ylabel('Speed (pu)')
    axs[0].legend()
    axs[0].grid(True)
    axs[0].set_title('Generator Speeds vs Time')

    # Generator Current Injections
    axs[1].plot(res['t'], res["Generator current injections"], label='Generator current injections')
    axs[1].set_ylabel('Current (pu)')
    axs[1].legend()
    axs[1].grid(True)
    axs[1].set_title('Generator Current Injections vs Time')

    # WT Injection
    axs[2].plot(res['t'], res["WT i_inj"], label='WT i_inj')
    axs[2].set_xlabel('Time (s)')
    axs[2].set_ylabel('Injected Current (pu)')
    axs[2].legend()
    axs[2].grid(True)
    axs[2].set_title('WT Injection vs Time')

    plt.tight_layout()
    plt.savefig("Figures/TOPS_WT_system_response.png")
    # plt.show()