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


vsc_params = {
    "name": "VSC1",
    "bus": "B1",
    "S_n": 50,
    "p_ref": 1,
    "q_ref": 0,
    "k_p": 1,
    "k_q": 1,
    "T_p": 0.1,
    "T_q": 0.1,
    "k_pll": 5,
    "T_pll": 1,
    "T_i": 0.01,
    "i_max": 1.2
}


# Create IPMSM instance
ipmsm = IPMSM(ipmsm_params, MSC_params = MSC_params, prime_mover_params = prime_mover_params)         # Initiate the IPMSM

# Create GridSideConverter instance
WindTurbine = GridSideConverter(DAEModel, ipmsm = ipmsm, vsc_params = vsc_params)               # Initiate the grid side converter / Wind turbine to be connected to TOPS


if __name__ == '__main__':

    # Load model
    import tops.ps_models.user_ps_models.k2a_vsc as model_data
    model = model_data.load()

    # Power system model
    ps = dps.PowerSystemModel(model=model)
    ps.init_dyn_sim()
    print(max(abs(ps.state_derivatives(0, ps.x_0, ps.v_0))))

    t_end = 15
    x_0 = ps.x_0.copy()

    dt = 5e-3

    # Solver
    sol = dps_sol.ModifiedEulerDAE(ps.state_derivatives, ps.solve_algebraic, 0, x_0, t_end, dt)

    # Initialize simulation
    t = 0
    res = defaultdict(list)
    t_0 = time.time()


    # Run simulation
    while t < t_end:
        sys.stdout.write("\r%d%%" % (t/(t_end)*100))



        if t > 1:
            ps.vsc1['VSC_PQ'].set_input('p_ref', 0.8)
        if t > 5:
            ps.vsc1['VSC_PQ'].set_input('q_ref', 0.1)

        # Simulate next step
        result = sol.step()
        x = sol.y
        v = sol.v
        t = sol.t

        dx = ps.ode_fun(0, ps.x_0)

        # Store result
        res['t'].append(t)

        res['gen_speed'].append(ps.gen['GEN'].speed(x, v).copy())

        res['VSC_P'].append(ps.vsc1["VSC_PQ"].p_e(x, v).copy())
        # res['VSC_P_setp'].append(ps.vsc["VSC_PQ"].Pref.copy())

        res['VSC_Q'].append(ps.vsc1["VSC_PQ"].q_e(x, v).copy())
        # res['VSC_Q_setp'].append(ps.vsc["VSC_PQ"].Q_setp(x, v).copy())

        res["VSC_abs(I)"].append(abs(ps.vsc1["VSC_PQ"].i_inj(x, v).copy()))

        # res["VSC_Id"].append(ps.vsc["VSC_PQ"].I_d(x, v).copy())
        # res["VSC_Iq"].append(ps.vsc["VSC_PQ"].I_q(x, v).copy())

    print('\n Simulation completed in {:.2f} seconds.'.format(time.time() - t_0))

    plt.figure()
    plt.plot(res['t'], res['gen_speed'], label = "Generator speed")
    plt.xlabel('Time [s]')
    plt.ylabel('Gen. speed')
    plt.legend(loc='upper right')
    # plt.show()

    fig, ax = plt.subplots(2)
    ax[0].plot(res['t'], res['VSC_P'], label = "VSC_P")
    # # ax[0].plot(res['t'], res['VSC_P_setp'], label = "VSC_P_setp")
    ax[0].legend(loc='upper right')
    
    ax[1].plot(res['t'], res['VSC_Q'], label = "VSC_Q")
    # # ax[1].plot(res['t'], res['VSC_Q_setp'], label = "VSC_Q_setp")
    ax[1].legend(loc="upper right")
    # # plt.show()

    # plt.figure()
    # plt.plot(res['t'], res['VSC_abs(I)'], label = "VSC_abs(I)")
    # plt.xlabel('Time [s]')
    # plt.ylabel('VSC current injection [A]')
    # plt.legend(loc='upper right')

    # fig, ax = plt.subplots(2)
    # ax[0].plot(res['t'], res['VSC_Id'], label = "VSC_Id")
    # ax[0].plot(res['t'], res['VSC_Iq'], label = "VSC_Iq")
    # ax[0].legend(loc='upper right')
    
    # ax[1].plot(res['t'], res['VSC_Q'], label = "VSC_Q")
    # # ax[1].plot(res['t'], res['VSC_Q_setp'], label = "VSC_Q_setp")
    # ax[1].legend(loc="upper right")
    # # plt.show()

    plt.show()