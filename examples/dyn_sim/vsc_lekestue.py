import sys
from collections import defaultdict
import matplotlib.pyplot as plt
import time
import tops.dynamic as dps
import tops.solvers as dps_sol

if __name__ == '__main__':

    # Load model
    import tops.ps_models.user_ps_models.k2a_vsc as model_data
    model = model_data.load()

    # Power system model
    ps = dps.PowerSystemModel(model=model)
    ps.init_dyn_sim()
    print(max(abs(ps.state_derivatives(0, ps.x_0, ps.v_0))))

    t_end = 10
    x_0 = ps.x_0.copy()

    # Solver
    sol = dps_sol.ModifiedEulerDAE(ps.state_derivatives, ps.solve_algebraic, 0, x_0, t_end, max_step=5e-3)

    # Initialize simulation
    t = 0
    res = defaultdict(list)
    t_0 = time.time()

    sc_bus_idx = ps.gen['GEN'].bus_idx_red['terminal'][0]

    # Run simulation
    while t < t_end:
        sys.stdout.write("\r%d%%" % (t/(t_end)*100))

        # Short circuit
        if t >= 1 and t <= 1.05:
            ps.y_bus_red_mod[4,4] = 1e6
        else:
            ps.y_bus_red_mod[4,4] = 0

        if t > 1:
            ps.vsc['VSC'].set_input('P_setp', 500)
        if t > 5:
            ps.vsc['VSC'].set_input('Q_setp', 200)

        # Simulate next step
        result = sol.step()
        x = sol.y
        v = sol.v
        t = sol.t

        dx = ps.ode_fun(0, ps.x_0)

        # Store result
        res['t'].append(t)

        res['gen_speed'].append(ps.gen['GEN'].speed(x, v).copy())

        res['VSC_P'].append(ps.vsc['VSC'].P(x, v).copy())
        res['VSC_P_setp'].append(ps.vsc['VSC'].P_setp(x, v).copy())

        res['VSC_Q'].append(ps.vsc['VSC'].Q(x, v).copy())
        res['VSC_Q_setp'].append(ps.vsc['VSC'].Q_setp(x, v).copy())

        res["VSC_I"].append(abs(ps.vsc["VSC"].I_inj(x, v).copy()))

    print('Simulation completed in {:.2f} seconds.'.format(time.time() - t_0))

    plt.figure()
    plt.plot(res['t'], res['gen_speed'], label = "Generator speed")
    plt.xlabel('Time [s]')
    plt.ylabel('Gen. speed')
    plt.legend(loc='upper right')
    # plt.show()

    fig, ax = plt.subplots(2)
    ax[0].plot(res['t'], res['VSC_P'], label = "VSC_P")
    ax[0].plot(res['t'], res['VSC_P_setp'], label = "VSC_P_setp")
    ax[0].legend(loc='upper right')
    
    ax[1].plot(res['t'], res['VSC_Q'], label = "VSC_Q")
    ax[1].plot(res['t'], res['VSC_Q_setp'], label = "VSC_Q_setp")
    ax[1].legend(loc="upper right")
    # plt.show()

    plt.figure()
    plt.plot(res['t'], res['VSC_I'], label = "VSC_I")
    plt.xlabel('Time [s]')
    plt.ylabel('VSC current injection [A]')
    plt.legend(loc='upper right')



    plt.show()