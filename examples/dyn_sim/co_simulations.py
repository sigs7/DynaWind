## Now lets try to combine all three simulation tools

# WindTurbine imports
import sys
from tops.cosim_models.windturbine import WindTurbine
from tops.cosim_models.results import Results

# TOPS imports
from collections import defaultdict
import matplotlib.pyplot as plt
import time
import tops.dynamic as dps
import tops.solvers as dps_sol
import importlib
importlib.reload(dps)

if __name__ == '__main__':

    import tops.ps_models.k2a_2WT as model_data
    model = model_data.load()

    model["vsc"] = {"GridSideConverter": [ 
    ['name',   'bus',    'S_n',      "p_ref_grid",      "q_ref_grid",       "p_ref_gen",     'Cdc',      'k_p',      'k_q',    'T_p',     'T_q',     'k_pll',   'T_pll',    'T_i',    'K_p_dc',   'T_i_dc',    'i_max',   'vdc_ref'],
    ['WT1',    'B1',      30,         0.0,               0,                  0.0,             0.2,          5,          1,        0.1,        0.1,        5,        1,         0.01,     1.0,          0.2,         1.5,       1.0],
    # ['WT2',    'B3',      50,         0.0,               0,                  0.0,             0.1,          5,          1,        0.1,        0.1,        5,        1,         0.01,     1.0,          0.2,         1.2,       1.0],
    ]}

    # Power system model
    ps = dps.PowerSystemModel(model=model)
    ps.init_dyn_sim()
    print(max(abs(ps.state_derivatives(0, ps.x_0, ps.v_0))))

    x_0 = ps.x_0.copy()

    ### SIMULATION SETTINGS ###
    simulation_name = "CO5"
    t = 0
    dt = 5e-3
    t_end = 45

    # Solver
    sol = dps_sol.ModifiedEulerDAE(ps.state_derivatives, ps.solve_algebraic, 0, x_0, t_end, max_step=dt)

    res = defaultdict(list)

    unique_timesteps = set()
    t_0 = time.time()

    ## Dict to store the results
    results = Results()

    # Create Wind Turbine instance
    # print("Chechpoint 1")
    WT1 = WindTurbine(name='WT1', index = 0)
    # print("Checkpoint 1.1")
    # WT2 = WindTurbine(name='WT2', index = 1)
    # print("Checkpoint 1.2")

    while t < t_end:
        sys.stdout.write("\r%d%%" % (t/(t_end)*100))

        # Step the Wind Turbine
        # print("Chechpoint 2")
        WT1.step_windturbine(t, dt)

        # Step TOPS
        result = sol.step()
        x = sol.y
        v = sol.v
        t = sol.t

        dx = ps.ode_fun(0, ps.x_0)

        Pref_gen1 = WT1.pmsm.get_p_e()
        ps.vsc["GridSideConverter"].set_pref_gen(ref=Pref_gen1, index=0)
        ps.vsc["GridSideConverter"].set_pref_grid(index=0, x=x, v=v)


        # Pref_gen2 = -WT2.pmsm.get_p_e()
        # ps.vsc["GridSideConverter"].set_pref_gen(ref=Pref_gen2, index=1)
        # ps.vsc["GridSideConverter"].set_pref_grid(index=1, x=x, v=v)
        
        # Store the results
        results.store_time(t)
        results.store_fmu_results(WT1)
        results.store_pmsm_results(WT1)
        results.store_vsc_results(WT1, ps, x, v)
        results.store_generator_results(ps, x, v)
        # results.store_vsc_results(WT2, ps, x, v, index=1)


        # Update time
        t += dt

    ## Simulation is finished ##

    # Terminate the FMU
    WT1.fast.terminate_fmu()
    # WT2.fast.terminate_fmu()


    results.plot_pmsm_overview(sim_name=simulation_name, WT = WT1)
    results.plot_fmu_overview(sim_name=simulation_name, WT = WT1)
    results.plot_vsc_overview(sim_name=simulation_name, WT = WT1)
    results.plot_tops_overview(sim_name=simulation_name, WT = WT1)


