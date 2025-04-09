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

    import tops.ps_models.simulink_validation as model_data
    model = model_data.load()

    # Power system model
    ps = dps.PowerSystemModel(model=model)
    ps.init_dyn_sim()
    print(max(abs(ps.state_derivatives(0, ps.x_0, ps.v_0))))
    # Print initial values

    x_0 = ps.x_0.copy()

    ### SIMULATION SETTINGS ###
    simulation_name = "Short circuit test"
    t = 0
    dt = 5e-3
    t_end = 30

    # Solver
    sol = dps_sol.ModifiedEulerDAE(ps.state_derivatives, ps.solve_algebraic, 0, x_0, t_end, max_step=dt)

    res = defaultdict(list)

    t_0 = time.time()

    ## Dict to store the results
    results = Results()

    # Create Wind Turbine instance
    WT1 = WindTurbine(name='WT1', index = 0)

    sc_bus_idx = ps.vsc['GridSideConverter'].bus_idx_red['terminal'][0]


    while t < t_end:
        sys.stdout.write("\r%d%%" % (t/(t_end)*100))

        # Short circuit
        if t >= 15 and t <= 15.30:
            ps.y_bus_red_mod[sc_bus_idx,sc_bus_idx] = 1e6
        else:
            ps.y_bus_red_mod[sc_bus_idx,sc_bus_idx] = 0

        # Step TOPS
        result = sol.step()
        x = sol.y
        v = sol.v
        t = sol.t

        # Step Wind Turbine
        WT1.step_windturbine(ps, t, dt, x, v)
        dx = ps.ode_fun(0, ps.x_0)
    
        # Store the results
        results.store_time(t)
        results.store_fmu_results(WT1)
        results.store_pmsm_results(WT1)
        results.store_msc_results(WT1)
        results.store_dclink_results(WT1, ps, x, v)
        results.store_gsc_results(WT1, ps, x, v)
        results.store_generator_results(ps, x, v)
        

        # Update time
        t += dt


    ## Simulation is finished #


    # Terminate the FMU
    WT1.fast.terminate_fmu()
    # WT2.fast.terminate_fmu()

    results.plot_fmu_overview(sim_name=simulation_name, WT = WT1)
    results.plot_pmsm_overview(sim_name=simulation_name, WT = WT1)
    results.plot_msc_overview(sim_name=simulation_name, WT = WT1)
    results.plot_dclink_overview(sim_name=simulation_name, WT = WT1)
    results.plot_tops_overview(sim_name=simulation_name, WT = WT1)
    results.plot_pmsm_overview_interactive(sim_name=simulation_name, WT = WT1)
    results.plot_controllers(sim_name=simulation_name, WT = WT1)



