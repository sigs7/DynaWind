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

    import tops.ps_models.k2a_highwind as model_data
    model = model_data.load()

    # Create Wind Turbine instance


    # Power system model
    ps = dps.PowerSystemModel(model=model)

    WT1 = WindTurbine(name='WT1', index = 0, gsc_control="PV")

    ps.init_dyn_sim()
    print(max(abs(ps.state_derivatives(0, ps.x_0, ps.v_0))))

    x_0 = ps.x_0.copy()

    ### SIMULATION SETTINGS ###
    # simulation_name = "HighWind-Multirate-SC"
    
    simulation_name = "Revisited_testing_PV"
    t = 0
    step_size_mech = 5e-3
    step_size_elec = 5e-6
    t_end = 7

    sc_time = 4

    # Solver
    sol = dps_sol.ModifiedEulerDAE(ps.state_derivatives, ps.solve_algebraic, 0, x_0, t_end, max_step=step_size_mech)

    t_0 = time.time()

    ## Dict to store the results
    results = Results()
    
    sc_bus_idx = ps.vsc['GridSideConverter_PV'].bus_idx_red['terminal'][0]

    while t < t_end:
        sys.stdout.write("\r%d%%" % (t/(t_end)*100))

        # Short circuit
        if t >= sc_time and t <= (sc_time + 0.150):
            ps.y_bus_red_mod[sc_bus_idx,sc_bus_idx] = 1e5
        else:
            ps.y_bus_red_mod[sc_bus_idx,sc_bus_idx] = 0

        # Step TOPS
        result = sol.step()
        x = sol.y
        v = sol.v
        t = sol.t

        # Step the Wind Turbine
        WT1.step_windturbine_multirate(ps, x, v, t, step_size_mech, step_size_elec, results)
        # WT1.step_windturbine(ps, t, step_size_mech, x, v)
        dx = ps.ode_fun(0, ps.x_0)

        # Store the results
        results.store_time(t)
        results.store_fmu_results(WT1)
        results.store_time_elec(t)
        results.store_pmsm_results(WT1)
        results.store_dclink_results(WT1, ps, x, v)


        if WT1.gsc_control == "PQ":
            results.store_gsc_results_PQ(WT1, ps, x, v)
        elif WT1.gsc_control == "PV":
            results.store_gsc_results_PV(WT1, ps, x, v)
        results.store_generator_results(ps, x, v)


        # Update time
        t += step_size_mech

    ## Simulation is finished ##

    # Terminate the FMU
    WT1.fast.terminate_fmu()
    # WT2.fast.terminate_fmu()

    results.plot_fmu_overview(sim_name=simulation_name, WT = WT1)
    results.plot_pmsm_overview(sim_name=simulation_name, WT = WT1)
    results.plot_dclink_overview(sim_name=simulation_name, WT = WT1)
    results.plot_tops_overview(sim_name=simulation_name, WT = WT1)
    results.plot_pmsm_overview_interactive(sim_name=simulation_name, WT = WT1)
    results.plot_controllers(sim_name=simulation_name, WT = WT1)
    results.plot_vdc_controller_terms(sim_name=simulation_name, WT = WT1)

