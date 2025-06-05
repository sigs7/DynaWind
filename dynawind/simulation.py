## Simulation script for DynaWind 

# WindTurbine imports
import sys
from dynawind.dynawind_models.windturbine import WindTurbine
from dynawind.dynawind_models.results import Results

# TOPS imports
import time
import tops.dynamic as dps
import tops.solvers as dps_sol
import importlib
importlib.reload(dps)

if __name__ == '__main__':

    import tops.ps_models.k2a_highwind as model_data
    model = model_data.load()

    # Power system model
    ps = dps.PowerSystemModel(model=model)

    # Create Wind Turbine instance
    WT1 = WindTurbine(name='WT1', index = 0, gsc_control="PV")

    # Initiate the power system model
    ps.init_dyn_sim()
    print(max(abs(ps.state_derivatives(0, ps.x_0, ps.v_0))))        # Checks for unstable initiation

    x_0 = ps.x_0.copy()

    ### SIMULATION SETTINGS ###
    simulation_name = "Paper_results_dynawindsolo"
    t = 0
    t_end = 0.1
    step_size_mech = 0.01
    step_size_elec = 5e-6

    # Short circuit settings
    sc_time = 30.0
    sc_duration = 0.100
    sc = True

    # Solver
    sol = dps_sol.ModifiedEulerDAE(ps.state_derivatives, ps.solve_algebraic, 0, x_0, t_end, max_step=step_size_mech)
    t_0 = time.time()

    ## Dict to store the results
    results = Results()
    
    sc_bus_idx = ps.vsc['GridSideConverter_PV'].bus_idx_red['terminal'][0]
    event_flag1 = True

    while t < t_end:
        sys.stdout.write("\r%d%%" % (t/(t_end)*100))

        # Short circuit
        if t >= sc_time and t <= (sc_time + sc_duration) and sc == True:
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

        # Calculate the derivatives of the power system
        dx = ps.ode_fun(0, ps.x_0)

        # Store the results
        results.store_time(t)
        results.store_time_elec(t)
        results.store_fmu_results(WT1)
        results.store_pmsm_results(WT1)
        results.store_dclink_results(WT1, ps, x, v)

        # Store the results of the GSC
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

    # Saves the result class to a file, see plotting.py for loading and plotting of the results
    results.save_to_file(simulation_name)
