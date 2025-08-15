## Simulation script for DynaWind 

# WindTurbine imports
import sys
sys.path.append(r"C:\git\DynaWind-1")
from dynawind.dynawind_models.windturbine import WindTurbine
from dynawind.dynawind_models.results import Results

# TOPS imports
import time
import tops.dynamic as dps
import tops.solvers as dps_sol
import importlib
importlib.reload(dps)
import matplotlib.pyplot as plt
from collections import defaultdict

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
    simulation_name = "dynawind_example"
    t = 0
    t_end = 60
    step_size_mech = 0.01
    step_size_elec = 5e-6

    # Short circuit settings
    sc_time = 15
    sc_duration = 0.100
    sc = False

    # Solver
    sol = dps_sol.ModifiedEulerDAE(ps.state_derivatives, ps.solve_algebraic, 0, x_0, t_end, max_step=step_size_mech)
    t_0 = time.time()

    ## Dict to store the results
    results = Results()

    res = defaultdict(list)
    
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


        res['t'].append(t)
        res['gen_speed'].append(ps.gen['GEN'].speed(x, v).copy())
        res['v'].append(v.copy())
        res['gen_p'].append(ps.gen['GEN'].p_e(x, v).copy())
        res['gen_q'].append(ps.gen['GEN'].q_e(x, v).copy())

        # Update time
        t += step_size_mech

    ## Simulation is finished ##

    # Terminate the FMU
    WT1.fast.terminate_fmu()

    
    plt.figure()
    plt.plot(res['t'], [abs(v_i) for v_i in res['v']])
    plt.xlabel('Time [s]')
    plt.ylabel('Voltage [p.u.]')
    plt.title(f"Voltages at buses k2a_highwind")
    #plt.title(f"Voltages at buses during short circuit at {model['generators']['GEN'][sc_bus_idx+1][0]} at t={sc_time}s to t={sc_time+sc_duration}s")
    plt.legend([buses[0] for buses in model['buses'][1:]])
    plt.ticklabel_format(useOffset=False)
    plt.grid()
    plt.show()

    plt.figure()
    plt.plot(res['t'], res['gen_p'])
    plt.xlabel('Time [s]')
    plt.ylabel('Active power [p.u.]')
    #plt.title(f"Active power of generators during short circuit at {model['generators']['GEN'][sc_bus_idx+1][0]} at t={sc_time}s to t={sc_time+sc_duration}s")
    plt.legend([gens[0] for gens in model['generators']['GEN'][1:]])
    plt.ticklabel_format(useOffset=False)   
    plt.grid()
    plt.show()

    plt.figure()
    plt.plot(res['t'], res['gen_q'])
    plt.xlabel('Time [s]')
    plt.ylabel('Reactive power [p.u.]')
    #plt.title(f"Reactive power of generators during short circuit at {model['generators']['GEN'][sc_bus_idx+1][0]} at t={sc_time}s to t={sc_time+sc_duration}s")
    plt.legend([gens[0] for gens in model['generators']['GEN'][1:]])
    plt.ticklabel_format(useOffset=False)
    plt.grid()
    plt.show()

    
    plt.figure()
    plt.plot(res['t'], res['gen_speed'])
    plt.xlabel('Time [s]')
    plt.ylabel('Gen. speed (relative)')
    #plt.title(f"Generator speed during short circuit at {model['generators']['GEN'][sc_bus_idx+1][0]} at t={sc_time}s to t={sc_time+sc_duration}s")
    plt.legend([f"Gen {i+1}" for i in range(len(res['gen_speed'][0]))])
    plt.grid()
    plt.show()

    # Saves the result class to a file, see plotting.py for loading and plotting of the results
    results.save_to_file(simulation_name)
