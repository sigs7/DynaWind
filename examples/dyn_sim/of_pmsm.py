# Coupled PMSM and OpenFAST lets goooo

import sys
from tops.cosim_models.windturbine import WindTurbine
from tops.cosim_models.results import Results

if __name__ == '__main__':

    ### SIMULATION SETTINGS ###
    t = 0
    step_size_mech = 5e-3
    step_size_elec = 5e-6
    t_end = 13

    sim_name = "Send_hjelp"

    ## Dict to store the results
    results = Results()

    # Create Wind Turbine instance
    WT1 = WindTurbine(name='WT1', index=0, gsc_control="PQ")

    # # Store the init results
    # results.store_time(t)
    # results.store_fmu_results(WT1)
    # results.store_pmsm_results(WT1)

    while t < t_end:
        sys.stdout.write("\r%d%%" % (t/(t_end)*100))

        # Step the Wind Turbine
        # WT1.step_windturbine_multirate(t, dt)
        WT1.fast.step_fmu(WT1.pmsm, t, step_size_mech)

        elec_steps_per_mech_step = int(step_size_mech / step_size_elec)

        for _ in range(elec_steps_per_mech_step):
            WT1.pmsm.step_pmsm(WT1.fast, time=t, step_size_elec=step_size_elec)
            t += step_size_elec
        
            # Store the results
            results.store_time(t)
            results.store_fmu_results(WT1)
            results.store_pmsm_results(WT1)


        # Update time
        t += step_size_mech

    ## Simulation is finished ##

    # Terminate the FMU
    WT1.fast.terminate_fmu()

    print("Simulation finished, plotting:")

    results.plot_pmsm_overview(sim_name=sim_name, WT=WT1)
    results.plot_pmsm_overview_interactive(sim_name=sim_name, WT=WT1)
    results.plot_fmu_overview(sim_name=sim_name, WT=WT1)
    results.plot_controllers(sim_name=sim_name, WT=WT1)





    
    