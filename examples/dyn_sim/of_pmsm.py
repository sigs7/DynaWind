# Coupled PMSM and OpenFAST lets goooo

import sys
from tops.cosim_models.windturbine import WindTurbine
from tops.cosim_models.results import Results

if __name__ == '__main__':

    ### SIMULATION SETTINGS ###
    t = 0
    dt = 5e-3
    t_end = 5

    ## Dict to store the results
    results = Results()

    # Create Wind Turbine instance
    WT1 = WindTurbine(name='WT1', index=0)

    while t < t_end:
        sys.stdout.write("\r%d%%" % (t/(t_end)*100))
        # Step the Wind Turbine
        WT1.step_windturbine(t, dt)
        
        # Store the results
        results.store_time(t)
        results.store_fmu_results(WT1)
        results.store_pmsm_results(WT1)

        # Update time
        t += dt

    ## Simulation is finished ##

    # Terminate the FMU
    WT1.fast.terminate_fmu()
    print("Successfully terminated FMU")


    results.plot_pmsm_overview(WT1)
    results.plot_fmu_overview(WT1)





    
    