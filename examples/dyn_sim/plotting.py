from tops.cosim_models.results import Results
from tops.cosim_models.windturbine import WindTurbine
import tops.dynamic as dps
import tops.ps_models.k2a_highwind as model_data


# This script is used to plot the results of a wind turbine simulation.
# Should be run after the simulation has been completed and the results have been saved.

# Very easy fix for now, need to initiate the power system model again amongst other things.
# Should be fixed or improved in the future, but works for now.


def main():
    simulation_name = "Paper_results_60_SC"       # Name of the simulation file
    start_time = 28                                  # Start time for plotting
    stop_time = 37.5                               # Stop time for plotting

    # Load the results
    results = Results()
    results.load_from_file(simulation_name)
    model = model_data.load()

    # Power system model
    ps = dps.PowerSystemModel(model=model)

    # Create Wind Turbine instance
    WT1 = WindTurbine(name='WT1', index=0, gsc_control="PV")
    WT1.get_converter_rating(ps)


    # Generate selected plots (comment/uncomment as needed)
    # results.plot_paper_dclink_3x1(sim_name=simulation_name, WT=WT1, start_time=start_time, stop_time=stop_time)
    # results.plot_paper_dclink_2x1(sim_name=simulation_name, WT=WT1, start_time=start_time, stop_time=stop_time)
    results.plot_paper_dclink_2x2(sim_name=simulation_name, WT=WT1, start_time=start_time, stop_time=stop_time)
    # results.plot_paper_fmugen(sim_name=simulation_name, WT=WT1, start_time=start_time, stop_time=stop_time)
    # results.plot_paper_gscgrid(sim_name=simulation_name, WT=WT1, ps=ps, start_time=start_time, stop_time=stop_time)
    # results.plot_paper_hsshfttq_vs_genspdortq(sim_name=simulation_name, WT=WT1, start_time=start_time)
    # results.plot_yawbrtaxp_and_y(sim_name=simulation_name, WT=WT1, start_time=start_time, stop_time=stop_time)
    # results.plot_fmu_overview(sim_name=simulation_name, WT=WT1)

if __name__ == "__main__":
    main()

