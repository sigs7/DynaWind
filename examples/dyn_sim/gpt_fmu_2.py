# GPT FMU aaah nr 2

# Custom input - adapted from https://github.com/CATIA-Systems/FMPy/blob/main/fmpy/examples/custom_input.py

""" This example demonstrates how to use the FMU.get*() and FMU.set*() functions
 to set custom input and control the simulation """

from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave
from fmpy.util import plot_result
import numpy as np
import shutil
import os
import matplotlib.pyplot as plt

def simulate_custom_input(show_plot=True):

    # define the model name and simulation parameters
    fmu_filename = 'C:/Users/larsi/OpenFAST/OpenFASTFMU2PF_Export/fast.fmu'
    start_time = 0.0
    stop_time = 10.0
    step_size = 0.1

    # read the model description
    model_description = read_model_description(fmu_filename)

    # collect the value references
    vrs = {}
    for variable in model_description.modelVariables:
        vrs[variable.name] = variable.valueReference

    # get the value references for the variables we want to get/set
    # vr_mode = vrs['GenTq']  # Replace 'Mode' with the actual variable name from your FMU
    vr_gen_pwr = vrs['GenPwr']  # kW
    vr_gen_speed = vrs['GenSpeed']  # rpm

    # extract the FMU
    unzipdir = extract(fmu_filename)

    fmu = FMU2Slave(guid=model_description.guid,
                    unzipDirectory=unzipdir,
                    modelIdentifier=model_description.coSimulation.modelIdentifier,
                    instanceName='instance1')

    # initialize
    fmu.instantiate()
    fmu.setupExperiment(startTime=start_time)
    fmu.enterInitializationMode()
    fmu.exitInitializationMode()

    time = start_time
    rows = []  # list to record the results

    # simulation loop
    while time < stop_time:

        # Example: Set inputs from the input data (replace with actual logic)
        # Assuming input_data contains time-series data for each input variable
        # if 'input1.txt' in input_data:
        #     input_value = np.interp(time, input_data['input1.txt'][:, 0], input_data['input1.txt'][:, 1])
        #     fmu.setReal([vr_mode], [input_value])

        # perform one step
        fmu.doStep(currentCommunicationPoint=time, communicationStepSize=step_size)

        # advance the time
        time += step_size

        # get the values for 'Mode'
        gen_pwr = fmu.getReal([vr_gen_pwr])[0]
        gen_speed = fmu.getReal([vr_gen_speed])[0]

        # append the results
        rows.append((time, gen_pwr, gen_speed))

    fmu.terminate()
    fmu.freeInstance()

    # clean up
    shutil.rmtree(unzipdir, ignore_errors=True)

    # convert the results to a structured NumPy array
    result = np.array(rows, dtype=np.dtype([('time', np.float64), ('mode', np.float64)]))

    # plot the results
    if show_plot:
        plt.figure()
        plt.plot(result['time'], result['mode'], label='Mode')
        plt.xlabel('Time [s]')
        plt.ylabel('Mode [-]')
        plt.title('Mode vs Time')
        plt.legend()
        plt.tight_layout()
        plt.show()

    return time

if __name__ == '__main__':
    simulate_custom_input(show_plot=True)