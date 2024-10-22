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
    # fmu_filename = 'C:/Users/larsi/OpenFAST/OpenFASTFMU2PF_Export/fast.fmu'
    
    #debugging fmu
    fmu_filename = 'CoupledClutches.fmu'
    start_time = 0.0
    stop_time = 2
    step_size = 1e-3

    # validate_fmu(fmu_filename)
    # read the model description
    model_description = read_model_description(fmu_filename)

    # collect the value references
    vrs = {}
    for variable in model_description.modelVariables:
        vrs[variable.name] = variable.valueReference


    # # Print the value references to verify them
    # for name, vr in vrs.items():
    #     print(f"Variable: {name}, Value Reference: {vr}")


    # get the value references for the variables we want to get/set
    vr_output1 = vrs['outputs[1]']  # Replace 'Output1' with the actual variable name
    vr_output2 = vrs['outputs[2]']  # Replace 'Output2' with the actual variable name
    vr_output3 = vrs['outputs[3]']  # Replace 'Output3' with the actual variable name
    vr_output4 = vrs['outputs[4]']

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
        # fmu.doStep(currentCommunicationPoint=time, communicationStepSize=step_size)

        try:
            fmu.doStep(currentCommunicationPoint=time, communicationStepSize=step_size)
        except Exception as e:
            print(f"Error during doStep at time {time}: {e}")
            break

        # Get the outputs
        try:
            output1 = fmu.getReal([vr_output1])[0]
            output2 = fmu.getReal([vr_output2])[0]
            output3 = fmu.getReal([vr_output3])[0]
            output4 = fmu.getReal([vr_output4])[0]
        except Exception as e:
            print(f"Error getting output at time {time}: {e}")
            break

        # Append the results
        rows.append((time, output1, output2, output3, output4))       #, output3, output4))


        # advance the time
        time += step_size


    fmu.terminate()
    fmu.freeInstance()

    # clean up
    shutil.rmtree(unzipdir, ignore_errors=True)

    # Convert results to a structured numpy array
    result = np.array(rows, dtype=np.dtype([
        ('time', np.float64),
        ('output1', np.float64),
        ('output2', np.float64),
        ('output3', np.float64),
        ('output4', np.float64)
    ]))

    # Save the results
    np.savetxt('results.csv', result, delimiter=',', header='time,output1,output2,output3,output4', comments='')

    if show_plot:
        times = result['time']
        outputs1 = result['output1']
        outputs2 = result['output2']
        outputs3 = result['output3']
        outputs4 = result['output4']
        plt.plot(times, outputs1, label='Output1')
        plt.plot(times, outputs2, label='Output2')
        plt.plot(times, outputs3, label='Output3')
        plt.plot(times, outputs4, label='Output4')
        plt.xlabel('Time (s)')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    simulate_custom_input(show_plot=True)
