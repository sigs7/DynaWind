
# hard copy of gptfmu2.py to adjust into OpenFAST FMU
# Custom input - adapted from https://github.com/CATIA-Systems/FMPy/blob/main/fmpy/examples/custom_input.py

""" This example demonstrates how to use the FMU.get*() and FMU.set*() functions
 to set custom input and control the simulation """

from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave
from fmpy.util import plot_result
import os
import numpy as np
import shutil
import matplotlib.pyplot as plt

def simulate_custom_input(show_plot=True):

    # define the model name and simulation parameters
    # fmu_filename = 'C:/Users/larsi/OpenFAST/OpenFASTFMU2PF_Export/fast.fmu'
    fmu_filename = 'fast.fmu'
    
    start_time = 0.0
    stop_time = 5
    step_size = 0.01

    # read the model description
    model_description = read_model_description(fmu_filename, validate=False)    #, validate=False

    # collect the value references
    vrs = {}
    for variable in model_description.modelVariables:
        vrs[variable.name] = variable.valueReference


    # Print the value references to verify them
    for name, vr in vrs.items():
        print(f"Variable: {name}, Value Reference: {vr}")

    # extract the FMU
    unzipdir = extract(fmu_filename, os.path.abspath('openfast_debug_fmu'))

    wd_file_path = 'openfast_fmu/resources/wd.txt'
    new_directory = 'C:/Users/larsi/Master/TOPS_LAH/TOPS_LAH'

    # Write the new directory to the wd.txt file
    with open(wd_file_path, 'w') as f:
        f.write(new_directory)

    print(f"Extracted to {unzipdir}")

    fmu = FMU2Slave(guid=model_description.guid,
                    unzipDirectory=unzipdir,
                    modelIdentifier=model_description.coSimulation.modelIdentifier,
                    instanceName='instance1')

    fmu.instantiate()
    fmu.setReal([vrs['testNr']], [1002])
    fmu.setupExperiment(startTime=start_time)
    fmu.enterInitializationMode()
    fmu.exitInitializationMode()

    time = start_time
    rows = []  # list to record the results

    print("Starting simulation")
    # Simulation loop
    while time < stop_time:

        # Perform one step
        # print(f" \n \nPerforming step at time {time} \n \n")
        fmu.doStep(currentCommunicationPoint=time, communicationStepSize=step_size)
        # print(f"Step performed at time {time}")

        # Get the outputs
        GenTrq = fmu.getReal([vrs['GenTq']])[0]
        GenSpeed = fmu.getReal([vrs['GenSpeed']])[0]
        RotSpeed = fmu.getReal([vrs['RotSpeed']])[0]
        BlPitchCom = fmu.getReal([vrs['HSShftTq']])[0]
        
        # Append the results
        rows.append((time, GenTrq, GenSpeed, RotSpeed, BlPitchCom))

        # advance the time
        time += step_size

    print(f"Simulation finished at time {time}")

    # clean up
    shutil.rmtree(unzipdir, ignore_errors=True)

    # Convert results to a structured numpy array
    result = np.array(rows, dtype=np.dtype([
        ('time', np.float64),
        ('GenTrq', np.float64),
        ('GenSpeed', np.float64),
        ('RotSpeed', np.float64),
        ('BlPitchCom', np.float64),
    ]))

    # Save the results
    np.savetxt('results.csv', result, delimiter=',', header='time, GenTrq, GenSpeed, RotSpeed, BlPitchCom', comments='')


    if show_plot:
        times = result['time']
        outputs1 = result['GenTrq']
        outputs2 = result['GenSpeed']
        outputs3 = result['RotSpeed']
        outputs4 = result['BlPitchCom']

        fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
    
        axs[0].plot(times, outputs1, label='GenTrq', color='blue')
        axs[0].set_ylabel('GenTrq')
        axs[0].legend()
        
        axs[1].plot(times, outputs2, label='GenSpeed', color='green')
        axs[1].set_ylabel('GenSpeed')
        axs[1].legend()
        
        axs[2].plot(times, outputs3, label='RotSpeed', color='red')
        axs[2].set_ylabel('RotSpeed')
        axs[2].legend()
        
        axs[3].plot(times, outputs4, label='BlPitchCom', color='purple')
        axs[3].set_ylabel('BlPitchCom')
        axs[3].legend()
        
        
        plt.tight_layout()
        plt.show()

    fmu.terminate()
    fmu.freeInstance()

if __name__ == '__main__':
    simulate_custom_input(show_plot=True)

### C:/Users/larsi/OpenFAST/OpenFASTFMU2PF_Export
### C:/Users/larsi/Master/TOPS_LAH/TOPS_LAH

