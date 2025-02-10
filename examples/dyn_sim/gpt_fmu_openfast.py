# GPT FMU aaah nr 2
# hard copy of gptfmu2.py to adjust into OpenFAST FMU

# Custom input - adapted from https://github.com/CATIA-Systems/FMPy/blob/main/fmpy/examples/custom_input.py

""" This example demonstrates how to use the FMU.get*() and FMU.set*() functions
 to set custom input and control the simulation """

from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave
from fmpy.util import plot_result
import numpy as np
import shutil
import matplotlib.pyplot as plt

def simulate_custom_input(show_plot=True):

    # define the model name and simulation parameters
    fmu_filename = 'C:/Users/larsi/OpenFAST/OpenFASTFMU4CTRL_Export/fast.fmu'
    

    start_time = 0.0
    stop_time = 3.0
    # step_size = 0.025
    step_size = 5e-3
    
    # traditional stepsize 0.025


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
    unzipdir = extract(fmu_filename)

    fmu = FMU2Slave(guid=model_description.guid,
                    unzipDirectory=unzipdir,
                    modelIdentifier=model_description.coSimulation.modelIdentifier,
                    instanceName='instance1')

    # initialize
    # fmu.setReal([0])=1002
    # fmu.setReal([0],[1002])  # Set the testNr value
    # fmu.setReal([vrs['testNr']], [1002])

    fmu.instantiate()
    fmu.setReal([vrs['testNr']], [1002])
    fmu.setupExperiment(startTime=start_time)
    fmu.enterInitializationMode()
    fmu.exitInitializationMode()

    time = start_time
    rows = []  # list to record the results

    # Simulation loop
    while time < stop_time:
        # Set initial values for the variables
        # fmu.setReal([3],[7.55])
        # fmu.setReal([4],[0.0])
        # fmu.setReal([5],[20_000])

        # Perform one step
        try:
            # print(f"Performing doStep at time {time}")
            fmu.doStep(currentCommunicationPoint=time, communicationStepSize=step_size)
            # status = fmu.doStep(currentCommunicationPoint=time, communicationStepSize=step_size)
            # if status != 0:
            #     print(f"doStep returned status {status} at time {time}")
            #     break

        except Exception as e:
            print(f"Error during doStep at time {time}: {e}")
            break

        # Get the outputs
        try:
            GenTrq = fmu.getReal([vrs['GenTrq']])[0]
            GenSpeed = fmu.getReal([vrs['GenSpeed']])[0]
            RotSpeed = fmu.getReal([vrs['RotSpeed']])[0]
            BlPitchCom = fmu.getReal([vrs['BlPitchCom']])[0]

            # print("GenTrq: ", GenTrq)
            # output1 = output2 * np.pi() * output4 / 30
            # print(f"Outputs at time {time}: {output1}, {output2}, {output3}, {output4}")
        except Exception as e:
            print(f"Error getting output at time {time}: {e}")
            break
        
        
        # Append the results
        rows.append((time, GenTrq, GenSpeed, RotSpeed, BlPitchCom))
        # print(f"rows : {rows}")

        # advance the time
        time += step_size

    print(f"Simulation finished at time {time}")
    # Terminate the FMU
    
    # fmu.terminate()

    # fmu.freeInstance()
    # print("FMU freed")

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

    # Print the results
    # print("results = " , result)

    # Save the results
    np.savetxt('results.csv', result, delimiter=',', header='time, GenTrq, GenSpeed, RotSpeed, BlPitchCom', comments='')

    # fmu.terminate()
    # fmu.freeInstance()


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

