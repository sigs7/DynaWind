
from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave
from fmpy.util import plot_result
import numpy as np
import shutil
import matplotlib.pyplot as plt


def initiate_FAST_FMU(fmu_filename, start_time : float = 0.0) ->  FMU2Slave:

    # define the model name and simulation parameters
    # fmu_filename = 'C:/Users/larsi/OpenFAST/OpenFASTFMU4CTRL_Export/fast.fmu'    

    # start_time = 0.0
    # stop_time = 3.0
    # step_size = 0.025
    
    # traditional stepsize 0.025


    # read the model description
    model_description = read_model_description(fmu_filename, validate=False)    #, validate=False

    # # collect the value references
    vrs = {}
    for variable in model_description.modelVariables:
        vrs[variable.name] = variable.valueReference


    # Print the value references to verify them
    for name, vr in vrs.items():
        print(f"Variable: {name}, Value Reference: {vr}")

    # # extract the FMU
    unzipdir = extract(fmu_filename)

    fmu = FMU2Slave(guid=model_description.guid,
                    unzipDirectory=unzipdir,
                    modelIdentifier=model_description.coSimulation.modelIdentifier,
                    instanceName='instance1')

    fmu.instantiate()
    fmu.setReal([vrs['testNr']], [1002])
    fmu.setupExperiment(startTime=start_time)
    fmu.enterInitializationMode()
    fmu.exitInitializationMode()

    return fmu, vrs