# GPT til unsetting med FMU hjelp

# import matplotlib.pyplot as plt
# import numpy as np
# from fmpy import read_model_description, extract
# import os

# fmu_path = 'C:\\Users\\larsi\\OpenFAST\\OpenFASTFMU2PF_Export\\fast.fmu'

# # Read the model description
# model_description = read_model_description(fmu_path)

# # Extract the FMU to a temporary directory
# unzipdir = extract(fmu_path)


# from fmpy import simulate_fmu

# # Run the simulation
# result = simulate_fmu(fmu_path)

# # If you want to plot the results


# time = result['time']
# output_var = result['your_output_variable']  # Replace with the name of the variable you want to plot

# plt.plot(time, output_var)
# plt.xlabel('Time')
# plt.ylabel('Output Variable')
# plt.title('FMU Simulation Result')
# plt.show()


# ########### Second itteration ##############

# import os
# import shutil
# import numpy as np
# import matplotlib.pyplot as plt
# from fmpy import read_model_description, extract
# from fmpy.fmi2 import FMU2Slave
# from fmpy.util import plot_result
# from fmpy.model_description import ValidationError

# def read_input_data(folder):
#     # Example function to read input data from files in the folder
#     # Replace this with actual logic to read your specific input files
#     input_data = {}
#     for filename in os.listdir(folder):
#         if filename.endswith('.txt'):  # Assuming input files are .txt
#             filepath = os.path.join(folder, filename)
#             data = np.loadtxt(filepath)
#             input_data[filename] = data
#     return input_data

# def simulate_custom_input(show_plot=True):

#     # define the model name and simulation parameters
#     fmu_filename = 'C:/Users/larsi/OpenFAST/OpenFASTFMU2PF_Export/fast.fmu'
#     input_folder = 'C:/Users/larsi/OpenFAST/OpenFASTFMU2PF_Export/test1001'  # Path to the input folder
#     start_time = 0.0
#     stop_time = 2.0
#     step_size = 1e-3

#     try:
#         # read the model description
#         model_description = read_model_description(fmu_filename)
#     except ValidationError as e:
#         print(f"Validation error: {e}")
#         # Proceed despite validation errors
#         model_description = read_model_description(fmu_filename, validate=False)

#     # collect the value references
#     vrs = {variable.name: variable.valueReference for variable in model_description.modelVariables}

#     # get the value references for the variables we want to get/set
#     vr_gen_pwr = vrs['GenPwr']  # kW
#     vr_gen_speed = vrs['GenSpeed']  # rpm
#     vr_mode = vrs['Mode']  # -

#     # extract the FMU
#     unzipdir = extract(fmu_filename)

#     fmu = FMU2Slave(guid=model_description.guid,
#                     unzipDirectory=unzipdir,
#                     modelIdentifier=model_description.coSimulation.modelIdentifier,
#                     instanceName='instance1')

#     # initialize
#     fmu.instantiate()
#     fmu.setupExperiment(startTime=start_time)
#     fmu.enterInitializationMode()
#     fmu.exitInitializationMode()

#     # Read input data from the folder
#     input_data = read_input_data(input_folder)

#     time = start_time
#     rows = []  # list to record the results

#     # simulation loop
#     while time < stop_time:

#         # Example: Set inputs from the input data (replace with actual logic)
#         # Assuming input_data contains time-series data for each input variable
#         if 'input1.txt' in input_data:
#             input_value = np.interp(time, input_data['input1.txt'][:, 0], input_data['input1.txt'][:, 1])
#             fmu.setReal([vr_mode], [input_value])

#         # perform one step
#         fmu.doStep(currentCommunicationPoint=time, communicationStepSize=step_size)

#         # advance the time
#         time += step_size

#         # get the values for 'GenPwr', 'GenSpeed', and 'Mode'
#         gen_pwr, gen_speed, mode = fmu.getReal([vr_gen_pwr, vr_gen_speed, vr_mode])

#         # append the results
#         rows.append((time, gen_pwr, gen_speed, mode))

#     fmu.terminate()
#     fmu.freeInstance()

#     # clean up
#     shutil.rmtree(unzipdir, ignore_errors=True)

#     # convert the results to a structured NumPy array
#     result = np.array(rows, dtype=np.dtype([('time', np.float64), ('GenPwr', np.float64), ('GenSpeed', np.float64), ('Mode', np.float64)]))

#     # plot the results
#     if show_plot:
#         plt.figure()
#         plt.subplot(3, 1, 1)
#         plt.plot(result['time'], result['GenPwr'])
#         plt.xlabel('Time [s]')
#         plt.ylabel('GenPwr [kW]')
#         plt.title('Generator Power')

#         plt.subplot(3, 1, 2)
#         plt.plot(result['time'], result['GenSpeed'])
#         plt.xlabel('Time [s]')
#         plt.ylabel('GenSpeed [rpm]')
#         plt.title('Generator Speed')

#         plt.subplot(3, 1, 3)
#         plt.plot(result['time'], result['Mode'])
#         plt.xlabel('Time [s]')
#         plt.ylabel('Mode [-]')
#         plt.title('Mode')

#         plt.tight_layout()
#         plt.show()

#     return time

# if __name__ == '__main__':
#     simulate_custom_input()


############ Third itteration ##############
import matplotlib.pyplot as plt
import numpy as np
from fmpy import read_model_description, extract, simulate_fmu
import os

##     fmu_filename = 'C:/Users/larsi/OpenFAST/OpenFASTFMU2PF_Export/fast.fmu'

fmu_path = 'C:/Users/larsi/OpenFAST/OpenFASTFMU2PF_Export/fast.fmu'

# Read the model description
model_description = read_model_description(fmu_path)

# Extract the FMU to a temporary directory
unzipdir = extract(fmu_path)

# Define the simulation parameters
start_time = 0.0
stop_time = 10.0
step_size = 0.1

# Run the simulation
result = simulate_fmu(fmu_path, start_time=start_time, stop_time=stop_time, step_size=step_size)

# Extract the time and mode variables from the result
time = result['time']
mode = result['mode']  # Replace 'mode' with the actual variable name from your FMU

# Plot the results
plt.figure()
plt.plot(time, mode, label='Mode')
plt.xlabel('Time [s]')
plt.ylabel('Mode [-]')
plt.title('Mode vs Time')
plt.legend()
plt.tight_layout()
plt.show()