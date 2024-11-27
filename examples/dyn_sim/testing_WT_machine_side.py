# Testing the machine side modelling

from tops.dyn_models.IPMSM_drives import *

# Define parameters
ipmsm_params = {
    "rs": 0.015,
    "x_d": 0.4,
    "x_q": 0.4,
    "Psi_m": 0.9,
    "Tm": 4,
    "w_n" : 2*np.pi*50  # nominal rad
}

MSC_params = {"T_conv" : 1e-2,
              "vq_0" : 0.5,
              "vd_0" : 0.0
}

prime_mover_params = {"T_pm" : 5e-1,
                      "speed_0" : 0.5,
                      "torque_0" : 0.5
}

# Create IPMSM instance
ipmsm = IPMSM(ipmsm_params = ipmsm_params, MSC_params = MSC_params, prime_mover_params = prime_mover_params)         # Initiate the IPMSM

# Lists to store the results
time_values = []

# States mechanical
speed_values = []
motor_torque_values = []
motor_speed_values = []
primemover_torque_values = []
primemover_speed_values = []

#States electrical
i_d_values = []
i_q_values = []
v_d_values = []
v_q_values = []

motor_Pe_values = []

# Control signals
i_d_ref_values = []
i_q_ref_values = []
v_d_ref_values = []
v_q_ref_values = []
torque_ref_values = []
speed_ref_values = []

# Error signals
error_speed_values = []
error_iq_values = []
error_id_values = []

event_flag1 = True
event_flag2 = True
event_flag3 = True



# Simulation parameters
t = 0
dt = 1e-3  # Time step
# tol = 1e-10  # Tolerance
simulation_time = 25  # Total simulation time

unique_timesteps = set()

# Run the simulation
while t < simulation_time:
    # Print the percentage of simulation completed
    if int(t / dt) % int(simulation_time / 100 / dt) == 0:
        print(f"\rSimulation {int(t / simulation_time * 100)}% complete", end='')

    ipmsm.update_states(t , dt)

    # Update the states
    if t > 5 and event_flag1:
        ipmsm.set_prime_mover_reference(speed_ref=0.8, torque_ref=0.3, ramp_time=3, dt=dt, current_time=t)
        event_flag1 = False

    if t > 15 and event_flag2:
        ipmsm.set_prime_mover_reference(speed_ref=0.4, torque_ref=0.7, ramp_time=3, dt=dt, current_time=t)
        event_flag2 = False

    # if t > 15 and event_flag3:
    #     ipmsm.set_prime_mover_reference(speed_ref=0.7, torque_ref=0.7, ramp_time=10, dt=dt, current_time=t)
    #     event_flag3 = False


    # if dt > 1e-4 and int(t / dt) % int(simulation_time / 10 / dt) == 0:
    #     print("")
    #     print(f"dt changed to {dt}")
    
    # Add the current timestep to the set
    unique_timesteps.add(dt)

    t += dt
    # Get mechanical states
    # motor_speed = ipmsm.speed

    # Get electrical states
    i_d = ipmsm.i_d
    i_q = ipmsm.i_q
    v_d = ipmsm.converter.vd
    v_q = ipmsm.converter.vq

    # Get controller signals and outputs
    error_speed = ipmsm.primemover.speed_ref - ipmsm.primemover.speed
    motor_torque_ref = ipmsm.torque_ref
    i_q_ref = ipmsm.i_q_ref
    i_d_ref = ipmsm.i_d_ref

    # Current control signals
    v_d_ref = ipmsm.converter.vd_ref
    v_q_ref = ipmsm.converter.vq_ref

    # Store the results
    time_values.append(t)

    # Store the mechanical states
    motor_torque_values.append(ipmsm.get_Te())
    motor_speed_values.append(ipmsm.speed)
    primemover_torque_values.append(ipmsm.primemover.torque)
    primemover_speed_values.append(ipmsm.primemover.speed)

    # Store the electrical states
    i_d_values.append(i_d)
    i_q_values.append(i_q)
    v_d_values.append(v_d)
    v_q_values.append(v_q)

    motor_Pe_values.append(ipmsm.get_Pe())

    # Store the control signals
    i_d_ref_values.append(ipmsm.i_d_ref)
    i_q_ref_values.append(i_q_ref)
    v_d_ref_values.append(v_d_ref)
    v_q_ref_values.append(v_q_ref)
    torque_ref_values.append(motor_torque_ref)
    speed_ref_values.append(ipmsm.primemover.speed_ref)

    # Store the error signals
    error_speed_values.append(error_speed)
    error_iq_values.append(ipmsm.i_q_ref - i_q)
    error_id_values.append(ipmsm.i_d_ref - i_d)

print("")
print("Unique timesteps: ", len(unique_timesteps))
print("")
print("Plotting in progress...")
# Plot the results of votlage and current
fig, axs = plt.subplots(4, 1, figsize=(12, 12), sharex=True)

# Speed plot
axs[0].plot(time_values, speed_ref_values, label='Speed Reference', linestyle='--')
axs[0].plot(time_values, motor_speed_values, label='Motor Speed')
axs[0].plot(time_values, primemover_speed_values, label='turbine speed')
axs[0].set_ylabel('Speed (pu)')
axs[0].legend()
axs[0].grid(True)

# Torque plots
axs[1].plot(time_values, motor_torque_values, label='motor torque')
axs[1].plot(time_values, primemover_torque_values, label='turbine torque')
axs[1].plot(time_values, torque_ref_values, label='Torque Reference', linestyle='--')
axs[1].set_ylabel('Torque (pu)')
axs[1].legend()
axs[1].grid(True)

# Current plots
axs[2].plot(time_values, i_d_values, label='i_d')
axs[2].plot(time_values, i_q_values, label='i_q')
axs[2].plot(time_values, motor_Pe_values, label='motor Pe')
axs[2].set_ylabel('Currents (pu)')
axs[2].legend()
axs[2].grid(True)

# Voltage plots

axs[3].plot(time_values, v_d_values, label='v_d')
axs[3].plot(time_values, v_q_values, label='v_q')
axs[3].set_xlabel('Time (s)')
axs[3].set_ylabel('Voltage (pu)')
axs[3].legend()
axs[3].grid(True)

plt.tight_layout()
plt.savefig("WT_machine_side_test.png")
# plt.show()

# Plotting the control signals
# Plot the reference currents, voltages, and speed error
fig, axs = plt.subplots(3, 1, figsize=(12, 12), sharex=True)

# Current reference plots
axs[0].plot(time_values, i_d_ref_values, label='i_d_ref')
axs[0].plot(time_values, i_d_values, label='i_d')
axs[0].plot(time_values, i_q_ref_values, label='i_q_ref')
axs[0].plot(time_values, i_q_values, label='i_q')
axs[2].set_xlabel('Time (s)')
axs[0].set_ylabel('Current References (pu)')
axs[0].legend()
axs[0].grid(True)

# Voltage reference plots
axs[1].plot(time_values, v_d_ref_values, label='v_d_ref')
axs[1].plot(time_values, v_q_ref_values, label='v_q_ref')
axs[1].plot(time_values, v_d_values, label='v_d')
axs[1].plot(time_values, v_q_values, label='v_q')
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Voltage References (pu)')
axs[1].legend()
axs[1].grid(True)

# Speed error plot

axs[2].plot(time_values, error_speed_values, label='Speed Error')
axs[2].plot(time_values, torque_ref_values, label='Torque Ref')
# axs[2].plot(time_values, i_q_ref_values, label='i_q_ref')
axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel('Speed Error (pu)')
axs[2].legend()
axs[2].grid(True)

plt.tight_layout()
plt.savefig("WT_machine_side_voltage_references.png")
# plt.show()

print("Simulation completed successfully")