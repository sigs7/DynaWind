# Testing the machine side modelling

from tops.dyn_models.IPMSM_drives import *

# Define parameters
params = {
    "rs": 0.03,
    "x_d": 0.3,
    "x_q": 0.8,
    "speed": 1,
    "Psi_m": 1.6,
    "Tm": 4,
    "wn" : 2*np.pi*50  # nominal rad
}

# Initiate the converter and the prime mover of the generator
converter = Converter(vd=0.0, vq=0.0, alpha=0.85)               # Initiate the converter
prime_mover = PrimeMover(torque=0.5, speed=0.8, alpha=0.5)      # Initiate the prime mover

# Create IPMSM instance
ipmsm = IPMSM(params, converter=converter, prime_mover=prime_mover, i_d0=0.0, i_q0=0.0)         # Initiate the IPMSM

# Simulation parameters
dt = 1e-4  # Time step 

ipmsm.set_prime_mover_reference(speed_ref=1, torque_ref=0.5)

#Multi timestep?
simulation_time = 30  # Total simulation time

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


# Run the simulation
for t in range(int(simulation_time / dt)):
    # Update the states
    ipmsm.update_states(dt)

    # Get mechanical states
    motor_speed = ipmsm.speed

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
    time_values.append(t * dt)

    # Store the mechanical states
    motor_torque_values.append(ipmsm.get_Te())
    motor_speed_values.append(motor_speed)
    primemover_torque_values.append(ipmsm.primemover.torque)
    primemover_speed_values.append(ipmsm.primemover.speed)

    # Store the electrical states
    i_d_values.append(i_d)
    i_q_values.append(i_q)
    v_d_values.append(v_d)
    v_q_values.append(v_q)

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



# Plot the results of votlage and current
fig, axs = plt.subplots(4, 1, figsize=(12, 12), sharex=True)

# Speed plot
axs[0].plot(time_values, motor_speed_values, label='Motor Speed')
axs[0].plot(time_values, primemover_speed_values, label='turbine speed')
axs[0].set_ylabel('Speed (pu)')
axs[0].legend()
axs[0].grid(True)

# Torque plots
axs[1].plot(time_values, motor_torque_values, label='motor torque')
axs[1].plot(time_values, primemover_torque_values, label='turbine torque')
axs[1].set_ylabel('Torque (pu)')
axs[1].legend()
axs[1].grid(True)

# Current plots
axs[2].plot(time_values, i_d_values, label='i_d')
axs[2].plot(time_values, i_q_values, label='i_q')
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
axs[2].plot(time_values, i_q_ref_values, label='i_q_ref')
axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel('Speed Error (pu)')
axs[2].legend()
axs[2].grid(True)

plt.tight_layout()
plt.savefig("WT_machine_side_voltage_references.png")
# plt.show()