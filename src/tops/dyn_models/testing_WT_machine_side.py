# Testing the machine side modelling

from tops.dyn_models.IPMSM_drives import *

# Define parameters
params = {
    "rs": 0.0173,
    "x_d": 0.543,
    "x_q": 1.086,
    "speed": 1,
    "Psi_m": 0.8,
    "Tm": 5,
    "wn" : 2*np.pi*50  # nominal rad
}

# Initiate the converter and the prime mover of the generator
converter = Converter(vd=0.0, vq=0.0, alpha=0.85)
prime_mover = PrimeMover(torque=1.0, speed=0.8, alpha=0.5)

# Create IPMSM instance
ipmsm = IPMSM(params, converter=converter, prime_mover=prime_mover, i_d0=0.0, i_q0=0.0)

# Set desired currents
# ipmsm.set_current_references(i_d_ref=0.1, i_q_ref=0.9)

# Set initial prime mover reference values
ipmsm.set_prime_mover_reference(torque_ref=0.2, speed_ref=1.0)

# Simulation parameters
dt = 1e-3  # Time step
simulation_time = 20  # Total simulation time

# Lists to store the results
time_values = []
speed_values = []
i_d_values = []
i_q_values = []
v_d_values = []
v_q_values = []
motor_torque = []

# Run the simulation
for t in range(int(simulation_time / dt)):
    # Update the states
    ipmsm.update_states(dt)

    # Get current states
    speed = ipmsm.speed
    i_d = ipmsm.i_d
    i_q = ipmsm.i_q
    v_d = ipmsm.converter.vd
    v_q = ipmsm.converter.vq

    # Store the results
    time_values.append(t * dt)
    speed_values.append(speed)
    i_d_values.append(i_d)
    i_q_values.append(i_q)
    v_d_values.append(v_d)
    v_q_values.append(v_q)
    motor_torque.append(ipmsm.get_Te())

# Plot the results of votlage and current
fig, axs = plt.subplots(4, 1, figsize=(12, 12), sharex=True)

# Speed plot
axs[0].plot(time_values, speed_values, label='Speed')
axs[0].set_ylabel('Speed (rad/s)')
axs[0].legend()
axs[0].grid(True)

# Current plots
axs[1].plot(time_values, i_d_values, label='i_d')
axs[1].plot(time_values, i_q_values, label='i_q')
axs[1].set_ylabel('Current (A)')
axs[1].legend()
axs[1].grid(True)

# Voltage plots
axs[2].plot(time_values, v_d_values, label='v_d')
axs[2].plot(time_values, v_q_values, label='v_q')
axs[2].set_ylabel('Voltage (V)')
axs[2].legend()
axs[2].grid(True)

# Combined plot for i_d, i_q, v_d, v_q
axs[3].plot(time_values, i_d_values, label='i_d')
axs[3].plot(time_values, i_q_values, label='i_q')
axs[3].plot(time_values, v_d_values, label='v_d')
axs[3].plot(time_values, v_q_values, label='v_q')
axs[3].set_xlabel('Time (s)')
axs[3].set_ylabel('Values')
axs[3].legend()
axs[3].grid(True)

plt.tight_layout()
plt.savefig("WT_machine_side_test.png")
# plt.show()


# Gammalt


#     # print(f"Time: {t*dt:.2f}, Speed: {speed:.2f}, i_d: {i_d:.2f}, i_q: {i_q:.2f}")


# # Plot the results
# fig, axs = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

# axs[0].plot(time_values, speed_values, label='Speed')
# axs[0].set_ylabel('Speed (rad/s)')
# axs[0].legend()
# axs[0].grid(True)

# axs[1].plot(time_values, i_d_values, label='i_d')
# axs[1].set_ylabel('i_d (A)')
# axs[1].legend()
# axs[1].grid(True)

# axs[2].plot(time_values, i_q_values, label='i_q')
# axs[2].set_ylabel('i_q (A)')
# axs[2].legend()
# axs[2].grid(True)

# axs[3].plot(time_values, motor_torque, label='Motor Torque')
# axs[3].set_xlabel('Time (s)')
# axs[3].set_ylabel('Torque (Nm)')
# axs[3].legend()
# axs[3].grid(True)

# plt.tight_layout()
# # plt.show()
