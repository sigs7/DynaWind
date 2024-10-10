# Testing the machine side modelling

from tops.dyn_models.IPMSM_drives import *

# Define parameters
params = {
    "rs": 0.0173,
    "x_d": 0.543,
    "x_q": 1.086,
    "speed": 1,
    "Psi_m": 0.8,
    "Tm": 3.7  # Mechanical time constant
}
# Create IPMSM instance
ipmsm = IPMSM(params)

# Set initial converter reference voltages
ipmsm.set_converter_reference_voltages(vd_ref=0.5, vq_ref=0.5)

# Set initial prime mover reference values
ipmsm.set_prime_mover_reference(torque_ref=1.0, speed_ref=1.0)

# Simulation parameters
dt = 0.01  # Time step
simulation_time = 10  # Total simulation time

# Lists to store the results
time_values = []
speed_values = []
i_d_values = []
i_q_values = []
motor_torque = []

# Run the simulation
for t in range(int(simulation_time / dt)):
    ipmsm.update_states(dt)

    # Get current states
    speed = ipmsm.speed
    i_d = ipmsm.i_d
    i_q = ipmsm.i_q

        # Store the results
    time_values.append(t * dt)
    speed_values.append(speed)
    i_d_values.append(i_d)
    i_q_values.append(i_q)
    motor_torque.append(ipmsm.get_Te())

    # print(f"Time: {t*dt:.2f}, Speed: {speed:.2f}, i_d: {i_d:.2f}, i_q: {i_q:.2f}")
# Plot the results
plt.figure(figsize=(12, 8))

plt.subplot(4, 1, 1)
plt.plot(time_values, speed_values, label='Speed')
plt.xlabel('Time (s)')
plt.ylabel('Speed (rad/s)')
plt.legend()
plt.grid(True)

plt.subplot(4, 1, 2)
plt.plot(time_values, i_d_values, label='i_d')
plt.xlabel('Time (s)')
plt.ylabel('i_d (A)')
plt.legend()
plt.grid(True)

plt.subplot(4, 1, 3)
plt.plot(time_values, i_q_values, label='i_q')
plt.xlabel('Time (s)')
plt.ylabel('i_q (A)')
plt.legend()
plt.grid(True)

plt.subplot(4, 1, 4)
plt.plot(time_values, motor_torque, label='motor torque Te')
plt.xlabel('Time (s)')
plt.ylabel('Motor torque')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()