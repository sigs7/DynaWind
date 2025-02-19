# Testing the machine side modelling

from tops.dyn_models.pmsm_1 import *

pmsm_params = {
        "rs": 0.03,
        "x_d": 0.4,
        "x_q": 0.4,
        "Psi_m": 0.9,
        "Tm": 4,
        "w_n" : 2*np.pi*50  # nominal rad
    }

MSC_params = {"T_conv" : 2/(300),
            "vq_0" : 0.5,
            "vd_0" : 0.0
}

prime_mover_params = {"T_pm" : 1e-1,
                    "speed_0" : 0.5,
                    "torque_0" : 0.5
}

# Create PMSM instance
pmsm = PMSM(pmsm_params = pmsm_params, MSC_params = MSC_params, prime_mover_params = prime_mover_params)         # Initiate the PMSM


res = {}
# Lists to store the results
time_values = []


event_flag1 = True
event_flag2 = True
event_flag3 = True


# Simulation parameters
t = 0
dt = 5e-3

simulation_time = 20  # Total simulation time

unique_timesteps = set()

# Run the simulation
while t < simulation_time:
    # Print the percentage of simulation completed
    if int(t / dt) % int(simulation_time / 100 / dt) == 0:
        print(f"\rSimulation {int(t / simulation_time * 100)}% complete", end='')

    pmsm.update_states(t , dt)

    # Update the states
    if t > 5 and event_flag1:
        pmsm.set_prime_mover_reference(speed_ref=0.8, torque_ref=0.3, ramp_time=3, dt=dt, current_time=t)
        event_flag1 = False

    if t > 15 and event_flag2:
        pmsm.set_prime_mover_reference(speed_ref=0.4, torque_ref=0.7, ramp_time=3, dt=dt, current_time=t)
        event_flag2 = False

    # if t > 15 and event_flag3:
    #     pmsm.set_prime_mover_reference(speed_ref=0.7, torque_ref=0.7, ramp_time=10, dt=dt, current_time=t)
    #     event_flag3 = False


    # if dt > 1e-4 and int(t / dt) % int(simulation_time / 10 / dt) == 0:
    #     print("")
    #     print(f"dt changed to {dt}")
    
    # Add the current timestep to the set
    unique_timesteps.add(dt)

    t += dt

    # Get mechanical states and controller signals
    error_speed = pmsm.primemover.speed_ref - pmsm.primemover.speed

    # Store the results
    time_values.append(t)

    # Store the mechanical states
    res.setdefault('motor_torque_values', []).append(pmsm.get_Te())
    res.setdefault('motor_speed_values', []).append(pmsm.speed)
    res.setdefault('primemover_torque_values', []).append(pmsm.primemover.torque)
    res.setdefault('primemover_speed_values', []).append(pmsm.primemover.speed)

    # Store the electrical states
    res.setdefault('i_d_values', []).append(pmsm.i_d)
    res.setdefault('i_q_values', []).append(pmsm.i_q)
    res.setdefault('v_d_values', []).append(pmsm.converter.vd)
    res.setdefault('v_q_values', []).append(pmsm.converter.vq)
    res.setdefault('motor_Pe_values', []).append(pmsm.get_Pe())

    # Store the control signals
    res.setdefault('i_d_ref_values', []).append(pmsm.i_d_ref)
    res.setdefault('i_q_ref_values', []).append(pmsm.i_q_ref)
    res.setdefault('v_d_ref_values', []).append(pmsm.converter.vd_ref)
    res.setdefault('v_q_ref_values', []).append(pmsm.converter.vq_ref)
    res.setdefault('torque_ref_values', []).append(pmsm.torque_ref)
    res.setdefault('speed_ref_values', []).append(pmsm.primemover.speed_ref)

    # Store the error signals
    res.setdefault('error_speed_values', []).append(error_speed)
    res.setdefault('error_iq_values', []).append(pmsm.i_q_ref - pmsm.i_q)
    res.setdefault('error_id_values', []).append(pmsm.i_d_ref - pmsm.i_d)

print("")
print("Unique timesteps: ", len(unique_timesteps))
print("")
print("Plotting in progress...")
# Plot the results of voltage and current
fig, axs = plt.subplots(4, 1, figsize=(12, 12), sharex=True)

# Speed plot
axs[0].plot(time_values, res['speed_ref_values'], label='Speed Reference', linestyle='--')
axs[0].plot(time_values, res['motor_speed_values'], label='Motor Speed')
axs[0].plot(time_values, res['primemover_speed_values'], label='Turbine Speed')
axs[0].set_ylabel('Speed (pu)')
axs[0].legend()
axs[0].grid(True)

# Torque plots
axs[1].plot(time_values, res['motor_torque_values'], label='Motor Torque')
axs[1].plot(time_values, res['primemover_torque_values'], label='Turbine Torque')
axs[1].plot(time_values, res['torque_ref_values'], label='Torque Reference', linestyle='--')
axs[1].set_ylabel('Torque (pu)')
axs[1].legend()
axs[1].grid(True)

# Current plots
axs[2].plot(time_values, res['i_d_values'], label='i_d')
axs[2].plot(time_values, res['i_q_values'], label='i_q')
axs[2].plot(time_values, res['motor_Pe_values'], label='Motor Pe')
axs[2].set_ylabel('Currents (pu)')
axs[2].legend()
axs[2].grid(True)

# Voltage plots
axs[3].plot(time_values, res['v_d_values'], label='v_d')
axs[3].plot(time_values, res['v_q_values'], label='v_q')
axs[3].set_xlabel('Time (s)')
axs[3].set_ylabel('Voltage (pu)')
axs[3].legend()
axs[3].grid(True)

plt.tight_layout()
plt.savefig("Figures/MST_motorperformance.png")
# plt.show()

# Plotting the control signals
# Plot the reference currents, voltages, and speed error
fig, axs = plt.subplots(3, 1, figsize=(12, 12), sharex=True)

# Current reference plots
axs[0].plot(time_values, res['i_d_ref_values'], label='i_d_ref')
axs[0].plot(time_values, res['i_d_values'], label='i_d')
axs[0].plot(time_values, res['i_q_ref_values'], label='i_q_ref')
axs[0].plot(time_values, res['i_q_values'], label='i_q')
axs[2].set_xlabel('Time (s)')
axs[0].set_ylabel('Current References (pu)')
axs[0].legend()
axs[0].grid(True)

# Voltage reference plots
axs[1].plot(time_values, res['v_d_ref_values'], label='v_d_ref')
axs[1].plot(time_values, res['v_q_ref_values'], label='v_q_ref')
axs[1].plot(time_values, res['v_d_values'], label='v_d')
axs[1].plot(time_values, res['v_q_values'], label='v_q')
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Voltage References (pu)')
axs[1].legend()
axs[1].grid(True)

# Speed error plot
axs[2].plot(time_values, res['error_speed_values'], label='Speed Error')
axs[2].plot(time_values, res['torque_ref_values'], label='Torque Ref')
# axs[2].plot(time_values, res['i_q_ref_values'], label='i_q_ref')
axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel('Speed Error (pu)')
axs[2].legend()
axs[2].grid(True)

plt.tight_layout()
plt.savefig("Figures/MST_controlsignals.png")
# plt.show()

print("Simulation completed successfully")

