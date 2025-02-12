# Lets try combining all three of the previous examples into one program

# TOPS imports
import sys
from collections import defaultdict
import matplotlib.pyplot as plt
import time
import tops.dynamic as dps
import tops.solvers as dps_sol
import importlib
importlib.reload(dps)
import importlib

# PMSM imports
from tops.dyn_models.pmsm import *

# FAST FMU imports
from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave
from fmpy.util import plot_result
import tops.dyn_models.FAST_FMU as fast

if __name__ == '__main__':
    #TODO
    # Prepare PMSM class for FMU input
    # Make functions for easier storage of results from PMSM
    # Take care of per unit considerations in the integration of the FMU

    ### SIMULATION SETTINGS ###
    t = 0
    dt = 5e-3
    t_end = 5
    ###

    ###  Initiate TOPS ###
    import tops.ps_models.k2a_2WT as model_data
    model = model_data.load()

    # Power system model
    ps = dps.PowerSystemModel(model=model)
    ps.init_dyn_sim()
    print(max(abs(ps.state_derivatives(0, ps.x_0, ps.v_0))))
    x_0 = ps.x_0.copy()

    # Solver
    sol = dps_sol.ModifiedEulerDAE(ps.state_derivatives, ps.solve_algebraic, 0, x_0, t_end, max_step=dt)

    res = defaultdict(list)
    t_0 = time.time()
    t = 0

    ### Initiate the PMSMs ###
    # Define parameters
    pmsm_params = {
        "rs": 0.03,
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

    prime_mover_params = {"T_pm" : 1e-1,
                        "speed_0" : 0.5,
                        "torque_0" : 0.5
    }

    # Create PMSM instance
    pmsm1 = PMSM(pmsm_params, MSC_params = MSC_params, prime_mover_params = prime_mover_params)         # Initiate the PMSM
    pmsm2 = PMSM(pmsm_params, MSC_params = MSC_params, prime_mover_params = prime_mover_params)         # Initiate the PMSM


    # Initiate the FAST FMU
    fmu_filename = 'C:/Users/larsi/OpenFAST/OpenFASTFMU2PF_Export/fast.fmu'

    fmu, vrs = fast.initiate_FAST_FMU(fmu_filename, start_time=t)

    ### Simulation events ###
    def event1(t, dt):
        pmsm1.set_prime_mover_reference(speed_ref=0.3, torque_ref=0.4, ramp_time=3, dt=dt, current_time=t)

    def event2(t, dt):
        pmsm1.set_prime_mover_reference(speed_ref=0.5, torque_ref=0.7, ramp_time=3, dt=dt, current_time=t)

    def event3(t, dt):
        pmsm2.set_prime_mover_reference(speed_ref=0.4, torque_ref=0.5, ramp_time=3, dt=dt, current_time=t)

    def event4(t, dt):
        pmsm2.set_prime_mover_reference(speed_ref=0.6, torque_ref=0.8, ramp_time=3, dt=dt, current_time=t)


    event_flag1 = True
    event_flag2 = True
    event_flag3 = True
    event_flag4 = True

    ######### START SIMULATION #########
    print("Starting simulation")
    while t < t_end:
        sys.stdout.write("\r%d%%" % (t/(t_end)*100))

        # Update the states
        if t > 5 and event_flag1:
            event1(t, dt)
            event_flag1 = False

        if t > 7 and event_flag2:
            event3(t, dt)
            event_flag2 = False

        if t > 11 and event_flag3:
            event2(t, dt)
            event_flag3 = False

        if t > 14 and event_flag4:
            event4(t, dt)
            event_flag4 = False

        pmsm1.update_states(t=t, dt=sol.dt)
        Pref1 = pmsm1.get_Pe()
        ps.vsc["GridSideConverter"].set_pref(Pref1, 0)

        pmsm2.update_states(t=t, dt=sol.dt)
        Pref2 = pmsm2.get_Pe()
        ps.vsc["GridSideConverter"].set_pref(Pref2, 1)

        # Set the torque input to the FMU
        fmu.setReal([vrs['GenPwr']], [pmsm1.get_Pe()])
        fmu.doStep(currentCommunicationPoint=t, communicationStepSize=dt)
        
        # fast.get_FMU_results(fmu, vrs, res, t)


        # Simulate next step
        result = sol.step()
        x = sol.y
        v = sol.v
        t = sol.t

        dx = ps.ode_fun(0, ps.x_0)

        # Store result
        res['t'].append(t)

        # Getting FMU results
        for name, value in vrs.items():
            res[name].append(fmu.getReal([value])[0])
        
        ### PMSM results ###
        res['pmsm1_speed'].append(pmsm1.speed)
        res['pmsm1_torque'].append(pmsm1.get_Te())
        res['pmsm1_speed_ref'].append(pmsm1.primemover.speed_ref)
        res['pmsm1_speed_error'].append(-(pmsm1.primemover.speed_ref - pmsm1.primemover.speed))
        res['pmsm1_i_d'].append(pmsm1.i_d)
        res['pmsm1_i_q'].append(pmsm1.i_q)
        res['pmsm1_v_d'].append(pmsm1.converter.vd)
        res['pmsm1_v_q'].append(pmsm1.converter.vq)

        res['pmsm2_speed'].append(pmsm2.speed)
        res['pmsm2_torque'].append(pmsm2.get_Te())
        res['pmsm2_speed_ref'].append(pmsm2.primemover.speed_ref)
        res['pmsm2_speed_error'].append(-(pmsm2.primemover.speed_ref - pmsm2.primemover.speed))
        res['pmsm2_i_d'].append(pmsm2.i_d)
        res['pmsm2_i_q'].append(pmsm2.i_q)
        res['pmsm2_v_d'].append(pmsm2.converter.vd)
        res['pmsm2_v_q'].append(pmsm2.converter.vq)

        res['primemover1_speed'].append(pmsm1.primemover.speed)
        res['primemover1_torque'].append(pmsm1.primemover.torque)

        res['primemover2_speed'].append(pmsm2.primemover.speed)
        res['primemover2_torque'].append(pmsm2.primemover.torque)

        res['electrical_power1'].append(pmsm1.get_Pe())
        res['mechanical_power1'].append(pmsm1.get_Pm())

        res['electrical_power2'].append(pmsm2.get_Pe())
        res['mechanical_power2'].append(pmsm2.get_Pm())

        # All generators
        res["Generator speeds"].append(ps.gen['GEN'].speed(x, v).copy())
        res["Generator current injections"].append(abs(ps.gen['GEN'].i(x, v).copy()))

        # GridSideConverter results
        # WT1
        res["WT1 p_e"].append(ps.vsc["GridSideConverter"].p_e(x, v)[0].copy())
        res["WT1 q_e"].append(ps.vsc["GridSideConverter"].q_e(x, v)[0].copy())
        res["WT1 i_inj"].append(abs(ps.vsc["GridSideConverter"].i_inj(x, v)[0].copy()))
        res["WT1 p_ref"].append(ps.vsc["GridSideConverter"].p_ref(x, v)[0].copy())

        # WT2
        res["WT2 p_e"].append(ps.vsc["GridSideConverter"].p_e(x, v)[1].copy())
        res["WT2 q_e"].append(ps.vsc["GridSideConverter"].q_e(x, v)[1].copy())
        res["WT2 i_inj"].append(abs(ps.vsc["GridSideConverter"].i_inj(x, v)[1].copy()))
        res["WT2 p_ref"].append(ps.vsc["GridSideConverter"].p_ref(x, v)[1].copy())

    

    # Plot the mechanical results
    fig1, axs1 = plt.subplots(3, 1, figsize=(15, 15), sharex=True)

    # Mechanical results
    axs1[0].plot(res['t'], res['pmsm1_speed'], label='PMSM1 Speed')
    axs1[0].plot(res['t'], res['pmsm2_speed'], label='PMSM2 Speed')
    axs1[0].plot(res['t'], res['primemover1_speed'], label='Prime Mover 1 Speed')
    axs1[0].plot(res['t'], res['primemover2_speed'], label='Prime Mover 2 Speed')
    axs1[0].plot(res['t'], res['pmsm1_speed_ref'], label='PMSM1 Speed Reference', linestyle='--')
    axs1[0].plot(res['t'], res['pmsm2_speed_ref'], label='PMSM2 Speed Reference', linestyle='--')
    axs1[0].set_ylabel('Speed (pu)')
    axs1[0].legend()
    axs1[0].grid(True)

    axs1[1].plot(res['t'], res['pmsm1_torque'], label='PMSM1 Torque')
    axs1[1].plot(res['t'], res['primemover1_torque'], label='Prime Mover 1 Torque')
    axs1[1].plot(res['t'], res['pmsm1_speed_error'], label='PMSM1 Speed Error', linestyle='--')
    axs1[1].set_ylabel('Torque (pu)')
    axs1[1].legend()
    axs1[1].grid(True)

    axs1[2].plot(res['t'], res['pmsm2_torque'], label='PMSM2 Torque')
    axs1[2].plot(res['t'], res['primemover2_torque'], label='Prime Mover 2 Torque')
    axs1[2].plot(res['t'], res['pmsm2_speed_error'], label='PMSM2 Speed Error', linestyle='--')
    axs1[2].set_ylabel('Torque (pu)')
    axs1[2].legend()
    axs1[2].grid(True)

    plt.xlabel('Time (s)')
    plt.tight_layout()
    plt.savefig("Figures/COSIM_Mechanical_response.png")

    fig = plt.figure(figsize=(12, 8))

    # Top plot (Electrical Power) spanning both columns
    ax1 = plt.subplot2grid((3, 2), (0, 0), colspan=2, fig=fig)
    ax1.plot(res['t'], res['electrical_power1'], label='PMSM1 Electrical Power')
    ax1.plot(res['t'], res['electrical_power2'], label='PMSM2 Electrical Power')
    ax1.set_ylabel('Power (pu)')
    ax1.legend()
    ax1.grid(True)
    ax1.set_title('Electrical Power')

    # Second row (PMSM1 Currents and Voltages)
    ax2 = plt.subplot2grid((3, 2), (1, 0), fig=fig)
    ax2.plot(res['t'], res['pmsm1_i_d'], label='PMSM1 i_d')
    ax2.plot(res['t'], res['pmsm1_i_q'], label='PMSM1 i_q')
    ax2.set_ylabel('Current (pu)')
    ax2.legend()
    ax2.grid(True)
    ax2.set_title('PMSM1 Currents')

    ax3 = plt.subplot2grid((3, 2), (1, 1), fig=fig)
    ax3.plot(res['t'], res['pmsm1_v_d'], label='PMSM1 v_d')
    ax3.plot(res['t'], res['pmsm1_v_q'], label='PMSM1 v_q')
    ax3.set_ylabel('Voltage (pu)')
    ax3.legend()
    ax3.grid(True)
    ax3.set_title('PMSM1 Voltages')

    # Third row (PMSM2 Currents and Voltages)
    ax4 = plt.subplot2grid((3, 2), (2, 0), fig=fig)
    ax4.plot(res['t'], res['pmsm2_i_d'], label='PMSM2 i_d')
    ax4.plot(res['t'], res['pmsm2_i_q'], label='PMSM2 i_q')
    ax4.set_ylabel('Current (pu)')
    ax4.legend()
    ax4.grid(True)
    ax4.set_title('PMSM2 Currents')

    ax5 = plt.subplot2grid((3, 2), (2, 1), fig=fig)
    ax5.plot(res['t'], res['pmsm2_v_d'], label='PMSM2 v_d')
    ax5.plot(res['t'], res['pmsm2_v_q'], label='PMSM2 v_q')
    ax5.set_ylabel('Voltage (pu)')
    ax5.legend()
    ax5.grid(True)
    ax5.set_title('PMSM2 Voltages')

    # Set the common xlabel and adjust spacing
    fig.text(0.5, 0.04, 'Time (s)', ha='center')
    plt.tight_layout(rect=[0, 0.05, 1, 1])  # Adjust bottom margin
    plt.savefig("Figures/COSIM_Electrical_response.png")

    # Generator Speeds
    fig, axs = plt.subplots(3, 1, figsize=(15, 10), sharex=True)

    axs[0].plot(res['t'], res["Generator speeds"], label='Generator speeds')
    axs[0].set_ylabel('Speed (pu)')
    axs[0].legend()
    axs[0].grid(True)
    axs[0].set_title('Generator Speeds vs Time')

    # Generator Current Injections
    axs[1].plot(res['t'], res["Generator current injections"], label='Generator current injections')
    axs[1].set_ylabel('Current (pu)')
    axs[1].legend()
    axs[1].grid(True)
    axs[1].set_title('Generator Current Injections vs Time')

    # Grid Reaction to Power Injections
    axs[2].plot(res['t'], res["WT1 p_e"], label='WT1 Power Injection')
    axs[2].plot(res['t'], res["WT2 p_e"], label='WT2 Power Injection')
    axs[2].set_xlabel('Time (s)')
    axs[2].set_ylabel('Power Injection (pu)')
    axs[2].legend()
    axs[2].grid(True)
    axs[2].set_title('Grid Reaction to Power Injections from WT1 and WT2')

    plt.tight_layout()
    plt.savefig("Figures/COSIM_WT_system_response.png")

    # FMU results
    fig, axs = plt.subplots(4, 1, figsize=(15, 12), sharex=True)

    axs[0].plot(res['t'], res['GenSpeed'], label='Generator Speed')
    axs[0].set_ylabel('Speed (pu)')
    axs[0].legend()
    axs[0].grid(True)
    axs[0].set_title('FMU Generator Speed vs Time')

    axs[1].plot(res['t'], res['GenSpdOrTrq'], label='Generator Speed or Torque')
    axs[1].set_ylabel('Speed or Torque (pu)')
    axs[1].legend()
    axs[1].grid(True)
    axs[1].set_title('FMU Generator Speed or Torque vs Time')

    axs[2].plot(res['t'], res['GenPwr'], label='Generator Power')
    axs[2].set_ylabel('Power (pu)')
    axs[2].legend()
    axs[2].grid(True)
    axs[2].set_title('FMU Generator Torque vs Time')

    axs[3].plot(res['t'], res['RotSpeed'], label='Generator Power')
    axs[3].set_xlabel('Time (s)')
    axs[3].set_ylabel('Power (pu)')
    axs[3].legend()
    axs[3].grid(True)
    axs[3].set_title('FMU Generator Power vs Time')

    axs[3].plot(res['t'], res['HSShftTq'], label='High Speed Shaft Torque')
    axs[3].set_xlabel('Time (s)')
    axs[3].set_ylabel('Torque (pu)')
    axs[3].legend()
    axs[3].grid(True)
    axs[3].set_title('FMU High Speed Shaft Torque vs Time')

    plt.tight_layout()
    plt.savefig("Figures/COSIM_FMU_response.png")


    fast.terminate_FAST_FMU(fmu)




