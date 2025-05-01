# Result class to store the result of the simulation
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from tops.cosim_models.windturbine import WindTurbine
from tops.dynamic import PowerSystemModel
import os
import plotly.graph_objects as go
import plotly.io as pio



class Results:
    def __init__(self):
        self.results = defaultdict(list)

    # region Store results
    def store_time(self, t : float):
        self.results["Time"].append(t)

    def store_time_elec(self, t : float):
        self.results["Time elec"].append(t)

    def store_fmu_results(self, WT : WindTurbine):
        for key, value in WT.fast.vrs.items():
            self.results[f"{WT.name}_{key}"].append(WT.fast.fmu.getReal([value])[0])

    # region store pmsm    
    def store_pmsm_results(self, WT : WindTurbine):
        self.results[f"{WT.name}_PMSM T_e"].append(WT.pmsm.T_e())
        self.results[f"{WT.name}_PMSM T_e_ref"].append(WT.pmsm.torque_ref * WT.pmsm.params["T_r"])
        self.results[f"{WT.name}_PMSM t_e_ref"].append(WT.pmsm.torque_ref)
        self.results[f"{WT.name}_PMSM speed"].append(WT.pmsm.speed)
        self.results[f"{WT.name}_PMSM SPEED"].append(WT.pmsm.SPEED)
        self.results[f"{WT.name}_PMSM p_e"].append(WT.pmsm.p_e())
        self.results[f"{WT.name}_PMSM P_e"].append(WT.pmsm.P_e())
        self.results[f"{WT.name}_PMSM q_e"].append(WT.pmsm.q_e())
        self.results[f"{WT.name}_PMSM Q_e"].append(WT.pmsm.Q_e())
        self.results[f"{WT.name}_PMSM i_d"].append(WT.pmsm.i_d)
        self.results[f"{WT.name}_PMSM i_q"].append(WT.pmsm.i_q)
        self.results[f"{WT.name}_PMSM v_d"].append(WT.pmsm.v_d)
        self.results[f"{WT.name}_PMSM v_q"].append(WT.pmsm.v_q)
        self.results[f"{WT.name}_PMSM i_q_ref"].append(WT.pmsm.i_q_ref)
        self.results[f"{WT.name}_PMSM i_d_ref"].append(WT.pmsm.i_d_ref)
        self.results[f"{WT.name}_PMSM v_qII"].append(WT.pmsm.v_qII)
        self.results[f"{WT.name}_PMSM v_dII"].append(WT.pmsm.v_dII)
        self.results[f"{WT.name}_PMSM theta elec"].append(WT.pmsm.theta_elec)
        self.results[f"{WT.name}_PMSM iq_controller P_term"].append(WT.pmsm.pi_controller_iq.P_term)
        self.results[f"{WT.name}_PMSM iq_controller I_term"].append(WT.pmsm.pi_controller_iq.I_term)
        self.results[f"{WT.name}_PMSM id_controller D_term"].append(WT.pmsm.pi_controller_id.D_term)
        self.results[f"{WT.name}_PMSM vdc_controller P_term"].append(WT.pmsm.pi_controller_vdc.P_term)
        self.results[f"{WT.name}_PMSM vdc_controller I_term"].append(WT.pmsm.pi_controller_vdc.I_term)
        self.results[f"{WT.name}_PMSM vdc_controller D_term"].append(WT.pmsm.pi_controller_vdc.D_term)
        self.results[f"{WT.name}_PMSM P_mech"].append(WT.pmsm.P_mech())
    
    # endregion

    # region store DClink
    def store_dclink_results(self, WT : WindTurbine, ps : PowerSystemModel , x, v):
        self.results[f"{WT.name}_DC_vdc"].append(WT.dclink.vdc)
        self.results[f"{WT.name}_DC_vdc_ref"].append(WT.dclink.vdc_ref)
        self.results[f"{WT.name}_DC_i_chopper"].append(WT.dclink.i_chopper())
        self.results[f"{WT.name}_DC_p_adjust"].append(WT.dclink.p_adjust())
        self.results[f"{WT.name}_DC_p_gsc"].append(WT.p_gsc)
        # self.results[f"{WT.name}_DC_pref_pu"].append(WT.pref_dclink_pu)
        self.results[f"{WT.name}_DC_gsc_p_ref"].append(WT.dclink.gsc_p_ref())
        self.results[f"{WT.name}_DC_duty"].append(WT.dclink.duty())
        # self.results[f"{WT.name}_MSC_GSC_pe_diff"].append(WT.msc.p_e_dq() - WT.dclink.gsc_p_ref())
        self.results[f"{WT.name}_MSC_GSC_pe_diff"].append(WT.pmsm.p_e() - WT.calculate_p_gsc(ps, x, v))
        self.results[f"{WT.name}_DC_x_pref_adj"].append(WT.dclink.x_pref_adj)
    # endregion

    # region store PQ
    def store_gsc_results_PQ(self, WT : WindTurbine , ps : PowerSystemModel, x, v):
        self.results[f"{WT.name}_GSC_p_e"].append(ps.vsc["GridSideConverter_PQ"].p_e(x, v)[WT.index].copy())
        self.results[f"{WT.name}_GSC_q_e"].append(ps.vsc["GridSideConverter_PQ"].q_e(x, v)[WT.index].copy())
        self.results[f"{WT.name}_GSC_i_inj"].append(abs(ps.vsc["GridSideConverter_PQ"].i_inj(x, v)[WT.index].copy()))
        self.results[f"{WT.name}_GSC_p_ref_grid"].append(ps.vsc["GridSideConverter_PQ"].par["p_ref_grid"][WT.index].copy())
        self.results[f"{WT.name}_GSC_q_ref_grid"].append(ps.vsc["GridSideConverter_PQ"].par["q_ref_grid"][WT.index].copy())
        self.results[f"{WT.name}_GSC_bus_voltage"].append(abs(ps.vsc["GridSideConverter_PQ"].v_t(x, v)[WT.index].copy()))
    # endregion

    # region store PV
    def store_gsc_results_PV(self, WT : WindTurbine , ps : PowerSystemModel, x, v):
        self.results[f"{WT.name}_GSC_p_e"].append(ps.vsc["GridSideConverter_PV"].p_e(x, v)[WT.index].copy())
        self.results[f"{WT.name}_GSC_q_e"].append(ps.vsc["GridSideConverter_PV"].q_e(x, v)[WT.index].copy())
        self.results[f"{WT.name}_GSC_i_inj"].append(abs(ps.vsc["GridSideConverter_PV"].i_inj(x, v)[WT.index].copy()))
        self.results[f"{WT.name}_GSC_p_ref_grid"].append(ps.vsc["GridSideConverter_PV"].par["p_ref_grid"][WT.index].copy())
        # self.results[f"{WT.name}_GSC_q_ref_grid"].append(ps.vsc["GridSideConverter_PV"].q_ref_grid[WT.index].copy())
        self.results[f"{WT.name}_GSC_v_ref_grid"].append(ps.vsc["GridSideConverter_PV"].par["v_ref_grid"][WT.index].copy())
        self.results[f"{WT.name}_GSC_bus_voltage"].append(abs(ps.vsc["GridSideConverter_PV"].v_t(x, v)[WT.index].copy()))
        self.results[f"{WT.name}_GSC_q_ref_grid"].append(ps.vsc["GridSideConverter_PV"].q_ref_grid(x, v)[WT.index].copy())
    # endregion


    # region store generator
    def store_generator_results(self, ps : PowerSystemModel, x, v, index = None):
        if index is not None:
            self.results[f"Generator_{index}_speed"].append(ps.gen['GEN'].speed(x, v)[index].copy())
            self.results[f"Generator_{index}_i_inj"].append(abs(ps.gen['GEN'].i(x, v)[index].copy()))
            self.results[f"Generator_{index}_p_e"].append(ps.gen['GEN'].p_e(x, v)[index].copy())

        else:
            self.results["Generators_speed"].append(ps.gen['GEN'].speed(x, v).copy())
            self.results["Generators_i_inj"].append(abs(ps.gen['GEN'].i(x, v).copy()))
            self.results["Generators_p_e"].append(ps.gen['GEN'].p_e(x, v).copy())


    # endregion


    # region Plot PMSM

    def plot_pmsm_overview(self, sim_name : str, WT : WindTurbine, t_init : float = None):
        fig, axs = plt.subplots(3, 2, figsize=(15, 10), sharex=True)


        # Add text on the top of the figure with the current controller tuning
        controller_tuning_text = f"iq_controller: Kp={WT.pmsm.pi_controller_iq.kp}, Ki={WT.pmsm.pi_controller_iq.ki}"
        fig.suptitle(controller_tuning_text, fontsize=12, fontweight='bold')
        # First column, first row: PMSM p_e and q_e
        axs[0, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM p_e'], label='p_e')
        axs[0, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM q_e'], label='q_e')
        axs[0, 0].set_ylabel('Power (pu)')
        axs[0, 0].legend()
        axs[0, 0].grid(True)
        axs[0, 0].set_title('PMSM Electrical Power')

        # First column, second row: PMSM currents
        axs[1, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM i_d'], label='i_d')
        axs[1, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM i_q'], label='i_q')
        axs[1, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM i_d_ref'], label='i_d_ref', linestyle='--')
        axs[1, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM i_q_ref'], label='i_q_ref', linestyle='--')
        
        axs[1, 0].set_ylabel('Current (pu)')
        axs[1, 0].legend()
        axs[1, 0].grid(True)
        axs[1, 0].set_title('PMSM Currents')

        # First column, third row: PMSM voltages
        axs[2, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM v_d'], label='v_d')
        axs[2, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM v_q'], label='v_q')
        axs[2, 0].set_ylabel('Voltage (pu)')
        axs[2, 0].legend()
        axs[2, 0].grid(True)
        axs[2, 0].set_title('PMSM Voltages')

        # Second column, first row: PMSM T_e and T_e_ref
        axs[0, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM T_e'], label='T_e')
        axs[0, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM T_e_ref'], label='T_e_ref')
        axs[0, 1].set_ylabel('Torque (Nm)')
        axs[0, 1].legend()
        axs[0, 1].grid(True)
        axs[0, 1].set_title('PMSM Torque')

        # Second column, second row: GenSpeed and PMSM speed
        # axs[1, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_GenSpeed'], label='GenSpeed')
        axs[1, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM SPEED'], label='PMSM SPEED')
        axs[1, 1].set_ylabel('Speed (rpm)')
        axs[1, 1].legend()
        axs[1, 1].grid(True)
        axs[1, 1].set_title('Generator and PMSM Speed')

        # Second column, third row: Current controller feedforward terms
        axs[2, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM v_qII'], label='v_qII')
        axs[2, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM v_dII'], label='v_dII')
        axs[2, 1].set_ylabel('Voltage (pu)')
        axs[2, 1].legend()
        axs[2, 1].grid(True)
        axs[2, 1].set_title('Current Controller Feedforward Terms')

        # Set the common xlabel and adjust spacing
        fig.text(0.5, 0.04, 'Time (s)', ha='center')
        plt.tight_layout(rect=[0, 0.05, 1, 1])  # Adjust bottom margin
        # plt.savefig(f"Figures/co_sim/PMSM_Overview_{sim_name}.png")

        # Set the x-axis limit if t_init is provided
        if t_init is not None:
            plt.xlim(left=t_init)

        # Set the common xlabel and adjust spacing
        fig.text(0.5, 0.04, 'Time (s)', ha='center')
        plt.tight_layout(rect=[0, 0.05, 1, 1])  # Adjust bottom margin

        # Create directory if it does not exist
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/PMSM_Overview_{sim_name}.pdf", format='pdf')


    def plot_pmsm_overview_interactive(self, sim_name: str, WT: WindTurbine):
        fig = go.Figure()

        # First column, first row: PMSM p_e and q_e
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_PMSM p_e'], mode='lines', name='p_e'))
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_PMSM q_e'], mode='lines', name='q_e'))

        # First column, second row: PMSM currents
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_PMSM i_d'], mode='lines', name='i_d'))
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_PMSM i_q'], mode='lines', name='i_q'))
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_PMSM i_d_ref'], mode='lines', name='i_d_ref', line=dict(dash='dash')))
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_PMSM i_q_ref'], mode='lines', name='i_q_ref', line=dict(dash='dash')))

        # First column, third row: PMSM voltages
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_PMSM v_d'], mode='lines', name='v_d'))
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_PMSM v_q'], mode='lines', name='v_q'))

        # Second column, first row: PMSM T_e and T_e_ref
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_PMSM T_e'], mode='lines', name='T_e'))
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_PMSM T_e_ref'], mode='lines', name='T_e_ref'))

        # Second column, second row: GenSpeed and PMSM speed
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_GenSpeed'], mode='lines', name='GenSpeed'))
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_PMSM SPEED'], mode='lines', name='PMSM SPEED'))

        # Second column, third row: PMSM DC-link voltage vs DC-link voltage
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_PMSM_vdc'], mode='lines', name='PMSM vdc'))
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_DC_vdc'], mode='lines', name='DC vdc'))
        fig.add_trace(go.Scatter(x=self.results['Time elec'], y=self.results[f'{WT.name}_DC_vdc_ref'], mode='lines', name='DC vdc_ref'))

        # Set the common xlabel and adjust spacing
        fig.update_layout(title='PMSM Overview', xaxis_title='Time (s)', yaxis_title='Values', legend_title='Legend')

        # Create directory if it does not exist
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Save the figure as an interactive HTML file
        pio.write_html(fig, file=f"{output_dir}/PMSM_Overview_{sim_name}.html", auto_open=False)
        

        # Example usage
        # plot_pmsm_overview_interactive('simulation_name', wind_turbine_instance)
    
    # endregion

    # region plot FMU

    def plot_fmu_overview(self, sim_name : str, WT : WindTurbine):
        fig, axs = plt.subplots(3, 2, figsize=(15, 10), sharex=True)

        # # First column, first row: HSShftTq
        # axs[0, 0].plot(self.results['Time'], self.results[f'{WT.name}_HSShftTq'], label='HSShftTq')
        # axs[0, 0].set_ylabel('Torque (Nm)')
        # axs[0, 0].legend()
        # axs[0, 0].grid(True)
        # axs[0, 0].set_title('High Speed Shaft Torque')

        # First column, first row: Difference between GenTq and PMSM T_e
        diff = np.array(self.results[f'{WT.name}_GenTq']) - np.array(self.results[f'{WT.name}_GenSpdOrTrq'])
        axs[0, 0].plot(self.results['Time'], diff, label='GenTq - PMSM T_e')
        axs[0, 0].set_ylabel('Torque Difference (Nm)')
        axs[0, 0].legend()
        axs[0, 0].grid(True)
        axs[0, 0].set_title('Torque Difference (GenTq - PMSM T_e)')

        # First column, second row: GenTq and GenSpdOrTrq
        axs[1, 0].plot(self.results['Time'], self.results[f'{WT.name}_GenTq'], label='GenTq')
        axs[1, 0].plot(self.results['Time'], self.results[f'{WT.name}_GenSpdOrTrq'], label='GenSpdOrTrq')
        axs[1, 0].set_ylabel('Torque (Nm)')
        axs[1, 0].legend()
        axs[1, 0].grid(True)
        axs[1, 0].set_title('Generator Torque')

        # First column, third row: Speeds
        axs[2, 0].plot(self.results['Time'], self.results[f'{WT.name}_GenSpeed'], label='GenSpeed')
        axs[2, 0].plot(self.results['Time'], self.results[f'{WT.name}_RefGenSpd'], label='RefGenSpd')
        axs[2, 0].plot(self.results['Time'], self.results[f'{WT.name}_RotSpeed'], label='RotSpeed')
        axs[2, 0].set_ylabel('Speed (rpm)')
        axs[2, 0].legend()
        axs[2, 0].grid(True)
        axs[2, 0].set_title('Generator Speed')

        # Second column, first row: Windspeed
        axs[0, 1].plot(self.results['Time'], self.results[f'{WT.name}_Wind1VelX'], label='Wind1VelX')
        axs[0, 1].plot(self.results['Time'], self.results[f'{WT.name}_RtVAvgxh'], label='RtVAvgxh')
        axs[0, 1].set_ylabel('Wind speed [m/s]')
        axs[0, 1].legend()
        axs[0, 1].grid(True)
        axs[0, 1].set_title('Wind speeds')

        # Second column, third row: Blade pitch
        axs[2, 1].plot(self.results['Time'], self.results[f'{WT.name}_BldPitch1'], label='BldPitch1')
        axs[2, 1].set_ylabel('Pitch (deg)')
        axs[2, 1].legend()
        axs[2, 1].grid(True)
        axs[2, 1].set_title('Blade pitch')

        # Second column, second row: GenPwr and ElecPwrCom
        axs[1, 1].plot(self.results['Time'], self.results[f'{WT.name}_GenPwr'], label='GenPwr')
        axs[1, 1].plot(self.results['Time'], self.results[f'{WT.name}_ElecPwrCom'], label='ElecPwrCom')
        axs[1, 1].plot(self.results['Time'], np.array(self.results[f'{WT.name}_GenTq']) * np.array(self.results[f'{WT.name}_GenSpeed']) * (2 * np.pi / 60), label='Mechanical Power')
        axs[1, 1].set_ylabel('Power [kW]')
        axs[1, 1].legend()
        axs[1, 1].grid(True)
        axs[1, 1].set_title('Generator and Electrical Power')

        # Set the common xlabel and adjust spacing
        fig.text(0.5, 0.04, 'Time (s)', ha='center')
        plt.tight_layout(rect=[0, 0.05, 1, 1])  # Adjust bottom margin
        
        # Create directory if it does not exist
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/FMU_Overview_{sim_name}.pdf", format='pdf')

    # endregion    
        
    # region plot MSC
    def plot_msc_overview(self, sim_name : str, WT : WindTurbine):
        fig, axs = plt.subplots(3, 2, figsize=(15, 10), sharex=True)

        # First column, first row: MSC v_q, v_d, v_q_ref, and v_d_ref
        axs[0, 0].plot(self.results['Time'], self.results[f'{WT.name}_MSC_v_q'], label='v_q')
        axs[0, 0].plot(self.results['Time'], self.results[f'{WT.name}_MSC_v_d'], label='v_d')
        axs[0, 0].plot(self.results['Time'], self.results[f'{WT.name}_MSC_v_q_ctrl'], label='v_q_ctrl')
        axs[0, 0].plot(self.results['Time'], self.results[f'{WT.name}_MSC_v_d_ctrl'], label='v_d_ctrl')
        axs[0, 0].set_ylabel('Voltage (pu)')
        axs[0, 0].legend()
        axs[0, 0].grid(True)
        axs[0, 0].set_title('MSC Voltages and References')

        # First column, second row: MSC v_qII and v_dII
        axs[1, 0].plot(self.results['Time'], self.results[f'{WT.name}_PMSM v_qII'], label='v_qII')
        axs[1, 0].plot(self.results['Time'], self.results[f'{WT.name}_PMSM v_dII'], label='v_dII')
        axs[1, 0].set_ylabel('Voltage (pu)')
        axs[1, 0].legend()
        axs[1, 0].grid(True)
        axs[1, 0].set_title('MSC Feedforward Terms')

        # First column, third row: MSC i_dc and v_dc
        axs[2, 0].plot(self.results['Time'], self.results[f'{WT.name}_MSC_i_dc'], label='i_dc')
        axs[2, 0].plot(self.results['Time'], self.results[f'{WT.name}_MSC_v_dc'], label='v_dc')
        axs[2, 0].set_ylabel('Current (pu)')
        axs[2, 0].legend()
        axs[2, 0].grid(True)
        axs[2, 0].set_title('MSC DC Current and Voltage')

        # Second column, first row: MSC p_e_dc and p_e_dq
        axs[0, 1].plot(self.results['Time'], self.results[f'{WT.name}_MSC_p_e_dc'], label='p_e_dc')
        axs[0, 1].plot(self.results['Time'], self.results[f'{WT.name}_MSC_p_e_dq'], label='p_e_dq')
        axs[0, 1].set_ylabel('Power (pu)')
        axs[0, 1].legend()
        axs[0, 1].grid(True)
        axs[0, 1].set_title('MSC Power')

        # Second column, second row: PMSM i_q and i_q_ref
        axs[1, 1].plot(self.results['Time'], self.results[f'{WT.name}_PMSM i_q'], label='i_q')
        axs[1, 1].plot(self.results['Time'], self.results[f'{WT.name}_PMSM i_q_ref'], label='i_q_ref')
        axs[1, 1].set_ylabel('Current (pu)')
        axs[1, 1].legend()
        axs[1, 1].grid(True)
        axs[1, 1].set_title('PMSM q-axis Current')

        # Second column, third row: PMSM i_d and i_d_ref
        axs[2, 1].plot(self.results['Time'], self.results[f'{WT.name}_PMSM i_d'], label='i_d')
        axs[2, 1].plot(self.results['Time'], self.results[f'{WT.name}_PMSM i_d_ref'], label='i_d_ref')
        axs[2, 1].set_ylabel('Current (pu)')
        axs[2, 1].legend()
        axs[2, 1].grid(True)
        axs[2, 1].set_title('PMSM d-axis Current')

        # Set the common xlabel and adjust spacing
        fig.text(0.5, 0.04, 'Time (s)', ha='center')
        plt.tight_layout(rect=[0, 0.05, 1, 1])  # Adjust bottom margin
        
        # Create directory if it does not exist
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/MSC_Overview_{sim_name}.pdf", format='pdf')

    # endregion

    # region Plot DC-link
    
    def plot_dclink_overview(self, sim_name : str, WT : WindTurbine):
        fig, axs = plt.subplots(3, 2, figsize=(15, 10), sharex=True)

        # First column, first row: DC-link voltage and reference voltage
        axs[0, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_vdc'], label='vdc')
        axs[0, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_vdc_ref'], label='vdc_ref')
        axs[0, 0].set_ylabel('Voltage (pu)')
        axs[0, 0].legend()
        axs[0, 0].grid(True)
        axs[0, 0].set_title('DC-link Voltage and Reference')

        # First column, second row: PMSM power
        axs[1, 0].plot(self.results['Time elec'], self.results[f"{WT.name}_PMSM p_e"], label='p_e_pmsm')
        axs[1, 0].set_ylabel('Active power (pu)')
        axs[1, 0].legend()
        axs[1, 0].grid(True)
        axs[1, 0].set_title('PMSM Active Power')


        # First column, third row: Chopper duty cycle
        axs[2, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_duty'], label='duty')
        axs[2, 0].set_ylabel('Duty cycle (pu)')
        axs[2, 0].legend()
        axs[2, 0].grid(True)
        axs[2, 0].set_title('Chopper Duty Cycle')

        ##############################################################################################################

        # Second column, first row: Power adjustment
        axs[0, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_p_adjust'], label='p_adjust')
        axs[0, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_x_pref_adj'], label='x_pref_adj')
        axs[0, 1].set_ylabel('Power (pu)')
        axs[0, 1].legend()
        axs[0, 1].grid(True)
        axs[0, 1].set_title('Power Adjustment and Power Difference')

        # Second column, second row: GSC power and reference power
        axs[1, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_gsc_p_ref'], label='gsc_p_ref')
        # axs[1, 1].plot(self.results['Time elec'], (WT.gsc_sn/(WT.pmsm.params["s_n"]*1e-3))*np.array(self.results[f'{WT.name}_GSC_p_e']), label='gsc p_e')
        axs[1, 1].set_ylabel('Power (pu)')
        axs[1, 1].legend()
        axs[1, 1].grid(True)
        axs[1, 1].set_title('GSC Power Reference')

        # Second column, third row: Chopper current
        axs[2, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_i_chopper'], label='i_chopper')
        axs[2, 1].set_ylabel('Current (pu)')
        axs[2, 1].legend()
        axs[2, 1].grid(True)
        axs[2, 1].set_title('Chopper Current')


        plt.xlabel('Time (s)')
        plt.tight_layout()

        # Create directory if it does not exist
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/DClink_Overview_{sim_name}.pdf", format='pdf')

    # region Plot GSC
    def plot_gsc_overview(self, sim_name : str, WT : WindTurbine):
        fig, axs = plt.subplots(4, 1, figsize=(15, 10), sharex=True)

        # GSC p_e and p_ref_grid
        axs[0].plot(self.results['Time'], self.results[f'{WT.name}_GSC_p_ref_grid'], label='p_ref_grid')
        axs[0].plot(self.results['Time'], self.results[f'{WT.name}_GSC_p_e'], label='p_e')
        axs[0].set_ylabel('Power (pu)')
        axs[0].legend()
        axs[0].grid(True)
        axs[0].set_title('GSC Active Power and Reference')

        # GSC q_e and q_ref_grid
        axs[1].plot(self.results['Time'], self.results[f'{WT.name}_GSC_q_ref_grid'], label='q_ref_grid')
        axs[1].plot(self.results['Time'], self.results[f'{WT.name}_GSC_q_e'], label='q_e')
        axs[1].set_ylabel('Reactive Power (pu)')
        axs[1].legend()
        axs[1].grid(True)
        axs[1].set_title('GSC Reactive Power and Reference')

        # GSC i_inj
        axs[2].plot(self.results['Time'], self.results[f'{WT.name}_GSC_i_inj'], label='i_inj')
        axs[2].set_ylabel('Current (pu)')
        axs[2].legend()
        axs[2].grid(True)
        axs[2].set_title('GSC Current Injection')

        # GSC v_q and v_d
        axs[3].plot(self.results['Time'], self.results[f'{WT.name}_MSC_v_q'], label='v_q')
        axs[3].plot(self.results['Time'], self.results[f'{WT.name}_MSC_v_d'], label='v_d')
        axs[3].set_ylabel('Voltage (pu)')
        axs[3].legend()
        axs[3].grid(True)
        axs[3].set_title('MSC Voltages')

        plt.xlabel('Time (s)')
        plt.tight_layout()

        # Create directory if it does not exist
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/GSC_Overview_{sim_name}.pdf", format='pdf')

    # endregion

    # region Plot TOPS

    def plot_tops_overview(self, sim_name : str, WT : WindTurbine):
        fig, axs = plt.subplots(3, 2, figsize=(15, 10), sharex=True)

        # Generator speeds
        axs[0, 0].plot(self.results['Time'], self.results["Generators_speed"], label='Generator speeds')
        axs[0, 0].set_ylabel('Speed (rpm)')
        axs[0, 0].legend()
        axs[0, 0].grid(True)
        axs[0, 0].set_title('Generator speeds')

        # Generator current injections
        axs[0, 1].plot(self.results['Time'], self.results["Generators_p_e"], label='P_e')
        axs[0, 1].set_ylabel('Active power (pu)')
        axs[0, 1].legend()
        axs[0, 1].grid(True)
        axs[0, 1].set_title('Generators p_e injections')

        # WT GSC p_e
        axs[1, 0].plot(self.results['Time'], self.results[f'{WT.name}_GSC_p_ref_grid'], label='p_ref_grid')
        axs[1, 0].plot(self.results['Time'], self.results[f'{WT.name}_GSC_p_e'], label=f'{WT.name} GSC p_e')
        axs[1, 0].set_ylabel('Power (pu)')
        axs[1, 0].legend()
        axs[1, 0].grid(True)
        axs[1, 0].set_title(f'{WT.name} GSC Active Power')

        # Terminal voltage of the bus to which the VSC is connected
        if WT.gsc_control == "PV":
            axs[1, 1].plot(self.results['Time'], self.results[f'{WT.name}_GSC_v_ref_grid'], label='Bus Voltage Reference')
        axs[1, 1].plot(self.results['Time'], self.results[f'{WT.name}_GSC_bus_voltage'], label='Bus Voltage')
        # axs[1, 1].plot(self.results['Time'], self.results[f'{WT.name}_GSC_v_d'], label='v_d')
        # axs[1, 1].plot(self.results['Time'], self.results[f'{WT.name}_GSC_v_q'], label='v_q')
        axs[1, 1].set_ylabel('Voltage (pu)')
        axs[1, 1].legend()
        axs[1, 1].grid(True)
        axs[1, 1].set_title('Bus Terminal Voltage')

        # Current injection
        axs[2, 0].plot(self.results['Time'], self.results[f'{WT.name}_GSC_i_inj'], label='i_inj')
        axs[2, 0].set_ylabel('Current (pu)')
        axs[2, 0].legend()
        axs[2, 0].grid(True)
        axs[2, 0].set_title('Current Injection')

        # Reactive power vs reference
        if WT.gsc_control == "PQ":
            axs[2, 1].plot(self.results['Time'], self.results[f'{WT.name}_GSC_q_ref_grid'], label='q_ref_grid')
            axs[2, 1].plot(self.results['Time'], self.results[f'{WT.name}_GSC_q_e'], label='q_e')
            axs[2, 1].set_ylabel('Reactive Power (pu)')
            axs[2, 1].legend()
            axs[2, 1].grid(True)
            axs[2, 1].set_title('Reactive Power vs Reference')
        elif WT.gsc_control == "PV":
            axs[2, 1].plot(self.results['Time'], self.results[f'{WT.name}_GSC_q_e'], label='q_e')
            axs[2, 1].plot(self.results['Time'], self.results[f'{WT.name}_GSC_q_ref_grid'], label='q_ref_grid')
            axs[2, 1].set_ylabel('Reactive Power (pu)')
            axs[2, 1].legend()
            axs[2, 1].grid(True)
            axs[2, 1].set_title('Reactive Power')

        plt.xlabel('Time (s)')
        plt.tight_layout()

        # Create directory if it does not exist
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/TOPS_Overview_{sim_name}.pdf", format='pdf')
    # endregion

    # region Plot Controllers
    def plot_controllers(self, sim_name : str, WT : WindTurbine):
        fig, axs = plt.subplots(3, 1, figsize=(15, 15), sharex=True)

        # iq_controller P_term
        axs[0].plot(self.results['Time'], self.results[f'{WT.name}_PMSM iq_controller P_term'], label='P_term')
        axs[0].set_ylabel('P_term')
        axs[0].legend()
        axs[0].grid(True)
        axs[0].set_title('iq_controller P_term')

        # iq_controller I_term
        axs[1].plot(self.results['Time'], self.results[f'{WT.name}_PMSM iq_controller I_term'], label='I_term')
        axs[1].set_ylabel('I_term')
        axs[1].legend()
        axs[1].grid(True)
        axs[1].set_title('iq_controller I_term')

        # id_controller D_term
        axs[2].plot(self.results['Time'], self.results[f'{WT.name}_PMSM id_controller D_term'], label='D_term')
        axs[2].set_ylabel('D_term')
        axs[2].legend()
        axs[2].grid(True)
        axs[2].set_title('id_controller D_term')

        plt.xlabel('Time (s)')
        plt.tight_layout()

        # Create directory if it does not exist
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/Controllers_{sim_name}.pdf", format='pdf')


    def plot_vdc_controller_terms(self, sim_name: str, WT: WindTurbine):
        fig, axs = plt.subplots(3, 1, figsize=(15, 15), sharex=True)

        # vdc_controller P_term
        axs[0].plot(self.results['Time'], self.results[f'{WT.name}_PMSM vdc_controller P_term'], label='P_term')
        axs[0].set_ylabel('P_term')
        axs[0].legend()
        axs[0].grid(True)
        axs[0].set_title('vdc_controller P_term')

        # vdc_controller I_term
        axs[1].plot(self.results['Time'], self.results[f'{WT.name}_PMSM vdc_controller I_term'], label='I_term')
        axs[1].set_ylabel('I_term')
        axs[1].legend()
        axs[1].grid(True)
        axs[1].set_title('vdc_controller I_term')

        # vdc_controller D_term
        axs[2].plot(self.results['Time'], self.results[f'{WT.name}_PMSM vdc_controller D_term'], label='D_term')
        axs[2].set_ylabel('D_term')
        axs[2].legend()
        axs[2].grid(True)
        axs[2].set_title('vdc_controller D_term')

        plt.xlabel('Time (s)')
        plt.tight_layout()

        # Create directory if it does not exist
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/VDC_Controller_{sim_name}.pdf", format='pdf')

    # endregion

    # region multirate
    def plot_pmsm_multirate(self, sim_name : str, WT : WindTurbine):
        fig, axs = plt.subplots(3, 1, figsize=(15, 15), sharex=True)

        # PMSM speed
        axs[0].plot(self.results['Time'], self.results[f'{WT.name}_PMSM SPEED'], label='PMSM SPEED')
        axs[0].plot(self.results['Time'], self.results[f'{WT.name}_GenSpeed'], label='GenSpeed')
        axs[0].set_ylabel('Speed (rpm)')
        axs[0].legend()
        axs[0].grid(True)
        axs[0].set_title('PMSM Speed')

        # PMSM Torque
        axs[1].plot(self.results['Time'], self.results[f'{WT.name}_PMSM T_e'], label='T_e')
        axs[1].plot(self.results['Time'], self.results[f'{WT.name}_PMSM T_e_ref'], label='T_e_ref')
        axs[1].set_ylabel('Torque (Nm)')
        axs[1].legend()
        axs[1].grid(True)
        axs[1].set_title('PMSM Torque')

        # PMSM Currents
        axs[2].plot(self.results['Time'], self.results[f'{WT.name}_PMSM i_d'], label='i_d')
        axs[2].plot(self.results['Time'], self.results[f'{WT.name}_PMSM i_q'], label='i_q')
        axs[2].plot(self.results['Time'], self.results[f'{WT.name}_PMSM i_d_ref'], label='i_d_ref')
        axs[2].plot(self.results['Time'], self.results[f'{WT.name}_PMSM i_q_ref'], label='i_q_ref')
        axs[2].set_ylabel('Current (pu)')
        axs[2].legend()
        axs[2].grid(True)
        axs[2].set_title('PMSM Currents')

        plt.xlabel('Time (s)')
        plt.tight_layout()

        # Create directory if it does not exist
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/PMSM_{sim_name}.pdf", format='pdf')


    # region generator torque
    def plot_generator_torque(self, sim_name: str, WT: WindTurbine):
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot generator torque and reference torque
        ax.plot(self.results['Time'], self.results[f'{WT.name}_GenTq'], label='Generator Torque', linewidth=2)
        ax.plot(self.results['Time'], self.results[f'{WT.name}_GenSpdOrTrq'], label='Reference Torque', linestyle='--', linewidth=2)

        # Set axis labels with larger font size
        ax.set_ylabel('Torque (kNm)', fontsize=14)
        ax.set_xlabel('Time (s)', fontsize=14)

        # Set title with larger font size
        ax.set_title('Generator Torque vs Reference', fontsize=16)

        # Increase tick label size
        ax.tick_params(axis='both', which='major', labelsize=12)

        # Enable both major and minor grids
        ax.grid(True, which='both', linestyle=':', linewidth=0.8)

        # Add a clearer legend
        ax.legend(fontsize=12, frameon=True, loc='upper left')

        # Adjust layout for publication quality
        plt.tight_layout()

        # Create directory if it does not exist
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/Generator_Torque_{sim_name}.pdf", format='pdf')
    # endregion

    # region fmugen
    def plot_paper_fmugen(self, sim_name: str, WT: WindTurbine):
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import numpy as np
        import os

        # Use Times New Roman globally
        mpl.rcParams['font.family'] = 'Times New Roman'

        fig, axs = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
        lw = 2.0  # Line width

        # (1,1) Generator torque vs reference
        axs[0, 0].plot(self.results['Time'], self.results[f'{WT.name}_GenTq'], label='Generator Torque', linewidth=lw)
        axs[0, 0].plot(self.results['Time'], self.results[f'{WT.name}_GenSpdOrTrq'], label='Reference Torque', linestyle='--', linewidth=lw)
        
        axs[0, 0].set_ylabel('Torque (kNm)')
        axs[0, 0].legend(fontsize=8)
        axs[0, 0].grid(True, linestyle=':')
        axs[0, 0].set_title('Generator Torque vs Reference', fontsize=10)

        # (1,2) Power overview
        axs[0, 1].plot(self.results['Time'], self.results[f'{WT.name}_GenPwr'], label='Electric Power GSC', linewidth=lw)
        axs[0, 1].plot(self.results['Time'], self.results[f'{WT.name}_ElecPwrCom'], label='Electric Power Command', linewidth=lw)
        axs[0, 1].plot(self.results['Time'], np.array(self.results[f'{WT.name}_GenTq']) * np.array(self.results[f'{WT.name}_GenSpeed']) * (2 * np.pi / 60), label='Mechanical Power OpenFAST', linestyle=':', linewidth=lw)
        # axs[0, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_PMSM P_mech'], label='Mechanical Power PMSM', linewidth=lw)
        # axs[0, 1].plot(self.results['Time elec'], np.array(self.results[f'{WT.name}_PMSM P_e']), label='PMSM Active Power', linewidth=lw)
        axs[0, 1].set_ylabel('Power (kW)')
        axs[0, 1].legend(fontsize=8)
        axs[0, 1].grid(True, linestyle=':')
        axs[0, 1].set_title('Power Overview', fontsize=10)

        # (2,1) Speed overview
        axs[1, 0].plot(self.results['Time'], self.results[f'{WT.name}_GenSpeed'], label='Generator Speed', linewidth=lw)
        axs[1, 0].plot(self.results['Time'], self.results[f'{WT.name}_RotSpeed'], label='Rotor Speed', linewidth=lw)
        axs[1, 0].plot(self.results['Time'], self.results[f'{WT.name}_RefGenSpd'], label='RefGenSpeed', linestyle='--', linewidth=lw)
        axs[1, 0].set_ylabel('Speed (rpm)')
        axs[1, 0].legend(fontsize=8)
        axs[1, 0].grid(True, linestyle=':')
        axs[1, 0].set_title('Speed Overview', fontsize=10)

        # (2,2) Blade pitch
        axs[1, 1].plot(self.results['Time'], self.results[f'{WT.name}_BldPitch1'], label='Blade Pitch', linewidth=lw)
        axs[1, 1].set_ylabel('Pitch (deg)')
        axs[1, 1].legend(fontsize=8)
        axs[1, 1].grid(True, linestyle=':')
        axs[1, 1].set_title('Blade Pitch', fontsize=10)

        # Common X label
        fig.text(0.5, 0.02, 'Time (s)', ha='center', fontsize=12)

        # Ticks and layout
        for ax in axs.flat:
            ax.tick_params(labelsize=8)

        plt.tight_layout(rect=[0, 0.05, 1, 1], h_pad=1.0, w_pad=1.2)

        # Save
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/FMUgen_paper_{sim_name}.pdf", format='pdf')

    # endregion

    # region gscgrid
    def plot_paper_gscgrid(self, sim_name: str, WT: WindTurbine, ps: PowerSystemModel):
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import numpy as np
        import os

        mpl.rcParams['font.family'] = 'Times New Roman'

        fig, axs = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
        lw = 2.0

        gsc_nominal_power_mw = WT.gsc_sn
        tops_nominal_power_mw = ps.s_n

        # (1,1) GSC active power
        axs[0, 0].plot(self.results['Time'], np.array(self.results[f'{WT.name}_GSC_p_e']) * gsc_nominal_power_mw, label='GSC Active Power', linewidth=lw)
        axs[0, 0].plot(self.results['Time'], np.array(self.results[f'{WT.name}_GSC_p_ref_grid']) * gsc_nominal_power_mw, label='GSC Active Power Reference', linestyle='--', linewidth=lw)
        axs[0, 0].set_ylabel('Active Power (MW)')
        axs[0, 0].legend(fontsize=8)
        axs[0, 0].grid(True, linestyle=':')
        axs[0, 0].set_title('GSC Active Power', fontsize=10)

        # (1,2) GSC bus voltage
        axs[0, 1].plot(self.results['Time'], self.results[f'{WT.name}_GSC_bus_voltage'], label='Bus Voltage', linewidth=lw)
        axs[0, 1].plot(self.results['Time'], self.results[f'{WT.name}_GSC_v_ref_grid'], label='Voltage Reference', linestyle='--', linewidth=lw)
        axs[0, 1].set_ylabel('Voltage (pu)')
        axs[0, 1].legend(fontsize=8)
        axs[0, 1].grid(True, linestyle=':')
        axs[0, 1].set_title('GSC Bus Voltage', fontsize=10)

        # (2,1) Generator p_e injections
        axs[1, 0].plot(self.results['Time'], np.array(self.results["Generators_p_e"]) * tops_nominal_power_mw, 
                    label='Generators Active Power', linewidth=lw)
        axs[1, 0].set_ylabel('Active Power (MW)')
        axs[1, 0].legend(fontsize=8)
        axs[1, 0].grid(True, linestyle=':')
        axs[1, 0].set_title('Generators Active Power Injections', fontsize=10)

        # (2,2) Reactive power
        axs[1, 1].plot(self.results['Time'], np.array(self.results[f'{WT.name}_GSC_q_e']) * gsc_nominal_power_mw, 
                    label='GSC Reactive Power', linewidth=lw)
        axs[1, 1].set_ylabel('Reactive Power (MW)')
        axs[1, 1].legend(fontsize=8)
        axs[1, 1].grid(True, linestyle=':')
        axs[1, 1].set_title('GSC Reactive Power', fontsize=10)

        fig.text(0.5, 0.02, 'Time (s)', ha='center', fontsize=12)

        for ax in axs.flat:
            ax.tick_params(labelsize=8)

        plt.tight_layout(rect=[0, 0.05, 1, 1], h_pad=1.0, w_pad=1.2)

        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/GSCgrid_paper_{sim_name}.pdf", format='pdf')


    # endregion

    # region dclink

    def plot_paper_dclink(self, sim_name: str, WT: WindTurbine):
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import os

        mpl.rcParams['font.family'] = 'Times New Roman'

        fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        lw = 2.0

        # (1,1) DC-link voltage and reference
        axs[0].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_vdc'], label='vdc', linewidth=lw)
        axs[0].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_vdc_ref'], label='vdc_ref', linestyle='--', linewidth=lw)
        axs[0].set_ylabel('Voltage (pu)')
        axs[0].legend(fontsize=8)
        axs[0].grid(True, linestyle=':')
        axs[0].set_title('DC-link Voltage vs Reference', fontsize=10)

        # (2,1) Power adjustment
        axs[1].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_p_adjust'], label='p_adjust', linewidth=lw)
        axs[1].set_ylabel('Power (pu)')
        axs[1].legend(fontsize=8)
        axs[1].grid(True, linestyle=':')
        axs[1].set_title('DC Voltage Controller Power Adjustment', fontsize=10)

        axs[1].set_xlabel('Time (s)', fontsize=12)
        for ax in axs.flat:
            ax.tick_params(labelsize=8)

        plt.tight_layout(rect=[0, 0.03, 1, 1], h_pad=1.0)

        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/DClink_paper_{sim_name}.pdf", format='pdf')


    # endregion
    def plot_paper_dclink_2x2(self, sim_name: str, WT: WindTurbine):
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import os

        # Use Times New Roman globally
        mpl.rcParams['font.family'] = 'Times New Roman'

        fig, axs = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
        lw = 2.0  # Line width

        # (1,1) DC-link voltage and reference
        axs[0, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_vdc'], label='vdc', linewidth=lw)
        axs[0, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_vdc_ref'], label='vdc_ref', linestyle='--', linewidth=lw)
        axs[0, 0].set_ylabel('Voltage (pu)')
        axs[0, 0].legend(fontsize=8)
        axs[0, 0].grid(True, linestyle=':')
        axs[0, 0].set_title('DC-link Voltage vs Reference', fontsize=10)

        # (1,2) Power adjustment
        axs[0, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_p_adjust'], label='p_adjust', linewidth=lw)
        axs[0, 1].set_ylabel('Power (pu)')
        axs[0, 1].legend(fontsize=8)
        axs[0, 1].grid(True, linestyle=':')
        axs[0, 1].set_title('DC Voltage Controller Power Adjustment', fontsize=10)

        # (2,1) Chopper duty cycle
        axs[1, 0].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_duty'], label='Duty Cycle', linewidth=lw)
        axs[1, 0].set_ylabel('Duty Cycle (pu)')
        axs[1, 0].legend(fontsize=8)
        axs[1, 0].grid(True, linestyle=':')
        axs[1, 0].set_title('Chopper Duty Cycle', fontsize=10)

        # (2,2) Chopper current
        axs[1, 1].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_i_chopper'], label='Chopper Current', linewidth=lw)
        axs[1, 1].set_ylabel('Current (pu)')
        axs[1, 1].legend(fontsize=8)
        axs[1, 1].grid(True, linestyle=':')
        axs[1, 1].set_title('Chopper Current', fontsize=10)

        # Common X label
        fig.text(0.5, 0.02, 'Time (s)', ha='center', fontsize=12)

        # Ticks and layout
        for ax in axs.flat:
            ax.tick_params(labelsize=8)

        plt.tight_layout(rect=[0, 0.05, 1, 1], h_pad=1.0, w_pad=1.2)

        # Save
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/DClink_2x2_paper_{sim_name}.pdf", format='pdf')


    def plot_paper_dclink_3x1(self, sim_name: str, WT: WindTurbine):
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import os
        # Use Times New Roman globally
        mpl.rcParams['font.family'] = 'Times New Roman'

        fig, axs = plt.subplots(3, 1, figsize=(5, 9), sharex=True)
        lw = 2.0  # Line width

        # (1,1) DC-link voltage and reference
        axs[0].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_vdc'], label='vdc', linewidth=lw)
        axs[0].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_vdc_ref'], label='vdc_ref', linestyle='--', linewidth=lw)
        axs[0].set_ylabel('Voltage (pu)')
        axs[0].legend(fontsize=8)
        axs[0].grid(True, linestyle=':')
        axs[0].set_title('DC-link Voltage vs Reference', fontsize=10)

        # (2,1) Power adjustment
        axs[1].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_p_adjust'], label='p_adjust', linewidth=lw)
        axs[1].set_ylabel('Power (pu)')
        axs[1].legend(fontsize=8)
        axs[1].grid(True, linestyle=':')
        axs[1].set_title('DC Voltage Controller Power Adjustment', fontsize=10)

        # (3,1) Chopper current
        axs[2].plot(self.results['Time elec'], self.results[f'{WT.name}_DC_i_chopper'], label='Chopper Current', linewidth=lw)
        axs[2].set_ylabel('Current (pu)')
        axs[2].legend(fontsize=8)
        axs[2].grid(True, linestyle=':')
        axs[2].set_title('Chopper Current', fontsize=10)

        # Common X label
        axs[2].set_xlabel('Time (s)', fontsize=12)

        # Ticks and layout
        for ax in axs.flat:
            ax.tick_params(labelsize=8)

        plt.tight_layout(rect=[0, 0.03, 1, 1], h_pad=1.0)

        # Save
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/DClink_3x1_paper_{sim_name}.pdf", format='pdf')


    # region OpenFAST vs Generator Torque
    def plot_openfast_vs_generator_torque(self, sim_name: str, WT: WindTurbine):
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import os

        # Use Times New Roman globally
        mpl.rcParams['font.family'] = 'Times New Roman'

        fig, ax = plt.subplots(figsize=(10, 6))
        lw = 2.0  # Line width

        # Plot OpenFAST turbine torque and generator torque
        ax.plot(self.results['Time'], self.results[f'{WT.name}_GenSpdOrTrq'], label='Generator Torque (kNm)', linewidth=lw)
        ax.plot(self.results['Time'], self.results[f'{WT.name}_HSShftTq'], label='OpenFAST Turbine Torque (kNm)', linestyle='--', linewidth=lw)

        # Set axis labels with larger font size
        ax.set_ylabel('Torque (kNm)', fontsize=14)
        ax.set_xlabel('Time (s)', fontsize=14)

        # Set title with larger font size
        ax.set_title('OpenFAST Turbine Torque vs Generator Torque', fontsize=16)

        # Increase tick label size
        ax.tick_params(axis='both', which='major', labelsize=12)

        # Enable both major and minor grids
        ax.grid(True, which='both', linestyle=':', linewidth=0.8)

        # Add a clearer legend
        ax.legend(fontsize=12, frameon=True, loc='upper left')

        # Adjust layout for publication quality
        plt.tight_layout()

        # Create directory if it does not exist
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/OpenFAST_vs_Generator_Torque_{sim_name}.pdf", format='pdf')
    # endregion



        # region Multi-rate visualization
        # region Multi-rate visualization
    def plot_multirate_torque(self, sim_name: str, WT: WindTurbine, step_size_mech: float, step_size_elec: float):
        import matplotlib.pyplot as plt
        import os
        import numpy as np
        import matplotlib as mpl

        # Load data
        time_elec = np.array(self.results["Time elec"])
        te = np.array(self.results[f"{WT.name}_PMSM T_e"])
        te_ref = np.array(self.results[f"{WT.name}_PMSM T_e_ref"])

        # Reconstruct mechanical-step reference time
        mech_interval = int(step_size_mech / step_size_elec)
        time_mech = time_elec[::mech_interval]
        te_ref_mech = te_ref[::mech_interval]

        # Force the reference to lead: shift one mech step earlier in time
        time_mech_leading = time_mech - step_size_mech
        time_mech_leading[0] = time_mech[0]  # avoid negative time at start

        # Plot
        plt.figure(figsize=(10, 5))
        plt.plot(time_elec, te, label="Tₑ (electrical step)", linewidth=2)
        plt.step(time_mech_leading, te_ref_mech, where="post", label="Tₑ,ref (mechanical step)", linestyle='--', linewidth=2)
        plt.xlabel("Time (s)")
        plt.ylabel("Torque (kNm)")
        plt.title("Multi-rate Simulation: PMSM Torque vs Reference")
        plt.grid(True, linestyle=":")
        plt.legend()
        plt.tight_layout()

        # Set x-axis limit to cut off the first 10% and the last 20% of time
        start_time = time_elec[0] + (time_elec[-1] - time_elec[0]) * 0.1
        cutoff_time = time_elec[-1] * 0.8
        plt.xlim(left=start_time, right=cutoff_time)

        # Save
        output_dir = f"Figures/co_sim/{sim_name}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/Multirate_Torque_{sim_name}.pdf", format='pdf')
        plt.close()

    # endregion

    # endregion
