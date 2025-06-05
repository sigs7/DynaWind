from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave
from dynawind.dynawind_models.pmsm import PMSM

class FAST:
    def __init__(self, params: dict):
        self.params = params
        self.fmu, self.vrs = self.initiate_FAST_FMU(params.get("fast_fmu_filename", "fast.fmu"), params.get("start_time", 0.0), params.get("mode", 3))

    def initiate_FAST_FMU(self, fmu_filename: str, start_time: float, mode: int) -> FMU2Slave:
        """
        Initializes the FMU for co-simulation.

        Args:
            fmu_filename (str): Path to the FMU file.
            start_time (float): Simulation start time.
            mode (int): FMU operational mode.

        Returns:
            FMU2Slave: Initialized FMU instance.
            dict: Variable name to value reference mapping.
        """



        model_description = read_model_description(fmu_filename, validate=False)

        # # collect the value references
        vrs = {}
        for variable in model_description.modelVariables:
            vrs[variable.name] = variable.valueReference

        print("Value References: \n")
        for name, vr in vrs.items():
            print(f"Variable: {name}, Value Reference: {vr}")

        unzipdir = extract(fmu_filename)        # , unzipdir=f'Turbine_{nr}'

        wd_file_path = 'openfast_fmu/resources/wd.txt'
        new_directory = 'C:/Users/larsi/Master/TOPS_LAH/TOPS_LAH'

        with open(wd_file_path, 'w') as f:
            f.write(new_directory)

        fmu = FMU2Slave(guid=model_description.guid,
                         unzipDirectory=unzipdir,
                           modelIdentifier=model_description.coSimulation.modelIdentifier,
                             instanceName='instance1')
        
        fmu.instantiate()
        fmu.setReal([vrs['testNr']], [1002])
        fmu.setReal([vrs['Mode']], [mode])
        fmu.setupExperiment(startTime=start_time)
        fmu.enterInitializationMode()
        fmu.exitInitializationMode()
        
        return fmu, vrs
    
    # Connects the FAST FMU to the wind turbine model
    def connect_fmu(self, wind_turbine):
        self.wind_turbine = wind_turbine

    # Stepping the OpenFAST FMU
    def step_fmu(self, pmsm, ps, v, x, time: float, step_size: float):

        from dynawind.dynawind_models.pmsm import PMSM  # Lazy import to avoid circular import

        self.fmu.setReal([self.vrs['GenSpdOrTrq']], [-pmsm.T_e()])
        # self.fmu.setReal([self.vrs['GenPwr']], [pmsm.P_e()])
        if self.wind_turbine.gsc_control == "PV":
            p_e_gsc = ps.vsc["GridSideConverter_PV"].p_e(x,v)[self.wind_turbine.index] * ps.vsc["GridSideConverter_PV"].par["S_n"][self.wind_turbine.index] * 1e3
            self.fmu.setReal([self.vrs['GenPwr']], [p_e_gsc])       # Power in kW
        if self.wind_turbine.gsc_control == "PQ":
            p_e_gsc = ps.vsc["GridSideConverter_PQ"].p_e(x,v)[self.wind_turbine.index] * ps.vsc["GridSideConverter_PQ"].par["S_n"][self.wind_turbine.index] * 1e3
            self.fmu.setReal([self.vrs['GenPwr']], [p_e_gsc])       # Power in kW
        

        # Set Torque reference to pmsm in mechanical timestep
        pmsm.torque_ref_raw = -1 * self.fmu.getReal([self.vrs['GenTq']])[0]/pmsm.params["T_r"]       # Torque reference from FAST FMU in pu

        # Limit torque reference
        pmsm.torque_ref_raw = pmsm.clamp_value(pmsm.torque_ref_raw, max_value=0, min_value=-1.5)

        # Update mechanical states to the OMSM from FAST FMU
        pmsm.speed = self.fmu.getReal([self.vrs['GenSpeed']])[0] / pmsm.params["rpm_n"]
        pmsm.SPEED = pmsm.speed * pmsm.params["rpm_n"]        # Mulig drit kodepraksis menmen     (rpm)


        # # # Maximum power
        # self.fmu.setReal([self.vrs['ElecPwrCom']], [20e3])        # Setting a high reference translates to delivering maximum power avaiable

        # Derating 
        if time < 120:
            self.fmu.setReal([self.vrs['ElecPwrCom']], [20e3])      # Setting a high reference translates to delivering maximum power avaiable
        elif 120.01 < time < 240:
            self.fmu.setReal([self.vrs['ElecPwrCom']], [10e3])      # Derating to 10 MW
        else:
            self.fmu.setReal([self.vrs['ElecPwrCom']], [20e3])      

        self.fmu.doStep(currentCommunicationPoint=time, communicationStepSize=step_size)

        
    # FMU termination
    def terminate_fmu(self):
        self.fmu.terminate()
        self.fmu.freeInstance()






    # region Old code
    # def step_fmu(self, pmsm, time: float, step_size: float):

    #     from tops.cosim_models.pmsm import PMSM  # Lazy import to avoid circular import

    #     self.fmu.setReal([self.vrs['GenSpdOrTrq']], [-pmsm.get_T_e()])
    #     self.fmu.setReal([self.vrs['GenPwr']], [pmsm.get_P_e()])

    #     # Maximum power
    #     self.fmu.setReal([self.vrs['ElecPwrCom']], [20e3])

    #     # # Derating 
    #     # if time < 15:
    #     #     self.fmu.setReal([self.vrs['ElecPwrCom']], [20e3])      # Setting a high reference translates to delivering maximum power avaiable
    #     # elif 35 < time < 55:
    #     #     self.fmu.setReal([self.vrs['ElecPwrCom']], [10e3])      # Derating to 10 MW
    #     # else:
    #     #     self.fmu.setReal([self.vrs['ElecPwrCom']], [20e3])      

    #     self.fmu.doStep(currentCommunicationPoint=time, communicationStepSize=step_size)