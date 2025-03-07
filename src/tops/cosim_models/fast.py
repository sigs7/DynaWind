from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave
from tops.cosim_models.pmsm import PMSM

class FAST:
    def __init__(self, params: dict):
        self.params = params
        self.fmu, self.vrs = self.initiate_FAST_FMU(params.get("fast_fmu_filename", "fast.fmu"), params.get("start_time", 0.0), params.get("mode", 3))

    def initiate_FAST_FMU(self, fmu_filename: str, start_time: float, mode: int) -> FMU2Slave:

        model_description = read_model_description(fmu_filename, validate=False)

        # # collect the value references
        vrs = {}
        for variable in model_description.modelVariables:
            vrs[variable.name] = variable.valueReference

        print("Value References: \n")
        for name, vr in vrs.items():
            print(f"Variable: {name}, Value Reference: {vr}")

        unzipdir = extract(fmu_filename)

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

    def step_fmu(self, pmsm, time: float, step_size: float):

        from tops.cosim_models.pmsm import PMSM  # Lazy import to avoid circular import

        self.fmu.setReal([self.vrs['GenSpdOrTrq']], [-pmsm.get_T_e()])
        self.fmu.setReal([self.vrs['GenPwr']], [pmsm.get_P_e()])

        if time < 25:
            self.fmu.setReal([self.vrs['ElecPwrCom']], [20e3])      # Setting a high reference translates to delivering maximum power avaiable
        else:
            self.fmu.setReal([self.vrs['ElecPwrCom']], [10e3])      # Derating to 10 MW

        self.fmu.doStep(currentCommunicationPoint=time, communicationStepSize=step_size)

    def terminate_fmu(self):
        self.fmu.terminate()
        self.fmu.freeInstance()