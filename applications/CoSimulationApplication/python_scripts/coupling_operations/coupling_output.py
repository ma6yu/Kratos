# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

def Create(*args):
    return CouplingOutput(*args)

class CouplingOutput(CoSimulationCouplingOperation):
    """This operation is used to output at different points in the coupling.
    TODO:
    - add support for json, hdf5, gid
    - add tests
    - more cleanup
    """
    def __init__(self, settings, solver_wrappers, process_info):
        super().__init__(settings, process_info)
        self.model = solver_wrappers[self.settings["solver"].GetString()].model
        self.execution_point = self.settings["execution_point"].GetString()
        model_part_name = self.settings["output_parameters"]["model_part_name"].GetString()
        self.base_output_file_name = "{}_{}_{}_".format(self.settings["solver"].GetString(), model_part_name, self.execution_point)

        available_execution_points = [
            "initialize_solution_step",
            "finalize_solution_step",
            "initialize_coupling_iteration",
            "finalize_coupling_iteration"
        ]

        if self.execution_point not in available_execution_points:
            err_msg  = 'Execution point "{}" is not available, only the following options are available:\n    '.format(self.execution_point)
            err_msg += "\n    ".join(available_execution_points)
            raise Exception(err_msg)

        self.step = 0 # this should come from self.process_info
        self.output_control = settings["output_parameters"]["output_control_type"].GetString()
        self.output_interval = settings["output_parameters"]["output_interval"].GetDouble()
        if self.output_control != 'step':
            raise Exception('Only step output interval currently implemented for CoSim coupling output')

        self.next_output = 0.0
        # TODO check if restarted. If not delete the folder => check self.process_info
        self.output = KM.VtkOutput(self.model[model_part_name], self.settings["output_parameters"]) # currently hardcoded to vtk

    def InitializeSolutionStep(self):
        self.step += 1
        self.coupling_iteration = 0

        if self.execution_point == "initialize_solution_step":
            output_file_name = self.base_output_file_name + str(self.step)
            self.PrintOutput(output_file_name)

    def FinalizeSolutionStep(self):
        if self.execution_point == "finalize_solution_step":
            output_file_name = self.base_output_file_name + str(self.step)
            self.PrintOutput(output_file_name)

    def InitializeCouplingIteration(self):
        self.coupling_iteration += 1

        if self.execution_point == "initialize_coupling_iteration":
            output_file_name = self.base_output_file_name + "{}_{}".format(self.step, self.coupling_iteration)
            self.PrintOutput(output_file_name)

    def FinalizeCouplingIteration(self):
        if self.execution_point == "finalize_coupling_iteration":
            output_file_name = self.base_output_file_name + "{}_{}".format(self.step, self.coupling_iteration)
            self.PrintOutput(output_file_name)

    def PrintOutput(self,output_file_name):
        if self.IsOutputStep():
            self.output.PrintOutput(output_file_name)

            # schedule next output
            if self.output_interval > 0.0: # Note: if == 0, we'll just always print
                while self.next_output <= self.step:
                    self.next_output += self.output_interval

    def IsOutputStep(self):
        return ( self.step >= self.next_output )


    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"            : "UNSPECIFIED",
            "execution_point"   : "UNSPECIFIED",
            "output_format"     : "vtk",
            "output_parameters" : {
                "model_part_name" : "UNSPECIFIED",
                "output_control_type" : "step",
                "output_interval" : 0}
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
