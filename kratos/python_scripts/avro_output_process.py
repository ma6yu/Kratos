import os
import json

try:
    import avro.schema
    from avro.datafile import DataFileReader, DataFileWriter
    from avro.io import DatumReader, DatumWriter
except:
    print("Unable to import AVRO. Please install avro drivers for python.")

import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils

from  KratosMultiphysics.deprecation_management import DeprecationManager

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AvroOutputProcess(model, settings["Parameters"])


class AvroOutputProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        model_part_name = settings["model_part_name"].GetString()
        schema_name = settings["schema_name"].GetString()
        
        self.model_part = model[model_part_name]
        self.schema = None
        self.codec = 'null'
        self.settings = settings

        force_generate_schema = False

        # Try to open a existing schema, generate it otherwise
        # if not os.path.exists(schema_name) or force_generate_schema:
        with open(schema_name, "w+") as schema_file:
            schema = self.GenerateSchemaFromProjectParameters()
            # schema_file.write(schema)
            json.dump(schema, schema_file)
        
        with open(schema_name, "rb") as schema_file:
            self.schema = avro.schema.parse(schema_file.read())
        
        # Replace deprecations
        self.TranslateLegacyVariablesAccordingToCurrentStandard(settings, {})

        with DataFileWriter(open(self.settings["output_name"].GetString(), "wb"), DatumWriter(), self.schema, codec = self.codec) as writer:
            pass

        unsued = {
            # In the future, this will load the avro drvier from CPP to improve efficiency
            # self.avro_io = KratosMultiphysics.AvroOutput(self.model_part, settings)

            # Charlie: ???
            # if settings["save_output_files_in_folder"].GetBool():
            #     if self.model_part.GetCommunicator().MyPID() == 0:
            #         folder_name = settings["folder_name"].GetString()
            #         if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            #             kratos_utils.DeleteDirectoryIfExisting(folder_name)
            #     self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

            # self.output_interval = settings["output_interval"].GetDouble()
            # self.output_control = settings["output_control_type"].GetString()
            # self.next_output = 0.0

            # Charlie: ???
            # self.__ScheduleNextOutput() # required here esp for restart
        }

    def GenerateSchemaFromProjectParameters(self):
        schema = {
            'namespace' : 'kratos',
            'type': 'record',
            'name': 'simulation_results',
            'fields': [
                {
                    'name':'time_step', 'type': 'double'
                },
                {
                    'name':'results',
                    'type':{
                        'name':'nodal',
                        'type':'record',
                        'fields': [{
                            'name':nodal_solution_step_value.GetString()+"_ARRAY",
                            'type':{
                                'type':'array',
                                'items':'double'
                            }
                        } for nodal_solution_step_value in self.settings['postprocess_parameters']['result_file_configuration']['nodal_results']]
                    }
                }
            ]
        }

        return schema

        unsued = {
            # Charlie: May be usefull in the future.
            # # Time step of the record
            # schema['fields'].append({
            #     'name':'solution_step_values',
            #     'type':{
            #         'type':'array',
            #         'items':{}
            #     }
            # })

            # # Nodal solution step values
            # nodal_solution_step_data_variables = {
            #     'name': 'nodal_solution_step_data_variables',
            #     'type' : {}
            # }
            # for nodal_solution_step_value in self.settings['nodal_solution_step_data_variables']:
            #     nodal_solution_step_data_variables.append("name":nodal_solution_step_value.GetString(), "type":{
            #         "type": "array", "items": "double"
            #     })

            # # Nodal data values
            # for nodal_data_value in self.settings['nodal_data_value_variables']:
            #     schema[fields].append("name":nodal_data_value.GetString(), "type":{
            #         "type": "array", "items": "double"
            #     })
        }

    def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings, deprecations):
        # Defining a string to help the user understand where the warnings come from (in case any is thrown)
        context_string = type(self).__name__

        for dep in deprecations:
            if DeprecationManager.HasDeprecatedVariable(context_string, settings, dep['old_name'], dep['new_name']):
                DeprecationManager.ReplaceDeprecatedVariableName(settings, dep['old_name'], dep['new_name'])

    def PrintOutput(self):
        with DataFileWriter(open(self.settings["output_name"].GetString(), "ab+"), DatumWriter()) as writer:
            object_results = {
                "time_step": float(self.model_part.ProcessInfo[KratosMultiphysics.STEP]), 
                "results": {
                    nodal_solution_step_value.GetString()+"_ARRAY":[node.GetSolutionStepValue(KratosMultiphysics.KratosGlobals.GetVariable(nodal_solution_step_value.GetString())) for node in self.model_part.Nodes] for nodal_solution_step_value in self.settings['postprocess_parameters']['result_file_configuration']['nodal_results']
                }
            }

            writer.append(object_results)
    
    def ReadPrint(self):
        with DataFileReader(open(self.settings["output_name"].GetString(), "rb"), DatumReader()) as reader:
            for result in reader:
                print(result)


    def IsOutputStep(self):
        if self.output_control == "time":
            return self.__GetTime() >= self.next_output
        else:
            return self.model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.next_output

    def __ScheduleNextOutput(self):
        if self.output_interval > 0.0: # Note: if == 0, we'll just always print
            if self.output_control == "time":
                while self.next_output <= self.__GetTime():
                    self.next_output += self.output_interval
            else:
                while self.next_output <= self.model_part.ProcessInfo[KratosMultiphysics.STEP]:
                    self.next_output += self.output_interval

    def __GetTime(self):
        # remove rounding errors that mess with the comparison
        # e.g. 1.99999999999999999 => 2.0
        return float("{0:.12g}".format(self.model_part.ProcessInfo[KratosMultiphysics.TIME]))
