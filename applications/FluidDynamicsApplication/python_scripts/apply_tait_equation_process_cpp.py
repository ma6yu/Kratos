import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTaitEquationProcess(Model, settings["Parameters"])

class ApplyTaitEquationProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        # Check the default values
        default_settings = KratosMultiphysics.Parameters("""
            {
                "help"                     : "This proces sets a given scalar value for the density in all the nodes of FluidModelPart",
                "interval"                 : [0.0,1e+30],
                "model_part_name"          : "",
                "rho_0"                    : 0,                  
                "p_0"                      : 0,                    
                "theta"                    : 0,
                "k_0"                      : 0 
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        # getting the ModelPart from the Model
        self.model_part = Model[settings["model_part_name"].GetString()]
        
        # Set the Tait Equation process
        self.TaitEquationProcess = KratosFluid.TaitEquationProcess(self.model_part, settings)

    def ExecuteInitializeSolutionStep(self):
        self.TaitEquationProcess.ExecuteInitializeSolutionStep()
    