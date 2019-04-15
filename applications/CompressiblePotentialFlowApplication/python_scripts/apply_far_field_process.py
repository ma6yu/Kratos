import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

def DotProduct(A,B):
    result = 0
    for i,j in zip(A,B):
        result += i*j
    return result

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFarFieldProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyFarFieldProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "inlet_phi": 1.0,
                "velocity_infinity": [1.0,0.0,0]
            }  """ )


        settings.ValidateAndAssignDefaults(default_parameters);

        self.far_field_model_part = Model[settings["model_part_name"].GetString()]
        self.velocity_infinity = KratosMultiphysics.Vector(3)#array('d', [1.0, 2.0, 3.14])#np.array([0,0,0])#np.zeros(3)#vector(3)
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        #self.density_infinity = settings["density_infinity"].GetDouble() #TODO: must read this from the properties
        self.inlet_phi_0 = settings["inlet_phi"].GetDouble()
        self.far_field_model_part.ProcessInfo.SetValue(CPFApp.VELOCITY_INFINITY,self.velocity_infinity)



    def Execute(self):
        #KratosMultiphysics.VariableUtils().SetVectorVar(CPFApp.VELOCITY_INFINITY, self.velocity_infinity, self.far_field_model_part.Conditions)
        for cond in self.far_field_model_part.Conditions:
            cond.SetValue(CPFApp.VELOCITY_INFINITY, self.velocity_infinity)

        # Select the first node in the mesh as reference
        for node in self.far_field_model_part.Nodes:
            x0 = node.X
            y0 = node.Y
            z0 = node.Z
            break

        # Find smallest distance_to_reference
        pos = 1e30
        for node in self.far_field_model_part.Nodes:
            # Computing distance to reference
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0

            distance_to_reference = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]

            if(distance_to_reference < pos):
                pos = distance_to_reference
                self.reference_inlet_node = node

        # Fix nodes in the inlet
        for cond in self.far_field_model_part.Conditions:
            normal = cond.GetValue(KratosMultiphysics.NORMAL)
            projection = DotProduct(normal,self.velocity_infinity)

            if( projection < 0):
                for node in cond.GetNodes():
                    # Computing distance to reference
                    dx = node.X - self.reference_inlet_node.X
                    dy = node.Y - self.reference_inlet_node.Y
                    dz = node.Z - self.reference_inlet_node.Z

                    inlet_phi = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]
                    node.Fix(CPFApp.VELOCITY_POTENTIAL)
                    node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,0,inlet_phi + self.inlet_phi_0)
                    if self.far_field_model_part.HasNodalSolutionStepVariable(CPFApp.ADJOINT_VELOCITY_POTENTIAL):
                        node.Fix(CPFApp.ADJOINT_VELOCITY_POTENTIAL)
                        node.SetSolutionStepValue(CPFApp.ADJOINT_VELOCITY_POTENTIAL,0,inlet_phi)

    def ExecuteInitializeSolutionStep(self):
        self.Execute()

