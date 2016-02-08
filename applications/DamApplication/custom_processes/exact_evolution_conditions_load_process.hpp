//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:     January 2016 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_EXACT_EVOLUTION_CONDITIONS_LOAD_PROCESS )
#define  KRATOS_EXACT_EVOLUTION_CONDITIONS_LOAD_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "dam_application.h"

namespace Kratos
{

class ExactEvolutionConditionsLoadProcess : public Process
{
    
public:

    typedef double argument_type; // To be STL conformance.
    typedef double result_type; // To be STL conformance.

    typedef boost::array<result_type, 1>  result_row_type;

    typedef std::pair<argument_type, result_row_type> RecordType;

    typedef std::vector<RecordType> TableContainerType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    ExactEvolutionConditionsLoadProcess(ModelPart& r_model_part, double time_unit_converter) : mr_model_part(r_model_part)
    {
        mlast_id = 0;
        mtime_unit_converter = time_unit_converter;
    }

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~ExactEvolutionConditionsLoadProcess(){}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        double time = mr_model_part.GetProcessInfo()[TIME];
        time = time/mtime_unit_converter;
        double delta_time = mr_model_part.GetProcessInfo()[DELTA_TIME];
        delta_time = delta_time/mtime_unit_converter;
                
        const TableContainerType& table_water_level = mr_model_part.pGetTable(2)->Data();    // Table with information about water level in each time
                
        if( (time + delta_time*1.0e-10) >= table_water_level[mlast_id+1].first )
        {
            mlast_id = mlast_id + 1;
           
            this->Hydrostatic_pressure(table_water_level); 
        }
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;
    std::size_t mlast_id;
    double mtime_unit_converter;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Hydrostatic_pressure(TableContainerType const& table_water_level)
    {
        // We have to pick up the provided values by the mesh
        std::string direction = (*(mr_model_part.pGetMesh(3)))[GRAVITY_DIRECTION];
        const double& coordinate_base = (*(mr_model_part.pGetMesh(3)))[COORDINATE_BASE_DAM];
        const double& spe_weight =  (*(mr_model_part.pGetMesh(3)))[SPECIFIC_WEIGHT];
        
        double ref_coord;
        
        const double& water_level = table_water_level[mlast_id].second[0];                          // Water Level
                
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(3)->NodesBegin(); i <mr_model_part.pGetMesh(3)->NodesEnd(); i++)
        {
            double& pressure  = (i)->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
            
            if(direction=="X")
            {
                ref_coord = coordinate_base + water_level;
                pressure = spe_weight*(ref_coord- (i)->X());
                    if(pressure < 0.0)
                        pressure = 0.0; 
            }
            else if(direction=="Y")
            {
                ref_coord = coordinate_base + water_level;
                pressure = spe_weight*(ref_coord- (i)->Y());
                    if(pressure < 0.0)
                        pressure = 0.0; 
            }
            else if(direction=="Z")
            {
                ref_coord = coordinate_base + water_level;
                pressure = spe_weight*(ref_coord- (i)->Z());
                    if(pressure < 0.0)
                        pressure = 0.0; 
            }
        }
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_EXACT_EVOLUTION_CONDITIONS_LOAD_PROCESS defined */
