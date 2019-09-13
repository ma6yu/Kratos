//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (based on vtk_output.h)
//
//

// System includes

// External includes
#include "inc/vtu11.hpp"

// Project includes
#include "vtu_output.h"

namespace Kratos
{
VtuOutput::VtuOutput(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : mrModelPart(rModelPart),
        mOutputSettings(ThisParameters)
{
    // The default parameters
    Parameters default_parameters = GetDefaultParameters();
    mOutputSettings.ValidateAndAssignDefaults(default_parameters);

    // Initialize other variables
    const std::string file_format = mOutputSettings["file_format"].GetString();
    if (file_format == "ascii") {
        mFileFormat = VtuOutput::FileFormat::VTU_ASCII;
    } else if (file_format == "binary_raw") {
        mFileFormat = VtuOutput::FileFormat::VTU_BINARY_RAW;
    } else if (file_format == "binary_raw_compressed") {
        mFileFormat = VtuOutput::FileFormat::VTU_BINARY_RAW_COMPRESSED;
    } else if (file_format == "binary_base64") {
        mFileFormat = VtuOutput::FileFormat::VTU_BINARY_BASE64;
    } else if (file_format == "binary_base64_appended") {
        mFileFormat = VtuOutput::FileFormat::VTU_BINARY_BASE64_APPENDED;
    } else {
        KRATOS_ERROR << "Option for \"file_format\": " << file_format << " not recognised!\n Possible output formats options are: \"ascii\", \"binary_raw\", \"binary_raw_compressed\", \"binary_base64\", \"binary_base64_appended\"" << std::endl;
    }

    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();

    const int num_elements = r_data_comm.SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));
    const int num_conditions = r_data_comm.SumAll(static_cast<int>(r_local_mesh.NumberOfConditions()));

    KRATOS_WARNING_IF("VtuOutput", num_elements > 0 && num_conditions > 0) << r_data_comm << "Modelpart \"" << rModelPart.Name() << "\" has both elements and conditions.\nGiving precedence to elements and writing only elements!" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::PrintOutput()
{
    std::vector<double> points
    {
        0.0, 0.0, 0.5,    0.0, 0.3, 0.5,    0.0, 0.7, 0.5,    0.0, 1.0, 0.5, // 0,  1,  2,  3
        0.5, 0.0, 0.5,    0.5, 0.3, 0.5,    0.5, 0.7, 0.5,    0.5, 1.0, 0.5, // 4,  5,  6,  7
        1.0, 0.0, 0.5,    1.0, 0.3, 0.5,    1.0, 0.7, 0.5,    1.0, 1.0, 0.5  // 8,  9, 10, 11
    };

    std::vector<std::size_t> connectivity
    {
        0,  4,  5,  1, // 0
        1,  5,  6,  2, // 1
        2,  6,  7,  3, // 2
        4,  8,  9,  5, // 3
        5,  9, 10,  6, // 4
        6, 10, 11,  7  // 5
    };

    std::vector<std::size_t> offsets { 4, 8, 12, 16, 20, 24 };
    std::vector<vtu11::VtkCellType> types { 9, 9, 9, 9, 9, 9 };

    vtu11::Vtu11UnstructuredMesh vtu_mesh{ points, connectivity, offsets, types };
}

/***********************************************************************************/
/***********************************************************************************/

std::string VtuOutput::GetOutputFileName(const ModelPart& rModelPart) const
{
    const int rank = rModelPart.GetCommunicator().MyPID();
    std::string model_part_name(rModelPart.Name());

    std::string label;
    std::stringstream ss;
    const std::string output_control = mOutputSettings["output_control_type"].GetString();
    if (output_control == "step") {
        ss << std::fixed << std::setfill('0')
           << rModelPart.GetProcessInfo()[STEP];
        label = ss.str();
    } else if(output_control == "time") {
        ss << std::fixed << std::setfill('0')
           << rModelPart.GetProcessInfo()[TIME];
        label = ss.str();
    } else {
        KRATOS_ERROR << "Option for \"output_control_type\": " << output_control
            <<" not recognised!\nPossible output_control_type options "
            << "are: \"step\", \"time\"" << std::endl;
    }

    // Putting everything together
    std::string output_file_name;
    if (mOutputSettings["save_output_files_in_folder"].GetBool()) {
        output_file_name += mOutputSettings["folder_name"].GetString() + "/";
    }
    const std::string& custom_name_prefix = mOutputSettings["custom_name_prefix"].GetString();
    output_file_name += custom_name_prefix + model_part_name + "_" + std::to_string(rank) + "_" + label + ".vtk";

    return output_file_name;
}

Parameters VtuOutput::GetDefaultParameters()
{
    // IMPORTANT: when "output_control_type" is "time", then paraview will not be able to group them
    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"                    : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "file_format"                        : "ascii",
        "output_precision"                   : 7,
        "output_control_type"                : "step",
        "output_frequency"                   : 1.0,
        "output_sub_model_parts"             : false,
        "folder_name"                        : "VTK_Output",
        "custom_name_prefix"                 : "",
        "save_output_files_in_folder"        : true,
        "write_deformed_configuration"       : false,
        "write_properties_id"                : false,
        "nodal_solution_step_data_variables" : [],
        "nodal_data_value_variables"         : [],
        "nodal_flags"                        : [],
        "element_data_value_variables"       : [],
        "element_flags"                      : [],
        "condition_data_value_variables"     : [],
        "condition_flags"                    : [],
        "gauss_point_variables"              : []
    })" );

    return default_parameters;
}

} // namespace Kratos
