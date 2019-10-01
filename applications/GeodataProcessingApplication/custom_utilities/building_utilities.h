//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Nicola Germano
//
//

#if !defined(KRATOS_BUILDING_UTILITIES_H_INCLUDED )
#define  KRATOS_BUILDING_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geodata_processing_application_variables.h"
#include "includes/checks.h"
#include "includes/model_part.h"


namespace Kratos
{
  ///@addtogroup GeodataProcessingApplication
  ///@{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// Auxiliary utility to maintain the quality of the model part

  class KRATOS_API(GEODATA_PROCESSING_APPLICATION) BuildingUtilities
  {

  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of BuildingUtilities
    KRATOS_CLASS_POINTER_DEFINITION(BuildingUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    BuildingUtilities( ModelPart& rModelPart ) : mrModelPart(rModelPart)
    { };

    /// Destructor.
    ~BuildingUtilities() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to check if there are overlapping elements
     *
     */
    void CheckOverlap();

    /**
     * @brief Function to split the buildings that are overlapping elements
     *
     */
    void SplitBuilding();

    /**
     * @brief Function to check if a node is in the map
     *
     */
    bool CheckNodeInMap(std::unordered_map<IndexType, IndexType> new_nodes_map, IndexType node_id);

    /**
     * @brief Function to find the maximum node id
     *
     */
    void FindMaxNodeId(IndexType new_node_id);


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    BuildingUtilities& operator=(BuildingUtilities const& rOther);

    /// Copy constructor.
    BuildingUtilities(BuildingUtilities const& rOther);

    ///@}

}; // Class BuildingUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const BuildingUtilities& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CLEANING_UTILITIES_H_INCLUDED  defined
