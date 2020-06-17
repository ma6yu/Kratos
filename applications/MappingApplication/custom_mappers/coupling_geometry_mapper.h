//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Tobias Teschemachen
//

#if !defined(KRATOS_COUPLING_GEOMETRY_MAPPER_H_INCLUDED)
#define  KRATOS_COUPLING_GEOMETRY_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper.h"
#include "custom_utilities/interface_vector_container.h"
#include "custom_utilities/mapper_flags.h"
#include "custom_utilities/mapper_local_system.h"


namespace Kratos
{
///@name Kratos Classes
///@{

class CouplingGeometryLocalSystem : public MapperLocalSystem
{
public:

    explicit CouplingGeometryLocalSystem(GeometryPointerType pGeom) : mpGeom(pGeom) {}

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;

    CoordinatesArrayType& Coordinates() const override
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpGeom) << "Members are not intitialized!" << std::endl;
        // return mpGeom->Center(); // check why not compiling...
        KRATOS_ERROR << "not implemented, needs checking" << std::endl;
    }

    /// Turn back information as a string.
    std::string PairingInfo(const int EchoLevel) const override;

    void SetValue(const Variable<bool>& rVariable, const bool rValue) override
    {
        if (rVariable == IS_PROJECTED_LOCAL_SYSTEM) mIsProjection = rValue;
        else MapperLocalSystem::SetValue(rVariable, rValue);
    }

    bool GetValue(const Variable<bool>& rVariable) override
    {
        if (rVariable == IS_PROJECTED_LOCAL_SYSTEM) return mIsProjection;
        else return MapperLocalSystem::GetValue(rVariable);
    }

private:
    GeometryPointerType mpGeom;
    bool mIsProjection; // Set to true is we are projecting the master onto the slave.
                        // Set to false if we are projecting the slave onto the slave.

};

template<class TSparseSpace, class TDenseSpace>
class CouplingGeometryMapper : public Mapper<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of CouplingGeometryMapper
    KRATOS_CLASS_POINTER_DEFINITION(CouplingGeometryMapper);

    typedef Mapper<TSparseSpace, TDenseSpace> BaseType;

    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
    typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;

    typedef InterfaceVectorContainer<TSparseSpace, TDenseSpace> InterfaceVectorContainerType;
    typedef Kratos::unique_ptr<InterfaceVectorContainerType> InterfaceVectorContainerPointerType;

    typedef std::size_t IndexType;

    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;
    typedef Kratos::unique_ptr<TMappingMatrixType> TMappingMatrixUniquePointerType;

    typedef Matrix DenseMappingMatrixType;
    typedef Kratos::unique_ptr<DenseMappingMatrixType> DenseMappingMatrixUniquePointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    CouplingGeometryMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                         : mrModelPartOrigin(rModelPartOrigin),
                           mrModelPartDestination(rModelPartDestination) {}

    CouplingGeometryMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination,
                         Parameters JsonParameters)
                         : mrModelPartOrigin(rModelPartOrigin),
                           mrModelPartDestination(rModelPartDestination),
                           mMapperSettings(JsonParameters)
    {
        mpInterfaceVectorContainerOrigin = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartOrigin.GetSubModelPart("interface"));
        mpInterfaceVectorContainerDestination = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartDestination.GetSubModelPart("interface"));

        mpCouplingMP = rModelPartOrigin.pGetSubModelPart("coupling");
        mpCouplingQuadraturePointsMP = rModelPartOrigin.pGetSubModelPart("coupling_quadrature_points");

        for (auto condition_itr = mpCouplingMP->ConditionsBegin();
            condition_itr != mpCouplingMP->ConditionsEnd();
            ++condition_itr)
            condition_itr->GetGeometry().SetValue(IS_DUAL_MORTAR, mMapperSettings["dual_mortar"].GetBool());

        this->InitializeInterface();
    }

    /// Destructor.
    ~CouplingGeometryMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    void Map(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        MapInternal(rOriginVariable, rDestinationVariable, MappingOptions,true);
    }

    void InverseMap(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        MapInternal(rOriginVariable, rDestinationVariable, MappingOptions, true);
    }

    ///@}
    ///@name Access
    ///@{

    TMappingMatrixType* pGetMappingMatrix() override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        return Kratos::make_unique<CouplingGeometryMapper<TSparseSpace, TDenseSpace>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CouplingGeometryMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CouplingGeometryMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

private:

    ///@name Private Operations
    ///@{

    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;
    ModelPart* mpCouplingMP;
    ModelPart* mpCouplingQuadraturePointsMP;

    Parameters mMapperSettings;

    MapperUniquePointerType mpInverseMapper = nullptr;

    DenseMappingMatrixUniquePointerType mpMappingMatrix;

    MapperLocalSystemPointerVector mMapperLocalSystems;

    InterfaceVectorContainerPointerType mpInterfaceVectorContainerOrigin;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerDestination;

    bool mIsInverseFlag = false;


    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    void AssignInterfaceEquationIds()
    {
        MapperUtilities::AssignInterfaceEquationIds(mrModelPartDestination.GetSubModelPart("interface").GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIds(mrModelPartOrigin.GetSubModelPart("interface").GetCommunicator());
    }

    void MapInternal(const Variable<double>& rOriginVariable,
                     const Variable<double>& rDestinationVariable,
                     Kratos::Flags MappingOptions, const bool IsInverse = false);

    void MapInternal(const Variable<array_1d<double, 3>>& rOriginVariable,
                     const Variable<array_1d<double, 3>>& rDestinationVariable,
                     Kratos::Flags MappingOptions, const bool IsInverse = false);

    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems)
    {
        MapperUtilities::CreateMapperLocalSystemsFromGeometries<CouplingGeometryLocalSystem>(
            rModelPartCommunicator,
            rLocalSystems);
    }

    void EnforceConsistencyWithScaling(const Matrix& rInterfaceMatrixSlave,
        Matrix& rInterfaceMatrixProjected,
        const double scalingLimit = 1.1);

    void CheckMappingMatrixConsistency()
    {
        for (size_t row = 0; row < mpMappingMatrix->size1(); ++row) {
            double row_sum = 0.0;
            for (size_t col = 0; col < mpMappingMatrix->size2(); ++col) row_sum += (*mpMappingMatrix)(row, col);
            if (std::abs(row_sum - 1.0) > 1e-6) {
                KRATOS_WATCH(*mpMappingMatrix)
                KRATOS_WATCH(row_sum)
                KRATOS_ERROR << "mapping matrix is not consistent\n";
            }
        }
    }

    Parameters GetMapperDefaultSettings() const
    {
        // @tobiasteschemachen
        return Parameters( R"({
            "echo_level" : 0
        })");
    }

    ///@}

}; // Class CouplingGeometryMapper

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_COUPLING_GEOMETRY_MAPPER_H_INCLUDED defined