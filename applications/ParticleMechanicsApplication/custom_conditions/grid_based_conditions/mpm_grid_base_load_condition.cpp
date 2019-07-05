//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes

// External includes

// Project includes
#include "custom_conditions/grid_based_conditions/mpm_grid_base_load_condition.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************

    void MPMGridBaseLoadCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.size();
        const unsigned int dim = r_geometry.WorkingSpaceDimension();
        if (rResult.size() != dim * number_of_nodes)
        {
            rResult.resize(dim*number_of_nodes,false);
        }

        const unsigned int pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                const unsigned int index = i * 2;
                rResult[index    ] = r_geometry[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
            }
        }
        else
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                const unsigned int index = i * 3;
                rResult[index    ] = r_geometry[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
                rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z,pos + 2).EquationId();
            }
        }
        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************
    void MPMGridBaseLoadCondition::GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY

        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.size();
        const unsigned int dim =  r_geometry.WorkingSpaceDimension();
        ElementalDofList.resize(0);
        ElementalDofList.reserve(dim * number_of_nodes);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                ElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
            }
        }
        else
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                ElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
                ElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Z));
            }
        }
        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************

    void MPMGridBaseLoadCondition::GetValuesVector(
        Vector& rValues,
        int Step
        )
    {
        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.size();
        const unsigned int dim = r_geometry.WorkingSpaceDimension();
        const unsigned int MatSize = number_of_nodes * dim;

        if (rValues.size() != MatSize)
        {
            rValues.resize(MatSize, false);
        }

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            const array_1d<double, 3 > & Displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            unsigned int index = i * dim;
            for(unsigned int k = 0; k < dim; ++k)
            {
                rValues[index + k] = Displacement[k];
            }
        }
    }

    //***********************************************************************
    //***********************************************************************

    void MPMGridBaseLoadCondition::GetFirstDerivativesVector(
        Vector& rValues,
        int Step
        )
    {
        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.size();
        const unsigned int dim = r_geometry.WorkingSpaceDimension();
        const unsigned int MatSize = number_of_nodes * dim;

        if (rValues.size() != MatSize)
        {
            rValues.resize(MatSize, false);
        }

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            const array_1d<double, 3 > & Velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
            const unsigned int index = i * dim;
            for(unsigned int k = 0; k<dim; ++k)
            {
                rValues[index + k] = Velocity[k];
            }
        }
    }

    //***********************************************************************
    //***********************************************************************

    void MPMGridBaseLoadCondition::GetSecondDerivativesVector(
        Vector& rValues,
        int Step
        )
    {
        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.size();
        const unsigned int dim = r_geometry.WorkingSpaceDimension();
        const unsigned int MatSize = number_of_nodes * dim;

        if (rValues.size() != MatSize)
        {
            rValues.resize(MatSize, false);
        }

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            const array_1d<double, 3 > & Acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const unsigned int index = i * dim;
            for(unsigned int k = 0; k < dim; ++k)
            {
                rValues[index + k] = Acceleration[k];
            }
        }
    }

    //************************************************************************************
    //************************************************************************************

    void MPMGridBaseLoadCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = false;
        const bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    //************************************************************************************
    //************************************************************************************
    void MPMGridBaseLoadCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = true;

        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    //***********************************************************************
    //***********************************************************************

    void MPMGridBaseLoadCondition::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        if(rMassMatrix.size1() != 0)
        {
            rMassMatrix.resize(0, 0, false);
        }
    }

    //***********************************************************************
    //***********************************************************************

    void MPMGridBaseLoadCondition::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        if(rDampingMatrix.size1() != 0)
        {
            rDampingMatrix.resize(0, 0, false);
        }
    }

    //***********************************************************************
    //***********************************************************************

    void MPMGridBaseLoadCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag
        )
    {
        KRATOS_ERROR << "You are calling the CalculateAll from the base class for loads" << std::endl;
    }

    //***********************************************************************
    //***********************************************************************

    int MPMGridBaseLoadCondition::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        if ( DISPLACEMENT.Key() == 0 )
        {
            KRATOS_ERROR <<  "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        }

        //verify that the dofs exist
        for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
        {
            if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            {
                KRATOS_ERROR << "missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
            }

            if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false ||
                 this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false ||
                 this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
            {
                KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << " of condition " << Id() << std::endl;
            }
        }

        return 0;
    }

    //***********************************************************************
    //***********************************************************************

    double MPMGridBaseLoadCondition::GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber,
        const double detJ
        )
    {
        return IntegrationPoints[PointNumber].Weight() * detJ;
    }

    //***********************************************************************************
    //***********************************************************************************

    void MPMGridBaseLoadCondition::AddExplicitContribution(const VectorType& rRHS,
        const Variable<VectorType>& rRHSVariable,
        Variable<array_1d<double,3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.PointsNumber();
        const unsigned int dimension       = r_geometry.WorkingSpaceDimension();

        if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
        {

            for(SizeType i=0; i< number_of_nodes; i++)
            {
                SizeType index = dimension * i;

                if (r_geometry[i].SolutionStepsDataHas(FORCE_RESIDUAL))
                {
                    array_1d<double, 3 > &ForceResidual = r_geometry[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
                    for(SizeType j=0; j<dimension; j++)
                    {
                        #pragma omp atomic
                        ForceResidual[j] += rRHS[index + j];
                    }
                }
            }
        }

    KRATOS_CATCH( "" )
    }

} // Namespace Kratos


