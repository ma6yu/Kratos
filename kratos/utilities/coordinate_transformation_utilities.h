//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//
//


#ifndef KRATOS_COORDINATE_TRANSFORMATION_UTILITIES_H
#define KRATOS_COORDINATE_TRANSFORMATION_UTILITIES_H

// system includes

// external includes
#include "boost/numeric/ublas/matrix_proxy.hpp"

// kratos includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/variable.h"
#include "geometries/geometry.h"

namespace Kratos {

///@addtogroup KratosCore
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

/// A utility to rotate the local contributions of certain nodes to the system matrix, which is required to apply slip conditions in arbitrary directions.
template<class TLocalMatrixType, class TLocalVectorType, class TValueType>
class CoordinateTransformationUtils {
public:
	///@name Type Definitions
	///@{

	/// Pointer definition of CoordinateTransformationUtils
	KRATOS_CLASS_POINTER_DEFINITION(CoordinateTransformationUtils);

	typedef Node<3> NodeType;

	typedef Geometry< Node<3> > GeometryType;

//     typedef boost::numeric::ublas::matrix_row<TLocalMatrixType>  LocalRowType;
//
//     typedef boost::numeric::ublas::matrix_range<TLocalMatrixType> MatrixBlockType;

	///@}
	///@name Life Cycle
	///@{

	/// Constructor.
	/** @param DomainSize Number of space dimensions (2 or 3)
	 * @param NumRowsPerNode Number of matrix or vector rows associated to each node. Velocity DOFs are assumed to be the first mDomainSize rows in each block of rows.
	 * @param rSelectionFlag All nodes where the flag given by this argument is set to true will be transformed to a rotated coordinate system.
	 */
	CoordinateTransformationUtils(const unsigned int DomainSize,
			const unsigned int NumRowsPerNode,
			const Kratos::Flags& rSelectionFlag = SLIP):
	mDomainSize(DomainSize),
	mBlockSize(NumRowsPerNode),
	mrFlag(rSelectionFlag)
	{}

	/// Destructor.
	virtual ~CoordinateTransformationUtils() {}

	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	/**
	 * @brief Calculates rotation operator for given point
	 *
	 * This metod calculates rotation matrix for a given point. Nodal NORMAL variable should be
	 * assigned properly since rotation is calculated based on it.
	 *
	 * @param rRotationMatrix 	Output rotation matrix
	 * @param rThisPoint 		Current node
	 */
    virtual void CalculateRotationOperatorPure(
		TLocalMatrixType& rRotationMatrix,
        const GeometryType::PointType& rThisPoint) const
    {
        if (mDomainSize == 2) {
            BoundedMatrix<double, 2, 2> local_matrix;
            this->LocalRotationOperatorPure(local_matrix, rThisPoint);
            if (rRotationMatrix.size1() != 2 || rRotationMatrix.size2() != 2) {
                rRotationMatrix.resize(2, 2, false);
            }
            noalias(rRotationMatrix) = local_matrix;
        } else {
            BoundedMatrix<double, 3, 3> local_matrix;
            this->LocalRotationOperatorPure(local_matrix, rThisPoint);
            if (rRotationMatrix.size1() != 3 || rRotationMatrix.size2() != 3) {
                rRotationMatrix.resize(3, 3, false);
            }
            noalias(rRotationMatrix) = local_matrix;
        }
    }

	/**
	 * @brief Calculates rotation nodal matrix shape sensitivities
	 *
	 * This method calculates shape sensitivities of rotation matrix for given node.
	 * Nodal NORMAL(historical data container) and NORMAL_SHAPE_SENSITIVITY(non-historical data contaienr) variables
	 * should be properly initialized.
	 *
	 * NORMAL_SHAPE_SENSITIVITY matrix should be properly sized and initialized with proper shape sensitivity values
	 * 		rows: number_of_nodes contributing to NORMAL * DOMAIN_SIZE, columns: DOMAIN_SIZE
	 *
	 * @param rRotationMatrixShapeDerivative 	Output shape sensitivities matrix w.r.t. NodeIndex and DerivativeIndex
	 * @param DerivativeNodeIndex 				NodeIndex for which shape sensitivity matrix is computed
	 * @param DerivativeDirectionIndex			Direction index of the node for which shape sensitivity matrix is computed
	 * @param rThisPoint 						Current node where rotation matrix shape sensitivities are required
	 */
    virtual void CalculateRotationOperatorPureShapeSensitivities(
		TLocalMatrixType& rRotationMatrixShapeDerivative,
        const std::size_t DerivativeNodeIndex,
        const std::size_t DerivativeDirectionIndex,
        const GeometryType::PointType& rThisPoint) const
    {
        if (mDomainSize == 2) {
            BoundedMatrix<double, 2, 2> local_matrix;
            this->CalculateRotationOperatorPureShapeSensitivities(
                local_matrix, DerivativeNodeIndex, DerivativeDirectionIndex, rThisPoint);
            if (rRotationMatrixShapeDerivative.size1() != 2 ||
                rRotationMatrixShapeDerivative.size2() != 2) {
                rRotationMatrixShapeDerivative.resize(2, 2, false);
            }
            noalias(rRotationMatrixShapeDerivative) = local_matrix;
        } else {
            BoundedMatrix<double, 3, 3> local_matrix;
            this->CalculateRotationOperatorPureShapeSensitivities(
                local_matrix, DerivativeNodeIndex, DerivativeDirectionIndex, rThisPoint);
            if (rRotationMatrixShapeDerivative.size1() != 3 ||
                rRotationMatrixShapeDerivative.size2() != 3) {
                rRotationMatrixShapeDerivative.resize(3, 3, false);
            }
            noalias(rRotationMatrixShapeDerivative) = local_matrix;
        }
    }

	/**
	 * @brief Calculate 2d rotation nodal matrix shape sensitivities
	 *
	 * This method calculates shape sensitivities of 2D rotation matrix for given node.
	 * Nodal NORMAL(historical data container) and NORMAL_SHAPE_SENSITIVITY(non-historical data contaienr) variables
	 * should be properly initialized.
	 *
	 * NORMAL_SHAPE_SENSITIVITY matrix should be properly sized and initialized with proper shape sensitivity values
	 * 		rows: (number_of_neighbour_nodes + 1) * 2
	 * 		cols: 2
	 *
	 * @param rOutput							Output shape sensitivities matrix w.r.t. NodeIndex and DerivativeIndex
	 * @param DerivativeNodeIndex 				NodeIndex for which shape sensitivity matrix is computed
	 * @param DerivativeDirectionIndex			Direction index of the node for which shape sensitivity matrix is computed
	 * @param rThisPoint 						Current node where rotation matrix shape sensitivities are required
	 */
    virtual void CalculateRotationOperatorPureShapeSensitivities(
        BoundedMatrix<double, 2, 2>& rOutput,
        const std::size_t DerivativeNodeIndex,
        const std::size_t DerivativeDirectionIndex,
        const GeometryType::PointType& rThisPoint) const
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!rThisPoint.SolutionStepsDataHas(NORMAL))
            << "NORMAL is not found in node at " << rThisPoint.Coordinates() << ".";
        KRATOS_ERROR_IF(!rThisPoint.Has(NORMAL_SHAPE_DERIVATIVE))
            << "NORMAL_SHAPE_DERIVATIVE is not found in node at "
            << rThisPoint.Coordinates() << ".";

        const array_1d<double, 3>& r_nodal_normal =
            rThisPoint.FastGetSolutionStepValue(NORMAL);
		const double nodal_normal_magnitude = norm_2(r_nodal_normal);

        KRATOS_ERROR_IF(nodal_normal_magnitude == 0.0)
            << "NORMAL at node " << rThisPoint.Coordinates()
            << " is not properly initialized.";

        const Matrix& r_sensitivity_values = rThisPoint.GetValue(NORMAL_SHAPE_DERIVATIVE);

        KRATOS_DEBUG_ERROR_IF(r_sensitivity_values.size2() != 2)
            << "NORMAL_SHAPE_DERIVATIVE is not properly initialized at node "
            << rThisPoint.Coordinates() << " to calculate 2D rotation operator shape sensitivities. [ required number of columns = 2, available number of columns = "
            << r_sensitivity_values.size2() << " ].";

        const std::size_t require_rows = (DerivativeNodeIndex + 1) * 2;
        KRATOS_DEBUG_ERROR_IF(r_sensitivity_values.size1() < require_rows)
            << "NORMAL_SHAPE_DERIVATIVE is not properly initialized at node "
            << rThisPoint.Coordinates() << " to calculate 2D rotation operator shape sensitivities. [ required number of rows >= "
            << require_rows
            << ", available number of rows = " << r_sensitivity_values.size1() << " ].";

        const Vector& r_nodal_normal_derivatives =
            row(r_sensitivity_values, DerivativeNodeIndex * 2 + DerivativeDirectionIndex);

        rOutput(0, 0) = r_nodal_normal_derivatives[0] / nodal_normal_magnitude;
        rOutput(0, 1) = r_nodal_normal_derivatives[1] / nodal_normal_magnitude;
        rOutput(1, 0) = -r_nodal_normal_derivatives[1] / nodal_normal_magnitude;
        rOutput(1, 1) = r_nodal_normal_derivatives[0] / nodal_normal_magnitude;

        const double nodal_normal_magnitude_derivative =
            (r_nodal_normal[0] * r_nodal_normal_derivatives[0] +
             r_nodal_normal[1] * r_nodal_normal_derivatives[1]) /
            nodal_normal_magnitude;

        const double coeff = nodal_normal_magnitude_derivative /
                             (std::pow(nodal_normal_magnitude, 2));

        rOutput(0, 0) -= r_nodal_normal[0] * coeff;
        rOutput(0, 1) -= r_nodal_normal[1] * coeff;
        rOutput(1, 0) -= -r_nodal_normal[1] * coeff;
        rOutput(1, 1) -= r_nodal_normal[0] * coeff;

        KRATOS_CATCH("");
    }

	/**
	 * @brief Calculate 3d rotation nodal matrix shape sensitivities
	 *
	 * This method calculates shape sensitivities of 3D rotation matrix for given node.
	 * Nodal NORMAL(historical data container) and NORMAL_SHAPE_SENSITIVITY(non-historical data contaienr) variables
	 * should be properly initialized.
	 *
	 * NORMAL_SHAPE_SENSITIVITY matrix should be properly sized and initialized with proper shape sensitivity values
	 * 		rows: (number_of_neighbour_nodes + 1) * 3
	 * 		cols: 3
	 *
	 * @param rOutput						 	Output shape sensitivities matrix w.r.t. NodeIndex and DerivativeIndex
	 * @param DerivativeNodeIndex 				NodeIndex for which shape sensitivity matrix is computed
	 * @param DerivativeDirectionIndex			Direction index of the node for which shape sensitivity matrix is computed
	 * @param rThisPoint 						Current node where rotation matrix shape sensitivities are required
	 */
    virtual void CalculateRotationOperatorPureShapeSensitivities(
		BoundedMatrix<double, 3, 3>& rOutput,
        const std::size_t DerivativeNodeIndex,
        const std::size_t DerivativeDirectionIndex,
		const GeometryType::PointType& rThisPoint) const
	{
		KRATOS_TRY

        KRATOS_ERROR_IF(!rThisPoint.SolutionStepsDataHas(NORMAL))
            << "NORMAL is not found in node at " << rThisPoint.Coordinates() << ".";
        KRATOS_ERROR_IF(!rThisPoint.Has(NORMAL_SHAPE_DERIVATIVE))
            << "NORMAL_SHAPE_DERIVATIVE is not found in node at "
            << rThisPoint.Coordinates() << ".";

        const array_1d<double, 3>& r_nodal_normal =
            rThisPoint.FastGetSolutionStepValue(NORMAL);
		const double nodal_normal_magnitude = norm_2(r_nodal_normal);

        KRATOS_ERROR_IF(nodal_normal_magnitude == 0.0)
            << "NORMAL at node " << rThisPoint.Coordinates()
            << " is not properly initialized.";

        const Matrix& r_sensitivity_values = rThisPoint.GetValue(NORMAL_SHAPE_DERIVATIVE);

        KRATOS_DEBUG_ERROR_IF(r_sensitivity_values.size2() != 3)
            << "NORMAL_SHAPE_DERIVATIVE is not properly initialized at node "
            << rThisPoint.Coordinates() << " to calculate 3D rotation operator shape sensitivities. [ required number of columns = 3, available number of columns = "
            << r_sensitivity_values.size2() << " ].";

		const std::size_t require_rows = (DerivativeNodeIndex + 1) * 3;
        KRATOS_DEBUG_ERROR_IF(r_sensitivity_values.size1() < require_rows)
            << "NORMAL_SHAPE_DERIVATIVE is not properly initialized at node "
            << rThisPoint.Coordinates() << " to calculate 3D rotation operator shape sensitivities. [ required number of rows >= "
            << require_rows
            << ", available number of rows = " << r_sensitivity_values.size1() << " ].";

        const Vector& r_nodal_normal_derivative =
            row(r_sensitivity_values, DerivativeNodeIndex * 3 + DerivativeDirectionIndex);

        const double x0 = 1.0/nodal_normal_magnitude;
        const double x1 = r_nodal_normal_derivative[0]*x0;
        const double x2 = std::pow(nodal_normal_magnitude, -3);
        const double x3 = r_nodal_normal[0]*r_nodal_normal_derivative[0];
        const double x4 = r_nodal_normal[1]*r_nodal_normal_derivative[1];
        const double x5 = r_nodal_normal[2]*r_nodal_normal_derivative[2];
        const double x6 = x2*(-x3 - x4 - x5);
        const double x7 = r_nodal_normal[0]*x6;
        const double x8 = x1 + x7;
        const double x9 = r_nodal_normal_derivative[1]*x0;
        const double x10 = r_nodal_normal[1]*x6;
        const double x11 = x10 + x9;
        const double x12 = r_nodal_normal_derivative[2]*x0;
        const double x13 = r_nodal_normal[2]*x6;
        const double x14 = x12 + x13;
        const double x15 = std::pow(nodal_normal_magnitude, -2);
        const double x16 = 2*x3;
        const double x17 = 2*x4;
        const double x18 = -x16 - x17 - 2*x5;
        const double x19 = std::pow(r_nodal_normal[0], 2);
        const double x20 = std::pow(nodal_normal_magnitude, -4);
        const double x21 = x19*x20;
        const double x22 = x18*x21;
        const double x23 = std::pow(r_nodal_normal[1], 2);
        const double x24 = x21*x23;
        const double x25 = std::pow(r_nodal_normal[2], 2);
        const double x26 = -x15*x19 + 1;
        const double x27 = x21*x25 + x24 + std::pow(x26, 2);
        const double x28 = std::pow(x27, -1.0/2.0);
        const double x29 = x28*(-x15*x16 - x22);
        const double x30 = 4*x3;
        const double x31 = x20*x25;
        const double x32 = std::pow(nodal_normal_magnitude, -6);
        const double x33 = 4*x4;
        const double x34 = -x30 - x33 - 4*x5;
        const double x35 = (1.0/2.0)*x19*x32*x34;
        const double x36 = x20*x23;
        const double x37 = -x21*x4 - x23*x35 - x3*x36;
        const double x38 = (-x21*x5 - x25*x35 - 1.0/2.0*x26*(-x15*x30 - 2*x22) - x3*x31 + x37)/std::pow(x27, 3.0/2.0);
        const double x39 = x26*x38;
        const double x40 = x15*x28;
        const double x41 = r_nodal_normal[0]*x40;
        const double x42 = r_nodal_normal_derivative[0]*x40;
        const double x43 = r_nodal_normal[0]*r_nodal_normal[1];
        const double x44 = x18*x20;
        const double x45 = x28*x44;
        const double x46 = x15*x38;
        const double x47 = r_nodal_normal[0]*r_nodal_normal[2];
        const double x48 = r_nodal_normal_derivative[2]*x2;
        const double x49 = x19*x28;
        const double x50 = r_nodal_normal[2]*x2;
        const double x51 = x16*x28;
        const double x52 = x26*x28;
        const double x53 = (-3*x3 - 3*x4 - 3*x5)/std::pow(nodal_normal_magnitude, 5);
        const double x54 = x49*x53;
        const double x55 = x0*x29;
        const double x56 = x19*x38;
        const double x57 = x0*x39;
        const double x58 = r_nodal_normal[1]*x2;
        const double x59 = -x15*x23 + 1;
        const double x60 = x24 + x25*x36 + std::pow(x59, 2);
        const double x61 = std::pow(x60, -1.0/2.0);
        const double x62 = x15*x61;
        const double x63 = r_nodal_normal_derivative[1]*x62;
        const double x64 = r_nodal_normal[1]*x62;
        const double x65 = x44*x61;
        const double x66 = x18*x36;
        const double x67 = (-1.0/2.0*x23*x25*x32*x34 - x31*x4 - x36*x5 + x37 - 1.0/2.0*x59*(-x15*x33 - 2*x66))/std::pow(x60, 3.0/2.0);
        const double x68 = x15*x67;
        const double x69 = x61*(-x15*x17 - x66);
        const double x70 = x59*x67;
        const double x71 = r_nodal_normal[1]*r_nodal_normal[2];
        const double x72 = x17*x61;
        const double x73 = x23*x61;
        const double x74 = x59*x61;
        const double x75 = x53*x73;
        const double x76 = x0*x69;
        const double x77 = x23*x67;
        const double x78 = x0*x70;
        const double x79 = r_nodal_normal[0]*x2;


        if ((r_nodal_normal[0] / nodal_normal_magnitude) > 0.99) {
            rOutput(0, 0) = x8;
            rOutput(0, 1) = x11;
            rOutput(0, 2) = x14;
            rOutput(1, 0) = x29 + x39;
            rOutput(1, 1) = -r_nodal_normal[1]*x42 - r_nodal_normal_derivative[1]*x41 - x43*x45 - x43*x46;
            rOutput(1, 2) = -r_nodal_normal[2]*x42 - r_nodal_normal_derivative[2]*x41 - x45*x47 - x46*x47;
            rOutput(2, 0) = 0;
            rOutput(2, 1) = r_nodal_normal[2]*x54 + r_nodal_normal[2]*x55 + r_nodal_normal[2]*x57 + x12*x52 + x13*x52 + x48*x49 + x50*x51 + x50*x56;
            rOutput(2, 2) = -r_nodal_normal[1]*x54 - r_nodal_normal[1]*x55 - r_nodal_normal[1]*x57 - r_nodal_normal_derivative[1]*x2*x49 - x10*x52 - x51*x58 - x52*x9 - x56*x58;
        } else {
            rOutput(0, 0) = x8;
            rOutput(0, 1) = x11;
            rOutput(0, 2) = x14;
            rOutput(1, 0) = -r_nodal_normal[0]*x63 - r_nodal_normal_derivative[0]*x64 - x43*x65 - x43*x68;
            rOutput(1, 1) = x69 + x70;
            rOutput(1, 2) = -r_nodal_normal[2]*x63 - r_nodal_normal_derivative[2]*x64 - x65*x71 - x68*x71;
            rOutput(2, 0) = -r_nodal_normal[2]*x75 - r_nodal_normal[2]*x76 - r_nodal_normal[2]*x78 - x12*x74 - x13*x74 - x48*x73 - x50*x72 - x50*x77;
            rOutput(2, 1) = 0;
            rOutput(2, 2) = r_nodal_normal[0]*x75 + r_nodal_normal[0]*x76 + r_nodal_normal[0]*x78 + r_nodal_normal_derivative[0]*x2*x73 + x1*x74 + x7*x74 + x72*x79 + x77*x79;
        }

        KRATOS_CATCH("");
	}

    /// Rotate the local system contributions so that they are oriented with each node's normal.
	/**
	 @param rLocalMatrix Local system matrix
	 @param rLocalVector Local RHS vector
	 @param rGeometry A reference to the element's (or condition's) geometry
	 */
	virtual void Rotate(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		if(mBlockSize != mDomainSize) //Monolithic case
		{
			if(mDomainSize == 2) RotateAux<2,3>(rLocalMatrix,rLocalVector,rGeometry);
			if(mDomainSize == 3) RotateAux<3,4>(rLocalMatrix,rLocalVector,rGeometry);
		}
		else //fractional step case
		{
			if(mDomainSize == 2) RotateAuxPure<2>(rLocalMatrix,rLocalVector,rGeometry);
			if(mDomainSize == 3) RotateAuxPure<3>(rLocalMatrix,rLocalVector,rGeometry);
		}

	}

	/// RHS only version of Rotate
	virtual void Rotate(TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		//const unsigned int LocalSize = rLocalVector.size(); // We expect this to work both with elements (4 nodes) and conditions (3 nodes)

		unsigned int Index = 0;

		if (rLocalVector.size() > 0)
		{
			if(mBlockSize != mDomainSize) //Monolithic case
			{
				for(unsigned int j = 0; j < rGeometry.PointsNumber(); ++j)
				{
					if( this->IsSlip(rGeometry[j]) )
					{
						if(mDomainSize == 3)
						{
							array_1d<double,4> aux,aux1;
							BoundedMatrix<double,4,4> rRot;
							LocalRotationOperator3D<4>(rRot,rGeometry[j]);

							for(unsigned int k=0; k<4; k++)
							aux[k] = rLocalVector[j*mBlockSize+k];

							noalias(aux1) = prod(rRot,aux);

							for(unsigned int k=0; k<4; k++)
							rLocalVector[j*mBlockSize+k] = aux1[k];
						}
						else
						{
							array_1d<double,3> aux,aux1;
							BoundedMatrix<double,3,3> rRot;
							LocalRotationOperator2D<3>(rRot,rGeometry[j]);

							for(unsigned int k=0; k<3; k++)
							{
								aux[k] = rLocalVector[j*mBlockSize+k];
							}

							noalias(aux1) = prod(rRot,aux);

							for(unsigned int k=0; k<3; k++)
							rLocalVector[j*mBlockSize+k] = aux1[k];
						}
					}
					Index += mBlockSize;
				}

			}
			else //fractional step case
			{
				for(unsigned int j = 0; j < rGeometry.PointsNumber(); ++j)
				{
					if( this->IsSlip(rGeometry[j]) )
					{
						if(mDomainSize == 3)
						{
							array_1d<double,3> aux,aux1;
							BoundedMatrix<double,3,3> rRot;
							LocalRotationOperatorPure(rRot,rGeometry[j]);

							for(unsigned int k=0; k<3; k++)
							aux[k] = rLocalVector[j*mBlockSize+k];

							noalias(aux1) = prod(rRot,aux);

							for(unsigned int k=0; k<3; k++)
							rLocalVector[j*mBlockSize+k] = aux1[k];
						}
						else
						{
							array_1d<double,2> aux,aux1;
							BoundedMatrix<double,2,2> rRot;
							LocalRotationOperatorPure(rRot,rGeometry[j]);

							for(unsigned int k=0; k<2; k++)
								aux[k] = rLocalVector[j*mBlockSize+k];

							noalias(aux1) = prod(rRot,aux);

							for(unsigned int k=0; k<2; k++)
							rLocalVector[j*mBlockSize+k] = aux1[k];
						}
					}
					Index += mBlockSize;
				}

			}

		}

	}

	/// Apply slip boundary conditions to the rotated local contributions.
	/** This function takes the local system contributions rotated so each
	 node's velocities are expressed using a base oriented with its normal
	 and imposes that the normal velocity is equal to the mesh velocity in
	 the normal direction.
	 */
	virtual void ApplySlipCondition(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		const unsigned int LocalSize = rLocalVector.size(); // We expect this to work both with elements (4 nodes) and conditions (3 nodes)

		if (LocalSize > 0)
		{
			for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
			{
				if( this->IsSlip(rGeometry[itNode]))
				{
					// We fix the first dof (normal velocity) for each rotated block
					unsigned int j = itNode * mBlockSize;
					//const double k = rLocalMatrix(j,j)+rLocalMatrix(j,j+1)+rLocalMatrix(j,j+2);

					// If the mesh is moving, we must impose v_normal = vmesh_normal
					array_1d<double,3> VMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
					VMesh -= rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
					array_1d<double,3> rN = rGeometry[itNode].FastGetSolutionStepValue(NORMAL);
					this->Normalize(rN);

					for( unsigned int i = 0; i < j; ++i)// Skip term (i,i)
					{
						rLocalMatrix(i,j) = 0.0;
						rLocalMatrix(j,i) = 0.0;
					}
					for( unsigned int i = j+1; i < LocalSize; ++i)
					{
						rLocalMatrix(i,j) = 0.0;
						rLocalMatrix(j,i) = 0.0;
					}

					rLocalVector(j) = inner_prod(rN,VMesh);
					rLocalMatrix(j,j) = 1.0;
				}
			}
		}
	}

	/// RHS only version of ApplySlipCondition
	virtual void ApplySlipCondition(TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		if (rLocalVector.size() > 0)
		{
			for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
			{
				if( this->IsSlip(rGeometry[itNode]) )
				{
					// We fix the first dof (normal velocity) for each rotated block
					unsigned int j = itNode * mBlockSize;

					// If the mesh is moving, we must impose v_normal = vmesh_normal
					array_1d<double,3> VMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
					VMesh -= rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
					array_1d<double,3> rN = rGeometry[itNode].FastGetSolutionStepValue(NORMAL);
					this->Normalize(rN);

					rLocalVector[j] = inner_prod(rN,VMesh);
				}
			}
		}
	}

	/// Transform nodal velocities to the rotated coordinates (aligned with each node's normal)
	virtual void RotateVelocities(ModelPart& rModelPart) const
	{
		TLocalVectorType Vel(mDomainSize);
		TLocalVectorType Tmp(mDomainSize);

		ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();
#pragma omp parallel for firstprivate(Vel,Tmp)
		for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
		{
			ModelPart::NodeIterator itNode = it_begin+iii;
			if( this->IsSlip(*itNode) )
			{
				//this->RotationOperator<TLocalMatrixType>(Rotation,);
				if(mDomainSize == 3)
				{
					BoundedMatrix<double,3,3> rRot;
					LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
					for(unsigned int i = 0; i < 3; i++) Vel[i] = rVelocity[i];
					noalias(Tmp) = prod(rRot,Vel);
					for(unsigned int i = 0; i < 3; i++) rVelocity[i] = Tmp[i];
				}
				else
				{
					BoundedMatrix<double,2,2> rRot;
					LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
					for(unsigned int i = 0; i < 2; i++) Vel[i] = rVelocity[i];
					noalias(Tmp) = prod(rRot,Vel);
					for(unsigned int i = 0; i < 2; i++) rVelocity[i] = Tmp[i];
				}
			}
		}
	}

	/// Transform nodal velocities from the rotated system to the original one
	virtual void RecoverVelocities(ModelPart& rModelPart) const
	{
		TLocalVectorType Vel(mDomainSize);
		TLocalVectorType Tmp(mDomainSize);

		ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();
#pragma omp parallel for firstprivate(Vel,Tmp)
		for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
		{
			ModelPart::NodeIterator itNode = it_begin+iii;
			if( this->IsSlip(*itNode) )
			{
				if(mDomainSize == 3)
				{
					BoundedMatrix<double,3,3> rRot;
					LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
					for(unsigned int i = 0; i < 3; i++) Vel[i] = rVelocity[i];
					noalias(Tmp) = prod(trans(rRot),Vel);
					for(unsigned int i = 0; i < 3; i++) rVelocity[i] = Tmp[i];
				}
				else
				{
					BoundedMatrix<double,2,2> rRot;
					LocalRotationOperatorPure(rRot,*itNode);

					array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
					for(unsigned int i = 0; i < 2; i++) Vel[i] = rVelocity[i];
					noalias(Tmp) = prod(trans(rRot),Vel);
					for(unsigned int i = 0; i < 2; i++) rVelocity[i] = Tmp[i];
				}
			}
		}
	}

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
	virtual std::string Info() const
	{
		std::stringstream buffer;
		buffer << "CoordinateTransformationUtils";
		return buffer.str();
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const
	{
		rOStream << "CoordinateTransformationUtils";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const {}

	///@}
	///@name Friends
	///@{

	///@}

protected:
	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

	template<unsigned int TDim, unsigned int TBlockSize, unsigned int TSkip = 0>
	void RotateAux(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		const unsigned int LocalSize = rLocalVector.size();

		unsigned int Index = 0;
		int rotations_needed = 0;
		const unsigned int NumBlocks = LocalSize / TBlockSize;
		DenseVector<bool> NeedRotation( NumBlocks, false);

		std::vector< BoundedMatrix<double,TBlockSize,TBlockSize> > rRot(NumBlocks);
		for(unsigned int j = 0; j < NumBlocks; ++j)
		{
			if( this->IsSlip(rGeometry[j]) )
			{
				NeedRotation[j] = true;
				rotations_needed++;

				if (TDim == 2) LocalRotationOperator2D<TBlockSize,TSkip>(rRot[j],rGeometry[j]);
				else LocalRotationOperator3D<TBlockSize,TSkip>(rRot[j],rGeometry[j]);
			}

			Index += TBlockSize;
		}

		if(rotations_needed > 0)
		{
			BoundedMatrix<double,TBlockSize,TBlockSize> mat_block, tmp;
			array_1d<double,TBlockSize> aux, aux1;

			for(unsigned int i=0; i<NumBlocks; i++)
			{
				if(NeedRotation[i] == true)
				{
					for(unsigned int j=0; j<NumBlocks; j++)
					{
						if(NeedRotation[j] == true)
						{
							ReadBlockMatrix<TBlockSize>(mat_block, rLocalMatrix, i*TBlockSize, j*TBlockSize);
							noalias(tmp) = prod(mat_block,trans(rRot[j]));
							noalias(mat_block) = prod(rRot[i],tmp);
							WriteBlockMatrix<TBlockSize>(mat_block, rLocalMatrix, i*TBlockSize, j*TBlockSize);
						}
						else
						{
							ReadBlockMatrix<TBlockSize>(mat_block, rLocalMatrix, i*TBlockSize, j*TBlockSize);
							noalias(tmp) = prod(rRot[i],mat_block);
							WriteBlockMatrix<TBlockSize>(tmp, rLocalMatrix, i*TBlockSize, j*TBlockSize);
						}
					}

					for(unsigned int k=0; k<TBlockSize; k++)
					aux[k] = rLocalVector[i*TBlockSize+k];

					noalias(aux1) = prod(rRot[i],aux);

					for(unsigned int k=0; k<TBlockSize; k++)
					rLocalVector[i*TBlockSize+k] = aux1[k];

				}
				else
				{
					for(unsigned int j=0; j<NumBlocks; j++)
					{
						if(NeedRotation[j] == true)
						{
							ReadBlockMatrix<TBlockSize>(mat_block, rLocalMatrix, i*TBlockSize, j*TBlockSize);
							noalias(tmp) = prod(mat_block,trans(rRot[j]));
							WriteBlockMatrix<TBlockSize>(tmp, rLocalMatrix, i*TBlockSize, j*TBlockSize);
						}
					}
				}

			}
		}
	}

	//to be used when there is only velocity (no additional pressure or other var block)
	template<unsigned int TDim>
	void RotateAuxPure(TLocalMatrixType& rLocalMatrix,
			TLocalVectorType& rLocalVector,
			GeometryType& rGeometry) const
	{
		const unsigned int LocalSize = rLocalVector.size();

		unsigned int Index = 0;
		int rotations_needed = 0;
		const unsigned int NumBlocks = LocalSize / mBlockSize;
		DenseVector<bool> NeedRotation( NumBlocks, false);

		std::vector< BoundedMatrix<double,TDim,TDim> > rRot(NumBlocks);
		for(unsigned int j = 0; j < NumBlocks; ++j)
		{
			if( this->IsSlip(rGeometry[j]) )
			{
				NeedRotation[j] = true;
				rotations_needed++;

				LocalRotationOperatorPure(rRot[j],rGeometry[j]);
			}

			Index += mBlockSize;
		}

		if(rotations_needed > 0)
		{
			BoundedMatrix<double,TDim,TDim> mat_block, tmp;
			array_1d<double,TDim> aux, aux1;

			for(unsigned int i=0; i<NumBlocks; i++)
			{
				if(NeedRotation[i] == true)
				{
					for(unsigned int j=0; j<NumBlocks; j++)
					{
						if(NeedRotation[j] == true)
						{
							ReadBlockMatrix<TDim>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
							noalias(tmp) = prod(mat_block,trans(rRot[j]));
							noalias(mat_block) = prod(rRot[i],tmp);
							WriteBlockMatrix<TDim>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
						}
						else
						{
							ReadBlockMatrix<TDim>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
							noalias(tmp) = prod(rRot[i],mat_block);
							WriteBlockMatrix<TDim>(tmp, rLocalMatrix, i*mBlockSize, j*mBlockSize);
						}
					}

					for(unsigned int k=0; k<TDim; k++)
					aux[k] = rLocalVector[i*mBlockSize+k];

					noalias(aux1) = prod(rRot[i],aux);

					for(unsigned int k=0; k<TDim; k++)
					rLocalVector[i*mBlockSize+k] = aux1[k];

				}
				else
				{
					for(unsigned int j=0; j<NumBlocks; j++)
					{
						if(NeedRotation[j] == true)
						{
							ReadBlockMatrix<TDim>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
							noalias(tmp) = prod(mat_block,trans(rRot[j]));
							WriteBlockMatrix<TDim>(tmp, rLocalMatrix, i*mBlockSize, j*mBlockSize);
						}
					}
				}

			}
		}
	}

	template<unsigned int TBlockSize, unsigned int TSkip = 0>
	void LocalRotationOperator2D(
		BoundedMatrix<double,TBlockSize,TBlockSize>& rRot,
		GeometryType::PointType& rThisPoint) const
	{
		noalias(rRot) = IdentityMatrix(TBlockSize);

		// Get the normal evaluated at the node
		const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

		double aux = rNormal[0]*rNormal[0] + rNormal[1]*rNormal[1];
		aux = sqrt(aux);

		rRot(TSkip  ,TSkip  ) = rNormal[0]/aux;
		rRot(TSkip  ,TSkip+1) = rNormal[1]/aux;
		rRot(TSkip+1,TSkip  ) = -rNormal[1]/aux;
		rRot(TSkip+1,TSkip+1) = rNormal[0]/aux;
	}

	template<unsigned int TBlockSize, unsigned int TSkip = 0>
	void LocalRotationOperator3D(
		BoundedMatrix<double,TBlockSize,TBlockSize>& rRot,
		GeometryType::PointType& rThisPoint) const
	{
		noalias(rRot) = IdentityMatrix(TBlockSize);

		// Get the normal evaluated at the node
		const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

		double aux = rNormal[0]*rNormal[0] + rNormal[1]*rNormal[1] + rNormal[2]*rNormal[2];
		aux = sqrt(aux);
		rRot(TSkip,TSkip  ) = rNormal[0]/aux;
		rRot(TSkip,TSkip+1) = rNormal[1]/aux;
		rRot(TSkip,TSkip+2) = rNormal[2]/aux;
		// Define the new coordinate system, where the first vector is aligned with the normal

		// To choose the remaining two vectors, we project the first component of the cartesian base to the tangent plane
		array_1d<double,3> rT1;
		rT1(0) = 1.0;
		rT1(1) = 0.0;
		rT1(2) = 0.0;
		double dot = rRot(TSkip,TSkip);//this->Dot(rN,rT1);

		// It is possible that the normal is aligned with (1,0,0), resulting in norm(rT1) = 0
		// If this is the case, repeat the procedure using (0,1,0)
		if ( fabs(dot) > 0.99 )
		{
			rT1(0) = 0.0;
			rT1(1) = 1.0;
			rT1(2) = 0.0;

			dot = rRot(TSkip,TSkip+1); //this->Dot(rN,rT1);
		}

		// calculate projection and normalize
		rT1[0] -= dot*rRot(TSkip,TSkip);
		rT1[1] -= dot*rRot(TSkip,TSkip+1);
		rT1[2] -= dot*rRot(TSkip,TSkip+2);
		this->Normalize(rT1);
		rRot(TSkip+1,TSkip  ) = rT1[0];
		rRot(TSkip+1,TSkip+1) = rT1[1];
		rRot(TSkip+1,TSkip+2) = rT1[2];

		// The third base component is choosen as N x T1, which is normalized by construction
		rRot(TSkip+2,TSkip  ) = rRot(TSkip,TSkip+1)*rT1[2] - rRot(TSkip,TSkip+2)*rT1[1];
		rRot(TSkip+2,TSkip+1) = rRot(TSkip,TSkip+2)*rT1[0] - rRot(TSkip,TSkip  )*rT1[2];
		rRot(TSkip+2,TSkip+2) = rRot(TSkip,TSkip  )*rT1[1] - rRot(TSkip,TSkip+1)*rT1[0];
	}


	void LocalRotationOperatorPure(BoundedMatrix<double,3,3>& rRot,
			const GeometryType::PointType& rThisPoint) const
	{

		// Get the normal evaluated at the node
		const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

		double aux = rNormal[0]*rNormal[0] + rNormal[1]*rNormal[1] + rNormal[2]*rNormal[2];
		aux = sqrt(aux);
		rRot(0,0) = rNormal[0]/aux;
		rRot(0,1) = rNormal[1]/aux;
		rRot(0,2) = rNormal[2]/aux;
		// Define the new coordinate system, where the first vector is aligned with the normal

		// To choose the remaining two vectors, we project the first component of the cartesian base to the tangent plane
		array_1d<double,3> rT1;
		rT1(0) = 1.0;
		rT1(1) = 0.0;
		rT1(2) = 0.0;
		double dot = rRot(0,0);//this->Dot(rN,rT1);

		// It is possible that the normal is aligned with (1,0,0), resulting in norm(rT1) = 0
		// If this is the case, repeat the procedure using (0,1,0)
		if ( fabs(dot) > 0.99 )
		{
			rT1(0) = 0.0;
			rT1(1) = 1.0;
			rT1(2) = 0.0;

			dot = rRot(0,1); //this->Dot(rN,rT1);
		}

		// calculate projection and normalize
		rT1[0] -= dot*rRot(0,0);
		rT1[1] -= dot*rRot(0,1);
		rT1[2] -= dot*rRot(0,2);
		this->Normalize(rT1);
		rRot(1,0) = rT1[0];
		rRot(1,1) = rT1[1];
		rRot(1,2) = rT1[2];

		// The third base component is choosen as N x T1, which is normalized by construction
		rRot(2,0) = rRot(0,1)*rT1[2] - rRot(0,2)*rT1[1];
		rRot(2,1) = rRot(0,2)*rT1[0] - rRot(0,0)*rT1[2];
		rRot(2,2) = rRot(0,0)*rT1[1] - rRot(0,1)*rT1[0];
	}

	void LocalRotationOperatorPure(BoundedMatrix<double,2,2>& rRot,
			const GeometryType::PointType& rThisPoint) const
	{
		// Get the normal evaluated at the node
		const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

		double aux = rNormal[0]*rNormal[0] + rNormal[1]*rNormal[1];
		aux = sqrt(aux);

		rRot(0,0) = rNormal[0]/aux;
		rRot(0,1) = rNormal[1]/aux;
		rRot(1,0) = -rNormal[1]/aux;
		rRot(1,1) = rNormal[0]/aux;

	}

	bool IsSlip(const Node<3>& rNode) const
	{
		return rNode.Is(mrFlag);
	}

	/// Normalize a vector.
	/**
	 * @param rThis the vector
	 * @return Original norm of the input vector
	 */
	template< class TVectorType >
	double Normalize(TVectorType& rThis) const
	{
		double Norm = 0;
		for(typename TVectorType::iterator iComponent = rThis.begin(); iComponent < rThis.end(); ++iComponent)
		Norm += (*iComponent)*(*iComponent);
		Norm = sqrt(Norm);
		for(typename TVectorType::iterator iComponent = rThis.begin(); iComponent < rThis.end(); ++iComponent)
		*iComponent /= Norm;
		return Norm;
	}

	///@}
	///@name Protected  Access
	///@{

	unsigned int GetDomainSize() const
	{
		return mDomainSize;
	}

	unsigned int GetBlockSize() const
	{
		return mBlockSize;
	}

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	///@}

private:
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	/// Number of spatial dimensions
	const unsigned int mDomainSize;

	/// Number of matrix or vector rows associated to each node.
	/** @note Velocity Dofs are assumed to be the first mDomainSize rows.
	 */
	const unsigned int mBlockSize;

	const Kratos::Flags& mrFlag;

	///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

//     /// Compute a rotation matrix to transform values from the cartesian base to one oriented with the node's normal
//     /**
//      * The normal is read from solution step data NORMAL. Use NormalCalculationUtils::CalculateOnSimplex to
//      * obtain and store the nodal normal from the normals of the model's conditons.
//      * @param rRot The rotation matrix (output)
//      * @param rThisPoint The point used to orient the new coordinate system.
//      * @see NormalCalculationUtils
//      */
//     template<class TMatrixType>
//     void RotationOperator(TMatrixType& rRot,
//                           GeometryType::PointType& rThisPoint) const
//     {
//         typedef boost::numeric::ublas::matrix_row<TMatrixType> ThisRowType;
//         // Get the normal evaluated at the node
//         const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);
//
//         if(mDomainSize == 3)
//         {
//             // Define the new coordinate system, where the first vector is aligned with the normal
//             ThisRowType rN(rRot,0);
//             for( unsigned int i = 0; i < 3; ++i)
//                 rN[i] = rNormal[i];
//             this->Normalize(rN);
//
//             // To choose the remaining two vectors, we project the first component of the cartesian base to the tangent plane
//             ThisRowType rT1(rRot,1);
//             rT1(0) = 1.0;
//             rT1(1) = 0.0;
//             rT1(2) = 0.0;
//
//             double dot = this->Dot(rN,rT1);
//
//             // It is possible that the normal is aligned with (1,0,0), resulting in norm(rT1) = 0
//             // If this is the case, repeat the procedure using (0,1,0)
//             if ( fabs(dot) > 0.99 )
//             {
//                 rT1(0) = 0.0;
//                 rT1(1) = 1.0;
//                 rT1(2) = 0.0;
//
//                 dot = this->Dot(rN,rT1);
//             }
//
//             // calculate projection and normalize
//             rT1 -= dot * rN;
//             this->Normalize(rT1);
//
//             // The third base component is choosen as N x T1, which is normalized by construction
//             ThisRowType rT2(rRot,2);
//             rT2(0) = rN(1)*rT1(2) - rN(2)*rT1(1);
//             rT2(1) = rN(2)*rT1(0) - rN(0)*rT1(2);
//             rT2(2) = rN(0)*rT1(1) - rN(1)*rT1(0);
//         }
//         else //if(mDomainSize == 2)
//         {
//             /* The basis for the new coordinate system is (normal,tangent)
//                Tangent vector is chosen (-normal_y, normal_x) so that the resulting base
//                is right-handed.
//              */
//             ThisRowType rN(rRot,0);
//             ThisRowType rT(rRot,1);
//
//             rN[0] = rNormal[0];
//             rN[1] = rNormal[1];
//             this->Normalize(rN);
//             rT[0] = -rN[1];
//             rT[1] = rN[0];
//         }
//
//     }

	template< class TVectorType >
	double Dot(const TVectorType& rV1,
			const TVectorType& rV2) const
	{
		double dot = 0.0;
		for( typename TVectorType::const_iterator iV1 = rV1.begin(),iV2 = rV2.begin(); iV1 != rV1.end(); ++iV1, ++iV2)
		{
			dot += (*iV1) * (*iV2);
		}
		return dot;
	}

	/// Transform a local contribution from cartesian coordinates to rotated ones
//     void ApplyRotation(TLocalMatrixType& rMatrix,
//                        const TLocalMatrixType& rRotation) const
//     {
//         // compute B = R*A*transpose(R)
//         const unsigned int LocalSize = rMatrix.size1();
//         const unsigned int NumBlocks = LocalSize / mBlockSize;
//         //TLocalMatrixType Tmp = ZeroMatrix(LocalSize,LocalSize);
// /*
//         for (unsigned int iBlock = 0; iBlock < NumBlocks; iBlock++)
//         {
//             for (unsigned int jBlock = 0; jBlock < NumBlocks; jBlock++)
//             {
//                 for (unsigned int i = iBlock*mBlockSize; i < (iBlock+1)*mBlockSize; i++)
//                 {
//                     for(unsigned int j = jBlock*mBlockSize; j < (jBlock+1)*mBlockSize; j++)
//                     {
//                         double& tij = Tmp(i,j);
//                         for(unsigned int k = iBlock*mBlockSize; k < (iBlock+1)*mBlockSize; k++)
//                         {
//                             for(unsigned int l = jBlock*mBlockSize; l < (jBlock+1)*mBlockSize; l++)
//                             {
//                                 tij += rRotation(i,k)*rMatrix(k,l)*rRotation(j,l);
//                             }
//                         }
//                     }
//                 }
//             }
//         }*/
//
// 	Matrix Tmp = prod(rMatrix,trans(rRotation));
// 	noalias(rMatrix) = prod(rRotation,Tmp);
//
// //         noalias(rMatrix) = Tmp;
//     }

	//auxiliary functions
	template< unsigned int TBlockSize >
	void ReadBlockMatrix( BoundedMatrix<double,TBlockSize, TBlockSize>& block, const Matrix& origin, const unsigned int Ibegin, const unsigned int Jbegin) const
	{
		for(unsigned int i=0; i<TBlockSize; i++)
		{
			for(unsigned int j=0; j<TBlockSize; j++)
			{
				block(i,j) = origin(Ibegin+i, Jbegin+j);
			}
		}
	}

	template< unsigned int TBlockSize >
	void WriteBlockMatrix( const BoundedMatrix<double,TBlockSize, TBlockSize>& block, Matrix& destination, const unsigned int Ibegin, const unsigned int Jbegin) const
	{
		for(unsigned int i=0; i<TBlockSize; i++)
		{
			for(unsigned int j=0; j<TBlockSize; j++)
			{
				destination(Ibegin+i, Jbegin+j) = block(i,j);
			}
		}
	}

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
	CoordinateTransformationUtils& operator=(CoordinateTransformationUtils const& rOther) {}

	/// Copy constructor.
	CoordinateTransformationUtils(CoordinateTransformationUtils const& rOther) {}

	///@}
};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TLocalMatrixType, class TLocalVectorType, class TValueType>
inline std::istream& operator >>(std::istream& rIStream,
		CoordinateTransformationUtils<TLocalMatrixType, TLocalVectorType,
				TValueType>& rThis) {
	return rIStream;
}

/// output stream function
template<class TLocalMatrixType, class TLocalVectorType, class TValueType>
inline std::ostream& operator <<(std::ostream& rOStream,
		const CoordinateTransformationUtils<TLocalMatrixType, TLocalVectorType,
				TValueType>& rThis) {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}

///@}

///@} addtogroup block

}

#endif // KRATOS_COORDINATE_TRANSFORMATION_UTILITIES_H
