//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// External includes

// Project includes

// Include base h
#include "custom_utilities/fluid_calculation_utilities.h"

namespace Kratos
{

template<>
KRATOS_API(FLUID_DYNAMICS_APPLICATION)
void FluidCalculationUtilities::AssignValue(
    const array_1d<double, 3>& rInput,
    array_1d<double, 2>& rOutput)
{
    rOutput[0] = rInput[0];
    rOutput[1] = rInput[1];
}

template<class TOutputDataType, class TInputDataType>
void FluidCalculationUtilities::AssignValue(
    const TInputDataType& rInput,
    TOutputDataType& rOutput)
{
    rOutput = rInput;
}

template<>
KRATOS_API(FLUID_DYNAMICS_APPLICATION)
void FluidCalculationUtilities::UpdateValue(
    const double& rInput,
    double& rOutput)
{
    rOutput += rInput;
}

template<>
KRATOS_API(FLUID_DYNAMICS_APPLICATION)
void FluidCalculationUtilities::UpdateValue(
    const array_1d<double, 3>& rInput,
    array_1d<double, 2>& rOutput)
{
    rOutput[0] += rInput[0];
    rOutput[1] += rInput[1];
}

template<class TOutputDataType, class TInputDataType>
void FluidCalculationUtilities::UpdateValue(
    const TInputDataType& rInput,
    TOutputDataType& rOutput)
{
    noalias(rOutput) += rInput;
}

template<>
void FluidCalculationUtilities::AssignGradientValue(
    const double& rInput,
    const double rShapeFunctionDerivative,
    const IndexType DirectionIndex,
    array_1d<double, 3>& rOutput)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DirectionIndex > 2)
        << "Unsupported Direction Index. DirectionIndex should be  less than 3 "
           "[ DirectionIndex = "
        << DirectionIndex << " ].\n";

    rOutput[DirectionIndex] = rInput * rShapeFunctionDerivative;

    KRATOS_CATCH("");
}

template<>
void FluidCalculationUtilities::AssignGradientValue(
    const double& rInput,
    const double rShapeFunctionDerivative,
    const IndexType DirectionIndex,
    array_1d<double, 2>& rOutput)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DirectionIndex > 1)
        << "Unsupported Direction Index. DirectionIndex should be  less than 2 "
           "[ DirectionIndex = "
        << DirectionIndex << " ].\n";

    rOutput[DirectionIndex] = rInput * rShapeFunctionDerivative;

    KRATOS_CATCH("");
}

template<>
void FluidCalculationUtilities::AssignGradientValue(
    const array_1d<double, 3>& rInput,
    const double rShapeFunctionDerivative,
    const IndexType DirectionIndex,
    BoundedMatrix<double, 3, 3>& rOutput)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DirectionIndex > 2)
        << "Unsupported Direction Index. DirectionIndex should be  less than 3 "
           "[ DirectionIndex = "
        << DirectionIndex << " ].\n";

    rOutput(0, DirectionIndex) = rInput[0] * rShapeFunctionDerivative;
    rOutput(1, DirectionIndex) = rInput[1] * rShapeFunctionDerivative;
    rOutput(2, DirectionIndex) = rInput[2] * rShapeFunctionDerivative;

    KRATOS_CATCH("");
}

template<>
void FluidCalculationUtilities::AssignGradientValue(
    const array_1d<double, 3>& rInput,
    const double rShapeFunctionDerivative,
    const IndexType DirectionIndex,
    BoundedMatrix<double, 2, 2>& rOutput)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DirectionIndex > 1)
        << "Unsupported Direction Index. DirectionIndex should be  less than 2 "
           "[ DirectionIndex = "
        << DirectionIndex << " ].\n";

    rOutput(0, DirectionIndex) = rInput[0] * rShapeFunctionDerivative;
    rOutput(1, DirectionIndex) = rInput[1] * rShapeFunctionDerivative;

    KRATOS_CATCH("");
}

template<>
void FluidCalculationUtilities::UpdateGradientValue(
    const double& rInput,
    const double rShapeFunctionDerivative,
    const IndexType DirectionIndex,
    array_1d<double, 3>& rOutput)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DirectionIndex > 2)
        << "Unsupported Direction Index. DirectionIndex should be  less than 3 "
           "[ DirectionIndex = "
        << DirectionIndex << " ].\n";

    rOutput[DirectionIndex] += rInput * rShapeFunctionDerivative;

    KRATOS_CATCH("");
}

template<>
void FluidCalculationUtilities::UpdateGradientValue(
    const double& rInput,
    const double rShapeFunctionDerivative,
    const IndexType DirectionIndex,
    array_1d<double, 2>& rOutput)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DirectionIndex > 1)
        << "Unsupported Direction Index. DirectionIndex should be  less than 2 "
           "[ DirectionIndex = "
        << DirectionIndex << " ].\n";

    rOutput[DirectionIndex] += rInput * rShapeFunctionDerivative;

    KRATOS_CATCH("");
}

template<>
void FluidCalculationUtilities::UpdateGradientValue(
    const array_1d<double, 3>& rInput,
    const double rShapeFunctionDerivative,
    const IndexType DirectionIndex,
    BoundedMatrix<double, 3, 3>& rOutput)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DirectionIndex > 2)
        << "Unsupported Direction Index. DirectionIndex should be  less than 3 "
           "[ DirectionIndex = "
        << DirectionIndex << " ].\n";

    rOutput(0, DirectionIndex) += rInput[0] * rShapeFunctionDerivative;
    rOutput(1, DirectionIndex) += rInput[1] * rShapeFunctionDerivative;
    rOutput(2, DirectionIndex) += rInput[2] * rShapeFunctionDerivative;

    KRATOS_CATCH("");
}

template<>
void FluidCalculationUtilities::UpdateGradientValue(
    const array_1d<double, 3>& rInput,
    const double rShapeFunctionDerivative,
    const IndexType DirectionIndex,
    BoundedMatrix<double, 2, 2>& rOutput)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DirectionIndex > 1)
        << "Unsupported Direction Index. DirectionIndex should be  less than 2 "
           "[ DirectionIndex = "
        << DirectionIndex << " ].\n";

    rOutput(0, DirectionIndex) += rInput[0] * rShapeFunctionDerivative;
    rOutput(1, DirectionIndex) += rInput[1] * rShapeFunctionDerivative;

    KRATOS_CATCH("");
}

// template instantiations
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::AssignValue<double>(const double&, double&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::AssignValue<array_1d<double, 3>>(const array_1d<double, 3>&, array_1d<double, 3>&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::AssignValue<Vector>(const Vector&, Vector&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::AssignValue<Matrix>(const Matrix&, Matrix&);

template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::UpdateValue<array_1d<double, 3>>(const array_1d<double, 3>&, array_1d<double, 3>&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::UpdateValue<Vector>(const Vector&, Vector&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::UpdateValue<Matrix>(const Matrix&, Matrix&);

} // namespace Kratos
