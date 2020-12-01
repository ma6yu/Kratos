/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

#include "iga_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE(double, CROSS_AREA)
KRATOS_CREATE_VARIABLE(double, PRESTRESS_CAUCHY)

KRATOS_CREATE_VARIABLE(double, RAYLEIGH_ALPHA)
KRATOS_CREATE_VARIABLE(double, RAYLEIGH_BETA)

// 5p Director Shell Variables
KRATOS_CREATE_VARIABLE(bool, DIRECTOR_COMPUTED)
KRATOS_CREATE_VARIABLE(Vector, DIRECTOR)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DIRECTORINC)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MOMENTDIRECTORINC)
KRATOS_CREATE_VARIABLE(Matrix, DIRECTORTANGENTSPACE)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MOMENT_LINE_LOAD)

//Load Condition Variables
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DEAD_LOAD)

//Penalty Variables
KRATOS_CREATE_VARIABLE(double, PENALTY_FACTOR)

} // namespace Kratos
