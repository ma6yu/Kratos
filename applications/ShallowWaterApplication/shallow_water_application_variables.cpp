//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#include "shallow_water_application_variables.h"

namespace Kratos
{
    // Primary variables
    KRATOS_CREATE_VARIABLE(double, HEIGHT)
    KRATOS_CREATE_VARIABLE(double, FREE_SURFACE_ELEVATION)
    KRATOS_CREATE_VARIABLE(double, VERTICAL_VELOCITY)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(FLOW_RATE)

    // Physical variables
    KRATOS_CREATE_VARIABLE(double, BATHYMETRY)
    KRATOS_CREATE_VARIABLE(double, TOPOGRAPHY)
    KRATOS_CREATE_VARIABLE(double, RAIN)
    KRATOS_CREATE_VARIABLE(double, MANNING)
    KRATOS_CREATE_VARIABLE(double, PERMEABILITY)
    KRATOS_CREATE_VARIABLE(double, ATMOSPHERIC_PRESSURE)

    // Auxiliary variables
    KRATOS_CREATE_VARIABLE(double, SHOCK_STABILIZATION_FACTOR)
    KRATOS_CREATE_VARIABLE(double, EQUIVALENT_MANNING)
    KRATOS_CREATE_VARIABLE(double, DRY_HEIGHT)
    KRATOS_CREATE_VARIABLE(double, DRY_DISCHARGE_PENALTY)

    // Post-process variables
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(TOPOGRAPHY_GRADIENT)

    // Specific variables for PFEM2
    KRATOS_CREATE_VARIABLE( int, NUMBER_OF_PARTICLES)
    KRATOS_CREATE_VARIABLE( double, SUM_AREAS)
    KRATOS_CREATE_VARIABLE( double, PARTICLE_AREA)
    KRATOS_CREATE_VARIABLE( double, PARTICLE_WEIGHT)
    KRATOS_CREATE_VARIABLE( double, SUM_PARTICLES_WEIGHTS)
    KRATOS_CREATE_VARIABLE( double, MASS_WEIGHT)
    KRATOS_CREATE_VARIABLE( double, MEAN_SIZE)
    KRATOS_CREATE_VARIABLE( double, MEAN_VEL_OVER_ELEM_SIZE)
    KRATOS_CREATE_VARIABLE( double, PROJECTED_SCALAR1)
    KRATOS_CREATE_VARIABLE( double, DELTA_SCALAR1)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VECTOR1)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DELTA_VECTOR1)
    KRATOS_CREATE_VARIABLE( Element::Pointer, CURRENT_ELEMENT)

    // Variables for Flux Corrected Transport algorithm
    KRATOS_CREATE_VARIABLE(Vector, POSITIVE_FLUX)
    KRATOS_CREATE_VARIABLE(Vector, NEGATIVE_FLUX)
    KRATOS_CREATE_VARIABLE(double, POSITIVE_RATIO)
    KRATOS_CREATE_VARIABLE(double, NEGATIVE_RATIO)

    // Benchark variables
    KRATOS_CREATE_VARIABLE( double, EXACT_HEIGHT)
    KRATOS_CREATE_VARIABLE( double, HEIGHT_ERROR)
    KRATOS_CREATE_VARIABLE( double, EXACT_FREE_SURFACE)
    KRATOS_CREATE_VARIABLE( double, FREE_SURFACE_ERROR)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(EXACT_VELOCITY)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_ERROR)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(EXACT_MOMENTUM)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MOMENTUM_ERROR)

    // Units conversion
    KRATOS_CREATE_VARIABLE( double, TIME_UNIT_CONVERTER)
    KRATOS_CREATE_VARIABLE( double, WATER_HEIGHT_UNIT_CONVERTER)
}
