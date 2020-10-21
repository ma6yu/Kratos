// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi
//

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"
#include "containers/model.h"

// Application includes
#include "structural_mechanics_application_variables.h"

// Constitutive law
#include "custom_advanced_constitutive/small_strain_j2_plasticity_3d.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos
{
namespace Testing
{
// We test the associated plasticity Constitutive laws...
typedef Node<3> NodeType;

// Check the correct calculation of the integrated stress with the CL's in small strain
KRATOS_TEST_CASE_IN_SUITE(_ConstitutiveLaw_SmallStrainJ2Plasticity3D, KratosStructuralMechanicsFastSuite)
{
    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    Vector stress_vector(6), strain_vector(6);
    Matrix const_matrix;
    // Create gauss point
    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("Main");
    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    Tetrahedra3D4<NodeType> geometry = Tetrahedra3D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);
    // Set material properties
    material_properties.SetValue(YOUNG_MODULUS, 21000);
    material_properties.SetValue(POISSON_RATIO, 0.3);
    material_properties.SetValue(YIELD_STRESS, 5.5);
    material_properties.SetValue(ISOTROPIC_HARDENING_MODULUS, 0.12924);
    material_properties.SetValue(EXPONENTIAL_SATURATION_YIELD_STRESS, 0.0);
    material_properties.SetValue(HARDENING_EXPONENT, 1.0);
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=cl_parameters.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    // Set required constitutive law parameters:
    cl_parameters.SetElementGeometry(geometry);
    cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetConstitutiveMatrix(const_matrix);
    // Create the CL
    SmallStrainJ2Plasticity3D cl = SmallStrainJ2Plasticity3D();

    // Set variables for the test
    const double tolerance = 1.0e-3;
    std::size_t nr_ts = 10;  // timesteps tested
    Vector ser(nr_ts), dvr(nr_ts);  // reference strain energy, damage variable
    Matrix epr(nr_ts, 6), str(nr_ts, 6);  // reference strain ("epsilon"), stress
    Matrix cmr(nr_ts, 6 * 6);  // reference constitutive matrix

    //
    // Test: check correct behavior of internal and calculated variables
    //
    KRATOS_CHECK_IS_FALSE(cl.Has(STRAIN_ENERGY));  // = False, in order to use CalculateValue())
    KRATOS_CHECK_IS_FALSE(cl.Has(STRAIN));  // = False, in order to use CalculateValue())
    KRATOS_CHECK(cl.Has(ACCUMULATED_PLASTIC_STRAIN));  // = True
    KRATOS_CHECK(cl.Has(INTERNAL_VARIABLES));  // = True
    Vector internal_variables_w(7);
    internal_variables_w[0] = 0.123;
    internal_variables_w[1] = 0.234;
    internal_variables_w[2] = 0.345;
    internal_variables_w[3] = 0.456;
    internal_variables_w[4] = 0.789;
    internal_variables_w[5] = 0.890;
    internal_variables_w[6] = 0.901;
    cl.SetValue(INTERNAL_VARIABLES, internal_variables_w, test_model_part.GetProcessInfo());
    Vector internal_variables_r;  // CL should internally resize it to 7
    cl.GetValue(INTERNAL_VARIABLES, internal_variables_r);
    KRATOS_CHECK_NEAR(internal_variables_r.size(), 7., tolerance);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[0], 0.123, tolerance);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[1], 0.234, tolerance);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[2], 0.345, tolerance);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[3], 0.456, tolerance);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[4], 0.789, tolerance);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[5], 0.890, tolerance);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[6], 0.901, tolerance);  // = True

    //
    // Test: load - unload in traction
    //
    epr(0,0)=0.0001; epr(0,1)=0.0001; epr(0,2)=0.0000; epr(0,3)=0.0001; epr(0,4)=0.0000; epr(0,5)=0.0001;
    epr(1,0)=0.0002; epr(1,1)=0.0002; epr(1,2)=0.0000; epr(1,3)=0.0002; epr(1,4)=0.0000; epr(1,5)=0.0002;
    epr(2,0)=0.0003; epr(2,1)=0.0003; epr(2,2)=0.0000; epr(2,3)=0.0003; epr(2,4)=0.0000; epr(2,5)=0.0003;
    epr(3,0)=0.0004; epr(3,1)=0.0004; epr(3,2)=0.0000; epr(3,3)=0.0004; epr(3,4)=0.0000; epr(3,5)=0.0004;
    epr(4,0)=0.0005; epr(4,1)=0.0005; epr(4,2)=0.0000; epr(4,3)=0.0005; epr(4,4)=0.0000; epr(4,5)=0.0005;
    epr(5,0)=0.0006; epr(5,1)=0.0006; epr(5,2)=0.0000; epr(5,3)=0.0006; epr(5,4)=0.0000; epr(5,5)=0.0006;
    epr(6,0)=0.0007; epr(6,1)=0.0007; epr(6,2)=0.0000; epr(6,3)=0.0007; epr(6,4)=0.0000; epr(6,5)=0.0007;
    epr(7,0)=0.0008; epr(7,1)=0.0008; epr(7,2)=0.0000; epr(7,3)=0.0008; epr(7,4)=0.0000; epr(7,5)=0.0008;
    epr(8,0)=0.0009; epr(8,1)=0.0009; epr(8,2)=0.0000; epr(8,3)=0.0009; epr(8,4)=0.0000; epr(8,5)=0.0009;
    epr(9,0)=0.0010; epr(9,1)=0.0010; epr(9,2)=0.0000; epr(9,3)=0.0010; epr(9,4)=0.0000; epr(9,5)=0.0010;

    str(0,0) = 4.03846; str(0,1) = 4.03846; str(0,2) = 2.42308; str(0,3) = 0.80769; str(0,4) = 0.0; str(0,5) = 0.80769;
    str(1,0) = 8.07692; str(1,1) = 8.07692; str(1,2) = 4.84615; str(1,3) = 1.61538; str(1,4) = 0.0; str(1,5) = 1.61538;
    str(2,0) = 11.6595; str(2,1) = 11.6595; str(2,2) = 8.18099; str(2,3) = 1.73926; str(2,4) = 0.0; str(2,5) = 1.73926;
    str(3,0) = 15.1595; str(3,1) = 15.1595; str(3,2) = 11.681 ; str(3,3) = 1.73926; str(3,4) = 0.0; str(3,5) = 1.73926;
    str(4,0) = 18.6595; str(4,1) = 18.6595; str(4,2) = 15.181 ; str(4,3) = 1.73926; str(4,4) = 0.0; str(4,5) = 1.73926;
    str(5,0) = 22.1595; str(5,1) = 22.1595; str(5,2) = 18.681 ; str(5,3) = 1.73927; str(5,4) = 0.0; str(5,5) = 1.73927;
    str(6,0) = 25.6595; str(6,1) = 25.6595; str(6,2) = 22.181 ; str(6,3) = 1.73927; str(6,4) = 0.0; str(6,5) = 1.73927;
    str(7,0) = 29.1595; str(7,1) = 29.1595; str(7,2) = 25.681 ; str(7,3) = 1.73928; str(7,4) = 0.0; str(7,5) = 1.73928;
    str(8,0) = 32.6595; str(8,1) = 32.6595; str(8,2) = 29.181 ; str(8,3) = 1.73928; str(8,4) = 0.0; str(8,5) = 1.73928;
    str(9,0) = 36.1595; str(9,1) = 36.1595; str(9,2) = 32.681 ; str(9,3) = 1.73929; str(9,4) = 0.0; str(9,5) = 1.73929;

    cmr(0, 0)=2.014743e+03; cmr(0, 1)=1.350944e+01; cmr(0, 2)=6.084757e+02; cmr(0, 3)=-2.974831e+02; cmr(0, 4)=0.000000e+00; cmr(0, 5)=0.000000e+00;
    cmr(0, 6)=1.350944e+01; cmr(0, 7)=2.014743e+03; cmr(0, 8)=6.084757e+02; cmr(0, 9)=-2.974831e+02; cmr(0,10)=0.000000e+00; cmr(0,11)=0.000000e+00;
    cmr(0,12)=6.084757e+02; cmr(0,13)=6.084757e+02; cmr(0,14)=2.966689e+03; cmr(0,15)=-1.784899e+02; cmr(0,16)=0.000000e+00; cmr(0,17)=0.000000e+00;
    cmr(0,18)=-2.974831e+02; cmr(0,19)=-2.974831e+02; cmr(0,20)=-1.784899e+02; cmr(0,21)=9.411201e+02; cmr(0,22)=0.000000e+00; cmr(0,23)=0.000000e+00;
    cmr(0,24)=0.000000e+00; cmr(0,25)=0.000000e+00; cmr(0,26)=0.000000e+00; cmr(0,27)=0.000000e+00; cmr(0,28)=1.000617e+03; cmr(0,29)=0.000000e+00;
    cmr(0,30)=0.000000e+00; cmr(0,31)=0.000000e+00; cmr(0,32)=0.000000e+00; cmr(0,33)=0.000000e+00; cmr(0,34)=0.000000e+00; cmr(0,35)=1.000617e+03;
    cmr(1, 0)=1.007371e+03; cmr(1, 1)=6.754721e+00; cmr(1, 2)=3.042379e+02; cmr(1, 3)=-1.487416e+02; cmr(1, 4)=0.000000e+00; cmr(1, 5)=0.000000e+00;
    cmr(1, 6)=6.754721e+00; cmr(1, 7)=1.007371e+03; cmr(1, 8)=3.042379e+02; cmr(1, 9)=-1.487416e+02; cmr(1,10)=0.000000e+00; cmr(1,11)=0.000000e+00;
    cmr(1,12)=3.042379e+02; cmr(1,13)=3.042379e+02; cmr(1,14)=1.483344e+03; cmr(1,15)=-8.924494e+01; cmr(1,16)=0.000000e+00; cmr(1,17)=0.000000e+00;
    cmr(1,18)=-1.487416e+02; cmr(1,19)=-1.487416e+02; cmr(1,20)=-8.924494e+01; cmr(1,21)=4.705601e+02; cmr(1,22)=0.000000e+00; cmr(1,23)=0.000000e+00;
    cmr(1,24)=0.000000e+00; cmr(1,25)=0.000000e+00; cmr(1,26)=0.000000e+00; cmr(1,27)=0.000000e+00; cmr(1,28)=5.003084e+02; cmr(1,29)=0.000000e+00;
    cmr(1,30)=0.000000e+00; cmr(1,31)=0.000000e+00; cmr(1,32)=0.000000e+00; cmr(1,33)=0.000000e+00; cmr(1,34)=0.000000e+00; cmr(1,35)=5.003084e+02;
    cmr(2, 0)=7.262499e+02; cmr(2, 1)=1.837572e+02; cmr(2, 2)=2.730021e+02; cmr(2, 3)=-4.462247e+01; cmr(2, 4)=0.000000e+00; cmr(2, 5)=0.000000e+00;
    cmr(2, 6)=1.837572e+02; cmr(2, 7)=7.262499e+02; cmr(2, 8)=2.730021e+02; cmr(2, 9)=-4.462247e+01; cmr(2,10)=0.000000e+00; cmr(2,11)=0.000000e+00;
    cmr(2,12)=2.730021e+02; cmr(2,13)=2.730021e+02; cmr(2,14)=8.690418e+02; cmr(2,15)=-2.677348e+01; cmr(2,16)=0.000000e+00; cmr(2,17)=0.000000e+00;
    cmr(2,18)=-4.462247e+01; cmr(2,19)=-4.462247e+01; cmr(2,20)=-2.677348e+01; cmr(2,21)=2.623219e+02; cmr(2,22)=0.000000e+00; cmr(2,23)=0.000000e+00;
    cmr(2,24)=0.000000e+00; cmr(2,25)=0.000000e+00; cmr(2,26)=0.000000e+00; cmr(2,27)=0.000000e+00; cmr(2,28)=2.712464e+02; cmr(2,29)=0.000000e+00;
    cmr(2,30)=0.000000e+00; cmr(2,31)=0.000000e+00; cmr(2,32)=0.000000e+00; cmr(2,33)=0.000000e+00; cmr(2,34)=0.000000e+00; cmr(2,35)=2.712464e+02;
    cmr(3, 0)=9.493622e+02; cmr(3, 1)=4.068695e+02; cmr(3, 2)=4.068695e+02; cmr(3, 3)=0.000000e+00; cmr(3, 4)=0.000000e+00; cmr(3, 5)=0.000000e+00;
    cmr(3, 6)=4.068695e+02; cmr(3, 7)=9.493622e+02; cmr(3, 8)=4.068695e+02; cmr(3, 9)=0.000000e+00; cmr(3,10)=0.000000e+00; cmr(3,11)=0.000000e+00;
    cmr(3,12)=4.068695e+02; cmr(3,13)=4.068695e+02; cmr(3,14)=9.493622e+02; cmr(3,15)=0.000000e+00; cmr(3,16)=0.000000e+00; cmr(3,17)=0.000000e+00;
    cmr(3,18)=0.000000e+00; cmr(3,19)=0.000000e+00; cmr(3,20)=0.000000e+00; cmr(3,21)=2.712464e+02; cmr(3,22)=0.000000e+00; cmr(3,23)=0.000000e+00;
    cmr(3,24)=0.000000e+00; cmr(3,25)=0.000000e+00; cmr(3,26)=0.000000e+00; cmr(3,27)=0.000000e+00; cmr(3,28)=2.712464e+02; cmr(3,29)=0.000000e+00;
    cmr(3,30)=0.000000e+00; cmr(3,31)=0.000000e+00; cmr(3,32)=0.000000e+00; cmr(3,33)=0.000000e+00; cmr(3,34)=0.000000e+00; cmr(3,35)=2.712464e+02;

    ser[0]=5.503392e-05; ser[1]=2.476526e-04; ser[2]=1.491855e-03; ser[3]=1.491855e-03;
    dvr[0]=1.327988e-01; dvr[1]=5.663994e-01; dvr[2]=7.649198e-01; dvr[3]=7.649198e-01;


      // Simulate the call sequence of the element
      Vector dummy;
      cl.InitializeMaterial(material_properties, geometry, dummy);
      for (std::size_t t = 0; t < nr_ts; ++t){
          for (std::size_t comp = 0; comp < 6; ++comp) {
              strain_vector[comp] = epr(t, comp);
          }
          if (cl.RequiresInitializeMaterialResponse()){
              cl.InitializeMaterialResponseCauchy(cl_parameters);
          }
          cl.CalculateMaterialResponseCauchy(cl_parameters);
          if (cl.RequiresFinalizeMaterialResponse()){
              cl.FinalizeMaterialResponseCauchy(cl_parameters);
          }
          double value;
  
//          // Check damage variable
//          cl.CalculateValue(cl_parameters, DAMAGE_VARIABLE, value);
//          // TODO(marandra): NAN values are not handled correctly by KRATOS CHECK functions
//          if (std::isnan(dvr[t]/value)){
//              KRATOS_CHECK_NEAR(dvr[t], value, tolerance);
//          } else {
//              KRATOS_CHECK_NEAR(dvr[t]/value, 1, tolerance);
//          }
//  
//          // Check strain energy
//          cl.CalculateValue(cl_parameters, STRAIN_ENERGY, value);
//          if (std::isnan(ser[t]/value)){
//              KRATOS_CHECK_NEAR(ser[t], value, tolerance);
//          } else {
//              KRATOS_CHECK_NEAR(ser[t]/value, 1, tolerance);
//          }
//  
//          // Check stress
//          for (std::size_t comp = 0; comp < 6; ++comp){
//              KRATOS_CHECK_IS_FALSE(std::isnan(stress_vector[comp]));
//              if (std::isnan(stress_vector[comp]/str(t, comp))){
//                  KRATOS_CHECK_NEAR(stress_vector[comp], str(t, comp), tolerance);
//              } else {
//                  KRATOS_CHECK_NEAR(stress_vector[comp]/str(t, comp), 1, tolerance);
//              }
//          }
  
          // Check constitutive tensor
  //        for (std::size_t i = 0; i < 6; ++i){
  //            for (std::size_t j = 0; j < 6; ++j){
  //                std::size_t idx = i * 6 + j;
  //                KRATOS_CHECK_IS_FALSE(std::isnan(const_matrix(i, j)));
  //                if (std::isnan(const_matrix(i, j)/cmr(t, idx))){
  //                    KRATOS_CHECK_NEAR(const_matrix(i, j), cmr(t, idx), tolerance);
  //                } else {
  //                    KRATOS_CHECK_NEAR(const_matrix(i, j)/cmr(t, idx), 1, tolerance);
  //                    }
  //            }
  //        }
      }
  

//    //
//    // Test: load - unload in compression
//    //
//    epr(0,0)=-1.000000e-04; epr(0,1)=-1.000000e-04; epr(0,2)=0.000000e+00; epr(0,3)=-1.000000e-04; epr(0,4)=0.000000e+00; epr(0,5)=0.000000e+00;
//    epr(1,0)=-3.000000e-04; epr(1,1)=-3.000000e-04; epr(1,2)=0.000000e+00; epr(1,3)=-3.000000e-04; epr(1,4)=0.000000e+00; epr(1,5)=0.000000e+00;
//    epr(2,0)=-1.000000e-03; epr(2,1)=-1.000000e-03; epr(2,2)=0.000000e+00; epr(2,3)=-1.000000e-03; epr(2,4)=0.000000e+00; epr(2,5)=0.000000e+00;
//    epr(3,0)=1.000000e-03; epr(3,1)=1.000000e-03; epr(3,2)=0.000000e+00; epr(3,3)=1.000000e-03; epr(3,4)=0.000000e+00; epr(3,5)=0.000000e+00;
//
//    str(0,0)=-5.003084e-01; str(0,1)=-5.003084e-01; str(0,2)=-3.001850e-01; str(0,3)=-1.000617e-01; str(0,4)=0.000000e+00; str(0,5)=0.000000e+00;
//    str(1,0)=-7.504626e-01; str(1,1)=-7.504626e-01; str(1,2)=-4.502775e-01; str(1,3)=-1.500925e-01; str(1,4)=0.000000e+00; str(1,5)=0.000000e+00;
//    str(2,0)=-1.356232e+00; str(2,1)=-1.356232e+00; str(2,2)=-8.137391e-01; str(2,3)=-2.712464e-01; str(2,4)=0.000000e+00; str(2,5)=0.000000e+00;
//    str(3,0)=1.356232e+00; str(3,1)=1.356232e+00; str(3,2)=8.137391e-01; str(3,3)=2.712464e-01; str(3,4)=0.000000e+00; str(3,5)=0.000000e+00;
//
//    cmr(0, 0)=2.014743e+03; cmr(0, 1)=1.350944e+01; cmr(0, 2)=6.084757e+02; cmr(0, 3)=-2.974831e+02; cmr(0, 4)=0.000000e+00; cmr(0, 5)=0.000000e+00;
//    cmr(0, 6)=1.350944e+01; cmr(0, 7)=2.014743e+03; cmr(0, 8)=6.084757e+02; cmr(0, 9)=-2.974831e+02; cmr(0,10)=0.000000e+00; cmr(0,11)=0.000000e+00;
//    cmr(0,12)=6.084757e+02; cmr(0,13)=6.084757e+02; cmr(0,14)=2.966689e+03; cmr(0,15)=-1.784899e+02; cmr(0,16)=0.000000e+00; cmr(0,17)=0.000000e+00;
//    cmr(0,18)=-2.974831e+02; cmr(0,19)=-2.974831e+02; cmr(0,20)=-1.784899e+02; cmr(0,21)=9.411201e+02; cmr(0,22)=0.000000e+00; cmr(0,23)=0.000000e+00;
//    cmr(0,24)=0.000000e+00; cmr(0,25)=0.000000e+00; cmr(0,26)=0.000000e+00; cmr(0,27)=0.000000e+00; cmr(0,28)=1.000617e+03; cmr(0,29)=0.000000e+00;
//    cmr(0,30)=0.000000e+00; cmr(0,31)=0.000000e+00; cmr(0,32)=0.000000e+00; cmr(0,33)=0.000000e+00; cmr(0,34)=0.000000e+00; cmr(0,35)=1.000617e+03;
//    cmr(1, 0)=1.007371e+03; cmr(1, 1)=6.754721e+00; cmr(1, 2)=3.042379e+02; cmr(1, 3)=-1.487416e+02; cmr(1, 4)=0.000000e+00; cmr(1, 5)=0.000000e+00;
//    cmr(1, 6)=6.754721e+00; cmr(1, 7)=1.007371e+03; cmr(1, 8)=3.042379e+02; cmr(1, 9)=-1.487416e+02; cmr(1,10)=0.000000e+00; cmr(1,11)=0.000000e+00;
//    cmr(1,12)=3.042379e+02; cmr(1,13)=3.042379e+02; cmr(1,14)=1.483344e+03; cmr(1,15)=-8.924494e+01; cmr(1,16)=0.000000e+00; cmr(1,17)=0.000000e+00;
//    cmr(1,18)=-1.487416e+02; cmr(1,19)=-1.487416e+02; cmr(1,20)=-8.924494e+01; cmr(1,21)=4.705601e+02; cmr(1,22)=0.000000e+00; cmr(1,23)=0.000000e+00;
//    cmr(1,24)=0.000000e+00; cmr(1,25)=0.000000e+00; cmr(1,26)=0.000000e+00; cmr(1,27)=0.000000e+00; cmr(1,28)=5.003084e+02; cmr(1,29)=0.000000e+00;
//    cmr(1,30)=0.000000e+00; cmr(1,31)=0.000000e+00; cmr(1,32)=0.000000e+00; cmr(1,33)=0.000000e+00; cmr(1,34)=0.000000e+00; cmr(1,35)=5.003084e+02;
//    cmr(2, 0)=7.262499e+02; cmr(2, 1)=1.837572e+02; cmr(2, 2)=2.730021e+02; cmr(2, 3)=-4.462247e+01; cmr(2, 4)=0.000000e+00; cmr(2, 5)=0.000000e+00;
//    cmr(2, 6)=1.837572e+02; cmr(2, 7)=7.262499e+02; cmr(2, 8)=2.730021e+02; cmr(2, 9)=-4.462247e+01; cmr(2,10)=0.000000e+00; cmr(2,11)=0.000000e+00;
//    cmr(2,12)=2.730021e+02; cmr(2,13)=2.730021e+02; cmr(2,14)=8.690418e+02; cmr(2,15)=-2.677348e+01; cmr(2,16)=0.000000e+00; cmr(2,17)=0.000000e+00;
//    cmr(2,18)=-4.462247e+01; cmr(2,19)=-4.462247e+01; cmr(2,20)=-2.677348e+01; cmr(2,21)=2.623219e+02; cmr(2,22)=0.000000e+00; cmr(2,23)=0.000000e+00;
//    cmr(2,24)=0.000000e+00; cmr(2,25)=0.000000e+00; cmr(2,26)=0.000000e+00; cmr(2,27)=0.000000e+00; cmr(2,28)=2.712464e+02; cmr(2,29)=0.000000e+00;
//    cmr(2,30)=0.000000e+00; cmr(2,31)=0.000000e+00; cmr(2,32)=0.000000e+00; cmr(2,33)=0.000000e+00; cmr(2,34)=0.000000e+00; cmr(2,35)=2.712464e+02;
//    cmr(3, 0)=9.493622e+02; cmr(3, 1)=4.068695e+02; cmr(3, 2)=4.068695e+02; cmr(3, 3)=0.000000e+00; cmr(3, 4)=0.000000e+00; cmr(3, 5)=0.000000e+00;
//    cmr(3, 6)=4.068695e+02; cmr(3, 7)=9.493622e+02; cmr(3, 8)=4.068695e+02; cmr(3, 9)=0.000000e+00; cmr(3,10)=0.000000e+00; cmr(3,11)=0.000000e+00;
//    cmr(3,12)=4.068695e+02; cmr(3,13)=4.068695e+02; cmr(3,14)=9.493622e+02; cmr(3,15)=0.000000e+00; cmr(3,16)=0.000000e+00; cmr(3,17)=0.000000e+00;
//    cmr(3,18)=0.000000e+00; cmr(3,19)=0.000000e+00; cmr(3,20)=0.000000e+00; cmr(3,21)=2.712464e+02; cmr(3,22)=0.000000e+00; cmr(3,23)=0.000000e+00;
//    cmr(3,24)=0.000000e+00; cmr(3,25)=0.000000e+00; cmr(3,26)=0.000000e+00; cmr(3,27)=0.000000e+00; cmr(3,28)=2.712464e+02; cmr(3,29)=0.000000e+00;
//    cmr(3,30)=0.000000e+00; cmr(3,31)=0.000000e+00; cmr(3,32)=0.000000e+00; cmr(3,33)=0.000000e+00; cmr(3,34)=0.000000e+00; cmr(3,35)=2.712464e+02;
//
//    ser[0]=5.503392e-05; ser[1]=2.476526e-04; ser[2]=1.491855e-03; ser[3]=1.491855e-03;
//    dvr[0]=1.327988e-01; dvr[1]=5.663994e-01; dvr[2]=7.649198e-01; dvr[3]=7.649198e-01;
//
//    // Here we must simulate the call sequence of the element
//    cl.InitializeMaterial(material_properties, geometry, dummy);
//    for (std::size_t t = 0; t < nr_ts; ++t){
//        for (std::size_t comp = 0; comp < 6; ++comp) {
//            strain_vector[comp] = epr(t, comp);
//        }
//        cl.InitializeMaterialResponseCauchy(cl_parameters);
//        cl.CalculateMaterialResponseCauchy(cl_parameters);
//        cl.FinalizeMaterialResponseCauchy(cl_parameters);
//        double value;
//
//        // Check damage variable
//        cl.CalculateValue(cl_parameters, DAMAGE_VARIABLE, value);
//        if (std::isnan(dvr[t]/value)){
//            KRATOS_CHECK_NEAR(dvr[t], value, tolerance);
//        } else {
//            KRATOS_CHECK_NEAR(dvr[t]/value, 1, tolerance);
//        }
//
//        // Check strain energy
//        cl.CalculateValue(cl_parameters, STRAIN_ENERGY, value);
//        if (std::isnan(ser[t]/value)){
//            KRATOS_CHECK_NEAR(ser[t], value, tolerance);
//        } else {
//            KRATOS_CHECK_NEAR(ser[t]/value, 1, tolerance);
//        }
//
//        // Check stress
//        for (std::size_t comp = 0; comp < 6; ++comp){
//            KRATOS_CHECK_IS_FALSE(std::isnan(stress_vector[comp]));
//            if (std::isnan(stress_vector[comp]/str(t, comp))){
//                KRATOS_CHECK_NEAR(stress_vector[comp], str(t, comp), tolerance);
//            } else {
//                KRATOS_CHECK_NEAR(stress_vector[comp]/str(t, comp), 1, tolerance);
//            }
//        }
//
//        // Check constitutive tensor
//        for (std::size_t i = 0; i < 6; ++i){
//            for (std::size_t j = 0; j < 6; ++j){
//                std::size_t idx = i * 6 + j;
//                KRATOS_CHECK_IS_FALSE(std::isnan(const_matrix(i, j)));
//                if (std::isnan(const_matrix(i, j)/cmr(t, idx))){
//                    KRATOS_CHECK_NEAR(const_matrix(i, j), cmr(t, idx), tolerance);
//                } else {
//                    KRATOS_CHECK_NEAR(const_matrix(i, j)/cmr(t, idx), 1, tolerance);
//                }
//            }
//        }
//    }


}

} // namespace Testing
} // namespace Kratos
