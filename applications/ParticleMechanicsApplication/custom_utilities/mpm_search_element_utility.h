//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#ifndef KRATOS_MPM_SEARCH_ELEMENT_UTILITY
#define KRATOS_MPM_SEARCH_ELEMENT_UTILITY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/quadrature_points_utility.h"

#include "particle_mechanics_application_variables.h"

namespace Kratos
{
namespace MPMSearchElementUtility
{

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef typename ModelPart::GeometryType GeometryType;

    void SearchNeighbourGeometry(ModelPart& rBackgroundGridModelPart)
    {
        std::cout << "SearchNeighbourGeometry" << std::endl;
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rBackgroundGridModelPart.NumberOfElements()); ++i)
        {
            auto& r_geometry = (rBackgroundGridModelPart.ElementsBegin() + i)->GetGeometry();
            //GeometryType& r_geometry = rBackgroundGridModelPart.Elements()[i].GetGeometry();
            std::vector<typename Geometry<Node<3>>::Pointer> geometry_neighbours;
            for (IndexType j = 0; j < rBackgroundGridModelPart.NumberOfElements(); j++)
            {
                auto p_geometry_neighbour = (rBackgroundGridModelPart.ElementsBegin() + j)->pGetGeometry();
                //auto p_geometry_neighbour = rBackgroundGridModelPart.Elements()[j].pGetGeometry();
                for (IndexType n = 0; n < p_geometry_neighbour->size(); n++)
                {
                    for (IndexType k = 0; k < r_geometry.size(); k++)
                    {
                        if (r_geometry[k].Id() == (*p_geometry_neighbour)[n].Id()) {
                            geometry_neighbours.push_back(p_geometry_neighbour);
                            //break;
                        }
                    }
                }
            }
            r_geometry.SetValue(GEOMETRY_NEIGHBOURS, geometry_neighbours);
        }
    }


    /**
     * @brief Search element connectivity for each particle
     * @details A search is performed to know in which grid element the material point falls.
     * If one or more material points fall in the grid element, the grid element is
     * set to be active and its connectivity is associated to the material point
     * element.
     * STEPS:
     * 1) All the elements are set to be INACTIVE
     * 2) A searching is performed and the grid elements which contain at least a MP are set to be ACTIVE
     *
     */
    template<std::size_t TDimension>
    void SearchElement(ModelPart& rBackgroundGridModelPart, ModelPart& rMPMModelPart, const std::size_t MaxNumberOfResults,
        const double Tolerance)
    {
        const ProcessInfo& r_process_info = rBackgroundGridModelPart.GetProcessInfo();
        const bool is_explicit = r_process_info.GetValue(IS_EXPLICIT);

        // Reset elements to inactive
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rBackgroundGridModelPart.Elements().size()); ++i) {
            auto element_itr = rBackgroundGridModelPart.Elements().begin() + i;
            auto& r_geometry = element_itr->GetGeometry();
            element_itr->Reset(ACTIVE);

            for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                r_geometry[j].Reset(ACTIVE);

        }

        if (!rBackgroundGridModelPart.Elements().begin()->GetGeometry().Has(GEOMETRY_NEIGHBOURS)) SearchNeighbourGeometry(rBackgroundGridModelPart);



        bool all_found = true;
        bool neighbour_search = true;
        if (neighbour_search) {
            //#pragma omp for
            for (int i = 0; i < static_cast<int>(rMPMModelPart.Elements().size()); ++i) {
                auto element_itr = rMPMModelPart.Elements().begin() + i;

                std::vector<array_1d<double, 3>> xg;
                element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rMPMModelPart.GetProcessInfo());

                GeometryType& r_parent_geometry = element_itr->GetGeometry().GetGeometryParent(0);

                array_1d<double, 3> local_coordinates;
                array_1d<double, 3> local_coordinates_2;
                bool is_found = r_parent_geometry.IsInside(xg[0], local_coordinates_2, Tolerance);
                r_parent_geometry.PointLocalCoordinates(local_coordinates, xg[0]);
                if (!is_found)
                {
                    auto& geometry_neighbours = r_parent_geometry.GetValue(GEOMETRY_NEIGHBOURS);
                    for (IndexType k = 0; k < geometry_neighbours.size(); k++)
                    {
                        if (geometry_neighbours[k]->IsInside(xg[0], local_coordinates, Tolerance))
                        {
                            is_found = true;
                            auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromLocalCoordinates(
                                *(geometry_neighbours[k].get()), local_coordinates,
                                element_itr->GetGeometry().IntegrationPoints()[0].Weight());

                            // Update geometry of particle element
                            element_itr->SetGeometry(p_new_geometry);

                            for (IndexType j = 0; j < geometry_neighbours[k]->PointsNumber(); ++j)
                                geometry_neighbours[k]->Points()[j].Set(ACTIVE);
                        }
                        else
                        {
                            all_found = false;
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Tobi search failed! " << std::endl;
                            // bin search
                        }
                    }
                }
                else {
                    //pelem->Set(ACTIVE);

                    auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromLocalCoordinates(
                        r_parent_geometry, local_coordinates,
                        element_itr->GetGeometry().IntegrationPoints()[0].Weight());

                    // Update geometry of particle element
                    element_itr->SetGeometry(p_new_geometry);



                    for (IndexType j = 0; j < r_parent_geometry.PointsNumber(); ++j)
                        r_parent_geometry[j].Set(ACTIVE);
                }
                if (!is_found) {
                    KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Search Element for Material Point: " << element_itr->Id()
                        << " is failed. Geometry is cleared." << std::endl;

                    element_itr->GetGeometry().clear();
                    element_itr->Reset(ACTIVE);
                    element_itr->Set(TO_ERASE);
                }
            }
        }






        // Search background grid and make element active
        Vector N;
        const int max_result = 1000;

        if (rMPMModelPart.NumberOfConditions() > 0 || !neighbour_search || !all_found)
        {


        //#pragma omp parallel
        {
            BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
            SearchStructure.UpdateSearchDatabase();

            typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(max_result);
            if (!neighbour_search || !all_found)
            {
                KRATOS_WATCH("FALLING BACK TO NORMAL SEARCH")
            // Element search and assign background grid
            //#pragma omp for
            for (int i = 0; i < static_cast<int>(rMPMModelPart.Elements().size()); ++i) {
                auto element_itr = rMPMModelPart.Elements().begin() + i;

                std::vector<array_1d<double, 3>> xg;
                element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rMPMModelPart.GetProcessInfo());
                typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();

                Element::Pointer pelem;

                // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);


                if (is_found && is_explicit) {
                    // check if MP is exactly on the edge of the element, this gives spurious strains in explicit
                    bool isOnEdge = false;
                    for (SizeType i = 0; i < N.size(); ++i) {
                        if (std::abs(N[i]) < std::numeric_limits<double>::epsilon()) {
                            isOnEdge = true;
                            break;
                        }
                    }
                    if (isOnEdge) {
                        // MP is exactly on the edge. Now we give it a little 'nudge'
                        array_1d<double, 3> xg_nudged = array_1d<double, 3>(xg[0]);
                        const double& delta_time = r_process_info[DELTA_TIME];
                        std::vector<array_1d<double, 3>> mp_vel;
                        element_itr->CalculateOnIntegrationPoints(MP_VELOCITY, mp_vel, rMPMModelPart.GetProcessInfo());
                        array_1d<double, 3> nudge_displacement = delta_time / 1000.0 * mp_vel[0];
                        xg_nudged += nudge_displacement;
                        is_found = SearchStructure.FindPointOnMesh(xg_nudged, N, pelem, result_begin, MaxNumberOfResults, Tolerance);
                        // check if the nudged point is found...
                        if (is_found) {
                            // store the nudged MP position
                            element_itr->SetValuesOnIntegrationPoints(MP_COORD, { xg_nudged }, rMPMModelPart.GetProcessInfo());
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: To prevent spurious explicit stresses, Material Point " << element_itr->Id()
                                << " was nudged by " << nudge_displacement << std::endl;
                        }
                        else {
                            // find the un-nudged MP again
                            is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Material Point " << element_itr->Id()
                                << " lies exactly on an element edge and may give spurious results." << std::endl;
                        }
                    }
                }


                if (is_found == true) {
                    pelem->Set(ACTIVE);

                    auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                        pelem->pGetGeometry(), xg[0],
                        element_itr->GetGeometry().IntegrationPoints()[0].Weight());

                    // Update geometry of particle element
                    element_itr->SetGeometry(p_new_geometry);

                    for (IndexType j = 0; j < p_new_geometry->PointsNumber(); ++j)
                        (*p_new_geometry)[j].Set(ACTIVE);
                }
                else {
                    KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Search Element for Material Point: " << element_itr->Id()
                        << " is failed. Geometry is cleared." << std::endl;

                    element_itr->GetGeometry().clear();
                    element_itr->Reset(ACTIVE);
                    element_itr->Set(TO_ERASE);
                }
            }
        }
            // Condition search and assign background grid
            #pragma omp for
            for (int i = 0; i < static_cast<int>(rMPMModelPart.Conditions().size()); ++i) {

                auto condition_itr = rMPMModelPart.Conditions().begin() + i;
                std::vector<array_1d<double, 3>> xg;
                condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());

                if (xg.size() == 1) {
                    // Only search for particle based BCs!
                    // Grid BCs are still applied on MP_model_part but we don't want to search for them.
                    typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();

                    Element::Pointer pelem;

                    // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                    bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                    if (is_found == true) {
                        //pelem->Set(ACTIVE);
                        condition_itr->GetGeometry() = pelem->GetGeometry();
                        auto& r_geometry = condition_itr->GetGeometry();

                        for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                            r_geometry[j].Set(ACTIVE);
                    }
                    else {
                        KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Search Element for Material Point Condition: " << condition_itr->Id()
                            << " is failed. Geometry is cleared." << std::endl;

                        condition_itr->GetGeometry().clear();
                        condition_itr->Reset(ACTIVE);
                        condition_itr->Set(TO_ERASE);
                    }
                }
            }
        }
                }
    }
} // end namespace MPMSearchElementUtility

} // end namespace Kratos

#endif // KRATOS_MPM_SEARCH_ELEMENT_UTILITY

