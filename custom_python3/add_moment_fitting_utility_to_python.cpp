// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes

// Project includes
#include "includes/element.h"
#include "custom_python3/add_moment_fitting_utility_to_python.h"
#include "brep_application/custom_algebra/brep.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/moment_fitting_utility.h"
#include "custom_utilities/moment_fitted_quad_tree_subcell.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

/// Python interface for FitQuadrature
/// REMARKS: if fit_quadrature_method = 0x0X, this function uses the standard Gauss quadrature rule on the element, defined by fit_quadrature_order
/// Otherwise it will fit a custom quadrature, see quadrature_utility.h for decoding of the integration_method
template<class TIntegratorType>
void MomentFittingUtility_FitQuadrature(MomentFittingUtility& rDummy,
                                        Element::Pointer p_elem,
                                        const pybind11::list& r_funcs,
                                        const BRep& r_brep,
                                        const TIntegratorType& r_integrator,
                                        const int& fit_quadrature_method,
                                        const int& integrator_integration_method,
                                        const std::string& solver_type,
                                        const int& echo_level,
                                        const double& small_weight)
{
    std::vector<FunctionR3R1::Pointer> funcs;
    for (auto f : r_funcs)
    {
        funcs.push_back(f.cast<FunctionR3R1::Pointer>());
    }

    int fit_quadrature_type = QuadratureUtility::GetQuadratureType(fit_quadrature_method);

    int fit_quadrature_order = QuadratureUtility::GetQuadratureOrder(fit_quadrature_method);

    if (echo_level > -1)
    {
        std::cout << "#############prepare moment-fitting for element " << p_elem->Id()
                  << " of type " << typeid(*p_elem).name()
                  << ", geometry type " << typeid(p_elem->GetGeometry()).name()
                  << std::endl;
        KRATOS_WATCH(fit_quadrature_type)
        KRATOS_WATCH(fit_quadrature_order)
    }

    GeometryData::IntegrationMethod ElementalIntegrationMethod
        = Function<double, double>::GetIntegrationMethod(fit_quadrature_order);

    if (fit_quadrature_type == 0)
    {
        const Element::GeometryType::IntegrationPointsArrayType& integration_points
            = p_elem->GetGeometry().IntegrationPoints( ElementalIntegrationMethod );

        Vector Weight = rDummy.FitQuadrature<FunctionR3R1, TIntegratorType>(p_elem->GetGeometry(),
                        funcs, r_brep, r_integrator, integration_points,
                        integrator_integration_method, solver_type, echo_level, small_weight);

        /* create new quadrature and assign to the geometry */
        FiniteCellGeometryUtility::AssignGeometryData(p_elem->GetGeometry(), ElementalIntegrationMethod, Weight);
    }
    else
    {
        Element::GeometryType::IntegrationPointsArrayType integration_points
            = r_integrator.ConstructCustomQuadrature(fit_quadrature_type, fit_quadrature_order);

        Vector Weight = rDummy.FitQuadrature<FunctionR3R1, TIntegratorType>(p_elem->GetGeometry(),
                        funcs, r_brep, r_integrator, integration_points,
                        integrator_integration_method, solver_type, echo_level, small_weight);

        for (std::size_t i = 0; i < integration_points.size(); ++i)
        {
            integration_points[i].Weight() = Weight(i);
        }

        /* create new quadrature and assign to the geometry */
        FiniteCellGeometryUtility::AssignGeometryData(p_elem->GetGeometry(), ElementalIntegrationMethod, integration_points);
    }
}


/// Python interface for FitQuadrature
/// REMARKS: if fit_quadrature_method = 0x0X, this function uses the standard Gauss quadrature rule on the element, defined by fit_quadrature_order
/// Otherwise it will fit a custom quadrature, see quadrature_utility.h for decoding of the integration_method
template<class TIntegratorType>
void MomentFittingUtility_MultithreadedFitQuadrature(MomentFittingUtility& rDummy,
        const pybind11::list& r_elements,
        const pybind11::list& r_funcs,
        const BRep& r_brep,
        const pybind11::list& r_integrators,
        const int& fit_quadrature_method,
        const int& integrator_integration_method,
        const std::string& solver_type,
        const int& echo_level,
        const double& small_weight)
{
    std::vector<Element::Pointer> elements;
    for (auto e : r_elements)
    {
        elements.push_back(e.cast<Element::Pointer>());
    }
//        KRATOS_WATCH(elements.size())

    std::vector<FunctionR3R1::Pointer> funcs;
    for (auto f : r_funcs)
    {
        funcs.push_back(f.cast<FunctionR3R1::Pointer>());
    }
//        KRATOS_WATCH(funcs.size())

    typedef typename TIntegratorType::Pointer TIntegratorPointerType;
    std::vector<TIntegratorPointerType> integrators;
    for (auto i : r_integrators)
    {
        integrators.push_back(i.cast<TIntegratorPointerType>());
    }
//        KRATOS_WATCH(integrators.size())

    int fit_quadrature_type = QuadratureUtility::GetQuadratureType(fit_quadrature_method);

    int fit_quadrature_order = QuadratureUtility::GetQuadratureOrder(fit_quadrature_method);

    unsigned int number_of_threads = 1;
    std::vector<unsigned int> element_partition;
#ifdef _OPENMP
    number_of_threads = omp_get_max_threads();
#endif
    KRATOS_WATCH(number_of_threads)
    FiniteCellAuxiliaryUtility::CreatePartition(number_of_threads, elements.size(), element_partition);
    std::cout << "element_partition:";
    for (std::size_t i = 0; i < element_partition.size(); ++i)
    {
        std::cout << " " << element_partition[i];
    }
    std::cout << std::endl;

    boost::progress_display* show_progress = NULL;

    if (echo_level == -2)
    {
        show_progress = new boost::progress_display( elements.size() );
    }

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int k = 0; k < number_of_threads; ++k)
    {
        std::vector<Element::Pointer>::iterator it_first_element = elements.begin() + element_partition[k];
        std::vector<Element::Pointer>::iterator it_last_element = elements.begin() + element_partition[k + 1];

        typename std::vector<TIntegratorPointerType>::iterator it_integrator = integrators.begin() + element_partition[k];

        for (std::vector<Element::Pointer>::iterator it = it_first_element; it != it_last_element; ++it, ++it_integrator)
        {
            if (fit_quadrature_type == 0)
            {
                GeometryData::IntegrationMethod ElementalIntegrationMethod
                    = Function<double, double>::GetIntegrationMethod(fit_quadrature_order);

                const Element::GeometryType::IntegrationPointsArrayType& integration_points
                    = (*it)->GetGeometry().IntegrationPoints( ElementalIntegrationMethod );

                Vector Weight = rDummy.FitQuadrature<FunctionR3R1, TIntegratorType>((*it)->GetGeometry(),
                                funcs, r_brep, *(*it_integrator), integration_points,
                                integrator_integration_method, solver_type, echo_level, small_weight);

                /* create new quadrature and assign to the geometry */
                FiniteCellGeometryUtility::AssignGeometryData((*it)->GetGeometry(), ElementalIntegrationMethod, Weight);
            }
            else
            {
                Element::GeometryType::IntegrationPointsArrayType integration_points
                    = (*it_integrator)->ConstructCustomQuadrature(fit_quadrature_type, fit_quadrature_order);
//KRATOS_WATCH(fit_quadrature_type)
//KRATOS_WATCH(fit_quadrature_order)
//KRATOS_WATCH(integration_points.size())
                Vector Weight = rDummy.FitQuadrature<FunctionR3R1, TIntegratorType>((*it)->GetGeometry(),
                                funcs, r_brep, *(*it_integrator), integration_points,
                                integrator_integration_method, solver_type, echo_level, small_weight);

                for (std::size_t i = 0; i < integration_points.size(); ++i)
                {
                    integration_points[i].Weight() = Weight(i);
                }

//                    GeometryData::IntegrationMethod ElementalIntegrationMethod
//                        = Function<double, double>::GetIntegrationMethod(fit_quadrature_order);

                // it is a hack here, since the integration method can be larger than Kratos can accommodate. We set to minimum value. In the element this information is not important anyway.
                GeometryData::IntegrationMethod ElementalIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;

                /* create new quadrature and assign to the geometry */
                FiniteCellGeometryUtility::AssignGeometryData((*it)->GetGeometry(), ElementalIntegrationMethod, integration_points);
            }

            if (echo_level == -2)
            {
                ++(*show_progress);
            }
        }
    }
    if (echo_level == -2)
    {
        delete show_progress;
        std::cout << std::endl;
    }
}


template<class TTreeType, class TFunctionType> // = MomentFittedQuadTreeSubCell
void MomentFittingutility_FitQuadratureSubCell(MomentFittingUtility& rDummy,
        typename TTreeType::Pointer p_tree,
        const pybind11::list& r_funcs,
        const BRep& r_brep,
        const int& integrator_integration_method,
        const std::string& solver_type,
        const int& echo_level,
        const double& small_weight,
        const ProcessInfo& rCurrentProcessInfo)
{
    /* extract the fitting functions */
    std::vector<typename TFunctionType::Pointer> funcs;
    for (auto f : r_funcs)
    {
        funcs.push_back(f.cast<typename TFunctionType::Pointer>());
    }

    rDummy.FitQuadratureSubCell<TTreeType, TFunctionType>(p_tree, funcs, r_brep, integrator_integration_method, solver_type, echo_level, small_weight, rCurrentProcessInfo);
}


template<class TTreeType, class TFunctionType> // = MomentFittedQuadTreeSubCell
void MomentFittingutility_MultithreadedFitQuadratureSubCell(MomentFittingUtility& rDummy,
        const pybind11::list& r_trees,
        const pybind11::list& r_funcs,
        const BRep& r_brep,
        const int& integrator_integration_method,
        const std::string& solver_type,
        const int& echo_level,
        const double& small_weight,
        const ProcessInfo& rCurrentProcessInfo)
{
    /* extract the tree with subcell */
    std::vector<typename TTreeType::Pointer> trees;
    for (auto t : r_trees)
    {
        trees.push_back(t.cast<typename TTreeType::Pointer>());
    }

    /* extract the fitting functions */
    std::vector<typename TFunctionType::Pointer> funcs;
    for (auto f : r_funcs)
    {
        funcs.push_back(f.cast<typename TFunctionType::Pointer>());
    }

    rDummy.MultithreadedFitQuadratureSubCell<TTreeType, TFunctionType>(trees, funcs, r_brep, integrator_integration_method, solver_type, echo_level, small_weight, rCurrentProcessInfo);
}

void FiniteCellApplication_AddMomentFittingUtilityToPython(pybind11::module& m)
{

    class_<MomentFittingUtility, MomentFittingUtility::Pointer, QuadratureUtility>
    (m, "MomentFittingUtility")
    .def(init<>())
    .def("FitQuadrature", &MomentFittingUtility_FitQuadrature<FunctionIntegrator>)
    .def("MultithreadedFitQuadrature", &MomentFittingUtility_MultithreadedFitQuadrature<FunctionIntegrator>)
    .def("FitQuadratureSubCell", &MomentFittingutility_FitQuadratureSubCell<BaseMomentFittedQuadTreeSubCell, FunctionR3R1>)
    .def("MultithreadedFitQuadratureSubCell", &MomentFittingutility_MultithreadedFitQuadratureSubCell<BaseMomentFittedQuadTreeSubCell, FunctionR3R1>)
    .def("FitQuadratureSubCellUnique", &MomentFittingUtility_FitQuadrature<FunctionIntegrator>)
    .def("MultithreadedFitQuadratureSubCellUnique", &MomentFittingUtility_MultithreadedFitQuadrature<FunctionIntegrator>)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

