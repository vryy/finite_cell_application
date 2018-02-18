// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//



// Project includes
#include "includes/element.h"
#include "custom_python/add_moment_fitting_utility_to_python.h"
#include "custom_algebra/brep.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/moment_fitting_utility.h"
#include "custom_utilities/moment_fitted_quad_tree_subcell.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

/// Python interface for FitQuadrature
/// REMARKS: if fit_quadrature_method = 0x0X, this function uses the standard Gauss quadrature rule on the element, defined by fit_quadrature_order
/// Otherwise it will fit a custom quadrature, see quadrature_utility.h for decoding of the integration_method
template<class TIntegratorType>
void MomentFittingUtility_FitQuadrature(MomentFittingUtility& rDummy,
        Element::Pointer p_elem,
        boost::python::list& r_funcs,
        const BRep& r_brep,
        const TIntegratorType& r_integrator,
        const int& fit_quadrature_method,
        const int& integrator_integration_method,
        const std::string& solver_type,
        const int& echo_level,
        const double& small_weight)
{
    std::vector<FunctionR3R1::Pointer> funcs;
    typedef boost::python::stl_input_iterator<FunctionR3R1::Pointer> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& f,
                  std::make_pair(iterator_value_type(r_funcs), // begin
                    iterator_value_type() ) ) // end
    {
        funcs.push_back(f);
    }

    int fit_quadrature_type = QuadratureUtility::GetQuadratureType(fit_quadrature_method);

    int fit_quadrature_order = QuadratureUtility::GetQuadratureOrder(fit_quadrature_method);

    GeometryData::IntegrationMethod ElementalIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(fit_quadrature_order);

    if(fit_quadrature_type == 0)
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
        // TODO this quadrature is not compatible with Bezier element; check this again
        Element::GeometryType::IntegrationPointsArrayType integration_points
                = r_integrator.Get().ConstructCustomQuadrature(fit_quadrature_type, fit_quadrature_order);

        Vector Weight = rDummy.FitQuadrature<FunctionR3R1, TIntegratorType>(p_elem->GetGeometry(),
                funcs, r_brep, r_integrator, integration_points,
                integrator_integration_method, solver_type, echo_level, small_weight);

        for(std::size_t i = 0; i < integration_points.size(); ++i)
            integration_points[i].Weight() = Weight(i);

        /* create new quadrature and assign to the geometry */
        FiniteCellGeometryUtility::AssignGeometryData(p_elem->GetGeometry(), ElementalIntegrationMethod, integration_points);
    }
}


/// Python interface for FitQuadrature
/// REMARKS: if fit_quadrature_method = 0x0X, this function uses the standard Gauss quadrature rule on the element, defined by fit_quadrature_order
/// Otherwise it will fit a custom quadrature, see quadrature_utility.h for decoding of the integration_method
template<class TIntegratorType>
void MomentFittingUtility_MultithreadedFitQuadrature(MomentFittingUtility& rDummy,
        boost::python::list& r_elements,
        boost::python::list& r_funcs,
        const BRep& r_brep,
        boost::python::list& r_integrators,
        const int& fit_quadrature_method,
        const int& integrator_integration_method,
        const std::string& solver_type,
        const int& echo_level,
        const double& small_weight)
{
    std::vector<Element::Pointer> elements;
    typedef boost::python::stl_input_iterator<Element::Pointer> iterator_element_type;
    BOOST_FOREACH(const iterator_element_type::value_type& e,
                  std::make_pair(iterator_element_type(r_elements), // begin
                  iterator_element_type() ) ) // end
    {
        elements.push_back(e);
    }
//        KRATOS_WATCH(elements.size())

    std::vector<FunctionR3R1::Pointer> funcs;
    typedef boost::python::stl_input_iterator<FunctionR3R1::Pointer> iterator_func_type;
    BOOST_FOREACH(const iterator_func_type::value_type& f,
                  std::make_pair(iterator_func_type(r_funcs), // begin
                  iterator_func_type() ) ) // end
    {
        funcs.push_back(f);
    }
//        KRATOS_WATCH(funcs.size())

    typedef typename TIntegratorType::Pointer TIntegratorPointerType;
    std::vector<TIntegratorPointerType> integrators;
    typedef boost::python::stl_input_iterator<TIntegratorPointerType> iterator_integrator_type;
    BOOST_FOREACH(const typename iterator_integrator_type::value_type& i,
                  std::make_pair(iterator_integrator_type(r_integrators), // begin
                  iterator_integrator_type() ) ) // end
    {
        integrators.push_back(i);
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
    FiniteCellAuxilliaryUtility::CreatePartition(number_of_threads, elements.size(), element_partition);
    std::cout << "element_partition:";
    for(std::size_t i = 0; i < element_partition.size(); ++i)
        std::cout << " " << element_partition[i];
    std::cout << std::endl;

    boost::progress_display* show_progress = NULL;

    if(echo_level == -2)
        show_progress = new boost::progress_display( elements.size() );

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int k = 0; k < number_of_threads; ++k)
    {
        std::vector<Element::Pointer>::iterator it_first_element = elements.begin() + element_partition[k];
        std::vector<Element::Pointer>::iterator it_last_element = elements.begin() + element_partition[k+1];

        typename std::vector<TIntegratorPointerType>::iterator it_integrator = integrators.begin() + element_partition[k];

        for(std::vector<Element::Pointer>::iterator it = it_first_element; it != it_last_element; ++it, ++it_integrator)
        {
            if(fit_quadrature_type == 0)
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
                        = (*it_integrator)->Get(0).ConstructCustomQuadrature(fit_quadrature_type, fit_quadrature_order);
//KRATOS_WATCH(fit_quadrature_type)
//KRATOS_WATCH(fit_quadrature_order)
//KRATOS_WATCH(integration_points.size())
                Vector Weight = rDummy.FitQuadrature<FunctionR3R1, TIntegratorType>((*it)->GetGeometry(),
                        funcs, r_brep, *(*it_integrator), integration_points,
                        integrator_integration_method, solver_type, echo_level, small_weight);

                for(std::size_t i = 0; i < integration_points.size(); ++i)
                    integration_points[i].Weight() = Weight(i);

//                    GeometryData::IntegrationMethod ElementalIntegrationMethod
//                        = Function<double, double>::GetIntegrationMethod(fit_quadrature_order);

                // it is a hack here, since the integration method can be larger than Kratos can accommodate. We set to minimum value. In the element this information is not important anyway.
                GeometryData::IntegrationMethod ElementalIntegrationMethod = GeometryData::GI_GAUSS_1;

                /* create new quadrature and assign to the geometry */
                FiniteCellGeometryUtility::AssignGeometryData((*it)->GetGeometry(), ElementalIntegrationMethod, integration_points);
            }

            if(echo_level == -2)
                ++(*show_progress);
        }
    }
    if(echo_level == -2)
    {
        delete show_progress;
        std::cout << std::endl;
    }
}


template<class TTreeType, class TFunctionType> // = MomentFittedQuadTreeSubCell
void MomentFittingutility_MultithreadedFitQuadratureSubCell(MomentFittingUtility& rDummy,
        boost::python::list& r_trees,
        boost::python::list& r_funcs,
        const BRep& r_brep,
        const int& integrator_integration_method,
        const std::string& solver_type,
        const int& echo_level,
        const double& small_weight)
{
    /* extract the tree with subcell */
    std::vector<typename TTreeType::Pointer> trees;
    typedef boost::python::stl_input_iterator<typename TTreeType::Pointer> iterator_tree_type;
    BOOST_FOREACH(const typename iterator_tree_type::value_type& t,
                  std::make_pair(iterator_tree_type(r_trees), // begin
                  iterator_tree_type() ) ) // end
    {
        trees.push_back(t);
    }

    /* extract the fitting functions */
    std::vector<typename TFunctionType::Pointer> funcs;
    typedef boost::python::stl_input_iterator<typename TFunctionType::Pointer> iterator_value_type;
    BOOST_FOREACH(const typename iterator_value_type::value_type& f,
                std::make_pair(iterator_value_type(r_funcs), // begin
                iterator_value_type() ) ) // end
    {
        funcs.push_back(f);
    }

    rDummy.MultithreadedFitQuadratureSubCell<TTreeType, TFunctionType>(trees, funcs, r_brep, integrator_integration_method, solver_type, echo_level, small_weight);
}

void FiniteCellApplication_AddMomentFittingUtilityToPython()
{

    class_<MomentFittingUtility, MomentFittingUtility::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("MomentFittingUtility", init<>())
    .def("FitQuadrature", &MomentFittingUtility_FitQuadrature<QuadTree<1> >)
    .def("MultithreadedFitQuadrature", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<1> >)
    .def("MultithreadedFitQuadratureSubCell", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<1>, FunctionR3R1>)
    .def("MultithreadedFitQuadratureSubCellUnique", &MomentFittingUtility_MultithreadedFitQuadrature<MomentFittedQuadTreeSubCell<1> >)
    /////
    .def("FitQuadrature2", &MomentFittingUtility_FitQuadrature<QuadTree<2> >)
    .def("MultithreadedFitQuadrature2", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<2> >)
    .def("MultithreadedFitQuadratureSubCell2", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<2>, FunctionR3R1>)
    .def("MultithreadedFitQuadratureSubCellUnique2", &MomentFittingUtility_MultithreadedFitQuadrature<MomentFittedQuadTreeSubCell<2> >)
    /////
    .def("FitQuadrature3", &MomentFittingUtility_FitQuadrature<QuadTree<3> >)
    .def("MultithreadedFitQuadrature3", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<3> >)
    .def("MultithreadedFitQuadratureSubCell3", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<3>, FunctionR3R1>)
    .def("MultithreadedFitQuadratureSubCellUnique3", &MomentFittingUtility_MultithreadedFitQuadrature<MomentFittedQuadTreeSubCell<3> >)
    /////
    .def("FitQuadrature4", &MomentFittingUtility_FitQuadrature<QuadTree<4> >)
    .def("MultithreadedFitQuadrature4", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<4> >)
    .def("MultithreadedFitQuadratureSubCell4", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<4>, FunctionR3R1>)
    .def("MultithreadedFitQuadratureSubCellUnique4", &MomentFittingUtility_MultithreadedFitQuadrature<MomentFittedQuadTreeSubCell<4> >)
    /////
    .def("FitQuadrature5", &MomentFittingUtility_FitQuadrature<QuadTree<5> >)
    .def("MultithreadedFitQuadrature5", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<5> >)
    .def("MultithreadedFitQuadratureSubCell5", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<5>, FunctionR3R1>)
    .def("MultithreadedFitQuadratureSubCellUnique5", &MomentFittingUtility_MultithreadedFitQuadrature<MomentFittedQuadTreeSubCell<5> >)
    /////
    .def("FitQuadrature6", &MomentFittingUtility_FitQuadrature<QuadTree<6> >)
    .def("MultithreadedFitQuadrature6", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<6> >)
    .def("MultithreadedFitQuadratureSubCell6", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<6>, FunctionR3R1>)
    .def("MultithreadedFitQuadratureSubCellUnique6", &MomentFittingUtility_MultithreadedFitQuadrature<MomentFittedQuadTreeSubCell<6> >)
    /////
    .def("FitQuadrature7", &MomentFittingUtility_FitQuadrature<QuadTree<7> >)
    .def("MultithreadedFitQuadrature7", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<7> >)
    .def("MultithreadedFitQuadratureSubCell7", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<7>, FunctionR3R1>)
    .def("MultithreadedFitQuadratureSubCellUnique7", &MomentFittingUtility_MultithreadedFitQuadrature<MomentFittedQuadTreeSubCell<7> >)
    /////
    .def("FitQuadrature8", &MomentFittingUtility_FitQuadrature<QuadTree<8> >)
    .def("MultithreadedFitQuadrature8", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<8> >)
    .def("MultithreadedFitQuadratureSubCell8", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<8>, FunctionR3R1>)
    .def("MultithreadedFitQuadratureSubCellUnique8", &MomentFittingUtility_MultithreadedFitQuadrature<MomentFittedQuadTreeSubCell<8> >)
    /////
    .def("FitQuadrature9", &MomentFittingUtility_FitQuadrature<QuadTree<9> >)
    .def("MultithreadedFitQuadrature9", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<9> >)
    .def("MultithreadedFitQuadratureSubCell9", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<9>, FunctionR3R1>)
    .def("MultithreadedFitQuadratureSubCellUnique9", &MomentFittingUtility_MultithreadedFitQuadrature<MomentFittedQuadTreeSubCell<9> >)
    /////
    .def("FitQuadrature10", &MomentFittingUtility_FitQuadrature<QuadTree<10> >)
    .def("MultithreadedFitQuadrature10", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<10> >)
    .def("MultithreadedFitQuadratureSubCell10", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<10>, FunctionR3R1>)
    .def("MultithreadedFitQuadratureSubCellUnique10", &MomentFittingUtility_MultithreadedFitQuadrature<MomentFittedQuadTreeSubCell<10> >)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

