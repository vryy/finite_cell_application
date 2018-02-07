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
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_algebra/brep.h"
#include "custom_utilities/quadrature_utility.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/div_free_basis_utility.h"
#include "custom_utilities/moment_fitting_utility.h"
#include "custom_utilities/moment_fitted_quad_tree_subcell.h"
#include "custom_utilities/finite_cell_auxilliary_utility.h"
#include "custom_utilities/ghost_penalty_utility.h"
#include "custom_utilities/skeleton_penalty_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

Matrix ComputeDivFreeBasis(DivFreeBasisUtility& dummy, const std::size_t& Dim, const std::size_t& Degree, const DivFreeBasisUtility::CoordinatesArrayType& rPoint)
{
    if(Dim == 2)
    {
        if(Degree == 0)
            return dummy.GetValues<2, 0>(rPoint);
        else if(Degree == 1)
            return dummy.GetValues<2, 1>(rPoint);
        else if(Degree == 2)
            return dummy.GetValues<2, 2>(rPoint);
        else if(Degree == 3)
            return dummy.GetValues<2, 3>(rPoint);
        else if(Degree == 4)
            return dummy.GetValues<2, 4>(rPoint);
    }
}

std::size_t FiniteCellAuxilliaryUtility_GetLastNodeId(FiniteCellAuxilliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastNodeId(r_model_part);
}

std::size_t FiniteCellAuxilliaryUtility_GetLastElementId(FiniteCellAuxilliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastElementId(r_model_part);
}

std::size_t FiniteCellAuxilliaryUtility_GetLastConditionId(FiniteCellAuxilliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastConditionId(r_model_part);
}

std::size_t FiniteCellAuxilliaryUtility_GetLastPropertiesId(FiniteCellAuxilliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastPropertiesId(r_model_part);
}

void FiniteCellAuxilliaryUtility_AddElement(FiniteCellAuxilliaryUtility& rDummy, ModelPart::ElementsContainerType& rpElements,
        Element::Pointer pElement)
{
    rDummy.AddElement(rpElements, pElement);
}

void FiniteCellAuxilliaryUtility_RemoveElement(FiniteCellAuxilliaryUtility& rDummy, ModelPart::ElementsContainerType& rpElements,
        Element::Pointer pElement)
{
    rDummy.RemoveElement(rpElements, pElement);
}

int QuadratureUtility_GetQuadratureType(QuadratureUtility& rDummy, const int& integration_method)
{
    return rDummy.GetQuadratureType(integration_method);
}

int QuadratureUtility_GetQuadratureOrder(QuadratureUtility& rDummy, const int& integration_method)
{
    return rDummy.GetQuadratureOrder(integration_method);
}

template<class TTreeType, class TBRepType>
void FiniteCellAuxilliaryUtility_MultithreadedRefineBy(FiniteCellAuxilliaryUtility& rDummy, boost::python::list& r_trees,
        const TBRepType& r_brep)
{
    typedef typename TTreeType::Pointer TTreePointerType;
    std::vector<TTreePointerType> trees;
    typedef boost::python::stl_input_iterator<TTreePointerType> iterator_tree_type;
    BOOST_FOREACH(const typename iterator_tree_type::value_type& t,
                  std::make_pair(iterator_tree_type(r_trees), // begin
                  iterator_tree_type() ) ) // end
    {
        trees.push_back(t);
    }

    rDummy.MultithreadedRefineBy<TTreeType, TBRepType>(trees, r_brep);
}

/// Extract the element from the list of ids
ModelPart::ElementsContainerType FiniteCellAuxilliaryUtility_GetElements(FiniteCellAuxilliaryUtility& rDummy,
    ModelPart& r_model_part, boost::python::list& element_list)
{
    std::set<std::size_t> element_ids;

    typedef boost::python::stl_input_iterator<int> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& id,
                  std::make_pair(iterator_value_type(element_list), // begin
                    iterator_value_type() ) ) // end
    {
        element_ids.insert(static_cast<std::size_t>(id));
    }

    return rDummy.GetElements(r_model_part, element_ids);
}

/// Extract the element from the list of ids
void FiniteCellAuxilliaryUtility_GetElements2(FiniteCellAuxilliaryUtility& rDummy,
    ModelPart::ElementsContainerType& rpElements,
    ModelPart& r_model_part, boost::python::list& element_list)
{
    std::set<std::size_t> element_ids;

    typedef boost::python::stl_input_iterator<int> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& id,
                  std::make_pair(iterator_value_type(element_list), // begin
                    iterator_value_type() ) ) // end
    {
        element_ids.insert(static_cast<std::size_t>(id));
    }

    rDummy.GetElements(rpElements, r_model_part, element_ids);
}

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
                        = (*it_integrator)->Get().ConstructCustomQuadrature(fit_quadrature_type, fit_quadrature_order);
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

                // it is a hack here, since the integration method can be larger than Kratos can accomodate. We set to minimum value. In the element this information is not important anyway.
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

void GhostPenaltyUtility_ProbeNeighbourElements(GhostPenaltyUtility& rDummy, Element::Pointer p_elem)
{
    rDummy.ProbeNeighbourElements(p_elem);
}

ModelPart::ConditionsContainerType GhostPenaltyUtility_SetUpSurfacePenaltyConditions1(GhostPenaltyUtility& rDummy,
        Element::Pointer p_elem, GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t lastCondId, Properties::Pointer pProperties)
{
    return rDummy.SetUpSurfacePenaltyConditions(p_elem, p_sample_condition, r_brep, lastCondId, pProperties);
}

ModelPart::ConditionsContainerType GhostPenaltyUtility_SetUpSurfacePenaltyConditions2(GhostPenaltyUtility& rDummy,
        ModelPart& r_model_part, GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t lastCondId, Properties::Pointer pProperties)
{
    return rDummy.SetUpSurfacePenaltyConditions(r_model_part, p_sample_condition, r_brep, lastCondId, pProperties);
}

ModelPart::ConditionsContainerType GhostPenaltyUtility_SetUpSurfacePenaltyConditions3(GhostPenaltyUtility& rDummy,
        ModelPart& r_model_part, ModelPart::ElementsContainerType& pElements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t lastCondId, Properties::Pointer pProperties)
{
    return rDummy.SetUpSurfacePenaltyConditions(r_model_part, pElements, p_sample_condition, r_brep, lastCondId, pProperties);
}

Condition::Pointer GhostPenaltyUtility_SetUpSurfacePenaltyCondition(GhostPenaltyUtility& rDummy,
        Element::Pointer p_element_1, Element::Pointer p_element_2, GhostPenaltyCondition::Pointer p_sample_condition,
        std::size_t lastCondId, Properties::Pointer pProperties)
{
    return rDummy.SetUpSurfacePenaltyCondition(p_element_1, p_element_2, p_sample_condition, lastCondId, pProperties);
}

void GhostPenaltyUtility_ProbeShapeFunctionSecondDerivatives(GhostPenaltyUtility& rDummy, Element::Pointer p_element)
{
    rDummy.ProbeShapeFunctionSecondDerivatives(p_element->GetGeometry());
}

void FiniteCellApplication_AddCustomUtilitiesToPython()
{
    void(QuadratureUtility::*pointer_to_CreateConditionFromQuadraturePoint)(ModelPart&, boost::python::list&,
            const std::string&, const double&, const double&) const = &QuadratureUtility::PyCreateConditionFromQuadraturePoint;
    ModelPart::ConditionsContainerType(QuadratureUtility::*pointer_to_CreateConditionFromPoint)(ModelPart&,
            boost::python::list&, const std::string&) const = &QuadratureUtility::PyCreateConditionFromPoint;
    ModelPart::ConditionsContainerType(QuadratureUtility::*pointer_to_CreateConditionFromPoint2)(ModelPart&,
            boost::python::list&, const std::string&, Properties::Pointer) const = &QuadratureUtility::PyCreateConditionFromPoint;

    class_<QuadratureUtility, QuadratureUtility::Pointer, boost::noncopyable>
    ("QuadratureUtility", init<>())
    .def("GetDefaultIntegrationMethod", &QuadratureUtility::GetDefaultIntegrationMethod<Element>)
    .def("GetDefaultIntegrationMethod", &QuadratureUtility::GetDefaultIntegrationMethod<Condition>)
    .def("GetQuadratureType", QuadratureUtility_GetQuadratureType)
    .def("GetQuadratureOrder", QuadratureUtility_GetQuadratureOrder)
    .def("ScaleQuadrature", &QuadratureUtility::PyScaleQuadrature)
    .def("SaveQuadrature", &QuadratureUtility::PySaveQuadrature)
    .def("SaveQuadrature", &QuadratureUtility::PySaveQuadratureAdvanced)
    .def("SaveQuadratureSubCell", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<1> >)
    .def("SaveQuadratureSubCell2", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<2> >)
    .def("SaveQuadratureSubCell3", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<3> >)
    .def("SaveQuadratureSubCell4", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<4> >)
    .def("SaveQuadratureSubCell5", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<5> >)
    .def("SetQuadrature", &QuadratureUtility::PySetQuadrature)
    .def("CreateConditionFromQuadraturePoint", pointer_to_CreateConditionFromQuadraturePoint)
    .def("CreateConditionFromPoint", pointer_to_CreateConditionFromPoint)
    .def("CreateConditionFromPoint", pointer_to_CreateConditionFromPoint2)
    .def("CreatePoint", &QuadratureUtility::CreatePoint)
    ;

    void(DivFreeBasisUtility::*pointer_to_AssignQuadrature2D)(Element::Pointer&, const LevelSet&, const unsigned int&, const unsigned int&) = &DivFreeBasisUtility::AssignQuadrature2D;

    class_<DivFreeBasisUtility, DivFreeBasisUtility::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("DivFreeBasisUtility", init<>() )
    .def("GetValues", &ComputeDivFreeBasis)
    .def("AssignQuadrature2D", pointer_to_AssignQuadrature2D)
    ;

    class_<MomentFittingUtility, MomentFittingUtility::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("MomentFittingUtility", init<>())
    .def("FitQuadrature", &MomentFittingUtility_FitQuadrature<QuadTree<1> >)
    .def("MultithreadedFitQuadrature", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<1> >)
    .def("MultithreadedFitQuadratureSubCell", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<1>, FunctionR3R1>)
    /////
    .def("FitQuadrature2", &MomentFittingUtility_FitQuadrature<QuadTree<2> >)
    .def("MultithreadedFitQuadrature2", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<2> >)
    .def("MultithreadedFitQuadratureSubCell2", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<2>, FunctionR3R1>)
    /////
    .def("FitQuadrature3", &MomentFittingUtility_FitQuadrature<QuadTree<3> >)
    .def("MultithreadedFitQuadrature3", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<3> >)
    .def("MultithreadedFitQuadratureSubCell3", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<3>, FunctionR3R1>)
    /////
    .def("FitQuadrature4", &MomentFittingUtility_FitQuadrature<QuadTree<4> >)
    .def("MultithreadedFitQuadrature4", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<4> >)
    .def("MultithreadedFitQuadratureSubCell4", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<4>, FunctionR3R1>)
    /////
    .def("FitQuadrature5", &MomentFittingUtility_FitQuadrature<QuadTree<5> >)
    .def("MultithreadedFitQuadrature5", &MomentFittingUtility_MultithreadedFitQuadrature<QuadTree<5> >)
    .def("MultithreadedFitQuadratureSubCell5", &MomentFittingutility_MultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<5>, FunctionR3R1>)
    ;

    Condition::Pointer(FiniteCellAuxilliaryUtility::*pointer_to_PyCreateCondition)(ModelPart&, const std::string&,
            const std::size_t&, Properties::Pointer, boost::python::list&) const = &FiniteCellAuxilliaryUtility::PyCreateCondition;

    void(FiniteCellAuxilliaryUtility::*pointer_to_Clean)(ModelPart&,
            ModelPart::ConditionsContainerType&, const int&) const = &FiniteCellAuxilliaryUtility::Clean;

    void(FiniteCellAuxilliaryUtility::*pointer_to_PrintGeometry)(Element::GeometryType::Pointer) const = &FiniteCellAuxilliaryUtility::Print;

    class_<FiniteCellAuxilliaryUtility, FiniteCellAuxilliaryUtility::Pointer, boost::noncopyable>
    ("FiniteCellAuxilliaryUtility", init<>())
    .def("CreateCondition", pointer_to_PyCreateCondition)
    .def("GetElements", &FiniteCellAuxilliaryUtility_GetElements)
    .def("GetElements", &FiniteCellAuxilliaryUtility_GetElements2)
    .def("Clean", pointer_to_Clean)
    .def("GetLastNodeId", &FiniteCellAuxilliaryUtility_GetLastNodeId)
    .def("GetLastElementId", &FiniteCellAuxilliaryUtility_GetLastElementId)
    .def("GetLastConditionId", &FiniteCellAuxilliaryUtility_GetLastConditionId)
    .def("GetLastPropertiesId", &FiniteCellAuxilliaryUtility_GetLastPropertiesId)
    .def("AddElement", &FiniteCellAuxilliaryUtility_AddElement)
    .def("RemoveElement", &FiniteCellAuxilliaryUtility_RemoveElement)
    .def("MultithreadedRefineBy", &FiniteCellAuxilliaryUtility_MultithreadedRefineBy<RefinableTree, BRep>)
    .def("Print", pointer_to_PrintGeometry)
    ;

    class_<GhostPenaltyUtility, GhostPenaltyUtility::Pointer, boost::noncopyable>
    ("GhostPenaltyUtility", init<>())
    .def("ProbeNeighbourElements", &GhostPenaltyUtility_ProbeNeighbourElements)
    .def("SetUpSurfacePenaltyConditions", &GhostPenaltyUtility_SetUpSurfacePenaltyConditions1)
    .def("SetUpSurfacePenaltyConditions", &GhostPenaltyUtility_SetUpSurfacePenaltyConditions2)
    .def("SetUpSurfacePenaltyConditions", &GhostPenaltyUtility_SetUpSurfacePenaltyConditions3)
    .def("SetUpSurfacePenaltyCondition", &GhostPenaltyUtility_SetUpSurfacePenaltyCondition)
    .def("ProbeShapeFunctionSecondDerivatives", &GhostPenaltyUtility_ProbeShapeFunctionSecondDerivatives)
    ;

    class_<SkeletonPenaltyUtility, SkeletonPenaltyUtility::Pointer, bases<GhostPenaltyUtility>, boost::noncopyable>
    ("SkeletonPenaltyUtility", init<>())
    ;
}

}  // namespace Python.

}  // namespace Kratos.

