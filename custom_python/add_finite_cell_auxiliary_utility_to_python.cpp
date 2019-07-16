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
#include "custom_python/add_finite_cell_auxiliary_utility_to_python.h"
#include "custom_algebra/brep.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/finite_cell_auxiliary_utility.h"
#include "custom_utilities/moment_fitted_quad_tree_subcell.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

std::size_t FiniteCellAuxiliaryUtility_GetLastNodeId(FiniteCellAuxiliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastNodeId(r_model_part);
}

std::size_t FiniteCellAuxiliaryUtility_GetLastElementId(FiniteCellAuxiliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastElementId(r_model_part);
}

std::size_t FiniteCellAuxiliaryUtility_GetLastConditionId(FiniteCellAuxiliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastConditionId(r_model_part);
}

std::size_t FiniteCellAuxiliaryUtility_GetLastPropertiesId(FiniteCellAuxiliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastPropertiesId(r_model_part);
}

void FiniteCellAuxiliaryUtility_AddElement(FiniteCellAuxiliaryUtility& rDummy, ModelPart::ElementsContainerType& rpElements,
        Element::Pointer pElement)
{
    rDummy.AddElement(rpElements, pElement);
}

void FiniteCellAuxiliaryUtility_RemoveElement(FiniteCellAuxiliaryUtility& rDummy, ModelPart::ElementsContainerType& rpElements,
        Element::Pointer pElement)
{
    rDummy.RemoveElement(rpElements, pElement);
}

/// Create an element from sample element name and from list of nodes
Element::Pointer FiniteCellAuxiliaryUtility_CreateElement(FiniteCellAuxiliaryUtility& rDummy,
    ModelPart& r_model_part, const std::string& sample_elem_name,
    const std::size_t& Id, Properties::Pointer pProperties, boost::python::list& node_ids)
{
    std::vector<std::size_t> node_list;
    typedef boost::python::stl_input_iterator<int> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& id,
                  std::make_pair(iterator_value_type(node_ids), // begin
                  iterator_value_type() ) ) // end
    {
        node_list.push_back(static_cast<std::size_t>(id));
    }

    return rDummy.CreateElement(r_model_part, sample_elem_name, Id, pProperties, node_list);
}

/// Create an element from sample element name and from other element (to provide the geometry)
Element::Pointer FiniteCellAuxiliaryUtility_CreateElement2(FiniteCellAuxiliaryUtility& rDummy,
    ModelPart& r_model_part, const std::string& sample_elem_name,
    const std::size_t& Id, Properties::Pointer pProperties, Element::Pointer pElement)
{
    return rDummy.CreateElement(r_model_part, sample_elem_name, Id, pProperties, pElement);
}

/// Create an element from sample element name and from other condition (to provide the geometry)
Element::Pointer FiniteCellAuxiliaryUtility_CreateElement3(FiniteCellAuxiliaryUtility& rDummy,
    ModelPart& r_model_part, const std::string& sample_elem_name,
    const std::size_t& Id, Properties::Pointer pProperties, Condition::Pointer pCond)
{
    return rDummy.CreateElement(r_model_part, sample_elem_name, Id, pProperties, pCond);
}

/// Create a condition from sample condition and from list of nodes
Condition::Pointer FiniteCellAuxiliaryUtility_CreateCondition(FiniteCellAuxiliaryUtility& rDummy,
    ModelPart& r_model_part, const std::string& sample_cond_name,
    const std::size_t& Id, Properties::Pointer pProperties, boost::python::list& node_ids)
{
    std::vector<std::size_t> node_list;
    typedef boost::python::stl_input_iterator<int> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& id,
                  std::make_pair(iterator_value_type(node_ids), // begin
                  iterator_value_type() ) ) // end
    {
        node_list.push_back(static_cast<std::size_t>(id));
    }

    return rDummy.CreateCondition(r_model_part, sample_cond_name, Id, pProperties, node_list);
}

/// Create a condition from sample condition and from other condition (to provide the geometry)
Condition::Pointer FiniteCellAuxiliaryUtility_CreateCondition2(FiniteCellAuxiliaryUtility& rDummy,
    ModelPart& r_model_part, const std::string& sample_cond_name,
    const std::size_t& Id, Properties::Pointer pProperties, Condition::Pointer pCond)
{
    return rDummy.CreateCondition(r_model_part, sample_cond_name, Id, pProperties, pCond);
}

/// Create a condition from sample condition and from other element (to provide the geometry)
Condition::Pointer FiniteCellAuxiliaryUtility_CreateCondition3(FiniteCellAuxiliaryUtility& rDummy,
    ModelPart& r_model_part, const std::string& sample_cond_name,
    const std::size_t& Id, Properties::Pointer pProperties, Element::Pointer pElement)
{
    return rDummy.CreateCondition(r_model_part, sample_cond_name, Id, pProperties, pElement);
}

template<class TTreeType, class TBRepType>
void FiniteCellAuxiliaryUtility_MultithreadedRefineBy(FiniteCellAuxiliaryUtility& rDummy, boost::python::list& r_trees,
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

template<class TCellType, class TBRepType>
void FiniteCellAuxiliaryUtility_MultithreadedGeneratePhysicalIntegrationPoints(FiniteCellAuxiliaryUtility& rDummy,
    boost::python::list& r_cells, typename TBRepType::Pointer p_brep, int integrator_integration_method)
{
    PointerVectorSet<TCellType> cells;
    typedef boost::python::stl_input_iterator<typename TCellType::Pointer> iterator_tree_type;
    BOOST_FOREACH(const typename iterator_tree_type::value_type& t,
                  std::make_pair(iterator_tree_type(r_cells), // begin
                  iterator_tree_type() ) ) // end
    {
        cells.push_back(t);
    }

    rDummy.MultithreadedGeneratePhysicalIntegrationPoints<TCellType, TBRepType>(cells, *p_brep, integrator_integration_method);
}

/// Extract the element from the list of ids
ModelPart::ElementsContainerType FiniteCellAuxiliaryUtility_GetElements(FiniteCellAuxiliaryUtility& rDummy,
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
void FiniteCellAuxiliaryUtility_GetElements2(FiniteCellAuxiliaryUtility& rDummy,
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

void FiniteCellApplication_AddFiniteCellAuxiliaryUtilityToPython()
{

    void(FiniteCellAuxiliaryUtility::*pointer_to_Clean)(ModelPart&,
            ModelPart::ConditionsContainerType&, const int&) const = &FiniteCellAuxiliaryUtility::Clean;

    void(FiniteCellAuxiliaryUtility::*pointer_to_PrintGeometry)(Element::GeometryType::Pointer) const = &FiniteCellAuxiliaryUtility::Print;

    class_<FiniteCellAuxiliaryUtility, FiniteCellAuxiliaryUtility::Pointer, boost::noncopyable>
    ("FiniteCellAuxiliaryUtility", init<>())
    .def("CreateElement", &FiniteCellAuxiliaryUtility_CreateElement)
    .def("CreateElement", &FiniteCellAuxiliaryUtility_CreateElement2)
    .def("CreateElement", &FiniteCellAuxiliaryUtility_CreateElement3)
    .def("CreateCondition", &FiniteCellAuxiliaryUtility_CreateCondition)
    .def("CreateCondition", &FiniteCellAuxiliaryUtility_CreateCondition2)
    .def("CreateCondition", &FiniteCellAuxiliaryUtility_CreateCondition3)
    .def("GetElements", &FiniteCellAuxiliaryUtility_GetElements)
    .def("GetElements", &FiniteCellAuxiliaryUtility_GetElements2)
    .def("Clean", pointer_to_Clean)
    .def("GetLastNodeId", &FiniteCellAuxiliaryUtility_GetLastNodeId)
    .def("GetLastElementId", &FiniteCellAuxiliaryUtility_GetLastElementId)
    .def("GetLastConditionId", &FiniteCellAuxiliaryUtility_GetLastConditionId)
    .def("GetLastPropertiesId", &FiniteCellAuxiliaryUtility_GetLastPropertiesId)
    .def("AddElement", &FiniteCellAuxiliaryUtility_AddElement)
    .def("RemoveElement", &FiniteCellAuxiliaryUtility_RemoveElement)
    .def("MultithreadedRefineBy", &FiniteCellAuxiliaryUtility_MultithreadedRefineBy<RefinableTree, BRep>)
    .def("MultithreadedGeneratePhysicalIntegrationPoints", &FiniteCellAuxiliaryUtility_MultithreadedGeneratePhysicalIntegrationPoints<BaseMomentFittedQuadTreeSubCell, BRep>)
    .def("Print", pointer_to_PrintGeometry)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

