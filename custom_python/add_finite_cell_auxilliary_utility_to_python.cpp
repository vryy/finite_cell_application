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
#include "custom_python/add_finite_cell_auxilliary_utility_to_python.h"
#include "custom_algebra/brep.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/finite_cell_auxilliary_utility.h"
#include "custom_utilities/moment_fitted_quad_tree_subcell.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

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

/// Create a line condition from sample_cond_name and from 2 nodes
Condition::Pointer FiniteCellAuxilliaryUtility_CreateCondition(FiniteCellAuxilliaryUtility& rDummy,
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

template<class TCellType, class TBRepType>
void FiniteCellAuxilliaryUtility_MultithreadedGeneratePhysicalIntegrationPoints(FiniteCellAuxilliaryUtility& rDummy,
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

void FiniteCellApplication_AddFiniteCellAuxilliaryUtilityToPython()
{

    void(FiniteCellAuxilliaryUtility::*pointer_to_Clean)(ModelPart&,
            ModelPart::ConditionsContainerType&, const int&) const = &FiniteCellAuxilliaryUtility::Clean;

    void(FiniteCellAuxilliaryUtility::*pointer_to_PrintGeometry)(Element::GeometryType::Pointer) const = &FiniteCellAuxilliaryUtility::Print;

    class_<FiniteCellAuxilliaryUtility, FiniteCellAuxilliaryUtility::Pointer, boost::noncopyable>
    ("FiniteCellAuxilliaryUtility", init<>())
    .def("CreateCondition", &FiniteCellAuxilliaryUtility_CreateCondition)
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
    .def("MultithreadedGeneratePhysicalIntegrationPoints", &FiniteCellAuxilliaryUtility_MultithreadedGeneratePhysicalIntegrationPoints<BaseMomentFittedQuadTreeSubCell, BRep>)
    .def("Print", pointer_to_PrintGeometry)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

