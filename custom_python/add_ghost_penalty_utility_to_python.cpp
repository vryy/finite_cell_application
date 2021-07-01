// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/element.h"
#include "custom_python/add_ghost_penalty_utility_to_python.h"
#include "custom_utilities/ghost_penalty_utility.h"
#include "custom_utilities/skeleton_penalty_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void GhostPenaltyUtility_ProbeNeighbourElements(GhostPenaltyUtility& rDummy, Element::Pointer p_elem)
{
    rDummy.ProbeNeighbourElements(p_elem);
}

ModelPart::ConditionsContainerType GhostPenaltyUtility_SetUpSurfacePenaltyConditions1(GhostPenaltyUtility& rDummy,
        Element::Pointer p_elem, GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t lastCondId, Properties::Pointer pProperties, const int& echo_level)
{
    return rDummy.SetUpSurfacePenaltyConditions(p_elem, p_sample_condition, r_brep, lastCondId, pProperties, echo_level);
}

ModelPart::ConditionsContainerType GhostPenaltyUtility_SetUpSurfacePenaltyConditions2(GhostPenaltyUtility& rDummy,
        ModelPart& r_model_part, GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t lastCondId, Properties::Pointer pProperties, const int& echo_level)
{
    return rDummy.SetUpSurfacePenaltyConditions(r_model_part, p_sample_condition, r_brep, lastCondId, pProperties, echo_level);
}

ModelPart::ConditionsContainerType GhostPenaltyUtility_SetUpSurfacePenaltyConditions3(GhostPenaltyUtility& rDummy,
        ModelPart& r_model_part, ModelPart::ElementsContainerType& pElements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t lastCondId, Properties::Pointer pProperties, const int& echo_level)
{
    return rDummy.SetUpSurfacePenaltyConditions(r_model_part, pElements, p_sample_condition, r_brep, lastCondId, pProperties, echo_level);
}

ModelPart::ConditionsContainerType GhostPenaltyUtility_SetUpSurfacePenaltyConditions4(GhostPenaltyUtility& rDummy,
        ModelPart& r_model_part, boost::python::list& list_elements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t lastCondId, Properties::Pointer pProperties, const int& echo_level)
{
    ModelPart::ElementsContainerType pElements;
    typedef boost::python::stl_input_iterator<Element::Pointer> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& p_elem,
        std::make_pair(iterator_value_type(list_elements), iterator_value_type() ) ) pElements.push_back(p_elem);
    return rDummy.SetUpSurfacePenaltyConditions(r_model_part, pElements, p_sample_condition, r_brep, lastCondId, pProperties, echo_level);
}

ModelPart::ConditionsContainerType GhostPenaltyUtility_SetUpSurfacePenaltyConditions5(GhostPenaltyUtility& rDummy,
        ModelPart& r_model_part, ModelPart::ElementsContainerType& pElements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t lastCondId, Properties::Pointer pProperties,
        const int& nsampling, const int& configuration, const int& echo_level)
{
    return rDummy.SetUpSurfacePenaltyConditions(r_model_part, pElements, p_sample_condition, r_brep, lastCondId, pProperties,
            nsampling, configuration, echo_level);
}

ModelPart::ConditionsContainerType GhostPenaltyUtility_SetUpSurfacePenaltyConditions6(GhostPenaltyUtility& rDummy,
        ModelPart& r_model_part, boost::python::list& list_elements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t lastCondId, Properties::Pointer pProperties,
        const int& nsampling, const int& configuration, const int& echo_level)
{
    ModelPart::ElementsContainerType pElements;
    typedef boost::python::stl_input_iterator<Element::Pointer> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& p_elem,
        std::make_pair(iterator_value_type(list_elements), iterator_value_type() ) ) pElements.push_back(p_elem);
    return rDummy.SetUpSurfacePenaltyConditions(r_model_part, pElements, p_sample_condition, r_brep, lastCondId, pProperties,
            nsampling, configuration, echo_level);
}

Condition::Pointer GhostPenaltyUtility_SetUpSurfacePenaltyCondition(GhostPenaltyUtility& rDummy,
        Element::Pointer p_element_1, Element::Pointer p_element_2, GhostPenaltyCondition::Pointer p_sample_condition,
        std::size_t lastCondId, Properties::Pointer pProperties, const int& echo_level)
{
    return rDummy.SetUpSurfacePenaltyCondition(p_element_1, p_element_2, p_sample_condition, lastCondId, pProperties, echo_level);
}

void GhostPenaltyUtility_ProbeShapeFunctionSecondDerivatives(GhostPenaltyUtility& rDummy, Element::Pointer p_element)
{
    rDummy.ProbeShapeFunctionSecondDerivatives(p_element->GetGeometry());
}

void FiniteCellApplication_AddGhostPenaltyUtilityToPython()
{

    class_<GhostPenaltyUtility, GhostPenaltyUtility::Pointer, boost::noncopyable>
    ("GhostPenaltyUtility", init<>())
    .def("ProbeNeighbourElements", &GhostPenaltyUtility_ProbeNeighbourElements)
    .def("SetUpSurfacePenaltyConditions", &GhostPenaltyUtility_SetUpSurfacePenaltyConditions1)
    .def("SetUpSurfacePenaltyConditions", &GhostPenaltyUtility_SetUpSurfacePenaltyConditions2)
    .def("SetUpSurfacePenaltyConditions", &GhostPenaltyUtility_SetUpSurfacePenaltyConditions3)
    .def("SetUpSurfacePenaltyConditions", &GhostPenaltyUtility_SetUpSurfacePenaltyConditions4)
    .def("SetUpSurfacePenaltyConditions", &GhostPenaltyUtility_SetUpSurfacePenaltyConditions5)
    .def("SetUpSurfacePenaltyConditions", &GhostPenaltyUtility_SetUpSurfacePenaltyConditions6)
    .def("SetUpSurfacePenaltyCondition", &GhostPenaltyUtility_SetUpSurfacePenaltyCondition)
    .def("ProbeShapeFunctionSecondDerivatives", &GhostPenaltyUtility_ProbeShapeFunctionSecondDerivatives)
    ;

    class_<SkeletonPenaltyUtility, SkeletonPenaltyUtility::Pointer, bases<GhostPenaltyUtility>, boost::noncopyable>
    ("SkeletonPenaltyUtility", init<>())
    ;

}

}  // namespace Python.

}  // namespace Kratos.

