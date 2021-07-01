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
#include "custom_python3/add_ghost_penalty_utility_to_python.h"
#include "custom_utilities/ghost_penalty_utility.h"
#include "custom_utilities/skeleton_penalty_utility.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

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
        ModelPart& r_model_part, pybind11::list& list_elements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t lastCondId, Properties::Pointer pProperties, const int& echo_level)
{
    ModelPart::ElementsContainerType pElements;
    for (auto it : list_elements)
    {
        Element::Pointer p_elem = it.cast<Element::Pointer>();
        pElements.push_back(p_elem);
    }
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
        ModelPart& r_model_part, pybind11::list& list_elements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t lastCondId, Properties::Pointer pProperties,
        const int& nsampling, const int& configuration, const int& echo_level)
{
    ModelPart::ElementsContainerType pElements;
    for (auto it : list_elements)
    {
        Element::Pointer p_elem = it.cast<Element::Pointer>();
        pElements.push_back(p_elem);
    }
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

void FiniteCellApplication_AddGhostPenaltyUtilityToPython(pybind11::module& m)
{

    class_<GhostPenaltyUtility, GhostPenaltyUtility::Pointer>
    (m, "GhostPenaltyUtility")
    .def(init<>())
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

    class_<SkeletonPenaltyUtility, SkeletonPenaltyUtility::Pointer, GhostPenaltyUtility>
    (m, "SkeletonPenaltyUtility")
    .def(init<>())
    ;

}

}  // namespace Python.

}  // namespace Kratos.

