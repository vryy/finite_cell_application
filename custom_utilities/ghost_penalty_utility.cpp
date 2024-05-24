//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         finite_cell_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            6 Jan 2018
//


#include "custom_utilities/ghost_penalty_utility.h"
#include "structural_application/custom_utilities/sd_math_utils.h"

// #define ENABLE_DEBUG_GHOST_PENALTY
// #define ENABLE_DEBUG_GHOST_PENALTY_2

namespace Kratos
{

ModelPart::ConditionsContainerType GhostPenaltyUtility::SetUpSurfacePenaltyConditions(ModelPart& r_model_part,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t& lastCondId, Properties::Pointer pProperties, const int& echo_level) const
{
    // setup all the element neighbour couples in the model
    std::set<std::pair<std::size_t, std::size_t> > edges;
    for (ModelPart::ElementsContainerType::iterator it = r_model_part.Elements().begin();
            it != r_model_part.Elements().end(); ++it)
    {
#ifdef SD_APP_FORWARD_COMPATIBILITY
        GlobalPointersVector<Element>& rNeighbours = it->GetValue(NEIGHBOUR_ELEMENTS);
#else
        WeakPointerVector<Element>& rNeighbours = it->GetValue(NEIGHBOUR_ELEMENTS);
#endif

        for (std::size_t i = 0; i < rNeighbours.size(); ++i)
        {
            const std::size_t& id1 = it->Id();
#ifdef SD_APP_FORWARD_COMPATIBILITY
            const std::size_t& id2 = rNeighbours(i)->Id();
#else
            const std::size_t& id2 = rNeighbours(i).lock()->Id();
#endif
            if (id1 < id2)
            {
                edges.insert(std::make_pair(id1, id2));
            }
            else if (id2 > id1)
            {
                edges.insert(std::make_pair(id2, id1));
            }
        }
    }

    // setup the ghost penalty condition from the edges
    ModelPart::ConditionsContainerType pNewConditions;

    for (std::set<std::pair<std::size_t, std::size_t> >::iterator it = edges.begin();
            it != edges.end(); ++it)
    {
        Element::Pointer p_element_1 = r_model_part.pGetElement(it->first);
        Element::Pointer p_element_2 = r_model_part.pGetElement(it->second);

        Condition::Pointer pNewCond = SetUpSurfacePenaltyCondition(p_element_1, p_element_2,
                                      p_sample_condition, r_brep, lastCondId, pProperties, echo_level);
        if (pNewCond != NULL)
        {
            pNewConditions.push_back(pNewCond);
        }
    }

    for (ModelPart::ConditionsContainerType::ptr_iterator it = pNewConditions.ptr_begin();
            it != pNewConditions.ptr_end(); ++it)
    {
        r_model_part.Conditions().push_back(*it);
    }

    if (echo_level > 0)
        std::cout << __FUNCTION__ << " completed: " << pNewConditions.size() << " new ghost conditions"
                  << " of type " << typeid(*p_sample_condition).name() << " are created and added to model_part " << r_model_part.Name()
                  << std::endl;

    return pNewConditions;
}

ModelPart::ConditionsContainerType GhostPenaltyUtility::SetUpSurfacePenaltyConditions(ModelPart& r_model_part,
        ModelPart::ElementsContainerType& pElements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t& lastCondId, Properties::Pointer pProperties, const int& echo_level) const
{
    // setup all the element neighbour couples in the list of elements
    std::set<std::pair<std::size_t, std::size_t> > edges;
    for (ModelPart::ElementsContainerType::iterator it = pElements.begin();
            it != pElements.end(); ++it)
    {
#ifdef SD_APP_FORWARD_COMPATIBILITY
        GlobalPointersVector<Element>& rNeighbours = it->GetValue(NEIGHBOUR_ELEMENTS);
#else
        WeakPointerVector<Element>& rNeighbours = it->GetValue(NEIGHBOUR_ELEMENTS);
#endif

        if (echo_level > 3)
        {
            std::cout << "neighbours of element " << it->Id() << ":";
            for (std::size_t i = 0; i < rNeighbours.size(); ++i)
            {
#ifdef SD_APP_FORWARD_COMPATIBILITY
                std::cout << " " << rNeighbours(i)->Id();
#else
                std::cout << " " << rNeighbours(i).lock()->Id();
#endif
            }
            std::cout << std::endl;
        }

        for (std::size_t i = 0; i < rNeighbours.size(); ++i)
        {
#ifdef SD_APP_FORWARD_COMPATIBILITY
            const std::size_t& id2 = rNeighbours(i)->Id();
#else
            const std::size_t& id2 = rNeighbours(i).lock()->Id();
#endif
            if (pElements.find(id2) == pElements.end())
            {
                continue;
            }

            const std::size_t& id1 = it->Id();
            if (id1 < id2)
            {
                edges.insert(std::make_pair(id1, id2));
            }
            else if (id2 > id1)
            {
                edges.insert(std::make_pair(id2, id1));
            }
        }
    }

    if (echo_level > 3)
    {
        std::cout << " neighbour couples:";
        for (std::set<std::pair<std::size_t, std::size_t> >::iterator it = edges.begin(); it != edges.end(); ++it)
        {
            std::cout << " (" << it->first << ", " << it->second << ")";
        }
        std::cout << std::endl;
    }

    // setup the ghost penalty condition from the edges
    ModelPart::ConditionsContainerType pNewConditions;

    for (std::set<std::pair<std::size_t, std::size_t> >::iterator it = edges.begin();
            it != edges.end(); ++it)
    {
        Element::Pointer p_element_1 = pElements(it->first);
        Element::Pointer p_element_2 = pElements(it->second);

        Condition::Pointer pNewCond = SetUpSurfacePenaltyCondition(p_element_1, p_element_2,
                                      p_sample_condition, r_brep, lastCondId, pProperties, echo_level);
        if (pNewCond != NULL)
        {
            if (echo_level > 2)
            {
                std::cout << "ghost penalty condition " << pNewCond->Info() << " is created" << std::endl;
            }
            pNewConditions.push_back(pNewCond);
        }
    }

    for (ModelPart::ConditionsContainerType::ptr_iterator it = pNewConditions.ptr_begin();
            it != pNewConditions.ptr_end(); ++it)
    {
        r_model_part.Conditions().push_back(*it);
    }

    if (echo_level > 0)
        std::cout << __FUNCTION__ << " on " << pElements.size() << " elements completed: " << pNewConditions.size() << " new ghost conditions"
                  << " of type " << typeid(*p_sample_condition).name() << " are created and added to model_part " << r_model_part.Name()
                  << std::endl;

    return pNewConditions;
}

ModelPart::ConditionsContainerType GhostPenaltyUtility::SetUpSurfacePenaltyConditions(ModelPart& r_model_part,
        ModelPart::ElementsContainerType& pElements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t& lastCondId, Properties::Pointer pProperties,
        const int& nsampling, const int& configuration, const int& echo_level) const
{
    // setup all the element neighbour couples in the list of elements
    std::set<std::pair<std::size_t, std::size_t> > edges;
    for (ModelPart::ElementsContainerType::iterator it = pElements.begin();
            it != pElements.end(); ++it)
    {
#ifdef SD_APP_FORWARD_COMPATIBILITY
        GlobalPointersVector<Element>& rNeighbours = it->GetValue(NEIGHBOUR_ELEMENTS);
#else
        WeakPointerVector<Element>& rNeighbours = it->GetValue(NEIGHBOUR_ELEMENTS);
#endif

        if (echo_level > 3)
        {
            std::cout << "neighbours of element " << it->Id() << ":";
            for (std::size_t i = 0; i < rNeighbours.size(); ++i)
            {
#ifdef SD_APP_FORWARD_COMPATIBILITY
                std::cout << " " << rNeighbours(i)->Id();
#else
                std::cout << " " << rNeighbours(i).lock()->Id();
#endif
            }
            std::cout << std::endl;
        }

        for (std::size_t i = 0; i < rNeighbours.size(); ++i)
        {
#ifdef SD_APP_FORWARD_COMPATIBILITY
            const std::size_t& id2 = rNeighbours(i)->Id();
#else
            const std::size_t& id2 = rNeighbours(i).lock()->Id();
#endif
            if (pElements.find(id2) == pElements.end())
            {
                continue;
            }

            const std::size_t& id1 = it->Id();
            if (id1 < id2)
            {
                edges.insert(std::make_pair(id1, id2));
            }
            else if (id2 > id1)
            {
                edges.insert(std::make_pair(id2, id1));
            }
        }
    }

    if (echo_level > 3)
    {
        std::cout << " neighbour couples:";
        for (std::set<std::pair<std::size_t, std::size_t> >::iterator it = edges.begin(); it != edges.end(); ++it)
        {
            std::cout << " (" << it->first << ", " << it->second << ")";
        }
        std::cout << std::endl;
    }

    // setup the ghost penalty condition from the edges
    ModelPart::ConditionsContainerType pNewConditions;

    for (std::set<std::pair<std::size_t, std::size_t> >::iterator it = edges.begin();
            it != edges.end(); ++it)
    {
        Element::Pointer p_element_1 = pElements(it->first);
        Element::Pointer p_element_2 = pElements(it->second);

        Condition::Pointer pNewCond = SetUpSurfacePenaltyCondition(p_element_1, p_element_2,
                                      p_sample_condition, r_brep, lastCondId, pProperties, nsampling, configuration, echo_level);
        if (pNewCond != NULL)
        {
            if (echo_level > 2)
            {
                std::cout << "ghost penalty condition " << pNewCond->Info() << " is created" << std::endl;
            }
            pNewConditions.push_back(pNewCond);
        }
    }

    for (ModelPart::ConditionsContainerType::ptr_iterator it = pNewConditions.ptr_begin();
            it != pNewConditions.ptr_end(); ++it)
    {
        r_model_part.Conditions().push_back(*it);
    }

    if (echo_level > 0)
        std::cout << __FUNCTION__ << " on " << pElements.size() << " elements completed: " << pNewConditions.size() << " new ghost conditions"
                  << " of type " << typeid(*p_sample_condition).name() << " are created and added to model_part " << r_model_part.Name()
                  << std::endl;

    return pNewConditions;
}

ModelPart::ConditionsContainerType GhostPenaltyUtility::SetUpSurfacePenaltyConditions(Element::Pointer p_element,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep, std::size_t& lastCondId,
        Properties::Pointer pProperties, const int& echo_level) const
{
    // firstly obtain all neighbour elements of the current element
#ifdef SD_APP_FORWARD_COMPATIBILITY
    GlobalPointersVector<Element>& rNeighbours = p_element->GetValue(NEIGHBOUR_ELEMENTS);
#else
    WeakPointerVector<Element>& rNeighbours = p_element->GetValue(NEIGHBOUR_ELEMENTS);
#endif

    ModelPart::ConditionsContainerType pNewConditions;

    // for each neighbour, find the common edge
    for (std::size_t i = 0; i < rNeighbours.size(); ++i)
    {
        if (rNeighbours[i].Id() != p_element->Id())
        {
            Condition::Pointer pNewCond = SetUpSurfacePenaltyCondition(p_element,
#ifdef SD_APP_FORWARD_COMPATIBILITY
                                          * rNeighbours(i),
#else
                                          rNeighbours(i).lock(),
#endif
                                          p_sample_condition, r_brep, lastCondId, pProperties, echo_level);
            if (pNewCond != NULL)
            {
                pNewConditions.push_back(pNewCond);
            }
        }
    }

    if (echo_level > 0)
        std::cout << __FUNCTION__ << " completed: " << pNewConditions.size() << " new ghost conditions"
                  << " of type " << typeid(*p_sample_condition).name() << " are created" << std::endl;

    return pNewConditions;
}


Condition::Pointer GhostPenaltyUtility::SetUpSurfacePenaltyCondition(Element::Pointer p_element_1, Element::Pointer p_element_2,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep, std::size_t& lastCondId,
        Properties::Pointer pProperties, const int& echo_level) const
{
    Condition::Pointer pNewCond = NULL;

    std::pair<GeometryData::KratosGeometryType, std::vector<std::size_t> > edge = FindCommonFace(p_element_1->GetGeometry(), p_element_2->GetGeometry());

    if (edge.first == GeometryData::KratosGeometryType::Kratos_generic_type)
    {
        return pNewCond;
    }

    if (edge.first == GeometryData::KratosGeometryType::Kratos_Line2D2)
    {
        if (edge.second.size () != 2)
        {
            return pNewCond;
        }
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Line2D3)
    {
        if (edge.second.size () != 3)
        {
            return pNewCond;
        }
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
    {
        if (edge.second.size () != 3)
        {
            return pNewCond;
        }
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Triangle3D6)
    {
        if (edge.second.size () != 6)
        {
            return pNewCond;
        }
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)
    {
        if (edge.second.size () != 4)
        {
            return pNewCond;
        }
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8)
    {
        if (edge.second.size () != 8)
        {
            return pNewCond;
        }
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9)
    {
        if (edge.second.size () != 9)
        {
            return pNewCond;
        }
    }

    // create the edge geometry
    typename Element::NodesArrayType temp_nodes;
    for (std::size_t j = 0; j < edge.second.size(); ++j)
    {
        temp_nodes.push_back(p_element_1->GetGeometry().pGetPoint(edge.second[j]));
    }

    // if (temp_nodes.size() == 0)
    // {
    //     KRATOS_WATCH(p_element_1->Id())
    //     KRATOS_WATCH(typeid(p_element_1->GetGeometry()).name())
    //     KRATOS_WATCH(p_element_2->Id())
    //     KRATOS_WATCH(typeid(p_element_2->GetGeometry()).name())
    //     std::cout << "Element " << p_element_1->Id() << " nodes:";
    //     for (std::size_t i = 0; i < p_element_1->GetGeometry().size(); ++i)
    //         std::cout << " " << p_element_1->GetGeometry()[i].Id();
    //     std::cout << std::endl;
    //     std::cout << "Element " << p_element_2->Id() << " nodes:";
    //     for (std::size_t i = 0; i < p_element_2->GetGeometry().size(); ++i)
    //         std::cout << " " << p_element_2->GetGeometry()[i].Id();
    //     std::cout << std::endl;
    //     KRATOS_THROW_ERROR(std::logic_error, "The two elements do not have common face", "")
    // }

    GeometryType::Pointer p_temp_geometry;

    if (edge.first == GeometryData::KratosGeometryType::Kratos_Line2D2)
    {
        p_temp_geometry = GeometryType::Pointer(new Line2D2<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Line2D3)
    {
        p_temp_geometry = GeometryType::Pointer(new Line2D3<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
    {
        p_temp_geometry = GeometryType::Pointer(new Triangle3D3<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Triangle3D6)
    {
        p_temp_geometry = GeometryType::Pointer(new Triangle3D6<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)
    {
        p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D4<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8)
    {
        p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D8<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9)
    {
        p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D9<NodeType>(temp_nodes));
    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "Unknown geometry type", static_cast<int>(edge.first))

        // check if this edge is cut by the brep or totally inside. If yes, then the new ghost condition is created.
        bool is_ghost = false;

    int configuration = 0; // check the cut status in the initial configuration
    int stat1 = r_brep.CutStatus(p_element_1->GetGeometry(), configuration);
    int stat2 = r_brep.CutStatus(p_element_2->GetGeometry(), configuration);

    if (stat1 == BRep::_CUT)
    {
        if (stat2 == BRep::_CUT || stat2 == BRep::_IN)
        {
            is_ghost = true;
        }
    }
    else if (stat1 == BRep::_IN)
    {
        if (stat2 == BRep::_CUT)
        {
            is_ghost = true;
        }
    }

    if (is_ghost)
    {
        // create the ghost penalty condition
        pNewCond = p_sample_condition->Create(++lastCondId, p_temp_geometry, p_element_1, p_element_2, pProperties);
#ifndef SD_APP_FORWARD_COMPATIBILITY
        pNewCond->SetValue(IS_INACTIVE, false);
#endif
        pNewCond->Set(ACTIVE, true);
    }

    return pNewCond;
}


Condition::Pointer GhostPenaltyUtility::SetUpSurfacePenaltyCondition(Element::Pointer p_element_1, Element::Pointer p_element_2,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep, std::size_t& lastCondId,
        Properties::Pointer pProperties, const int& nsampling, const int& configuration, const int& echo_level) const
{
    Condition::Pointer pNewCond = NULL;

    std::pair<GeometryData::KratosGeometryType, std::vector<std::size_t> > edge = FindCommonFace(p_element_1->GetGeometry(), p_element_2->GetGeometry());

    if (edge.first == GeometryData::KratosGeometryType::Kratos_generic_type)
    {
        return pNewCond;
    }

    if (edge.first == GeometryData::KratosGeometryType::Kratos_Line2D2)
    {
        if (edge.second.size () != 2)
        {
            return pNewCond;
        }
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Line2D3)
    {
        if (edge.second.size () != 3)
        {
            return pNewCond;
        }
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
    {
        if (edge.second.size () != 3)
        {
            return pNewCond;
        }
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Triangle3D6)
    {
        if (edge.second.size () != 6)
        {
            return pNewCond;
        }
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)
    {
        if (edge.second.size () != 4)
        {
            return pNewCond;
        }
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8)
    {
        if (edge.second.size () != 8)
        {
            return pNewCond;
        }
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9)
    {
        if (edge.second.size () != 9)
        {
            return pNewCond;
        }
    }

    // create the edge geometry
    typename Element::NodesArrayType temp_nodes;
    for (std::size_t j = 0; j < edge.second.size(); ++j)
    {
        temp_nodes.push_back(p_element_1->GetGeometry().pGetPoint(edge.second[j]));
    }

    // if (temp_nodes.size() == 0)
    // {
    //     KRATOS_WATCH(p_element_1->Id())
    //     KRATOS_WATCH(typeid(p_element_1->GetGeometry()).name())
    //     KRATOS_WATCH(p_element_2->Id())
    //     KRATOS_WATCH(typeid(p_element_2->GetGeometry()).name())
    //     std::cout << "Element " << p_element_1->Id() << " nodes:";
    //     for (std::size_t i = 0; i < p_element_1->GetGeometry().size(); ++i)
    //         std::cout << " " << p_element_1->GetGeometry()[i].Id();
    //     std::cout << std::endl;
    //     std::cout << "Element " << p_element_2->Id() << " nodes:";
    //     for (std::size_t i = 0; i < p_element_2->GetGeometry().size(); ++i)
    //         std::cout << " " << p_element_2->GetGeometry()[i].Id();
    //     std::cout << std::endl;
    //     KRATOS_THROW_ERROR(std::logic_error, "The two elements do not have common face", "")
    // }

    GeometryType::Pointer p_temp_geometry;

    if (edge.first == GeometryData::KratosGeometryType::Kratos_Line2D2)
    {
        p_temp_geometry = GeometryType::Pointer(new Line2D2<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Line2D3)
    {
        p_temp_geometry = GeometryType::Pointer(new Line2D3<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
    {
        p_temp_geometry = GeometryType::Pointer(new Triangle3D3<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Triangle3D6)
    {
        p_temp_geometry = GeometryType::Pointer(new Triangle3D6<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)
    {
        p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D4<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8)
    {
        p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D8<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9)
    {
        p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D9<NodeType>(temp_nodes));
    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "Unknown geometry type", static_cast<int>(edge.first))

        // check if this edge is cut by the brep or totally inside. If yes, then the new ghost condition is created.
        bool is_ghost = false;

    int stat1 = r_brep.CutStatusBySampling(p_element_1->GetGeometry(), nsampling, configuration);
    int stat2 = r_brep.CutStatusBySampling(p_element_2->GetGeometry(), nsampling, configuration);

    if (stat1 == BRep::_CUT)
    {
        if (stat2 == BRep::_CUT || stat2 == BRep::_IN)
        {
            is_ghost = true;
        }
    }
    else if (stat1 == BRep::_IN)
    {
        if (stat2 == BRep::_CUT)
        {
            is_ghost = true;
        }
    }

    if (is_ghost)
    {
        // create the ghost penalty condition
        pNewCond = p_sample_condition->Create(++lastCondId, p_temp_geometry, p_element_1, p_element_2, pProperties);
#ifndef SD_APP_FORWARD_COMPATIBILITY
        pNewCond->SetValue(IS_INACTIVE, false);
#endif
        pNewCond->Set(ACTIVE, true);
    }

    return pNewCond;
}


Condition::Pointer GhostPenaltyUtility::SetUpSurfacePenaltyCondition(Element::Pointer p_element_1, Element::Pointer p_element_2,
        GhostPenaltyCondition::Pointer p_sample_condition, std::size_t& lastCondId,
        Properties::Pointer pProperties, const int& echo_level) const
{
    Condition::Pointer pNewCond = NULL;

    std::pair<GeometryData::KratosGeometryType, std::vector<std::size_t> > edge = FindCommonFace(p_element_1->GetGeometry(), p_element_2->GetGeometry());

    if (edge.first == GeometryData::KratosGeometryType::Kratos_generic_type)
    {
        return pNewCond;
    }

    // create the edge geometry
    typename Element::NodesArrayType temp_nodes;
    for (std::size_t j = 0; j < edge.second.size(); ++j)
    {
        temp_nodes.push_back(p_element_1->GetGeometry().pGetPoint(edge.second[j]));
    }

    GeometryType::Pointer p_temp_geometry;

    if (edge.first == GeometryData::KratosGeometryType::Kratos_Line2D2)
    {
        p_temp_geometry = GeometryType::Pointer(new Line2D2<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Line2D3)
    {
        p_temp_geometry = GeometryType::Pointer(new Line2D3<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
    {
        p_temp_geometry = GeometryType::Pointer(new Triangle3D3<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Triangle3D6)
    {
        p_temp_geometry = GeometryType::Pointer(new Triangle3D6<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)
    {
        p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D4<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8)
    {
        p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D8<NodeType>(temp_nodes));
    }
    else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9)
    {
        p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D9<NodeType>(temp_nodes));
    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "Unknown geometry type", static_cast<int>(edge.first))

        // create the ghost penalty condition
        pNewCond = p_sample_condition->Create(++lastCondId, p_temp_geometry, p_element_1, p_element_2, pProperties);
#ifndef SD_APP_FORWARD_COMPATIBILITY
    pNewCond->SetValue(IS_INACTIVE, false);
#endif
    pNewCond->Set(ACTIVE, true);

    if (echo_level > 0)
        std::cout << __FUNCTION__ << " completed, new ghost condition of type "
                  << typeid(*p_sample_condition).name() << " is created" << std::endl;

    return pNewCond;
}

void GhostPenaltyUtility::ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{
    if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::ComputeShapeFunction(N, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D6)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::ComputeShapeFunction(N, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::ComputeShapeFunction(N, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::ComputeShapeFunction(N, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::ComputeShapeFunction(N, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::ComputeShapeFunction(N, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::ComputeShapeFunction(N, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::ComputeShapeFunction(N, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::ComputeShapeFunction(N, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D27)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::ComputeShapeFunction(N, r_element_geometry, r_edge_geometry, integration_points);
    }
    else
    {
        std::stringstream ss;
        ss << __FUNCTION__ << "does not work on geometry " << static_cast<int>(r_element_geometry.GetGeometryType());
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }
}


void GhostPenaltyUtility::ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{
    if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::ComputeShapeFunctionGradient(DN_DX, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D6)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::ComputeShapeFunctionGradient(DN_DX, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::ComputeShapeFunctionGradient(DN_DX, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::ComputeShapeFunctionGradient(DN_DX, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::ComputeShapeFunctionGradient(DN_DX, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::ComputeShapeFunctionGradient(DN_DX, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::ComputeShapeFunctionGradient(DN_DX, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::ComputeShapeFunctionGradient(DN_DX, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::ComputeShapeFunctionGradient(DN_DX, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D27)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::ComputeShapeFunctionGradient(DN_DX, r_element_geometry, r_edge_geometry, integration_points);
    }
    else
    {
        std::stringstream ss;
        ss << __FUNCTION__ << "does not work on geometry " << static_cast<int>(r_element_geometry.GetGeometryType());
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }
}

void GhostPenaltyUtility::ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{
    if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::ComputeShapeFunctionNormalGradient(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D6)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::ComputeShapeFunctionNormalGradient(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::ComputeShapeFunctionNormalGradient(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::ComputeShapeFunctionNormalGradient(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::ComputeShapeFunctionNormalGradient(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::ComputeShapeFunctionNormalGradient(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::ComputeShapeFunctionNormalGradient(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::ComputeShapeFunctionNormalGradient(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::ComputeShapeFunctionNormalGradient(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D27)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::ComputeShapeFunctionNormalGradient(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else
    {
        std::stringstream ss;
        ss << __FUNCTION__ << "does not work on geometry " << static_cast<int>(r_element_geometry.GetGeometryType());
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }
}

void GhostPenaltyUtility::ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{
    if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::ComputeShapeFunctionNormalSecondDerivatives(d2Ndn2, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D6)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::ComputeShapeFunctionNormalSecondDerivatives(d2Ndn2, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::ComputeShapeFunctionNormalSecondDerivatives(d2Ndn2, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::ComputeShapeFunctionNormalSecondDerivatives(d2Ndn2, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::ComputeShapeFunctionNormalSecondDerivatives(d2Ndn2, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::ComputeShapeFunctionNormalSecondDerivatives(d2Ndn2, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::ComputeShapeFunctionNormalSecondDerivatives(d2Ndn2, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::ComputeShapeFunctionNormalSecondDerivatives(d2Ndn2, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::ComputeShapeFunctionNormalSecondDerivatives(d2Ndn2, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D27)
    {
        GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::ComputeShapeFunctionNormalSecondDerivatives(d2Ndn2, r_element_geometry, r_edge_geometry, integration_points);
    }
    else
    {
        std::stringstream ss;
        ss << __FUNCTION__ << "does not work on geometry " << static_cast<int>(r_element_geometry.GetGeometryType());
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }
}

int GhostPenaltyUtility::FindSide(GeometryType& r_element_geometry, GeometryType& r_edge_geometry)
{
    if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
    {
        return FindSide_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3> >::Execute(r_element_geometry, r_edge_geometry);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D6)
    {
        return FindSide_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6> >::Execute(r_element_geometry, r_edge_geometry);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
    {
        return FindSide_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4> >::Execute(r_element_geometry, r_edge_geometry);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8)
    {
        return FindSide_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8> >::Execute(r_element_geometry, r_edge_geometry);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9)
    {
        return FindSide_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9> >::Execute(r_element_geometry, r_edge_geometry);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
    {
        return FindSide_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4> >::Execute(r_element_geometry, r_edge_geometry);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10)
    {
        return FindSide_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10> >::Execute(r_element_geometry, r_edge_geometry);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8)
    {
        return FindSide_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8> >::Execute(r_element_geometry, r_edge_geometry);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20)
    {
        return FindSide_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20> >::Execute(r_element_geometry, r_edge_geometry);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D27)
    {
        return FindSide_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27> >::Execute(r_element_geometry, r_edge_geometry);
    }
    else
    {
        std::stringstream ss;
        ss << __FUNCTION__ << "does not work on geometry " << static_cast<int>(r_element_geometry.GetGeometryType());
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    return -1;
}

const int* GhostPenaltyUtility::Faces(GeometryType& r_element_geometry, const std::size_t& side)
{
    if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::Faces(side);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D6)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::Faces(side);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::Faces(side);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::Faces(side);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::Faces(side);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::Faces(side);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::Faces(side);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::Faces(side);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::Faces(side);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D27)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::Faces(side);
    }
    else
    {
        std::stringstream ss;
        ss << __FUNCTION__ << "does not work on geometry " << static_cast<int>(r_element_geometry.GetGeometryType());
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    return NULL;
}

GhostPenaltyUtility::IntegrationPointType GhostPenaltyUtility::ComputeIntegrationPoint(GhostPenaltyUtility::GeometryType& r_element_geometry, const int& side, const GhostPenaltyUtility::IntegrationPointType& edge_integration_point)
{
    if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::ComputeIntegrationPoint(side, edge_integration_point);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D6)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::ComputeIntegrationPoint(side, edge_integration_point);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::ComputeIntegrationPoint(side, edge_integration_point);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::ComputeIntegrationPoint(side, edge_integration_point);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::ComputeIntegrationPoint(side, edge_integration_point);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::ComputeIntegrationPoint(side, edge_integration_point);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::ComputeIntegrationPoint(side, edge_integration_point);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::ComputeIntegrationPoint(side, edge_integration_point);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::ComputeIntegrationPoint(side, edge_integration_point);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D27)
    {
        return GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::ComputeIntegrationPoint(side, edge_integration_point);
    }
    else
    {
        std::stringstream ss;
        ss << __FUNCTION__ << "does not work on geometry " << static_cast<int>(r_element_geometry.GetGeometryType());
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }
}

void GhostPenaltyUtility::ProbeNeighbourElements(Element::Pointer p_element)
{
#ifdef SD_APP_FORWARD_COMPATIBILITY
    GlobalPointersVector<Element>& rNeighbours = p_element->GetValue(NEIGHBOUR_ELEMENTS);
#else
    WeakPointerVector<Element>& rNeighbours = p_element->GetValue(NEIGHBOUR_ELEMENTS);
#endif
    std::cout << "Neighbour elements of element " << p_element->Id() << ":";
    for (std::size_t i = 0; i < rNeighbours.size(); ++i)
    {
        std::cout << " " << rNeighbours[i].Id();
    }
    std::cout << std::endl;
}

void GhostPenaltyUtility::ProbeShapeFunctionSecondDerivatives(GeometryType& r_element_geometry)
{

    const IntegrationPointsArrayType& integration_points = r_element_geometry.IntegrationPoints(r_element_geometry.GetDefaultIntegrationMethod());

    for (std::size_t PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        std::cout << "At point " << integration_points[PointNumber] << std::endl;
        std::vector<Vector> D2N_DX2;
        D2N_DX2 = GhostPenalty_Helper::ComputeShapeFunctionSecondDerivatives(D2N_DX2, r_element_geometry, integration_points[PointNumber]);
        for (std::size_t i = 0; i < D2N_DX2.size(); ++i)
        {
            std::cout << " shape function " << i << " 2nd-der: " << D2N_DX2[i] << std::endl;
        }
    }
}

std::pair<GeometryData::KratosGeometryType, std::vector<std::size_t> > GhostPenaltyUtility::FindCommonFace(Element::GeometryType& r_geom_1, Element::GeometryType& r_geom_2)
{
    if (r_geom_1.GetGeometryType() != r_geom_2.GetGeometryType())
        KRATOS_THROW_ERROR(std::logic_error, "The geometry type is not the same", "")

        if (r_geom_1.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
        {
            return std::make_pair(GeometryData::KratosGeometryType::Kratos_Line2D2, GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::FindCommonFace(r_geom_1, r_geom_2));
        }
        else if (r_geom_1.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D6)
        {
            return std::make_pair(GeometryData::KratosGeometryType::Kratos_Line2D3, GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::FindCommonFace(r_geom_1, r_geom_2));
        }
        else if (r_geom_1.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
        {
            return std::make_pair(GeometryData::KratosGeometryType::Kratos_Line2D2, GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::FindCommonFace(r_geom_1, r_geom_2));
        }
        else if (r_geom_1.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8)
        {
            return std::make_pair(GeometryData::KratosGeometryType::Kratos_Line2D3, GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::FindCommonFace(r_geom_1, r_geom_2));
        }
        else if (r_geom_1.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9)
        {
            return std::make_pair(GeometryData::KratosGeometryType::Kratos_Line2D3, GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::FindCommonFace(r_geom_1, r_geom_2));
        }
        else if (r_geom_1.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
        {
            return std::make_pair(GeometryData::KratosGeometryType::Kratos_Triangle3D3, GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::FindCommonFace(r_geom_1, r_geom_2));
        }
        else if (r_geom_1.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10)
        {
            return std::make_pair(GeometryData::KratosGeometryType::Kratos_Triangle3D6, GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::FindCommonFace(r_geom_1, r_geom_2));
        }
        else if (r_geom_1.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8)
        {
            return std::make_pair(GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::FindCommonFace(r_geom_1, r_geom_2));
        }
        else if (r_geom_1.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20)
        {
            return std::make_pair(GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::FindCommonFace(r_geom_1, r_geom_2));
        }
        else if (r_geom_1.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D27)
        {
            return std::make_pair(GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9, GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::FindCommonFace(r_geom_1, r_geom_2));
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")

            return std::make_pair(GeometryData::KratosGeometryType::Kratos_generic_type, std::vector<std::size_t> {});
}

///////////////////////////////////////////////////////////////////////////////

bool GhostPenalty_Helper::IsBelongedToGeometry(const std::size_t& node_id, GeometryType& r_geom)
{
    for (std::size_t i = 0; i < r_geom.size(); ++i)
        if (r_geom[i].Id() == node_id)
        {
            return true;
        }
    return false;
}

bool GhostPenalty_Helper::IsSame(GeometryType& r_geom, const std::vector<std::size_t>& nodes)
{
    std::set<std::size_t> nodes_set(nodes.begin(), nodes.end());
    for (std::size_t i = 0; i < r_geom.size(); ++i)
    {
        nodes_set.erase(r_geom[i].Id());
    }
    return nodes_set.size() == 0;
}

bool GhostPenalty_Helper::IsSame(GeometryType& r_edge_geometry,
                                 GeometryType& r_element_geometry, const int* nodes_on_edge)
{
    for (std::size_t i = 0; i < r_edge_geometry.size(); ++i)
    {
        bool found = false;
        for (std::size_t j = 0; j < r_edge_geometry.size(); ++j)
        {
            if (r_edge_geometry[i].Id() == r_element_geometry[nodes_on_edge[j]].Id())
            {
                found = true;
                break;
            }
        }
        if (!found) { return false; }
    }
    return true;
}

bool GhostPenalty_Helper::BuildMapEdgeNodeIndexToElementNodeIndex(std::map<std::size_t, std::size_t>& map_edge_node_index_to_element_node_index,
        GeometryType& r_edge_geometry, GeometryType& r_element_geometry, const int* nodes_on_edge)
{
    for (std::size_t i = 0; i < r_edge_geometry.size(); ++i)
    {
        bool found = false;
        for (std::size_t j = 0; j < r_edge_geometry.size(); ++j)
        {
            if (r_edge_geometry[i].Id() == r_element_geometry[nodes_on_edge[j]].Id())
            {
                map_edge_node_index_to_element_node_index[i] = nodes_on_edge[j];
                found = true;
                break;
            }
        }
        if (!found) { return false; }
    }
    return true;
}

void GhostPenalty_Helper::ComputeLocalFrame(Vector& t, Vector& n, GeometryType& r_edge_geometry, const IntegrationPointType& integration_point)
{
    if (r_edge_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D2)
    {
        double length = norm_2(r_edge_geometry[0] - r_edge_geometry[1]);
        t(0) = (r_edge_geometry[1].X0() - r_edge_geometry[0].X0()) / length;
        t(1) = (r_edge_geometry[1].Y0() - r_edge_geometry[0].Y0()) / length;
        n(0) = -t(1);
        n(1) = t(0);
    }
    else if (r_edge_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D3)
    {
        Matrix aux;
        aux = r_edge_geometry.ShapeFunctionsLocalGradients(aux, integration_point);
        noalias(t) = ZeroVector(2);
        noalias(n) = ZeroVector(2);
        for (std::size_t i = 0; i < r_edge_geometry.size(); ++i)
        {
            t(0) += aux(i, 0) * r_edge_geometry[i].X0();
            t(1) += aux(i, 0) * r_edge_geometry[i].Y0();
        }
        double length = norm_2(t);
        t(0) = t(0) / length;
        t(1) = t(1) / length;
        n(0) = -t(1);
        n(1) = t(0);
    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "The geometry type is invalid:", static_cast<int>(r_edge_geometry.GetGeometryType()))
    }

void GhostPenalty_Helper::ComputeLocalFrame(Vector& t1, Vector& t2, Vector& n, GeometryType& r_face_geometry, const IntegrationPointType& integration_point)
{
    if (r_face_geometry.WorkingSpaceDimension() == 3)
    {
        Matrix aux;
        aux = r_face_geometry.ShapeFunctionsLocalGradients(aux, integration_point);
        noalias(t1) = ZeroVector(3);
        noalias(t2) = ZeroVector(3);
        noalias(n) = ZeroVector(3);
        for (std::size_t i = 0; i < r_face_geometry.size(); ++i)
        {
            t1(0) += aux(i, 0) * r_face_geometry[i].X0();
            t1(1) += aux(i, 0) * r_face_geometry[i].Y0();
            t1(2) += aux(i, 0) * r_face_geometry[i].Z0();
            t2(0) += aux(i, 1) * r_face_geometry[i].X0();
            t2(1) += aux(i, 1) * r_face_geometry[i].Y0();
            t2(2) += aux(i, 1) * r_face_geometry[i].Z0();
        }
        // Note that t1 and t2 is not orthogonal
        t1 *= (1.0 / norm_2(t1));
        t2 *= (1.0 / norm_2(t2));
        noalias(n) = MathUtils<double>::CrossProduct(t1, t2);
        n *= (1.0 / norm_2(n));
    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "The geometry type is invalid:", static_cast<int>(r_face_geometry.GetGeometryType()))
    }

Element::GeometryType::IntegrationPointType GhostPenalty_Helper::ComputeIntegrationPointOnTriSide(const int& side, const Element::GeometryType::IntegrationPointType& edge_integration_point)
{
    // TODO
    KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "not yet implemented")
}

Element::GeometryType::IntegrationPointType GhostPenalty_Helper::ComputeIntegrationPointOnQuadSide(const int& side, const Element::GeometryType::IntegrationPointType& edge_integration_point)
{
    IntegrationPointType integration_point;

    if (side == 0) // xi: n, eta: -1
    {
        integration_point.X() = edge_integration_point.X();
        integration_point.Y() = -1.0;
    }
    else if (side == 1) // xi: 1, eta: n
    {
        integration_point.X() = 1.0;
        integration_point.Y() = edge_integration_point.X();
    }
    else if (side == 2) // xi: n, eta: 1
    {
        integration_point.X() = edge_integration_point.X();
        integration_point.Y() = 1.0;
    }
    else if (side == 3) // xi: -1, eta: n
    {
        integration_point.X() = -1.0;
        integration_point.Y() = edge_integration_point.X();
    }

    return integration_point;
}

Element::GeometryType::IntegrationPointType GhostPenalty_Helper::ComputeIntegrationPointOnTetSide(const int& side, const Element::GeometryType::IntegrationPointType& face_integration_point)
{
    // TODO
    KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "not yet implemented")
}

Element::GeometryType::IntegrationPointType GhostPenalty_Helper::ComputeIntegrationPointOnHexSide(const int& side, const Element::GeometryType::IntegrationPointType& face_integration_point)
{
    IntegrationPointType integration_point;

    if (side == 0) // zeta: -1
    {
        integration_point.X() = face_integration_point.X();
        integration_point.Y() = face_integration_point.Y();
        integration_point.Z() = -1.0;
    }
    else if (side == 1) // eta: -1
    {
        integration_point.X() = face_integration_point.X();
        integration_point.Y() = -1.0;
        integration_point.Z() = face_integration_point.Y();
    }
    else if (side == 2) // xi: 1
    {
        integration_point.X() = 1.0;
        integration_point.Y() = face_integration_point.X();
        integration_point.Z() = face_integration_point.Y();
    }
    else if (side == 3) // eta: 1
    {
        integration_point.X() = face_integration_point.X();
        integration_point.Y() = 1.0;
        integration_point.Z() = face_integration_point.Y();
    }
    else if (side == 4) // xi: -1
    {
        integration_point.X() = -1.0;
        integration_point.Y() = face_integration_point.X();
        integration_point.Z() = face_integration_point.Y();
    }
    else if (side == 5) // zeta: 1
    {
        integration_point.X() = face_integration_point.X();
        integration_point.Y() = face_integration_point.Y();
        integration_point.Z() = 1.0;
    }

    return integration_point;
}

Vector& GhostPenalty_Helper::ComputeShapeFunction(Vector& N, GeometryType& r_element_geometry, const IntegrationPointType& integration_point)
{
    N = r_element_geometry.ShapeFunctionsValues(N, integration_point);
    return N;
}

Matrix& GhostPenalty_Helper::ComputeShapeFunctionGradient(Matrix& DN_DX, GeometryType& r_element_geometry, const IntegrationPointType& integration_point)
{
    Matrix DeltaPosition(r_element_geometry.size(), 3);
    for ( unsigned int node = 0; node < r_element_geometry.size(); ++node )
    {
        noalias( row( DeltaPosition, node ) ) = r_element_geometry[node].Coordinates() - r_element_geometry[node].GetInitialPosition();
    }

    Matrix DN_De;
    DN_De = r_element_geometry.ShapeFunctionsLocalGradients(DN_De, integration_point);

    Matrix InvJ;
    InvJ = r_element_geometry.InverseOfJacobian(InvJ, integration_point, DeltaPosition);

    noalias( DN_DX ) = prod( DN_De, InvJ );

    return DN_DX;
}

std::vector<Vector>& GhostPenalty_Helper::ComputeShapeFunctionSecondDerivatives(std::vector<Vector>& D2N_DX2, GeometryType& r_element_geometry, const IntegrationPointType& integration_point)
{
    unsigned int dim = r_element_geometry.WorkingSpaceDimension();
    unsigned int number_of_nodes = r_element_geometry.size();
    unsigned int strain_size = dim * (dim + 1) / 2;

    // compute the gradient and inversed of Jacobian
    Matrix DeltaPosition(r_element_geometry.size(), 3);
    for ( unsigned int node = 0; node < r_element_geometry.size(); ++node )
    {
        noalias( row( DeltaPosition, node ) ) = r_element_geometry[node].Coordinates() - r_element_geometry[node].GetInitialPosition();
    }

    Matrix DN_De;
    DN_De = r_element_geometry.ShapeFunctionsLocalGradients(DN_De, integration_point);

    Matrix J, InvJ;
    double DetJ;
    J = r_element_geometry.Jacobian(J, integration_point, DeltaPosition);
    MathUtils<double>::InvertMatrix(J, InvJ, DetJ);

    Matrix DN_DX(r_element_geometry.size(), dim);
    noalias( DN_DX ) = prod( DN_De, InvJ );

    // compute the second derivatives
    ShapeFunctionsSecondDerivativesType D2N_De2;
    D2N_De2 = r_element_geometry.ShapeFunctionsSecondDerivatives(D2N_De2, integration_point);
#ifdef ENABLE_DEBUG_GHOST_PENALTY_2
    KRATOS_WATCH(D2N_De2)
    KRATOS_WATCH(J)
#endif

    Matrix D2X_De2(dim, dim);
    Matrix D2Y_De2(dim, dim);
    Matrix D2Z_De2(dim, dim);
    for (unsigned int i = 0; i < dim; ++i)
    {
        for (unsigned int j = 0; j < dim; ++j)
        {
            D2X_De2(i, j) = 0.0;
            for (unsigned int k = 0; k < number_of_nodes; ++k)
            {
                D2X_De2(i, j) += r_element_geometry[k].X0() * D2N_De2[k](i, j);
            }

            D2Y_De2(i, j) = 0.0;
            for (unsigned int k = 0; k < number_of_nodes; ++k)
            {
                D2Y_De2(i, j) += r_element_geometry[k].Y0() * D2N_De2[k](i, j);
            }

            if (dim > 2)
            {
                D2Z_De2(i, j) = 0.0;
                for (unsigned int k = 0; k < number_of_nodes; ++k)
                {
                    D2Z_De2(i, j) += r_element_geometry[k].Z0() * D2N_De2[k](i, j);
                }
            }
        }
    }

    Matrix D2(strain_size, strain_size);
    if (dim == 2)
    {
        D2(0, 0) = pow(J(0, 0), 2);
        D2(0, 1) = pow(J(0, 1), 2);
        D2(0, 2) = J(0, 0) * J(0, 1);

        D2(1, 0) = pow(J(1, 0), 2);
        D2(1, 1) = pow(J(1, 1), 2);
        D2(1, 2) = J(1, 0) * J(1, 1);

        D2(2, 0) = 2.0 * J(0, 0) * J(1, 0);
        D2(2, 1) = 2.0 * J(0, 1) * J(1, 1);
        D2(2, 2) = J(0, 0) * J(1, 1) + J(1, 0) * J(0, 1);
    }
    else if (dim == 3)
    {
        //TODO there are still some bugs in computing D2, must be checked again
        D2(0, 0) = pow(J(0, 0), 2);
        D2(0, 1) = pow(J(0, 1), 2);
        D2(0, 2) = pow(J(0, 2), 2);
        D2(0, 3) = J(0, 0) * J(0, 1);
        D2(0, 4) = J(0, 1) * J(0, 2);
        D2(0, 5) = J(0, 0) * J(0, 2);

        D2(1, 0) = pow(J(1, 0), 2);
        D2(1, 1) = pow(J(1, 1), 2);
        D2(1, 2) = pow(J(1, 2), 2);
        D2(1, 3) = J(1, 0) * J(1, 1);
        D2(1, 4) = J(1, 1) * J(1, 2);
        D2(1, 5) = J(1, 0) * J(1, 2);

        D2(2, 0) = pow(J(2, 0), 2);
        D2(2, 1) = pow(J(2, 1), 2);
        D2(2, 2) = pow(J(2, 2), 2);
        D2(2, 3) = J(2, 0) * J(2, 1);
        D2(2, 4) = J(2, 1) * J(2, 2);
        D2(2, 5) = J(2, 0) * J(2, 2);

        D2(3, 0) = 2.0 * J(0, 0) * J(1, 0);
        D2(3, 1) = 2.0 * J(0, 1) * J(1, 1);
        D2(3, 2) = 2.0 * J(0, 2) * J(1, 2);
        D2(3, 3) = J(0, 0) * J(1, 1) + J(0, 1) * J(1, 0);
        D2(3, 4) = J(0, 1) * J(1, 2) + J(0, 2) * J(1, 1);
        D2(3, 5) = J(0, 0) * J(1, 2) + J(0, 2) * J(1, 0);

        D2(4, 0) = 2.0 * J(1, 0) * J(2, 0);
        D2(4, 1) = 2.0 * J(1, 1) * J(2, 1);
        D2(4, 2) = 2.0 * J(1, 2) * J(2, 2);
        D2(4, 3) = J(1, 0) * J(2, 1) + J(1, 1) * J(2, 0);
        D2(4, 4) = J(1, 1) * J(2, 2) + J(1, 2) * J(2, 1);
        D2(4, 5) = J(1, 0) * J(2, 2) + J(1, 2) * J(2, 0);

        D2(5, 0) = 2.0 * J(0, 0) * J(2, 0);
        D2(5, 1) = 2.0 * J(0, 1) * J(2, 1);
        D2(5, 2) = 2.0 * J(0, 2) * J(2, 2);
        D2(5, 3) = J(0, 0) * J(2, 1) + J(0, 1) * J(2, 0);
        D2(5, 4) = J(0, 1) * J(2, 2) + J(0, 2) * J(2, 1);
        D2(5, 5) = J(0, 0) * J(2, 2) + J(0, 2) * J(2, 0);
    }

    Matrix InvD2(strain_size, strain_size);
    SD_MathUtils<double>::InvertMatrix(D2, InvD2);

    if (D2N_DX2.size() != number_of_nodes)
    {
        D2N_DX2.resize(number_of_nodes);
    }

    Vector R2(strain_size);
#ifdef ENABLE_DEBUG_GHOST_PENALTY_2
    KRATOS_WATCH(D2)
    KRATOS_WATCH(InvD2)
    KRATOS_WATCH(D2X_De2)
    KRATOS_WATCH(D2Y_De2)
    KRATOS_WATCH(D2Z_De2)
#endif
    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        if (D2N_DX2[i].size() != strain_size)
        {
            D2N_DX2[i].resize(strain_size, false);
        }

        if (dim == 2)
        {
            R2(0) = D2N_De2[i](0, 0) - DN_DX(i, 0) * D2X_De2(0, 0) - DN_DX(i, 1) * D2Y_De2(0, 0);
            R2(1) = D2N_De2[i](1, 1) - DN_DX(i, 0) * D2X_De2(1, 1) - DN_DX(i, 1) * D2Y_De2(1, 1);
            R2(2) = D2N_De2[i](0, 1) - DN_DX(i, 0) * D2X_De2(0, 1) - DN_DX(i, 1) * D2Y_De2(0, 1);
        }
        else if (dim == 3)
        {
            R2(0) = D2N_De2[i](0, 0) - DN_DX(i, 0) * D2X_De2(0, 0) - DN_DX(i, 1) * D2Y_De2(0, 0) - DN_DX(i, 2) * D2Z_De2(0, 0);
            R2(1) = D2N_De2[i](1, 1) - DN_DX(i, 0) * D2X_De2(1, 1) - DN_DX(i, 1) * D2Y_De2(1, 1) - DN_DX(i, 2) * D2Z_De2(1, 1);
            R2(2) = D2N_De2[i](2, 2) - DN_DX(i, 0) * D2X_De2(2, 2) - DN_DX(i, 1) * D2Y_De2(2, 2) - DN_DX(i, 2) * D2Z_De2(2, 2);
            R2(3) = D2N_De2[i](0, 1) - DN_DX(i, 0) * D2X_De2(0, 1) - DN_DX(i, 1) * D2Y_De2(0, 1) - DN_DX(i, 2) * D2Z_De2(0, 1);
            R2(4) = D2N_De2[i](1, 2) - DN_DX(i, 0) * D2X_De2(1, 2) - DN_DX(i, 1) * D2Y_De2(1, 2) - DN_DX(i, 2) * D2Z_De2(1, 2);
            R2(5) = D2N_De2[i](0, 2) - DN_DX(i, 0) * D2X_De2(0, 2) - DN_DX(i, 1) * D2Y_De2(0, 2) - DN_DX(i, 2) * D2Z_De2(0, 2);
        }

        noalias(D2N_DX2[i]) = prod(R2, InvD2);
    }

    return D2N_DX2;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

const int GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::msEdges[][2] =
{
    {0, 1},
    {1, 2},
    {2, 0}
};

std::vector<std::size_t> GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::FindCommonFace(Element::GeometryType& r_geom_1, Element::GeometryType& r_geom_2)
{
    for (std::size_t i = 0; i < NumberOfSides(); ++i)
    {
        const std::size_t& i1 = Faces(i)[0];
        const std::size_t& i2 = Faces(i)[1];
        const std::size_t& n1 = r_geom_1[i1].Id();
        const std::size_t& n2 = r_geom_1[i2].Id();

        if (GhostPenalty_Helper::IsBelongedToGeometry(n1, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n2, r_geom_2))
            return std::vector<std::size_t> {i1, i2};
    }
    return std::vector<std::size_t> {};
}

Element::GeometryType::IntegrationPointType GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::ComputeIntegrationPoint(const int& side, const Element::GeometryType::IntegrationPointType& edge_integration_point)
{
    return GhostPenalty_Helper::ComputeIntegrationPointOnTriSide(side, edge_integration_point);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunction_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3> >::Execute(N, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunctionGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3> >::Execute(DN_DX, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::ComputeShapeFunctionNormalGradient(Matrix& dNdn, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_edge_geometry, const Element::GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunctionNormalGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3> >::Execute2D(dNdn, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D3>::ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_edge_geometry, const Element::GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "for T3 not yet implemented")
}

///////////////////////////////////////////////////////////////////////////////

const int GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::msEdges[][3] =
{
    {0, 1, 3},
    {1, 2, 4},
    {2, 0, 5}
};

std::vector<std::size_t> GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::FindCommonFace(Element::GeometryType& r_geom_1, Element::GeometryType& r_geom_2)
{
    for (std::size_t i = 0; i < NumberOfSides(); ++i)
    {
        const std::size_t& i1 = Faces(i)[0];
        const std::size_t& i2 = Faces(i)[1];
        const std::size_t& i3 = Faces(i)[2];
        const std::size_t& n1 = r_geom_1[i1].Id();
        const std::size_t& n2 = r_geom_1[i2].Id();
        const std::size_t& n3 = r_geom_1[i3].Id();

        if (GhostPenalty_Helper::IsBelongedToGeometry(n1, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n2, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n3, r_geom_2))
            return std::vector<std::size_t> {i1, i2, i3};
    }
    return std::vector<std::size_t> {};
}

Element::GeometryType::IntegrationPointType GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::ComputeIntegrationPoint(const int& side, const Element::GeometryType::IntegrationPointType& edge_integration_point)
{
    return GhostPenalty_Helper::ComputeIntegrationPointOnTriSide(side, edge_integration_point);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunction_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6> >::Execute(N, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunctionGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6> >::Execute(DN_DX, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::ComputeShapeFunctionNormalGradient(Matrix& dNdn, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_edge_geometry, const Element::GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunctionNormalGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6> >::Execute2D(dNdn, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Triangle2D6>::ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_edge_geometry, const Element::GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "for T6 not yet implemented")
}

///////////////////////////////////////////////////////////////////////////////

const int GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::msEdges[][2] =
{
    {0, 1},
    {1, 2},
    {2, 3},
    {3, 0}
};

std::vector<std::size_t> GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::FindCommonFace(Element::GeometryType& r_geom_1, Element::GeometryType& r_geom_2)
{
    for (std::size_t i = 0; i < NumberOfSides(); ++i)
    {
        const std::size_t& i1 = Faces(i)[0];
        const std::size_t& i2 = Faces(i)[1];
        const std::size_t& n1 = r_geom_1[i1].Id();
        const std::size_t& n2 = r_geom_2[i2].Id();

        if (   GhostPenalty_Helper::IsBelongedToGeometry(n1, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n2, r_geom_2))
            return std::vector<std::size_t> {i1, i2};
    }
    return std::vector<std::size_t> {};
}

Element::GeometryType::IntegrationPointType GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::ComputeIntegrationPoint(const int& side, const Element::GeometryType::IntegrationPointType& edge_integration_point)
{
    return GhostPenalty_Helper::ComputeIntegrationPointOnQuadSide(side, edge_integration_point);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunction_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4> >::Execute(N, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunctionGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4> >::Execute(DN_DX, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::ComputeShapeFunctionNormalGradient(Matrix& dNdn, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_edge_geometry, const Element::GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunctionNormalGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4> >::Execute2D(dNdn, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4>::ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_edge_geometry, const Element::GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "for Q4 not yet implemented")
}

///////////////////////////////////////////////////////////////////////////////

const int GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::msEdges[][3] =
{
    {0, 1, 4},
    {1, 2, 5},
    {2, 3, 6},
    {3, 0, 7}
};

std::vector<std::size_t> GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::FindCommonFace(Element::GeometryType& r_geom_1, Element::GeometryType& r_geom_2)
{
    for (std::size_t i = 0; i < NumberOfSides(); ++i)
    {
        const std::size_t& i1 = Faces(i)[0];
        const std::size_t& i2 = Faces(i)[1];
        const std::size_t& i3 = Faces(i)[2];
        const std::size_t& n1 = r_geom_1[i1].Id();
        const std::size_t& n2 = r_geom_1[i2].Id();
        const std::size_t& n3 = r_geom_1[i3].Id();

        if (   GhostPenalty_Helper::IsBelongedToGeometry(n1, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n2, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n3, r_geom_2))
            return std::vector<std::size_t> {i1, i2, i3};
    }
    return std::vector<std::size_t> {};
}

Element::GeometryType::IntegrationPointType GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::ComputeIntegrationPoint(const int& side, const Element::GeometryType::IntegrationPointType& edge_integration_point)
{
    return GhostPenalty_Helper::ComputeIntegrationPointOnQuadSide(side, edge_integration_point);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunction_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8> >::Execute(N, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunctionGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8> >::Execute(DN_DX, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::ComputeShapeFunctionNormalGradient(Matrix& dNdn, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_edge_geometry, const Element::GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunctionNormalGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8> >::Execute2D(dNdn, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8>::ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_edge_geometry, const Element::GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "for Q8 not yet implemented")
}

///////////////////////////////////////////////////////////////////////////////

const int GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::msEdges[][3] =
{
    {0, 1, 4},
    {1, 2, 5},
    {2, 3, 6},
    {3, 0, 7}
};

std::vector<std::size_t> GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::FindCommonFace(Element::GeometryType& r_geom_1, Element::GeometryType& r_geom_2)
{
    for (std::size_t i = 0; i < NumberOfSides(); ++i)
    {
        const std::size_t& i1 = Faces(i)[0];
        const std::size_t& i2 = Faces(i)[1];
        const std::size_t& i3 = Faces(i)[2];
        const std::size_t& n1 = r_geom_1[i1].Id();
        const std::size_t& n2 = r_geom_1[i2].Id();
        const std::size_t& n3 = r_geom_1[i3].Id();

        if (   GhostPenalty_Helper::IsBelongedToGeometry(n1, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n2, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n3, r_geom_2))
            return std::vector<std::size_t> {i1, i2, i3};
    }
    return std::vector<std::size_t> {};
}

Element::GeometryType::IntegrationPointType GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::ComputeIntegrationPoint(const int& side, const Element::GeometryType::IntegrationPointType& edge_integration_point)
{
    return GhostPenalty_Helper::ComputeIntegrationPointOnQuadSide(side, edge_integration_point);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunction_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9> >::Execute(N, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunctionGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9> >::Execute(DN_DX, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::ComputeShapeFunctionNormalGradient(Matrix& dNdn, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_edge_geometry, const Element::GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    ComputeShapeFunctionNormalGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9> >::Execute2D(dNdn, r_element_geometry, r_edge_geometry, edge_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9>::ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_edge_geometry, const Element::GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "for Q9 not yet implemented")
}

///////////////////////////////////////////////////////////////////////////////

const int GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::msFaces[][3] =
{
    {0, 2, 1},
    {0, 3, 2},
    {0, 1, 3},
    {2, 3, 1}
};

std::vector<std::size_t> GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::FindCommonFace(Element::GeometryType& r_geom_1, Element::GeometryType& r_geom_2)
{
    for (std::size_t i = 0; i < NumberOfSides(); ++i)
    {
        const std::size_t& i1 = Faces(i)[0];
        const std::size_t& i2 = Faces(i)[1];
        const std::size_t& i3 = Faces(i)[2];
        const std::size_t& n1 = r_geom_1[i1].Id();
        const std::size_t& n2 = r_geom_1[i2].Id();
        const std::size_t& n3 = r_geom_1[i3].Id();

        if (   GhostPenalty_Helper::IsBelongedToGeometry(n1, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n2, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n3, r_geom_2))
            return std::vector<std::size_t> {i1, i2, i3};
    }
    return std::vector<std::size_t> {};
}

Element::GeometryType::IntegrationPointType GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::ComputeIntegrationPoint(const int& side, const Element::GeometryType::IntegrationPointType& edge_integration_point)
{
    return GhostPenalty_Helper::ComputeIntegrationPointOnTetSide(side, edge_integration_point);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunction_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4> >::Execute(N, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunctionGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4> >::Execute(DN_DX, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::ComputeShapeFunctionNormalGradient(Matrix& dNdn, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_face_geometry, const Element::GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunctionNormalGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4> >::Execute3D(dNdn, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4>::ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_face_geometry, const Element::GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "for T4 not yet implemented")
}

///////////////////////////////////////////////////////////////////////////////

const int GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::msFaces[][6] =
{
    {0, 2, 1, 6, 5, 4},
    {0, 3, 2, 7, 9, 6},
    {0, 1, 3, 4, 8, 7},
    {2, 3, 1, 9, 8, 5}
};

std::vector<std::size_t> GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::FindCommonFace(Element::GeometryType& r_geom_1, Element::GeometryType& r_geom_2)
{
    for (std::size_t i = 0; i < NumberOfSides(); ++i)
    {
        const std::size_t& i1 = Faces(i)[0];
        const std::size_t& i2 = Faces(i)[1];
        const std::size_t& i3 = Faces(i)[2];
        const std::size_t& i4 = Faces(i)[3];
        const std::size_t& i5 = Faces(i)[4];
        const std::size_t& i6 = Faces(i)[5];
        const std::size_t& n1 = r_geom_1[i1].Id();
        const std::size_t& n2 = r_geom_1[i2].Id();
        const std::size_t& n3 = r_geom_1[i3].Id();
        const std::size_t& n4 = r_geom_1[i4].Id();
        const std::size_t& n5 = r_geom_1[i5].Id();
        const std::size_t& n6 = r_geom_1[i6].Id();

        if (   GhostPenalty_Helper::IsBelongedToGeometry(n1, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n2, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n3, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n4, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n5, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n6, r_geom_2))
            return std::vector<std::size_t> {i1, i2, i3, i4, i5, i6};
    }
    return std::vector<std::size_t> {};
}

Element::GeometryType::IntegrationPointType GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::ComputeIntegrationPoint(const int& side, const Element::GeometryType::IntegrationPointType& edge_integration_point)
{
    return GhostPenalty_Helper::ComputeIntegrationPointOnTetSide(side, edge_integration_point);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunction_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10> >::Execute(N, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunctionGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10> >::Execute(DN_DX, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::ComputeShapeFunctionNormalGradient(Matrix& dNdn, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_face_geometry, const Element::GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunctionNormalGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10> >::Execute3D(dNdn, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10>::ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_face_geometry, const Element::GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "for T10 not yet implemented")
}

///////////////////////////////////////////////////////////////////////////////

const int GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::msFaces[][4] =
{
    {3, 2, 1, 0},
    {0, 1, 5, 4},
    {2, 6, 5, 1},
    {7, 6, 2, 3},
    {7, 3, 0, 4},
    {4, 5, 6, 7}
};

std::vector<std::size_t> GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::FindCommonFace(Element::GeometryType& r_geom_1, Element::GeometryType& r_geom_2)
{
    for (std::size_t i = 0; i < NumberOfSides(); ++i)
    {
        const std::size_t& i1 = Faces(i)[0];
        const std::size_t& i2 = Faces(i)[1];
        const std::size_t& i3 = Faces(i)[2];
        const std::size_t& i4 = Faces(i)[3];
        const std::size_t& n1 = r_geom_1[i1].Id();
        const std::size_t& n2 = r_geom_1[i2].Id();
        const std::size_t& n3 = r_geom_1[i3].Id();
        const std::size_t& n4 = r_geom_1[i4].Id();

        if (   GhostPenalty_Helper::IsBelongedToGeometry(n1, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n2, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n3, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n4, r_geom_2))
            return std::vector<std::size_t> {i1, i2, i3, i4};
    }
    return std::vector<std::size_t> {};
}

Element::GeometryType::IntegrationPointType GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::ComputeIntegrationPoint(const int& side, const Element::GeometryType::IntegrationPointType& edge_integration_point)
{
    return GhostPenalty_Helper::ComputeIntegrationPointOnHexSide(side, edge_integration_point);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunction_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8> >::Execute(N, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunctionGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8> >::Execute(DN_DX, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::ComputeShapeFunctionNormalGradient(Matrix& dNdn, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_face_geometry, const Element::GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunctionNormalGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8> >::Execute3D(dNdn, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D8>::ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_face_geometry, const Element::GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "for H8 not yet implemented")
}

///////////////////////////////////////////////////////////////////////////////

const int GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::msFaces[][8] =
{
    {3, 2, 1, 0, 10, 9, 8, 11},
    {0, 1, 5, 4, 8, 13, 16, 12},
    {2, 6, 5, 1, 14, 17, 13, 9},
    {7, 6, 2, 3, 18, 14, 10, 15},
    {7, 3, 0, 4, 15, 11, 12, 19},
    {4, 5, 6, 7, 16, 17, 18, 19}
};

std::vector<std::size_t> GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::FindCommonFace(Element::GeometryType& r_geom_1, Element::GeometryType& r_geom_2)
{
    for (std::size_t i = 0; i < NumberOfSides(); ++i)
    {
        const std::size_t& i1 = Faces(i)[0];
        const std::size_t& i2 = Faces(i)[1];
        const std::size_t& i3 = Faces(i)[2];
        const std::size_t& i4 = Faces(i)[3];
        const std::size_t& i5 = Faces(i)[4];
        const std::size_t& i6 = Faces(i)[5];
        const std::size_t& i7 = Faces(i)[6];
        const std::size_t& i8 = Faces(i)[7];
        const std::size_t& n1 = r_geom_1[i1].Id();
        const std::size_t& n2 = r_geom_1[i2].Id();
        const std::size_t& n3 = r_geom_1[i3].Id();
        const std::size_t& n4 = r_geom_1[i4].Id();
        const std::size_t& n5 = r_geom_1[i5].Id();
        const std::size_t& n6 = r_geom_1[i6].Id();
        const std::size_t& n7 = r_geom_1[i7].Id();
        const std::size_t& n8 = r_geom_1[i8].Id();

        if (   GhostPenalty_Helper::IsBelongedToGeometry(n1, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n2, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n3, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n4, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n5, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n6, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n7, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n8, r_geom_2))
            return std::vector<std::size_t> {i1, i2, i3, i4, i5, i6, i7, i8};
    }
    return std::vector<std::size_t> {};
}

Element::GeometryType::IntegrationPointType GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::ComputeIntegrationPoint(const int& side, const Element::GeometryType::IntegrationPointType& edge_integration_point)
{
    return GhostPenalty_Helper::ComputeIntegrationPointOnHexSide(side, edge_integration_point);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunction_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20> >::Execute(N, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunctionGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20> >::Execute(DN_DX, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::ComputeShapeFunctionNormalGradient(Matrix& dNdn, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_face_geometry, const Element::GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunctionNormalGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20> >::Execute3D(dNdn, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D20>::ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_face_geometry, const Element::GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "for H20 not yet implemented")
}

///////////////////////////////////////////////////////////////////////////////

const int GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::msFaces[][9] =
{
    {3, 2, 1, 0, 10, 9, 8, 11, 20},
    {0, 1, 5, 4, 8, 13, 16, 12, 21},
    {2, 6, 5, 1, 14, 17, 13, 9, 22},
    {7, 6, 2, 3, 18, 14, 10, 15, 23},
    {7, 3, 0, 4, 15, 11, 12, 19, 24},
    {4, 5, 6, 7, 16, 17, 18, 19, 25}
};

std::vector<std::size_t> GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::FindCommonFace(Element::GeometryType& r_geom_1, Element::GeometryType& r_geom_2)
{
    for (std::size_t i = 0; i < NumberOfSides(); ++i)
    {
        const std::size_t& i1 = Faces(i)[0];
        const std::size_t& i2 = Faces(i)[1];
        const std::size_t& i3 = Faces(i)[2];
        const std::size_t& i4 = Faces(i)[3];
        const std::size_t& i5 = Faces(i)[4];
        const std::size_t& i6 = Faces(i)[5];
        const std::size_t& i7 = Faces(i)[6];
        const std::size_t& i8 = Faces(i)[7];
        const std::size_t& i9 = Faces(i)[8];
        const std::size_t& n1 = r_geom_1[i1].Id();
        const std::size_t& n2 = r_geom_1[i2].Id();
        const std::size_t& n3 = r_geom_1[i3].Id();
        const std::size_t& n4 = r_geom_1[i4].Id();
        const std::size_t& n5 = r_geom_1[i5].Id();
        const std::size_t& n6 = r_geom_1[i6].Id();
        const std::size_t& n7 = r_geom_1[i7].Id();
        const std::size_t& n8 = r_geom_1[i8].Id();
        const std::size_t& n9 = r_geom_1[i9].Id();

        if (   GhostPenalty_Helper::IsBelongedToGeometry(n1, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n2, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n3, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n4, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n5, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n6, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n7, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n8, r_geom_2)
                && GhostPenalty_Helper::IsBelongedToGeometry(n9, r_geom_2))
            return std::vector<std::size_t> {i1, i2, i3, i4, i5, i6, i7, i8, i9};
    }
    return std::vector<std::size_t> {};
}

Element::GeometryType::IntegrationPointType GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::ComputeIntegrationPoint(const int& side, const Element::GeometryType::IntegrationPointType& edge_integration_point)
{
    return GhostPenalty_Helper::ComputeIntegrationPointOnHexSide(side, edge_integration_point);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::ComputeShapeFunction(Matrix& N, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_face_geometry, const Element::GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunction_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27> >::Execute(N, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunctionGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27> >::Execute(DN_DX, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::ComputeShapeFunctionNormalGradient(Matrix& dNdn, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_face_geometry, const Element::GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    ComputeShapeFunctionNormalGradient_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27> >::Execute3D(dNdn, r_element_geometry, r_face_geometry, face_integration_points);
}

void GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27>::ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, Element::GeometryType& r_element_geometry,
        Element::GeometryType& r_face_geometry, const Element::GeometryType::IntegrationPointsArrayType& face_integration_points)
{
    std::map<std::size_t, std::size_t> map_edge_node_index_to_element_node_index;

#ifdef ENABLE_DEBUG_GHOST_PENALTY_2
    std::cout << "element geometry:";
    for (std::size_t i = 0; i < r_element_geometry.size(); ++i)
    {
        std::cout << " " << r_element_geometry[i].Id();
    }
    std::cout << std::endl;

    std::cout << "face geometry:";
    for (std::size_t i = 0; i < r_face_geometry.size(); ++i)
    {
        std::cout << " " << r_face_geometry[i].Id();
    }
    std::cout << std::endl;
#endif

    int side = FindSide_Helper<GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType::Kratos_Hexahedra3D27> >::Execute(r_element_geometry, r_face_geometry);
    if (side != -1)
    {
        GhostPenalty_Helper::BuildMapEdgeNodeIndexToElementNodeIndex(map_edge_node_index_to_element_node_index,
                r_face_geometry, r_element_geometry, Faces(side));

#ifdef ENABLE_DEBUG_GHOST_PENALTY_2
        KRATOS_WATCH(side)
#endif

        std::vector<Vector> D2N_DX2;

        if (d2Ndn2.size1() != r_face_geometry.size() || d2Ndn2.size2() != face_integration_points.size())
        {
            d2Ndn2.resize(r_face_geometry.size(), face_integration_points.size(), false);
        }

        // construct the local coordinates system associated with the face
        Vector t1(3), t2(3), n(3);

        for (std::size_t i = 0; i < face_integration_points.size(); ++i)
        {
            GhostPenalty_Helper::ComputeLocalFrame(t1, t2, n, r_face_geometry, face_integration_points[i]);
#ifdef ENABLE_DEBUG_GHOST_PENALTY_2
            KRATOS_WATCH(t1)
            KRATOS_WATCH(t2)
            KRATOS_WATCH(n)
#endif

            IntegrationPointType integration_point = ComputeIntegrationPoint(side, face_integration_points[i]);

            D2N_DX2 = GhostPenalty_Helper::ComputeShapeFunctionSecondDerivatives(D2N_DX2, r_element_geometry, integration_point);
#ifdef ENABLE_DEBUG_GHOST_PENALTY_2
            std::cout << "D2N_DX2:" << std::endl;
            for (std::size_t i = 0; i < D2N_DX2.size(); ++i)
            {
                std::cout << " " << D2N_DX2[i] << std::endl;
            }
            std::cout << std::endl;
#endif
            for (std::size_t j = 0; j < r_face_geometry.size(); ++j)
            {
                d2Ndn2(j, i) = pow(n(0), 2) * D2N_DX2[j](0) + pow(n(1), 2) * D2N_DX2[j](1) + pow(n(2), 2) * D2N_DX2[j](2)
                               + 2.0 * n(0) * n(1) * D2N_DX2[j](3) + 2.0 * n(1) * n(2) * D2N_DX2[j](4) + 2.0 * n(0) * n(2) * D2N_DX2[j](5);
            }
        }
    }

    if (side == -1)
        KRATOS_THROW_ERROR(std::logic_error, "The face geometry is not found on the element geometry. Check the input", "")
    }

///////////////////////////////////////////////////////////////////////////////

}  // namespace Kratos.

#undef ENABLE_DEBUG_GHOST_PENALTY
#undef ENABLE_DEBUG_GHOST_PENALTY_2

