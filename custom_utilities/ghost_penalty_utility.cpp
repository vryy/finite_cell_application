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


namespace Kratos
{

const int GhostPenaltyUtility::msEdgesT3[][2] = { {0, 1}, {1, 2}, {2, 0} };
const int GhostPenaltyUtility::msEdgesT6[][3] = { {0, 1, 3}, {1, 2, 4}, {2, 0, 5} };
const int GhostPenaltyUtility::msEdgesQ4[][2] = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
const int GhostPenaltyUtility::msEdgesQ8[][3] = { {0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7} };
const int GhostPenaltyUtility::msEdgesQ9[][3] = { {0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7} };
const int GhostPenaltyUtility::msFacesT4[][3] = { {0, 2, 1}, {0, 3, 2}, {0, 1, 3}, {2, 3, 1} };
const int GhostPenaltyUtility::msFacesT10[][6] = { {0, 2, 1, 6, 5, 4}, {0, 3, 2, 7, 9, 6}, {0, 1, 3, 4, 8, 7}, {2, 3, 1, 9, 8, 5} };
const int GhostPenaltyUtility::msFacesH8[][4] = { {3, 2, 1, 0}, {0, 1, 5, 4}, {2, 6, 5, 1}, {7, 6, 2, 3}, {7, 3, 0, 4}, {4, 5, 6, 7} };
const int GhostPenaltyUtility::msFacesH20[][8] = { {3, 2, 1, 0, 10, 9, 8, 11}, {0, 1, 5, 4, 8, 13, 16, 12}, {2, 6, 5, 1, 14, 17, 13, 9}, {7, 6, 2, 3, 14, 18, 10, 15}, {7, 3, 0, 4, 15, 11, 12, 19}, {4, 5, 6, 7, 16, 17, 18, 19} };
const int GhostPenaltyUtility::msFacesH27[][9] = { {3, 2, 1, 0, 10, 9, 8, 11, 20}, {0, 1, 5, 4, 8, 13, 16, 12, 21}, {2, 6, 5, 1, 14, 17, 13, 9, 22}, {7, 6, 2, 3, 14, 18, 10, 15, 23}, {7, 3, 0, 4, 15, 11, 12, 19, 24}, {4, 5, 6, 7, 16, 17, 18, 19, 25} };

ModelPart::ConditionsContainerType GhostPenaltyUtility::SetUpSurfacePenaltyConditions(ModelPart& r_model_part,
    GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
    std::size_t& lastCondId, Properties::Pointer pProperties)
{
    // setup all the element neighbour couples in the model
    std::set<std::pair<std::size_t, std::size_t> > edges;
    for (ModelPart::ElementsContainerType::iterator it = r_model_part.Elements().begin();
        it != r_model_part.Elements().end(); ++it)
    {
        WeakPointerVector<Element>& rNeighbours = it->GetValue(NEIGHBOUR_ELEMENTS);

        for (std::size_t i = 0; i < rNeighbours.size(); ++i)
        {
            const std::size_t& id1 = it->Id();
            const std::size_t& id2 = rNeighbours(i).lock()->Id();
            if (id1 < id2)
                edges.insert(std::make_pair(id1, id2));
            else if (id2 > id1)
                edges.insert(std::make_pair(id2, id1));
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
                p_sample_condition, r_brep, lastCondId, pProperties);
        if (pNewCond != NULL)
            pNewConditions.push_back(pNewCond);
    }

    for (ModelPart::ConditionsContainerType::ptr_iterator it = pNewConditions.ptr_begin();
        it != pNewConditions.ptr_end(); ++it)
    {
        r_model_part.Conditions().push_back(*it);
    }

    std::cout << __FUNCTION__ << " completed: " << pNewConditions.size() << " new ghost conditions"
              << " of type " << typeid(*p_sample_condition).name() << " are created and added to model_part" << std::endl;

    return pNewConditions;
}

ModelPart::ConditionsContainerType GhostPenaltyUtility::SetUpSurfacePenaltyConditions(ModelPart& r_model_part,
        ModelPart::ElementsContainerType& pElements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t& lastCondId, Properties::Pointer pProperties)
{
    // setup all the element neighbour couples in the list of elements
    std::set<std::pair<std::size_t, std::size_t> > edges;
    for (ModelPart::ElementsContainerType::iterator it = pElements.begin();
        it != pElements.end(); ++it)
    {
        WeakPointerVector<Element>& rNeighbours = it->GetValue(NEIGHBOUR_ELEMENTS);

        for (std::size_t i = 0; i < rNeighbours.size(); ++i)
        {
            const std::size_t& id2 = rNeighbours(i).lock()->Id();
            if (pElements.find(id2) == pElements.end())
                continue;

            const std::size_t& id1 = it->Id();
            if (id1 < id2)
                edges.insert(std::make_pair(id1, id2));
            else if (id2 > id1)
                edges.insert(std::make_pair(id2, id1));
        }
    }

    // setup the ghost penalty condition from the edges
    ModelPart::ConditionsContainerType pNewConditions;

    for (std::set<std::pair<std::size_t, std::size_t> >::iterator it = edges.begin();
        it != edges.end(); ++it)
    {
        Element::Pointer p_element_1 = pElements(it->first);
        Element::Pointer p_element_2 = pElements(it->second);

        Condition::Pointer pNewCond = SetUpSurfacePenaltyCondition(p_element_1, p_element_2,
                p_sample_condition, r_brep, lastCondId, pProperties);
        if (pNewCond != NULL)
            pNewConditions.push_back(pNewCond);
    }

    for (ModelPart::ConditionsContainerType::ptr_iterator it = pNewConditions.ptr_begin();
        it != pNewConditions.ptr_end(); ++it)
    {
        r_model_part.Conditions().push_back(*it);
    }

    std::cout << __FUNCTION__ << " on " << pElements.size() << " elements completed: " << pNewConditions.size() << " new ghost conditions"
              << " of type " << typeid(*p_sample_condition).name() << " are created and added to model_part" << std::endl;

    return pNewConditions;
}

ModelPart::ConditionsContainerType GhostPenaltyUtility::SetUpSurfacePenaltyConditions(Element::Pointer p_element,
        GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t& lastCondId, Properties::Pointer pProperties)
{
    // firstly obtain all neighbour elements of the current element
    WeakPointerVector<Element>& rNeighbours = p_element->GetValue(NEIGHBOUR_ELEMENTS);

    ModelPart::ConditionsContainerType pNewConditions;

    // for each neighbour, find the common edge
    for (std::size_t i = 0; i < rNeighbours.size(); ++i)
    {
        if (rNeighbours[i].Id() != p_element->Id())
        {
            Condition::Pointer pNewCond = SetUpSurfacePenaltyCondition(p_element, rNeighbours(i).lock(),
                p_sample_condition, r_brep, lastCondId, pProperties);
            if (pNewCond != NULL)
                pNewConditions.push_back(pNewCond);
        }
    }

    std::cout << __FUNCTION__ << " completed: " << pNewConditions.size() << " new ghost conditions"
              << " of type " << typeid(*p_sample_condition).name() << " are created" << std::endl;

    return pNewConditions;
}


Condition::Pointer GhostPenaltyUtility::SetUpSurfacePenaltyCondition(Element::Pointer p_element_1,
        Element::Pointer p_element_2, GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t& lastCondId, Properties::Pointer pProperties)
{
    Condition::Pointer pNewCond = NULL;

    std::pair<GeometryData::KratosGeometryType, std::vector<std::size_t> > edge = FindCommonEdge(*p_element_1, *p_element_2);

    if (edge.first == GeometryData::Kratos_generic_type)
        return pNewCond;

    // create the edge geometry
    typename Element::NodesArrayType temp_nodes;
    for (std::size_t j = 0; j < edge.second.size(); ++j)
        temp_nodes.push_back(p_element_1->GetGeometry().pGetPoint(edge.second[j]));

    GeometryType::Pointer p_temp_geometry;

    if (edge.first == GeometryData::Kratos_Line2D2)
        p_temp_geometry = GeometryType::Pointer(new Line2D2<NodeType>(temp_nodes));
    else if (edge.first == GeometryData::Kratos_Line2D3)
        p_temp_geometry = GeometryType::Pointer(new Line2D3<NodeType>(temp_nodes));
    else if (edge.first == GeometryData::Kratos_Triangle3D3)
        p_temp_geometry = GeometryType::Pointer(new Triangle3D3<NodeType>(temp_nodes));
    else if (edge.first == GeometryData::Kratos_Triangle3D6)
        p_temp_geometry = GeometryType::Pointer(new Triangle3D6<NodeType>(temp_nodes));
    else if (edge.first == GeometryData::Kratos_Quadrilateral3D4)
        p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D4<NodeType>(temp_nodes));
    else if (edge.first == GeometryData::Kratos_Quadrilateral3D8)
        p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D8<NodeType>(temp_nodes));
    else if (edge.first == GeometryData::Kratos_Quadrilateral3D9)
        p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D9<NodeType>(temp_nodes));
    else
        KRATOS_THROW_ERROR(std::logic_error, "Unknown geometry type", edge.first)

    // check if this edge is cut by the brep or totally inside. If yes, then the new ghost condition is created.
    bool is_ghost = false;

    int stat1 = r_brep.CutStatus(p_element_1->GetGeometry());
    int stat2 = r_brep.CutStatus(p_element_2->GetGeometry());

    if (stat1 == BRep::_CUT)
    {
        if (stat2 == BRep::_CUT || stat2 == BRep::_IN)
            is_ghost = true;
    }
    else if (stat1 == BRep::_IN)
    {
        if (stat2 == BRep::_CUT)
            is_ghost = true;
    }

    if (is_ghost)
    {
        // create the ghost penalty condition
        pNewCond = p_sample_condition->Create(++lastCondId, p_temp_geometry, p_element_1, p_element_2, pProperties);
        pNewCond->SetValue(IS_INACTIVE, false);
        pNewCond->Set(ACTIVE, true);
    }

    return pNewCond;
}


void GhostPenaltyUtility::ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{
    if (r_element_geometry.GetGeometryType() == GeometryData::Kratos_Triangle2D3)
    {
        ComputeShapeFunctionNormalGradientT3(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::Kratos_Triangle2D6)
    {
        ComputeShapeFunctionNormalGradientT6(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4)
    {
        ComputeShapeFunctionNormalGradientQ4(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8
          || r_element_geometry.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9)
    {
        ComputeShapeFunctionNormalGradientQ8(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::Kratos_Tetrahedra3D4)
    {
        ComputeShapeFunctionNormalGradientT4(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::Kratos_Tetrahedra3D10)
    {
        ComputeShapeFunctionNormalGradientT10(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::Kratos_Hexahedra3D8)
    {
        ComputeShapeFunctionNormalGradientH8(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::Kratos_Hexahedra3D20)
    {
        ComputeShapeFunctionNormalGradientH20(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
    else if (r_element_geometry.GetGeometryType() == GeometryData::Kratos_Hexahedra3D27)
    {
        ComputeShapeFunctionNormalGradientH27(dNdn, r_element_geometry, r_edge_geometry, integration_points);
    }
}


void GhostPenaltyUtility::ProbeNeighbourElements(Element::Pointer p_element)
{
    WeakPointerVector<Element>& rNeighbours = p_element->GetValue(NEIGHBOUR_ELEMENTS);
    std::cout << "Neighbour elements of element " << p_element->Id() << ":";
    for (std::size_t i = 0; i < rNeighbours.size(); ++i)
    {
        std::cout << " " << rNeighbours[i].Id();
    }
    std::cout << std::endl;
}


std::pair<GeometryData::KratosGeometryType, std::vector<std::size_t> > GhostPenaltyUtility::FindCommonEdge(Element& r_elem_1, Element& r_elem_2)
{
    if (r_elem_1.GetGeometry().GetGeometryType() != r_elem_2.GetGeometry().GetGeometryType())
        KRATOS_THROW_ERROR(std::logic_error, "The geometry type is not the same", "")

    if (r_elem_1.GetGeometry().GetGeometryType() == GeometryData::Kratos_Triangle2D3)
    {
        return std::make_pair(GeometryData::Kratos_Line2D2, FindCommonEdgeT3(r_elem_1, r_elem_2));
    }
    else if (r_elem_1.GetGeometry().GetGeometryType() == GeometryData::Kratos_Triangle2D6)
    {
        return std::make_pair(GeometryData::Kratos_Line2D3, FindCommonEdgeT6(r_elem_1, r_elem_2));
    }
    else if (r_elem_1.GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4)
    {
        return std::make_pair(GeometryData::Kratos_Line2D2, FindCommonEdgeQ4(r_elem_1, r_elem_2));
    }
    else if (r_elem_1.GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8)
    {
        return std::make_pair(GeometryData::Kratos_Line2D3, FindCommonEdgeQ8(r_elem_1, r_elem_2));
    }
    else if (r_elem_1.GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9)
    {
        return std::make_pair(GeometryData::Kratos_Line2D3, FindCommonEdgeQ8(r_elem_1, r_elem_2));
    }
    else if (r_elem_1.GetGeometry().GetGeometryType() == GeometryData::Kratos_Tetrahedra3D4)
    {
        return std::make_pair(GeometryData::Kratos_Triangle3D3, FindCommonFaceT4(r_elem_1, r_elem_2));
    }
    else if (r_elem_1.GetGeometry().GetGeometryType() == GeometryData::Kratos_Tetrahedra3D10)
    {
        return std::make_pair(GeometryData::Kratos_Triangle3D6, FindCommonFaceT10(r_elem_1, r_elem_2));
    }
    else if (r_elem_1.GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D8)
    {
        return std::make_pair(GeometryData::Kratos_Quadrilateral3D4, FindCommonFaceH8(r_elem_1, r_elem_2));
    }
    else if (r_elem_1.GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D20)
    {
        return std::make_pair(GeometryData::Kratos_Quadrilateral3D8, FindCommonFaceH20(r_elem_1, r_elem_2));
    }
    else if (r_elem_1.GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D27)
    {
        return std::make_pair(GeometryData::Kratos_Quadrilateral3D9, FindCommonFaceH27(r_elem_1, r_elem_2));
    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")

    return std::make_pair(GeometryData::Kratos_generic_type, std::vector<std::size_t>{});
}

std::vector<std::size_t> GhostPenaltyUtility::FindCommonEdgeT3(Element& r_elem_1, Element& r_elem_2)
{
    for (std::size_t i = 0; i < 3; ++i)
    {
        const std::size_t& i1 = msEdgesT3[i][0];
        const std::size_t& i2 = msEdgesT3[i][1];
        const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
        const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();

        if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n2, r_elem_2.GetGeometry()))
            return std::vector<std::size_t>{i1, i2};
    }
    return std::vector<std::size_t>{};
}

std::vector<std::size_t> GhostPenaltyUtility::FindCommonEdgeT6(Element& r_elem_1, Element& r_elem_2)
{
    for (std::size_t i = 0; i < 3; ++i)
    {
        const std::size_t& i1 = msEdgesT6[i][0];
        const std::size_t& i2 = msEdgesT6[i][1];
        const std::size_t& i3 = msEdgesT6[i][2];
        const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
        const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
        const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();

        if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n3, r_elem_2.GetGeometry()))
            return std::vector<std::size_t>{i1, i2, i3};
    }
    return std::vector<std::size_t>{};
}

std::vector<std::size_t> GhostPenaltyUtility::FindCommonEdgeQ4(Element& r_elem_1, Element& r_elem_2)
{
    for (std::size_t i = 0; i < 4; ++i)
    {
        const std::size_t& i1 = msEdgesQ4[i][0];
        const std::size_t& i2 = msEdgesQ4[i][1];
        const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
        const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();

        if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n2, r_elem_2.GetGeometry()))
            return std::vector<std::size_t>{i1, i2};
    }
    return std::vector<std::size_t>{};
}

std::vector<std::size_t> GhostPenaltyUtility::FindCommonEdgeQ8(Element& r_elem_1, Element& r_elem_2)
{
    for (std::size_t i = 0; i < 4; ++i)
    {
        const std::size_t& i1 = msEdgesQ8[i][0];
        const std::size_t& i2 = msEdgesQ8[i][1];
        const std::size_t& i3 = msEdgesQ8[i][2];
        const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
        const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
        const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();

        if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n3, r_elem_2.GetGeometry()))
            return std::vector<std::size_t>{i1, i2, i3};
    }
    return std::vector<std::size_t>{};
}

std::vector<std::size_t> GhostPenaltyUtility::FindCommonFaceT4(Element& r_elem_1, Element& r_elem_2)
{
    for (std::size_t i = 0; i < 4; ++i)
    {
        const std::size_t& i1 = msFacesT4[i][0];
        const std::size_t& i2 = msFacesT4[i][1];
        const std::size_t& i3 = msFacesT4[i][2];
        const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
        const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
        const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();

        if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n3, r_elem_2.GetGeometry()))
            return std::vector<std::size_t>{i1, i2, i3};
    }
    return std::vector<std::size_t>{};
}

std::vector<std::size_t> GhostPenaltyUtility::FindCommonFaceT10(Element& r_elem_1, Element& r_elem_2)
{
    for (std::size_t i = 0; i < 4; ++i)
    {
        const std::size_t& i1 = msFacesT10[i][0];
        const std::size_t& i2 = msFacesT10[i][1];
        const std::size_t& i3 = msFacesT10[i][2];
        const std::size_t& i4 = msFacesT10[i][3];
        const std::size_t& i5 = msFacesT10[i][4];
        const std::size_t& i6 = msFacesT10[i][5];
        const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
        const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
        const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();
        const std::size_t& n4 = r_elem_1.GetGeometry()[i4].Id();
        const std::size_t& n5 = r_elem_1.GetGeometry()[i5].Id();
        const std::size_t& n6 = r_elem_1.GetGeometry()[i6].Id();

        if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n3, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n4, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n5, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n6, r_elem_2.GetGeometry()))
            return std::vector<std::size_t>{i1, i2, i3, i4, i5, i6};
    }
    return std::vector<std::size_t>{};
}

std::vector<std::size_t> GhostPenaltyUtility::FindCommonFaceH8(Element& r_elem_1, Element& r_elem_2)
{
    for (std::size_t i = 0; i < 8; ++i)
    {
        const std::size_t& i1 = msFacesH8[i][0];
        const std::size_t& i2 = msFacesH8[i][1];
        const std::size_t& i3 = msFacesH8[i][2];
        const std::size_t& i4 = msFacesH8[i][3];
        const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
        const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
        const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();
        const std::size_t& n4 = r_elem_1.GetGeometry()[i4].Id();

        if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n3, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n4, r_elem_2.GetGeometry()))
            return std::vector<std::size_t>{i1, i2, i3, i4};
    }
    return std::vector<std::size_t>{};
}

std::vector<std::size_t> GhostPenaltyUtility::FindCommonFaceH20(Element& r_elem_1, Element& r_elem_2)
{
    for (std::size_t i = 0; i < 8; ++i)
    {
        const std::size_t& i1 = msFacesH20[i][0];
        const std::size_t& i2 = msFacesH20[i][1];
        const std::size_t& i3 = msFacesH20[i][2];
        const std::size_t& i4 = msFacesH20[i][3];
        const std::size_t& i5 = msFacesH20[i][4];
        const std::size_t& i6 = msFacesH20[i][5];
        const std::size_t& i7 = msFacesH20[i][6];
        const std::size_t& i8 = msFacesH20[i][7];
        const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
        const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
        const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();
        const std::size_t& n4 = r_elem_1.GetGeometry()[i4].Id();
        const std::size_t& n5 = r_elem_1.GetGeometry()[i5].Id();
        const std::size_t& n6 = r_elem_1.GetGeometry()[i6].Id();
        const std::size_t& n7 = r_elem_1.GetGeometry()[i7].Id();
        const std::size_t& n8 = r_elem_1.GetGeometry()[i8].Id();

        if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n3, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n4, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n5, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n6, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n7, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n8, r_elem_2.GetGeometry()))
            return std::vector<std::size_t>{i1, i2, i3, i4, i5, i6, i7, i8};
    }
    return std::vector<std::size_t>{};
}

std::vector<std::size_t> GhostPenaltyUtility::FindCommonFaceH27(Element& r_elem_1, Element& r_elem_2)
{
    for (std::size_t i = 0; i < 8; ++i)
    {
        const std::size_t& i1 = msFacesH27[i][0];
        const std::size_t& i2 = msFacesH27[i][1];
        const std::size_t& i3 = msFacesH27[i][2];
        const std::size_t& i4 = msFacesH27[i][3];
        const std::size_t& i5 = msFacesH27[i][4];
        const std::size_t& i6 = msFacesH27[i][5];
        const std::size_t& i7 = msFacesH27[i][6];
        const std::size_t& i8 = msFacesH27[i][7];
        const std::size_t& i9 = msFacesH27[i][8];
        const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
        const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
        const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();
        const std::size_t& n4 = r_elem_1.GetGeometry()[i4].Id();
        const std::size_t& n5 = r_elem_1.GetGeometry()[i5].Id();
        const std::size_t& n6 = r_elem_1.GetGeometry()[i6].Id();
        const std::size_t& n7 = r_elem_1.GetGeometry()[i7].Id();
        const std::size_t& n8 = r_elem_1.GetGeometry()[i8].Id();
        const std::size_t& n9 = r_elem_1.GetGeometry()[i9].Id();

        if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n3, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n4, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n5, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n6, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n7, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n8, r_elem_2.GetGeometry())
            && IsBelongedToGeometry(n9, r_elem_2.GetGeometry()))
            return std::vector<std::size_t>{i1, i2, i3, i4, i5, i6, i7, i8, i9};
    }
    return std::vector<std::size_t>{};
}

bool GhostPenaltyUtility::IsBelongedToGeometry(const std::size_t& node_id, GeometryType& r_geom)
{
    for (std::size_t i = 0; i < r_geom.size(); ++i)
        if (r_geom[i].Id() == node_id)
            return true;
    return false;
}

bool GhostPenaltyUtility::IsSame(GeometryType& r_geom, const std::vector<std::size_t>& nodes)
{
    std::set<std::size_t> nodes_set(nodes.begin(), nodes.end());
    for (std::size_t i = 0; i < r_geom.size(); ++i)
        nodes_set.erase(r_geom[i].Id());
    return nodes_set.size() == 0;
}

bool GhostPenaltyUtility::BuildMapEdgeNodeIndexToElementNodeIndex(std::map<std::size_t, std::size_t>& map_edge_node_index_to_element_node_index,
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
        if (!found) return false;
    }
    return true;
}

void GhostPenaltyUtility::ComputeShapeFunctionNormalGradientT3(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{}

void GhostPenaltyUtility::ComputeShapeFunctionNormalGradientT6(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{}

void GhostPenaltyUtility::ComputeShapeFunctionNormalGradientQ4(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
{
    std::map<std::size_t, std::size_t> map_edge_node_index_to_element_node_index;

    std::cout << "element geometry:";
    for (std::size_t i = 0; i < r_element_geometry.size(); ++i)
        std::cout << " " << r_element_geometry[i].Id();
    std::cout << std::endl;

    std::cout << "edge geometry:";
    for (std::size_t i = 0; i < r_edge_geometry.size(); ++i)
        std::cout << " " << r_edge_geometry[i].Id();
    std::cout << std::endl;

    bool found = false;
    for (std::size_t side = 0; side < 4; ++side)
    {
        map_edge_node_index_to_element_node_index.clear();

        if (BuildMapEdgeNodeIndexToElementNodeIndex(map_edge_node_index_to_element_node_index,
                r_edge_geometry, r_element_geometry, msEdgesQ4[side]))
        {
            found = true;
            KRATOS_WATCH(side)

            Matrix DN_DX;
            GeometryType::IntegrationPointsArrayType integration_points(edge_integration_points.size());

            // construct the local coordinates system associated with the edge
            Vector t(2), n(2);
            double length = norm_2(r_element_geometry[msEdgesQ4[side][0]] - r_element_geometry[msEdgesQ4[side][1]]);
            t(0) = (r_element_geometry[msEdgesQ4[side][1]].X0() - r_element_geometry[msEdgesQ4[side][0]].X0()) / length;
            t(1) = (r_element_geometry[msEdgesQ4[side][1]].Y0() - r_element_geometry[msEdgesQ4[side][0]].Y0()) / length;
            n(0) = -t(1);
            n(1) = t(0);
            KRATOS_WATCH(t)
            KRATOS_WATCH(n)

            if (side == 0) // xi: n, eta: -1
            {
                for (std::size_t i = 0; i < edge_integration_points.size(); ++i)
                {
                    integration_points[i].X() = edge_integration_points[i].X();
                    integration_points[i].Y() = -1.0;
                    DN_DX = r_element_geometry.ShapeFunctionsLocalGradients(DN_DX, integration_points[i]);
//                    dNdn(0, i) = DN_DX(map_edge_node_index_to_element_node_index[0], 0);
//                    dNdn(1, i) = DN_DX(map_edge_node_index_to_element_node_index[1], 0);
                    dNdn(0, i) = DN_DX(map_edge_node_index_to_element_node_index[0], 0) * n(0) + DN_DX(map_edge_node_index_to_element_node_index[0], 1) * n(1);
                    dNdn(1, i) = DN_DX(map_edge_node_index_to_element_node_index[1], 0) * n(0) + DN_DX(map_edge_node_index_to_element_node_index[1], 1) * n(1);
                }
            }
            else if (side == 1) // xi: 1, eta: n
            {
                for (std::size_t i = 0; i < edge_integration_points.size(); ++i)
                {
                    integration_points[i].X() = 1.0;
                    integration_points[i].Y() = edge_integration_points[i].X();
                    DN_DX = r_element_geometry.ShapeFunctionsLocalGradients(DN_DX, integration_points[i]);
//                    dNdn(0, i) = DN_DX(map_edge_node_index_to_element_node_index[0], 1);
//                    dNdn(1, i) = DN_DX(map_edge_node_index_to_element_node_index[1], 1);
                    dNdn(0, i) = DN_DX(map_edge_node_index_to_element_node_index[0], 0) * n(0) + DN_DX(map_edge_node_index_to_element_node_index[0], 1) * n(1);
                    dNdn(1, i) = DN_DX(map_edge_node_index_to_element_node_index[1], 0) * n(0) + DN_DX(map_edge_node_index_to_element_node_index[1], 1) * n(1);
                }
            }
            else if (side == 2) // xi: n, eta: 1
            {
                for (std::size_t i = 0; i < edge_integration_points.size(); ++i)
                {
                    integration_points[i].X() = edge_integration_points[i].X();
                    integration_points[i].Y() = 1.0;
                    DN_DX = r_element_geometry.ShapeFunctionsLocalGradients(DN_DX, integration_points[i]);
//                    dNdn(0, i) = DN_DX(map_edge_node_index_to_element_node_index[0], 0);
//                    dNdn(1, i) = DN_DX(map_edge_node_index_to_element_node_index[1], 0);
                    dNdn(0, i) = DN_DX(map_edge_node_index_to_element_node_index[0], 0) * n(0) + DN_DX(map_edge_node_index_to_element_node_index[0], 1) * n(1);
                    dNdn(1, i) = DN_DX(map_edge_node_index_to_element_node_index[1], 0) * n(0) + DN_DX(map_edge_node_index_to_element_node_index[1], 1) * n(1);
                }
            }
            else if (side == 3) // xi: -1, eta: n
            {
                for (std::size_t i = 0; i < edge_integration_points.size(); ++i)
                {
                    integration_points[i].X() = -1.0;
                    integration_points[i].Y() = edge_integration_points[i].X();
                    DN_DX = r_element_geometry.ShapeFunctionsLocalGradients(DN_DX, integration_points[i]);
//                    dNdn(0, i) = DN_DX(map_edge_node_index_to_element_node_index[0], 1);
//                    dNdn(1, i) = DN_DX(map_edge_node_index_to_element_node_index[1], 1);
                    dNdn(0, i) = DN_DX(map_edge_node_index_to_element_node_index[0], 0) * n(0) + DN_DX(map_edge_node_index_to_element_node_index[0], 1) * n(1);
                    dNdn(1, i) = DN_DX(map_edge_node_index_to_element_node_index[1], 0) * n(0) + DN_DX(map_edge_node_index_to_element_node_index[1], 1) * n(1);
                }
            }
        }
    }

    if (!found)
        KRATOS_THROW_ERROR(std::logic_error, "The edge geometry is not found on the element geometry. Check the input", "")
}

void GhostPenaltyUtility::ComputeShapeFunctionNormalGradientQ8(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{}

void GhostPenaltyUtility::ComputeShapeFunctionNormalGradientQ9(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{}

void GhostPenaltyUtility::ComputeShapeFunctionNormalGradientT4(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{}

void GhostPenaltyUtility::ComputeShapeFunctionNormalGradientT10(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{}

void GhostPenaltyUtility::ComputeShapeFunctionNormalGradientH8(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{}

void GhostPenaltyUtility::ComputeShapeFunctionNormalGradientH20(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{}

void GhostPenaltyUtility::ComputeShapeFunctionNormalGradientH27(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points)
{}


}  // namespace Kratos.

