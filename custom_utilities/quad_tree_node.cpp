// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//


// Project includes
#include "includes/define.h"
#include "custom_utilities/quad_tree_node.h"


namespace Kratos
{

template<>
struct QuadTreeNode_AddToModelPart_Helper<false, Element>
{
    template<class TTreeType>
    static void Execute(const TTreeType& r_tree,
        Element::GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Element const& r_clone_element,
        std::size_t& lastNodeId,
        std::size_t& lastElementId,
        std::vector<std::size_t>& new_node_ids,
        std::vector<std::size_t>& new_entity_ids,
        const int level)
    {
        // std::cout << __FUNCTION__ << " is called" << std::endl;
        Properties::Pointer p_properties = r_model_part.pGetProperties(level);
        Element::GeometryType::Pointer pThisGeometry = r_tree.pCreateGeometry(pParentGeometry);
        Element::Pointer pNewElement = r_clone_element.Create(++lastElementId, pThisGeometry, p_properties);
        new_entity_ids.push_back(lastElementId);
        r_model_part.AddElement(pNewElement);

        for(std::size_t i = 0; i < pThisGeometry->size(); ++i)
        {
            (*pThisGeometry)[i].SetId(++lastNodeId);
            new_node_ids.push_back(lastNodeId);
            (*pThisGeometry)[i].SetSolutionStepVariablesList(&r_model_part.GetNodalSolutionStepVariablesList());
            (*pThisGeometry)[i].SetBufferSize(r_model_part.GetBufferSize());
            r_model_part.AddNode((*pThisGeometry)(i));
        }
    }
};


template<>
struct QuadTreeNode_AddToModelPart_Helper<false, Condition>
{
    template<class TTreeType>
    static void Execute(const TTreeType& r_tree,
        Condition::GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Condition const& r_clone_condition,
        std::size_t& lastNodeId,
        std::size_t& lastCondId,
        std::vector<std::size_t>& new_node_ids,
        std::vector<std::size_t>& new_entity_ids,
        const int level)
    {
        // std::cout << __FUNCTION__ << " is called" << std::endl;
        Properties::Pointer p_properties = r_model_part.pGetProperties(level);
        Condition::GeometryType::Pointer pThisGeometry = r_tree.pCreateGeometry(pParentGeometry);
        Condition::Pointer pNewCondition = r_clone_condition.Create(++lastCondId, pThisGeometry, p_properties);
        new_entity_ids.push_back(lastCondId);
        r_model_part.AddCondition(pNewCondition);

        for(std::size_t i = 0; i < pThisGeometry->size(); ++i)
        {
            (*pThisGeometry)[i].SetId(++lastNodeId);
            new_node_ids.push_back(lastNodeId);
            (*pThisGeometry)[i].SetSolutionStepVariablesList(&r_model_part.GetNodalSolutionStepVariablesList());
            (*pThisGeometry)[i].SetBufferSize(r_model_part.GetBufferSize());
            r_model_part.AddNode((*pThisGeometry)(i));
        }
    }
};

template<>
struct QuadTreeNode_AddToModelPart_Helper<true, Element>
{
    template<class TTreeType>
    static void Execute(const TTreeType& r_tree,
        Element::GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Element const& r_clone_element,
        std::size_t& lastNodeId,
        std::size_t& lastElementId,
        std::vector<std::size_t>& new_node_ids,
        std::vector<std::size_t>& new_entity_ids,
        const int level)
    {
        // std::cout << __FUNCTION__ << " is called" << std::endl;
        if(r_tree.IsLeaf())
        {
            if (level > 1) // this is to avoid adding duplicated element with the original mesh
                QuadTreeNode_AddToModelPart_Helper<false, Element>::Execute(r_tree, pParentGeometry, r_model_part, r_clone_element, lastNodeId, lastElementId, new_node_ids, new_entity_ids, level);
        }
        else
        {
            for(std::size_t i = 0; i < r_tree.Size(); ++i)
                QuadTreeNode_AddToModelPart_Helper<true, Element>::Execute(*(r_tree.pChild(i)), pParentGeometry, r_model_part, r_clone_element, lastNodeId, lastElementId, new_node_ids, new_entity_ids, level + 1);
        }
    }
};

template<>
struct QuadTreeNode_AddToModelPart_Helper<true, Condition>
{
    template<class TTreeType>
    static void Execute(const TTreeType& r_tree,
        Condition::GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Condition const& r_clone_condition,
        std::size_t& lastNodeId,
        std::size_t& lastCondId,
        std::vector<std::size_t>& new_node_ids,
        std::vector<std::size_t>& new_entity_ids,
        const int level)
    {
        // std::cout << __FUNCTION__ << " is called" << std::endl;
        if(r_tree.IsLeaf())
        {
            if (level > 1) // this is to avoid adding duplicated element with the original mesh
                QuadTreeNode_AddToModelPart_Helper<false, Condition>::Execute(r_tree, pParentGeometry, r_model_part, r_clone_condition, lastNodeId, lastCondId, new_node_ids, new_entity_ids, level);
        }
        else
        {
            for(std::size_t i = 0; i < r_tree.Size(); ++i)
                QuadTreeNode_AddToModelPart_Helper<true, Condition>::Execute(*(r_tree.pChild(i)), pParentGeometry, r_model_part, r_clone_condition, lastNodeId, lastCondId, new_node_ids, new_entity_ids, level + 1);
        }
    }
};

template<int TFrameType>
bool QuadTreeNode<TFrameType>::CenterOfGravity(PointType& COG, GeometryType::Pointer pParentGeometry, const BRep& r_brep, const int& integration_method) const
{
    FunctionR3R1::Pointer FX = FunctionR3R1::Pointer(new MonomialFunctionR3R1<1, 0, 0>());
    FunctionR3R1::Pointer FY = FunctionR3R1::Pointer(new MonomialFunctionR3R1<0, 1, 0>());
    FunctionR3R1::Pointer FZ = FunctionR3R1::Pointer(new MonomialFunctionR3R1<0, 0, 1>());
    FunctionR3R1::Pointer FH = FunctionR3R1::Pointer(new HeavisideFunction<FunctionR3R1>(r_brep));

    GeometryType::IntegrationPointsArrayType integration_points;
    this->ConstructQuadrature(pParentGeometry, integration_points, integration_method);

    double X = 0.0, Y = 0.0, Z = 0.0, A = 0.0;
    this->IntegrateImpl<double, GeometryType::IntegrationPointsArrayType, GLOBAL>(*pParentGeometry, ProductFunction<FunctionR3R1>(FX, FH), X, integration_points);
    this->IntegrateImpl<double, GeometryType::IntegrationPointsArrayType, GLOBAL>(*pParentGeometry, ProductFunction<FunctionR3R1>(FY, FH), Y, integration_points);
    this->IntegrateImpl<double, GeometryType::IntegrationPointsArrayType, GLOBAL>(*pParentGeometry, ProductFunction<FunctionR3R1>(FZ, FH), Z, integration_points);
    this->IntegrateImpl<double, GeometryType::IntegrationPointsArrayType, GLOBAL>(*pParentGeometry, *FH, A, integration_points);
//        std::cout << "X: " << X << ", Y: " << Y << ", Z: " << Z << ", A: " << A << std::endl;
    if(A != 0.0)
    {
        // the quad-tree is able to integrate the COG
        COG[0] = X/A;
        COG[1] = Y/A;
        COG[2] = Z/A;
        return true;
    }
    else
        return false;
}

template class QuadTreeNode<GLOBAL_REFERENCE>;
template class QuadTreeNode<GLOBAL_CURRENT>;
template class QuadTreeNodeQ4<GLOBAL_REFERENCE>;
template class QuadTreeNodeQ4<GLOBAL_CURRENT>;
#ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
template class QuadTreeNodeBezier2D<GLOBAL_REFERENCE>;
template class QuadTreeNodeBezier2D<GLOBAL_CURRENT>;
template class QuadTreeNodeBezier3D<GLOBAL_REFERENCE>;
template class QuadTreeNodeBezier3D<GLOBAL_CURRENT>;
#endif
template class QuadTreeNodeH8<GLOBAL_REFERENCE>;
template class QuadTreeNodeH8<GLOBAL_CURRENT>;
template class QuadTreeNodeT3<GLOBAL_REFERENCE>;
template class QuadTreeNodeT3<GLOBAL_CURRENT>;
template class QuadTreeNodeT4<GLOBAL_REFERENCE>;
template class QuadTreeNodeT4<GLOBAL_CURRENT>;

}

