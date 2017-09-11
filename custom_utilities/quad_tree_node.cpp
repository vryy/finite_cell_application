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
void QuadTreeNode::AddToModelPart<false, Element>(GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Element const& r_clone_element,
        std::size_t& lastNodeId,
        std::size_t& lastElementId,
        const int level) const
{
    // std::cout << __FUNCTION__ << " is called" << std::endl;
    Properties::Pointer p_properties = r_model_part.pGetProperties(level);
    GeometryType::Pointer pThisGeometry = this->pCreateGeometry(pParentGeometry);
    Element::Pointer pNewElement = r_clone_element.Create(++lastElementId, pThisGeometry, p_properties);
    r_model_part.AddElement(pNewElement);

    for(std::size_t i = 0; i < pThisGeometry->size(); ++i)
    {
        (*pThisGeometry)[i].SetId(++lastNodeId);
        (*pThisGeometry)[i].SetSolutionStepVariablesList(&r_model_part.GetNodalSolutionStepVariablesList());
        (*pThisGeometry)[i].SetBufferSize(r_model_part.GetBufferSize());
        r_model_part.AddNode((*pThisGeometry)(i));
    }
}


template<>
void QuadTreeNode::AddToModelPart<false, Condition>(GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Condition const& r_clone_condition,
        std::size_t& lastNodeId,
        std::size_t& lastCondId,
        const int level) const
{
    // std::cout << __FUNCTION__ << " is called" << std::endl;
    Properties::Pointer p_properties = r_model_part.pGetProperties(level);
    GeometryType::Pointer pThisGeometry = this->pCreateGeometry(pParentGeometry);
    Condition::Pointer pNewCondition = r_clone_condition.Create(++lastCondId, pThisGeometry, p_properties);
    r_model_part.AddCondition(pNewCondition);

    for(std::size_t i = 0; i < pThisGeometry->size(); ++i)
    {
        (*pThisGeometry)[i].SetId(++lastNodeId);
        (*pThisGeometry)[i].SetSolutionStepVariablesList(&r_model_part.GetNodalSolutionStepVariablesList());
        (*pThisGeometry)[i].SetBufferSize(r_model_part.GetBufferSize());
        r_model_part.AddNode((*pThisGeometry)(i));
    }
}


template<>
void QuadTreeNode::AddToModelPart<true, Element>(GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Element const& r_clone_element,
        std::size_t& lastNodeId,
        std::size_t& lastElementId,
        const int level) const
{
    // std::cout << __FUNCTION__ << " is called" << std::endl;
    if(this->IsLeaf())
    {
        this->AddToModelPart<false, Element>(pParentGeometry, r_model_part, r_clone_element, lastNodeId, lastElementId, level);
    }
    else
    {
        for(std::size_t i = 0; i < mpChildren.size(); ++i)
            mpChildren[i]->AddToModelPart<true, Element>(pParentGeometry, r_model_part, r_clone_element, lastNodeId, lastElementId, level + 1);
    }
}

template<>
void QuadTreeNode::AddToModelPart<true, Condition>(GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Condition const& r_clone_condition,
        std::size_t& lastNodeId,
        std::size_t& lastCondId,
        const int level) const
{
    // std::cout << __FUNCTION__ << " is called" << std::endl;
    if(this->IsLeaf())
    {
        this->AddToModelPart<false, Condition>(pParentGeometry, r_model_part, r_clone_condition, lastNodeId, lastCondId, level);
    }
    else
    {
        for(std::size_t i = 0; i < mpChildren.size(); ++i)
            mpChildren[i]->AddToModelPart<true, Condition>(pParentGeometry, r_model_part, r_clone_condition, lastNodeId, lastCondId, level + 1);
    }
}


}

