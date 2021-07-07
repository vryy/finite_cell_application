// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_FINITE_CELL_APPLICATION_ADD_QUADTREE_TO_PYTHON_HPP_INCLUDED )
#define  KRATOS_FINITE_CELL_APPLICATION_ADD_QUADTREE_TO_PYTHON_HPP_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/moment_fitted_quad_tree_subcell.h"
#include "custom_utilities/finite_cell_geometry_utility.h"
#include "custom_utilities/finite_cell_auxiliary_utility.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

void FiniteCellApplication_AddRefinableTreeToPython(pybind11::module& m)
{
    class_<RefinableTree, RefinableTree::Pointer>
    (m, "RefinableTree")
    .def(init<>())
    .def("Refine", &RefinableTree::Refine)
    .def("RefineBy", &RefinableTree::RefineBy)
    .def("RefineWithGeometryBy", &RefinableTree::RefineWithGeometryBy)
    .def("RefineUpToLevelBy", &RefinableTree::RefineUpToLevelBy)
    .def("RefineWithGeometryUpToLevelBy", &RefinableTree::RefineWithGeometryUpToLevelBy)
    .def("Clear", &RefinableTree::Clear)
    ;
}

void FiniteCellApplication_AddFunctionIntegratorToPython(pybind11::module& m)
{
    double(FunctionIntegrator::*pointer_to_IntegrateLocal_double)(const FunctionR3R1&, const int&) const = &FunctionIntegrator::IntegrateLocal;
    double(FunctionIntegrator::*pointer_to_IntegrateGlobal_double)(const FunctionR3R1&, const int&) const = &FunctionIntegrator::IntegrateGlobal;

    class_<FunctionIntegrator, FunctionIntegrator::Pointer>
    (m, "FunctionIntegrator")
    .def(init<>())
    .def("IntegrateLocal", pointer_to_IntegrateLocal_double)
    .def("IntegrateGlobal", pointer_to_IntegrateGlobal_double)
    ;
}

void FiniteCellApplication_AddBaseMomentFittedQuadTreeSubCellToPython(pybind11::module& m)
{
    class_<BaseMomentFittedQuadTreeSubCell, BaseMomentFittedQuadTreeSubCell::Pointer>
    (m, "BaseMomentFittedQuadTreeSubCell")
    .def(init<Element::Pointer>())
    .def("GeneratePhysicalIntegrationPoints", &BaseMomentFittedQuadTreeSubCell::GeneratePhysicalIntegrationPoints)
    ;
}

/// construct the element out from quad-tree and add to model_part.
/// This is mainly for post-processing
template<class TTreeType>
pybind11::list QuadTree_AddToModelPart_WithLevel(TTreeType& rDummy,
    ModelPart& r_model_part, const std::string sample_entity_name,
    std::size_t lastNodeId, std::size_t lastEntityId, int start_level)
{
    std::vector<std::size_t> new_node_ids;
    std::vector<std::size_t> new_entity_ids;

    if( KratosComponents<Element>::Has(sample_entity_name) )
    {
        Element const& r_clone_element = KratosComponents<Element>::Get(sample_entity_name);
        rDummy.Get().template AddToModelPart<true, Element>(rDummy.pGetGeometry(), r_model_part, r_clone_element, lastNodeId, lastEntityId, new_node_ids, new_entity_ids, start_level);
    }
    else if( KratosComponents<Condition>::Has(sample_entity_name) )
    {
        Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_entity_name);
        rDummy.Get().template AddToModelPart<true, Condition>(rDummy.pGetGeometry(), r_model_part, r_clone_condition, lastNodeId, lastEntityId, new_node_ids, new_entity_ids, start_level);
    }
    else
    {
        std::stringstream ss;
        ss << sample_entity_name << " is not registerred to the database";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    pybind11::list list;
    list.append(lastNodeId);
    list.append(lastEntityId);

    pybind11::list list_new_nodes;
    pybind11::list list_new_entities;

    for (std::size_t i = 0; i < new_node_ids.size(); ++i)
        list_new_nodes.append(new_node_ids[i]);
    list.append(list_new_nodes);

    for (std::size_t i = 0; i < new_entity_ids.size(); ++i)
        list_new_entities.append(new_entity_ids[i]);
    list.append(list_new_entities);

    return list;
}

/// construct the element out from quad-tree and add to model_part.
/// This is mainly for post-processing
template<class TTreeType>
pybind11::list QuadTree_AddToModelPart(TTreeType& rDummy,
    ModelPart& r_model_part, const std::string sample_entity_name,
    std::size_t lastNodeId, std::size_t lastEntityId)
{
    return QuadTree_AddToModelPart_WithLevel(rDummy, r_model_part, sample_entity_name, lastNodeId, lastEntityId, 1);
}

template<class TTreeType>
pybind11::list QuadTree_ConstructQuadrature1(TTreeType& rDummy,
    const int& integration_method)
{
    std::size_t output = rDummy.ConstructQuadrature(integration_method);
    pybind11::list list_output;
    list_output.append(output);
    return list_output;
}

template<class TTreeType>
pybind11::list QuadTree_ConstructQuadrature2(TTreeType& rDummy,
    const BRep& r_brep, const int& integration_method, const double& small_weight)
{
    std::vector<std::size_t> output = rDummy.ConstructQuadrature(r_brep, integration_method, small_weight);
    pybind11::list list_output;
    for (std::size_t i = 0; i < output.size(); ++i)
        list_output.append(output[i]);
    return list_output;
}

template<class TTreeType>
pybind11::list QuadTree_ConstructQuadratureNoSkip(TTreeType& rDummy,
    const BRep& r_brep, const int& integration_method, const double& small_weight)
{
    std::vector<std::size_t> output = rDummy.ConstructQuadratureNoSkip(r_brep, integration_method, small_weight);
    pybind11::list list_output;
    for (std::size_t i = 0; i < output.size(); ++i)
        list_output.append(output[i]);
    return list_output;
}

/// construct the element out from quad-tree and add to model_part
/// This is mainly for post-processing
template<class TTreeType, bool TShallow = true>
pybind11::list QuadTreeSubCell_AddToModelPart(TTreeType& rDummy,
    ModelPart& r_model_part, const std::string sample_entity_name,
    std::size_t lastNodeId, std::size_t lastEntityId, std::size_t starting_prop_id)
{
    std::vector<std::size_t> new_node_ids;
    std::vector<std::size_t> new_entity_ids;

    if( KratosComponents<Element>::Has(sample_entity_name) )
    {
        Element const& r_clone_element = KratosComponents<Element>::Get(sample_entity_name);
        for(std::size_t i = 0; i < rDummy.NumberOfSubCells(); ++i)
        {
            rDummy.Get(i).template AddToModelPart<false, Element>(rDummy.pGetGeometry(), r_model_part, r_clone_element, lastNodeId, lastEntityId, new_node_ids, new_entity_ids, starting_prop_id);
            if(!TShallow)
                rDummy.Get(i).template AddToModelPart<true, Element>(rDummy.pGetGeometry(), r_model_part, r_clone_element, lastNodeId, lastEntityId, new_node_ids, new_entity_ids, starting_prop_id+1);
        }
    }
    else if( KratosComponents<Condition>::Has(sample_entity_name) )
    {
        Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_entity_name);
        for(std::size_t i = 0; i < rDummy.NumberOfSubCells(); ++i)
        {
            rDummy.Get(i).template AddToModelPart<false, Condition>(rDummy.pGetGeometry(), r_model_part, r_clone_condition, lastNodeId, lastEntityId, new_node_ids, new_entity_ids, starting_prop_id);
            if(!TShallow)
                rDummy.Get(i).template AddToModelPart<true, Condition>(rDummy.pGetGeometry(), r_model_part, r_clone_condition, lastNodeId, lastEntityId, new_node_ids, new_entity_ids, starting_prop_id+1);
        }
    }
    else
    {
        std::stringstream ss;
        ss << sample_entity_name << " is not registerred to the Kratos database";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    pybind11::list list;
    list.append(lastNodeId);
    list.append(lastEntityId);

    pybind11::list list_new_nodes;
    pybind11::list list_new_entities;

    for (std::size_t i = 0; i < new_node_ids.size(); ++i)
        list_new_nodes.append(new_node_ids[i]);
    list.append(list_new_nodes);

    for (std::size_t i = 0; i < new_entity_ids.size(); ++i)
        list_new_entities.append(new_entity_ids[i]);
    list.append(list_new_entities);

    return list;
}

/// construct the element out from quad-tree and add to model_part
/// This is mainly for post-processing
template<class TTreeType, bool TShallow = true>
pybind11::list QuadTreeSubCell_AddToModelPart2(TTreeType& rDummy,
    ModelPart& r_model_part, const std::string sample_entity_name,
    std::size_t lastNodeId, std::size_t lastEntityId)
{
    return QuadTreeSubCell_AddToModelPart<TTreeType, TShallow>(rDummy, r_model_part, sample_entity_name, lastNodeId, lastEntityId, 1);
}

/// Get the list of physical integration points
template<class TTreeType>
pybind11::list MomentFittedQuadTreeSubCell_GetPhysicalIntegrationPoints(TTreeType& rDummy)
{
    const Element::GeometryType::IntegrationPointsArrayType& physical_integration_points = rDummy.GetPhysicalIntegrationPoints();
//KRATOS_WATCH(physical_integration_points.size())
    pybind11::list Output;

    for (std::size_t i = 0; i < physical_integration_points.size(); ++i)
    {
        pybind11::list point;
        point.append(physical_integration_points[i].X());
        point.append(physical_integration_points[i].Y());
        point.append(physical_integration_points[i].Z());
        point.append(physical_integration_points[i].Weight());
        Output.append(point);
    }

    return Output;
}

/// Get the list of fictitious integration points
template<class TTreeType>
pybind11::list MomentFittedQuadTreeSubCell_GetFictitiousIntegrationPoints(TTreeType& rDummy)
{
    const Element::GeometryType::IntegrationPointsArrayType& fictitious_integration_points = rDummy.GetFictitiousIntegrationPoints();

    pybind11::list Output;

    for (std::size_t i = 0; i < fictitious_integration_points.size(); ++i)
    {
        pybind11::list point;
        point.append(fictitious_integration_points[i].X());
        point.append(fictitious_integration_points[i].Y());
        point.append(fictitious_integration_points[i].Z());
        point.append(fictitious_integration_points[i].Weight());
        Output.append(point);
    }

    return Output;
}

/// Create the element out from sub-cells. The element takes the same geometry of the original element, but the weight is given.
/// This is the python interface
template<class TTreeType>
pybind11::list MomentFittedQuadTreeSubCell_CreateSubCellElements(TTreeType& rDummy,
        ModelPart& r_model_part,
        const std::string& subcell_element_type,
        const int& cut_cell_quadrature_order,
        const pybind11::list& cut_cell_full_quadrature,
        const pybind11::list& subcell_weights,
        const std::size_t& lastElementId,
        const std::size_t& lastCondId,
        Properties::Pointer pProperties)
{
    typedef Element::GeometryType GeometryType;

    ModelPart::ElementsContainerType NewElements;

//    if(!pybind11::list::is_empty(subcell_weights))
//    {
    // TODO find a way to check if a list is empty
        GeometryType::IntegrationPointsArrayType integration_points;
        Matrix Weights;

        std::size_t num_physical_point = pybind11::len(subcell_weights);
        std::size_t weight_length = pybind11::len(subcell_weights[0]);
        Weights.resize(num_physical_point, weight_length, false);

        for(std::size_t i = 0; i < num_physical_point; ++i)
        {
            pybind11::list weights = subcell_weights[i].cast<pybind11::list>();
            for(std::size_t j = 0; j < weight_length; ++j)
                Weights(i, j) = weights[j].cast<double>();
        }
//        KRATOS_WATCH(Weights)

        for(std::size_t i = 0; i < pybind11::len(cut_cell_full_quadrature); ++i)
        {
            pybind11::list point = cut_cell_full_quadrature[i].cast<pybind11::list>();
            GeometryType::IntegrationPointType integration_point;
            integration_point.X() = point[0].cast<double>();
            integration_point.Y() = point[1].cast<double>();
            integration_point.Z() = point[2].cast<double>();
            integration_point.Weight() = point[3].cast<double>();
//            KRATOS_WATCH(integration_point)
            integration_points.push_back(integration_point);
        }

        std::size_t new_lastElementId = lastElementId;
        std::size_t new_lastCondId = lastCondId;
        NewElements = rDummy.CreateSubCellElements(r_model_part,
            rDummy.pGetElement(),
            subcell_element_type,
            cut_cell_quadrature_order,
            integration_points,
            Weights,
            new_lastElementId,
            new_lastCondId,
            pProperties);
//        std::cout << "----------------------" << std::endl;
//    }

    pybind11::list Output;
    Output.append(new_lastElementId);
    Output.append(new_lastCondId);
    Output.append(NewElements);

    return Output;
}

/// Create the element out from sub-cells. The element takes the same geometry of the original element.
/// This is the python interface
template<class TTreeType>
pybind11::list MomentFittedQuadTreeSubCell_CreateSubCellElements2(TTreeType& rDummy,
        ModelPart& r_model_part,
        const std::string& subcell_element_type,
        const int& cut_cell_quadrature_order,
        const std::size_t& lastElementId,
        const std::size_t& lastCondId,
        Properties::Pointer pProperties)
{
    typedef Element::GeometryType GeometryType;

    ModelPart::ElementsContainerType NewElements;

    const Matrix& Weights = rDummy.pGetElement()->GetValue(SUBCELL_WEIGHTS);
//    KRATOS_WATCH(Weights)

    const GeometryType::IntegrationPointsArrayType& integration_points = rDummy.GetRepresentativeIntegrationPoints();

    std::size_t new_lastElementId = lastElementId;
    std::size_t new_lastCondId = lastCondId;
    NewElements = rDummy.CreateSubCellElements(r_model_part,
        rDummy.pGetElement(),
        subcell_element_type,
        cut_cell_quadrature_order,
        integration_points,
        Weights,
        new_lastElementId,
        new_lastCondId,
        pProperties);

    pybind11::list Output;
    Output.append(new_lastElementId);
    Output.append(new_lastCondId);
    Output.append(NewElements);

    return Output;
}

/// Create the element out from sub-cells. The element takes the same geometry of the original element, but the weight is fitted by sub-cell.
template<class TTreeType>
ModelPart::ElementsContainerType MomentFittedQuadTreeSubCell_FitAndCreateSubCellElements(TTreeType& rDummy,
    ModelPart& r_model_part,
    const std::string& sample_element_name,
    const pybind11::list& r_funcs,
    const BRep& r_brep,
    const int& integrator_integration_method,
    const std::string& solver_type,
    const int& echo_level,
    const double& small_weight,
    const bool& compute_subcell_domain_size)
{
    /* firstly compute the physical integration point */
    typedef Element::GeometryType GeometryType;

    rDummy.GeneratePhysicalIntegrationPoints(r_brep, integrator_integration_method);
    const std::vector<std::size_t>& subcell_index = rDummy.SubCellIndices();
    const GeometryType::IntegrationPointsArrayType& physical_integration_points = rDummy.GetPhysicalIntegrationPoints();

    /* secondly assign the quadrature for parent element based on physical integration_points of the previous step */
    GeometryData::IntegrationMethod RepresentativeIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(rDummy.GetRepresentativeIntegrationOrder());
//        KRATOS_WATCH(RepresentativeIntegrationMethod)
    FiniteCellGeometryUtility::AssignGeometryData(*(rDummy.pGetGeometry()), RepresentativeIntegrationMethod, physical_integration_points);
    Variable<int>& INTEGRATION_ORDER_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("INTEGRATION_ORDER"));
    ProcessInfo& CurrentProcessInfo =  r_model_part.GetProcessInfo();
    rDummy.pGetElement()->SetValue(INTEGRATION_ORDER_var, rDummy.GetRepresentativeIntegrationOrder());
    rDummy.pGetElement()->Initialize(CurrentProcessInfo);

    /* thirdly fit the sub-cell */
    std::vector<FunctionR3R1::Pointer> funcs;
    for (auto f : r_funcs)
        funcs.push_back(f.cast<FunctionR3R1::Pointer>());

    Matrix Weights = rDummy.FitQuadratureSubCell(subcell_index, funcs, r_brep, integrator_integration_method, solver_type, echo_level, small_weight);

    /* next create the sub-elements */
    // find the last element id
    std::size_t lastElementId = FiniteCellAuxiliaryUtility::GetLastElementId(r_model_part);

    // find the last condition id
    std::size_t lastCondId = FiniteCellAuxiliaryUtility::GetLastConditionId(r_model_part);

    ModelPart::ElementsContainerType NewSubCellElements = rDummy.CreateSubCellElements(r_model_part, sample_element_name, Weights, lastElementId, lastCondId);

    // set the domain size for the sub-cell
    if(compute_subcell_domain_size)
    {
        std::size_t cnt = 0;
        for(typename ModelPart::ElementsContainerType::ptr_iterator it = NewSubCellElements.ptr_begin();
                it != NewSubCellElements.ptr_end(); ++it)
        {
            (*it)->SetValue(SUBCELL_DOMAIN_SIZE, rDummy.DomainSize(subcell_index[cnt++], r_brep, integrator_integration_method));
        }
    }

    return NewSubCellElements;
}

template<std::size_t TNsampling, int TFrameType>
void FiniteCellApplication_AddQuadTreeToPython(pybind11::module& m)
{
    typedef QuadTree<TNsampling, TFrameType> QuadTreeType;

    std::stringstream QuadTreeName;
    QuadTreeName << "QuadTreeLocal";
    if(TFrameType == GLOBAL_REFERENCE)
        QuadTreeName << "Reference";
    else if(TFrameType == GLOBAL_CURRENT)
        QuadTreeName << "Current";
    if(TNsampling > 1)
        QuadTreeName << TNsampling;

    class_<QuadTreeType, typename QuadTreeType::Pointer, RefinableTree, FunctionIntegrator>
    (m, QuadTreeName.str().c_str())
    .def(init<Element::Pointer>())
    .def(init<Condition::Pointer>())
    .def(init<Element::Pointer, const int&>())
    .def(init<Condition::Pointer, const int&>())
    .def("pGetGeometry", &QuadTreeType::pGetGeometry)
    .def("Reset", &QuadTreeType::Reset)
    .def("LastLevel", &QuadTreeType::LastLevel)
    .def("NumberOfCells", &QuadTreeType::NumberOfCells)
    .def("NumberOfPhysicalPoints", &QuadTreeType::NumberOfPhysicalPoints)
    .def("NumberOfFictitiousPoints", &QuadTreeType::NumberOfFictitiousPoints)
    .def("DomainSize", &QuadTreeType::DomainSize)
    .def("CenterOfGravity", &QuadTreeType::CenterOfGravity)
    .def("AddToModelPartWithLevel", &QuadTree_AddToModelPart_WithLevel<QuadTreeType>)
    .def("AddToModelPart", &QuadTree_AddToModelPart<QuadTreeType>)
    .def("ConstructQuadrature", &QuadTree_ConstructQuadrature1<QuadTreeType>)
    .def("ConstructQuadrature", &QuadTree_ConstructQuadrature2<QuadTreeType>)
    .def("ConstructQuadratureNoSkip", &QuadTree_ConstructQuadratureNoSkip<QuadTreeType>)
    .def("__str__", &PrintObject<QuadTreeType>)
    ;
}

template<std::size_t TNsampling, int TFrameType>
void FiniteCellApplication_AddQuadTreeSubCellToPython(pybind11::module& m)
{
    typedef QuadTreeSubCell<TNsampling, TFrameType> QuadTreeSubCellType;
    typedef MomentFittedQuadTreeSubCell<TNsampling, TFrameType> MomentFittedQuadTreeSubCellType;

    double(QuadTreeSubCellType::*pointer_to_DomainSize1)(const std::size_t&, const BRep&, const int&) const = &QuadTreeSubCellType::DomainSize;
    double(QuadTreeSubCellType::*pointer_to_DomainSize2)(const BRep&, const int&) const = &QuadTreeSubCellType::DomainSize;
    std::size_t(QuadTreeSubCellType::*pointer_to_ConstructQuadrature1)(const int&) const = &QuadTreeSubCellType::ConstructQuadrature;
    std::size_t(QuadTreeSubCellType::*pointer_to_ConstructQuadrature2)(const BRep&, const int&, const double&) const = &QuadTreeSubCellType::ConstructQuadrature;

    std::stringstream QuadTreeSubCellName;
    QuadTreeSubCellName << "QuadTreeSubCell";
    if(TFrameType == GLOBAL_REFERENCE)
        QuadTreeSubCellName << "Reference";
    else if(TFrameType == GLOBAL_CURRENT)
        QuadTreeSubCellName << "Current";
    if(TNsampling > 1)
        QuadTreeSubCellName << TNsampling;

    class_<QuadTreeSubCellType, typename QuadTreeSubCellType::Pointer, RefinableTree>
    (m, QuadTreeSubCellName.str().c_str())
    .def(init<Element::Pointer>())
    .def(init<Condition::Pointer>())
    .def("NumberOfSubCells", &QuadTreeSubCellType::NumberOfSubCells)
    .def("DomainSize", pointer_to_DomainSize1)
    .def("DomainSize", pointer_to_DomainSize2)
    .def("CreateQuadTree", &QuadTreeSubCellType::CreateQuadTree)
    .def("ConstructQuadrature", pointer_to_ConstructQuadrature1)
    .def("ConstructQuadrature", pointer_to_ConstructQuadrature2)
    .def("ShallowAddToModelPart", &QuadTreeSubCell_AddToModelPart<QuadTreeSubCellType, true>) // only add the sub-cell
    .def("DeepAddToModelPart", &QuadTreeSubCell_AddToModelPart<QuadTreeSubCellType, false>) // add the sub-cell and all the quad-trees
    .def("ShallowAddToModelPart", &QuadTreeSubCell_AddToModelPart2<QuadTreeSubCellType, true>) // only add the sub-cell, the properties id starts from 1
    .def("DeepAddToModelPart", &QuadTreeSubCell_AddToModelPart2<QuadTreeSubCellType, false>) // add the sub-cell and all the quad-trees, the properties id starts from 1
    ;

    std::stringstream MomentFittedQuadTreeSubCellName;
    MomentFittedQuadTreeSubCellName << "MomentFittedQuadTreeSubCell";
    if(TFrameType == GLOBAL_REFERENCE)
        MomentFittedQuadTreeSubCellName << "Reference";
    else if(TFrameType == GLOBAL_CURRENT)
        MomentFittedQuadTreeSubCellName << "Current";
    if(TNsampling > 1)
        MomentFittedQuadTreeSubCellName << TNsampling;

    void(MomentFittedQuadTreeSubCellType::*pointer_to_ConstructSubCellsBasedOnEqualDistribution1)(const int&) = &MomentFittedQuadTreeSubCellType::ConstructSubCellsBasedOnEqualDistribution;
    void(MomentFittedQuadTreeSubCellType::*pointer_to_ConstructSubCellsBasedOnEqualDistribution2)(const int&) = &MomentFittedQuadTreeSubCellType::ConstructSubCellsBasedOnEqualDistribution;

    class_<MomentFittedQuadTreeSubCellType, typename MomentFittedQuadTreeSubCellType::Pointer, QuadTreeSubCellType, BaseMomentFittedQuadTreeSubCell>
    (m, MomentFittedQuadTreeSubCellName.str().c_str())
    .def(init<Element::Pointer>())
    .def("ConstructSubCellsBasedOnGaussQuadrature", &MomentFittedQuadTreeSubCellType::ConstructSubCellsBasedOnGaussQuadrature)
    .def("ConstructSubCellsBasedOnEqualDistribution", pointer_to_ConstructSubCellsBasedOnEqualDistribution1)
    .def("ConstructSubCellsBasedOnEqualDistribution", pointer_to_ConstructSubCellsBasedOnEqualDistribution2)
    .def("GetRepresentativeIntegrationOrder", &MomentFittedQuadTreeSubCellType::GetRepresentativeIntegrationOrder)
    .def("GetElement", &MomentFittedQuadTreeSubCellType::pGetElement)
    .def("SetFictitiousWeight", &MomentFittedQuadTreeSubCellType::SetFictitiousWeight)
    .def("GetPhysicalIntegrationPoints", &MomentFittedQuadTreeSubCell_GetPhysicalIntegrationPoints<MomentFittedQuadTreeSubCellType>)
    .def("GetFictitiousIntegrationPoints", &MomentFittedQuadTreeSubCell_GetFictitiousIntegrationPoints<MomentFittedQuadTreeSubCellType>)
    .def("FitAndCreateSubCellElements", &MomentFittedQuadTreeSubCell_FitAndCreateSubCellElements<MomentFittedQuadTreeSubCellType>)
    .def("CreateSubCellElements", &MomentFittedQuadTreeSubCell_CreateSubCellElements<MomentFittedQuadTreeSubCellType>)
    .def("CreateSubCellElements", &MomentFittedQuadTreeSubCell_CreateSubCellElements2<MomentFittedQuadTreeSubCellType>)
    .def("CreateFictitiousElement", &MomentFittedQuadTreeSubCellType::CreateFictitiousElement)
    ;
}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_APPLICATION_ADD_QUADTREE_TO_PYTHON_HPP_INCLUDED  defined
