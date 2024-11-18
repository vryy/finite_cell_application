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
//  Date:            3 Feb 2018
//



// Project includes
#include "custom_utilities/finite_cell_mesh_utility.h"
#include "custom_utilities/finite_cell_geometry_utility.h"
#include "custom_utilities/finite_cell_auxiliary_utility.h"
#include "brep_application/custom_utilities/brep_math_utility.h"
#include "finite_cell_application_variables.h"


namespace Kratos
{

void FiniteCellMeshUtility::GenerateSampling(std::vector<double>& sampling,
        const double& s_min, const double& s_max, const std::size_t& nsampling)
{
    if (sampling.size() != nsampling + 1)
    {
        sampling.resize(nsampling + 1);
    }

    double ds = (s_max - s_min) / nsampling;
    for (std::size_t i = 0; i < nsampling + 1; ++i)
    {
        sampling[i] = s_min + i * ds;
    }
}


void FiniteCellMeshUtility::GenerateSampling(std::vector<double>& sampling,
        const double& s_min, const double& s_max,
        const double& w1, const double& w2,
        const std::size_t& nsampling)
{
    if (sampling.size() != nsampling + 1)
    {
        sampling.resize(nsampling + 1);
    }

    double t, dt = 1.0 / nsampling, N1, N2;
    for (std::size_t i = 0; i < nsampling + 1; ++i)
    {
        t = i * dt;
        N1 = 1.0 - t;
        N2 = t;
        sampling[i] = w1 * N1 / (w1 * N1 + w2 * N2) * s_min + w2 * N2 / (w1 * N1 + w2 * N2) * s_max;
    }
}


void FiniteCellMeshUtility::GenerateSamplingPoints(std::vector<PointType>& SamplingPoints,
        FiniteCellMeshUtility::GeometryType& r_geom, const std::size_t& nsampling)
{
    BRepMeshUtility::GenerateSamplingPoints<1>(SamplingPoints, r_geom, nsampling);
}


void FiniteCellMeshUtility::GenerateSamplingPoints(std::vector<PointType>& SamplingPoints,
        const PointType& StartPoint, const PointType& EndPoint, const std::size_t& nsampling)
{
    BRepMeshUtility::GenerateSamplingPoints(SamplingPoints, StartPoint, EndPoint, nsampling);
}


void FiniteCellMeshUtility::GenerateStructuredPoints2D(std::vector<std::vector<PointType> >& sampling_points,
        const int& type,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::size_t>& nsampling)
{
    if (type == 1)
    {
        GenerateStructuredPoints_Helper<2, 1>::Execute(sampling_points, StartPoint, EndPoint, nsampling);
    }
    else if (type == 2)
    {
        GenerateStructuredPoints_Helper<2, 2>::Execute(sampling_points, StartPoint, EndPoint, nsampling);
    }
    else if (type == 3)
    {
        GenerateStructuredPoints_Helper<2, 3>::Execute(sampling_points, StartPoint, EndPoint, nsampling);
    }
}


void FiniteCellMeshUtility::GenerateStructuredPoints2D(std::vector<std::vector<PointType> >& sampling_points,
        const int& type,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::vector<double> >& sampling)
{
    if (type == 1)
    {
        GenerateStructuredPoints_Helper<2, 1>::Execute(sampling_points, StartPoint, EndPoint, sampling);
    }
    else if (type == 2)
    {
        GenerateStructuredPoints_Helper<2, 2>::Execute(sampling_points, StartPoint, EndPoint, sampling);
    }
    else if (type == 3)
    {
        GenerateStructuredPoints_Helper<2, 3>::Execute(sampling_points, StartPoint, EndPoint, sampling);
    }
}


void FiniteCellMeshUtility::GenerateStructuredPoints2D(std::vector<std::vector<PointType> >& sampling_points,
        const int& type,
        const PointType& StartPoint,
        const std::vector<PointType>& Axis,
        const std::vector<std::size_t>& nsampling)
{
    if (type == 1)
    {
        GenerateStructuredPoints_Helper<2, 1>::Execute(sampling_points, StartPoint, Axis, nsampling);
    }
    else if (type == 2)
    {
        GenerateStructuredPoints_Helper<2, 2>::Execute(sampling_points, StartPoint, Axis, nsampling);
    }
    else if (type == 3)
    {
        GenerateStructuredPoints_Helper<2, 3>::Execute(sampling_points, StartPoint, Axis, nsampling);
    }
}


void FiniteCellMeshUtility::GenerateStructuredPoints3D(std::vector<std::vector<std::vector<PointType> > >& sampling_points,
        const int& type,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::size_t>& nsampling)
{
    if (type == 1)
    {
        GenerateStructuredPoints_Helper<3, 1>::Execute(sampling_points, StartPoint, EndPoint, nsampling);
    }
    else if (type == 2)
    {
        GenerateStructuredPoints_Helper<3, 2>::Execute(sampling_points, StartPoint, EndPoint, nsampling);
    }
    else if (type == 3)
    {
        GenerateStructuredPoints_Helper<3, 3>::Execute(sampling_points, StartPoint, EndPoint, nsampling);
    }
}


void FiniteCellMeshUtility::GenerateStructuredPoints3D(std::vector<std::vector<std::vector<PointType> > >& sampling_points,
        const int& type,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::vector<double> >& sampling)
{
    if (type == 1)
    {
        GenerateStructuredPoints_Helper<3, 1>::Execute(sampling_points, StartPoint, EndPoint, sampling);
    }
    else if (type == 2)
    {
        GenerateStructuredPoints_Helper<3, 2>::Execute(sampling_points, StartPoint, EndPoint, sampling);
    }
    else if (type == 3)
    {
        GenerateStructuredPoints_Helper<3, 3>::Execute(sampling_points, StartPoint, EndPoint, sampling);
    }
}


FiniteCellMeshUtility::ElementMeshInfoType FiniteCellMeshUtility::CreateLineElements(ModelPart& r_model_part,
        const std::vector<PointType>& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate L2 elements; 2: L3 elements;
        const bool& close, // if false: open loop; true: close loop
        Properties::Pointer pProperties)
{
    return BRepMeshUtility::CreateLineElements(r_model_part, sampling_points, sample_element_name, type, close, pProperties);
}


FiniteCellMeshUtility::ElementMeshInfoType FiniteCellMeshUtility::CreateQuadElements(ModelPart& r_model_part,
        const std::vector<std::vector<PointType> >& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
        const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
        const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
        Properties::Pointer pProperties)
{
    return BRepMeshUtility::CreateQuadElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);
}


FiniteCellMeshUtility::ElementMeshInfoType FiniteCellMeshUtility::CreateHexElements(ModelPart& r_model_part,
        const std::vector<std::vector<std::vector<PointType> > >& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate H8 elements; 2: H20 elements; 3: H27 elements
        Properties::Pointer pProperties)
{
    return BRepMeshUtility::CreateHexElements(r_model_part, sampling_points, sample_element_name, type, pProperties);
}


Element::Pointer FiniteCellMeshUtility::CreateParasiteElement(const std::string& sample_element_name,
        std::size_t& lastElementId,
        Element::Pointer pElement, Properties::Pointer pProperties )
{
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);

    // REMARK: when creating the element here, the integration rule is not passed. Instead the default integration rule of this element_type is applied, which is not the same as the original element.
    Element::Pointer pNewElement = r_clone_element.Create(++lastElementId, pElement->pGetGeometry(), pProperties);

    std::cout << "1 element of type " << sample_element_name << " is created" << std::endl;

    return pNewElement;
}


Element::Pointer FiniteCellMeshUtility::CreateParasiteElement(Element::Pointer pElement, // the parent element keep the geometry
        const std::string& sample_element_name,
        const int& RepresentativeIntegrationOrder,
        const GeometryType::IntegrationPointsArrayType& integration_points,
        std::size_t& lastElementId,
        Properties::Pointer pProperties)
{
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);
    const GeometryType& r_geom = *(pElement->pGetGeometry());

    GeometryData::IntegrationMethod RepresentativeIntegrationMethod
        = BRepMathUtility<>::GetIntegrationMethod(RepresentativeIntegrationOrder);
    Variable<int>& INTEGRATION_ORDER_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("INTEGRATION_ORDER"));

    // create the new elements from sub-cell
    // here we make a clone of the geometry because we want to assign different geometry data later on
    // this also works with Bezier element, because Bezier geometry has implemented the Create method
    Element::Pointer pNewElement;
    pNewElement = r_clone_element.Create(++lastElementId, r_geom.Create(r_geom.Points()), pProperties);

    FiniteCellGeometryUtility::AssignGeometryData(pNewElement->GetGeometry(), RepresentativeIntegrationMethod, integration_points);
    pNewElement->SetValue(INTEGRATION_ORDER_var, RepresentativeIntegrationOrder);

    return pNewElement;
}


ModelPart::NodesContainerType FiniteCellMeshUtility::ImportNodes(ModelPart& rThisModelPart, ModelPart& rOtherModelPart, const int& echo_level)
{
    std::size_t last_node_id = FiniteCellAuxiliaryUtility::GetLastNodeId(rThisModelPart);

    // create nodes and add to model_part
    ModelPart::NodesContainerType NewNodes;
    for (ModelPart::NodeIterator it = rOtherModelPart.NodesBegin(); it != rOtherModelPart.NodesEnd(); ++it)
    {
        std::size_t new_node_id = ++last_node_id;
        NodeType::Pointer pNewNode = rThisModelPart.CreateNewNode(new_node_id, it->X(), it->Y(), it->Z());
        it->SetValue(OTHER_NODE_ID, new_node_id);
        NewNodes.push_back(pNewNode);
    }

    if (echo_level > 0)
        std::cout << NewNodes.size() << " nodes from " << rOtherModelPart.Name()
                  << " are added to the model_part " << rThisModelPart.Name() << std::endl;

    return NewNodes;
}


ModelPart::NodesContainerType FiniteCellMeshUtility::ImportNodes(ModelPart& rThisModelPart, ModelPart& rOtherModelPart,
        const double& offset_x, const double& offset_y, const double& offset_z,
        const double& cx, const double& cy, const double& theta, const int& echo_level)
{
    std::size_t last_node_id = FiniteCellAuxiliaryUtility::GetLastNodeId(rThisModelPart);

    // create nodes and add to model_part
    ModelPart::NodesContainerType NewNodes;
    double x, y, z, xnew, ynew;
    double c = std::cos(theta);
    double s = std::sin(theta);
    for (ModelPart::NodeIterator it = rOtherModelPart.NodesBegin(); it != rOtherModelPart.NodesEnd(); ++it)
    {
        std::size_t new_node_id = ++last_node_id;

        x = it->X() - cx;
        y = it->Y() - cy;

        xnew = x * c - y * s;
        ynew = x * s + y * c;

        x = xnew + cx + offset_x;
        y = ynew + cy + offset_y;
        z = it->Z() + offset_z;

        NodeType::Pointer pNewNode = rThisModelPart.CreateNewNode(new_node_id, x, y, z);
        it->SetValue(OTHER_NODE_ID, new_node_id);
        NewNodes.push_back(pNewNode);
    }

    if (echo_level > 0)
        std::cout << NewNodes.size() << " nodes from " << rOtherModelPart.Name()
                  << " are added to the model_part " << rThisModelPart.Name() << std::endl;

    return NewNodes;
}


ModelPart::NodesContainerType FiniteCellMeshUtility::ImportNodes(ModelPart& rThisModelPart, ModelPart& rOtherModelPart,
        const Transformation<double>& rTrans, const int& echo_level)
{
    std::size_t last_node_id = FiniteCellAuxiliaryUtility::GetLastNodeId(rThisModelPart);

    // create nodes and add to model_part
    ModelPart::NodesContainerType NewNodes;
    array_1d<double, 3> point;
    for (ModelPart::NodeIterator it = rOtherModelPart.NodesBegin(); it != rOtherModelPart.NodesEnd(); ++it)
    {
        std::size_t new_node_id = ++last_node_id;

        point[0] = it->X();
        point[1] = it->Y();
        point[2] = it->Z();
        rTrans.ApplyTransformation(point);

        NodeType::Pointer pNewNode = rThisModelPart.CreateNewNode(new_node_id, point[0], point[1], point[2]);
        it->SetValue(OTHER_NODE_ID, new_node_id);
        NewNodes.push_back(pNewNode);
    }

    if (echo_level > 0)
        std::cout << NewNodes.size() << " nodes from " << rOtherModelPart.Name()
                  << " are added to the model_part " << rThisModelPart.Name() << std::endl;

    return NewNodes;
}

ModelPart::ElementsContainerType FiniteCellMeshUtility::ImportElements(ModelPart& rThisModelPart,
        ModelPart::ElementsContainerType& rOtherElements,
        const std::string& sample_element_name, Properties::Pointer pProperties, const int& echo_level)
{
    std::size_t last_element_id = FiniteCellAuxiliaryUtility::GetLastElementId(rThisModelPart);
    if (!KratosComponents<Element>::Has(sample_element_name))
        KRATOS_ERROR << sample_element_name << " is not registerred in Kratos database";
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);

    ModelPart::ElementsContainerType NewElements;
    FiniteCellMeshUtility_Helper<Element, ModelPart::ElementsContainerType>::ImportEntities(rThisModelPart,
            rThisModelPart.Elements(), NewElements, rOtherElements,
            last_element_id, r_clone_element, pProperties, echo_level);

    if (echo_level > 0)
        std::cout << NewElements.size() << " " << sample_element_name
                  << " elements are created and added to the model_part " << rThisModelPart.Name() << std::endl;

    return NewElements;
}

ModelPart::ElementsContainerType FiniteCellMeshUtility::ImportElements(ModelPart& rThisModelPart,
        ModelPart::ElementsContainerType& rOtherElements,
        const int& properties_id_offset, const int& echo_level)
{
    std::size_t last_element_id = FiniteCellAuxiliaryUtility::GetLastElementId(rThisModelPart);

    ModelPart::ElementsContainerType NewElements;
    FiniteCellMeshUtility_Helper<Element, ModelPart::ElementsContainerType>::ImportEntities(rThisModelPart,
            rThisModelPart.Elements(), NewElements, rOtherElements,
            last_element_id, properties_id_offset, echo_level);

    if (echo_level > 0)
        std::cout << NewElements.size()
                  << " elements are imported to the model_part " << rThisModelPart.Name() << std::endl;

    return NewElements;
}

ModelPart::ConditionsContainerType FiniteCellMeshUtility::ImportConditions(ModelPart& rThisModelPart,
        ModelPart::ConditionsContainerType& rOtherConditions,
        const std::string& sample_cond_name, Properties::Pointer pProperties, const int& echo_level)
{
    std::size_t last_cond_id = FiniteCellAuxiliaryUtility::GetLastConditionId(rThisModelPart);
    if (!KratosComponents<Condition>::Has(sample_cond_name))
        KRATOS_ERROR << sample_cond_name << " is not registerred in Kratos database";
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_cond_name);

    ModelPart::ConditionsContainerType NewConditions;
    FiniteCellMeshUtility_Helper<Condition, ModelPart::ConditionsContainerType>::ImportEntities(rThisModelPart,
            rThisModelPart.Conditions(), NewConditions, rOtherConditions,
            last_cond_id, r_clone_condition, pProperties, echo_level);

    if (echo_level > 0)
        std::cout << NewConditions.size() << " " << sample_cond_name
                  << " conditions are created and added to the model_part " << rThisModelPart.Name() << std::endl;

    return NewConditions;
}

ModelPart::ConditionsContainerType FiniteCellMeshUtility::ImportConditions(ModelPart& rThisModelPart,
        ModelPart::ConditionsContainerType& rOtherConditions,
        const int& properties_id_offset, const int& echo_level)
{
    std::size_t last_cond_id = FiniteCellAuxiliaryUtility::GetLastConditionId(rThisModelPart);

    ModelPart::ConditionsContainerType NewConditions;
    FiniteCellMeshUtility_Helper<Condition, ModelPart::ConditionsContainerType>::ImportEntities(rThisModelPart,
            rThisModelPart.Conditions(), NewConditions, rOtherConditions,
            last_cond_id, properties_id_offset, echo_level);

    if (echo_level > 0)
        std::cout << NewConditions.size()
                  << " conditions are imported to the model_part " << rThisModelPart.Name() << std::endl;

    return NewConditions;
}

template<class TEntityType, class TEntityContainerType>
void FiniteCellMeshUtility_Helper<TEntityType, TEntityContainerType>::ImportEntities(ModelPart& rThisModelPart,
        TEntityContainerType& rThisElements, // rThisModelPart.Elements() or rThisModelPart.Conditions()
        TEntityContainerType& rNewElements, // the added elements to rThisElements
        TEntityContainerType& rOtherElements, // other elements to be imported from
        std::size_t& last_element_id,
        TEntityType const& r_clone_element,
        Properties::Pointer pProperties,
        const int& echo_level)
{
    typedef typename TEntityType::GeometryType GeometryType;

    typename TEntityType::NodesArrayType temp_element_nodes;
    for (typename TEntityContainerType::ptr_iterator it = rOtherElements.ptr_begin(); it != rOtherElements.ptr_end(); ++it)
    {
        GeometryType& r_geom = (*it)->GetGeometry();

        temp_element_nodes.clear();

        for (std::size_t i = 0; i < r_geom.size(); ++i)
        {
            std::size_t other_node_id = static_cast<std::size_t>(r_geom[i].GetValue(OTHER_NODE_ID));
            temp_element_nodes.push_back(*(BRepUtility::FindKey(rThisModelPart.Nodes(), other_node_id, "Node").base()));
        }

        typename GeometryType::Pointer pNewGeometry = r_geom.Create(temp_element_nodes);

        typename TEntityType::Pointer pNewElement;
        pNewElement = r_clone_element.Create(++last_element_id, pNewGeometry, pProperties);
        pNewElement->SetValue(OTHER_ID, (*it)->Id());
        rNewElements.push_back(pNewElement);
    }

    for (typename TEntityContainerType::ptr_iterator it = rNewElements.ptr_begin(); it != rNewElements.ptr_end(); ++it)
    {
        rThisElements.push_back(*it);
    }

    rThisElements.Unique();
}

template<class TEntityType, class TEntityContainerType>
void FiniteCellMeshUtility_Helper<TEntityType, TEntityContainerType>::ImportEntities(ModelPart& rThisModelPart,
        TEntityContainerType& rThisElements, // rThisModelPart.Elements() or rThisModelPart.Conditions()
        TEntityContainerType& rNewElements, // the added elements to rThisElements
        TEntityContainerType& rOtherElements, // other elements to be imported from
        std::size_t& last_element_id,
        const int& properties_id_offset,
        const int& echo_level)
{
    typedef typename TEntityType::GeometryType GeometryType;

    typename TEntityType::NodesArrayType temp_element_nodes;
    std::size_t old_prop_id, new_prop_id;
    for (typename TEntityContainerType::ptr_iterator it = rOtherElements.ptr_begin(); it != rOtherElements.ptr_end(); ++it)
    {
        GeometryType& r_geom = (*it)->GetGeometry();

        temp_element_nodes.clear();

        for (std::size_t i = 0; i < r_geom.size(); ++i)
        {
            std::size_t other_node_id = static_cast<std::size_t>(r_geom[i].GetValue(OTHER_NODE_ID));
            temp_element_nodes.push_back(*(BRepUtility::FindKey(rThisModelPart.Nodes(), other_node_id, "Node").base()));
        }

        old_prop_id = (*it)->GetProperties().Id();
        new_prop_id = old_prop_id + properties_id_offset;

        typename TEntityType::Pointer pNewElement;
        Properties::Pointer pProperties = rThisModelPart.rProperties()(new_prop_id);
        pNewElement = (*it)->Create(++last_element_id, temp_element_nodes, pProperties);
        pNewElement->SetValue(OTHER_ID, (*it)->Id());
        rNewElements.push_back(pNewElement);
    }

    for (typename TEntityContainerType::ptr_iterator it = rNewElements.ptr_begin(); it != rNewElements.ptr_end(); ++it)
    {
        rThisElements.push_back(*it);
    }

    rThisElements.Unique();
}

template struct FiniteCellMeshUtility_Helper<Element, ModelPart::ElementsContainerType>;
template struct FiniteCellMeshUtility_Helper<Condition, ModelPart::ConditionsContainerType>;

}  // namespace Kratos.

