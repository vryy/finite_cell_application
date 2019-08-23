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


namespace Kratos
{

void FiniteCellMeshUtility::GenerateSampling(std::vector<double>& sampling,
    const double& s_min, const double& s_max, const std::size_t& nsampling)
{
    if(sampling.size() != nsampling + 1)
        sampling.resize(nsampling + 1);

    double ds = (s_max - s_min) / nsampling;
    for (std::size_t i = 0; i < nsampling + 1; ++i)
        sampling[i] = s_min + i*ds;
}


void FiniteCellMeshUtility::GenerateSampling(std::vector<double>& sampling,
        const double& s_min, const double& s_max,
        const double& w1, const double& w2,
        const std::size_t& nsampling)
{
    if(sampling.size() != nsampling + 1)
        sampling.resize(nsampling + 1);

    double t, dt = 1.0 / nsampling, N1, N2;
    for (std::size_t i = 0; i < nsampling + 1; ++i)
    {
        t = i*dt;
        N1 = 1.0 - t;
        N2 = t;
        sampling[i] = w1*N1/(w1*N1+w2*N2)*s_min + w2*N2/(w1*N1+w2*N2)*s_max;
    }
}


void FiniteCellMeshUtility::GenerateSamplingPoints(std::vector<PointType>& SamplingPoints,
        FiniteCellMeshUtility::GeometryType& r_geom, const std::size_t& nsampling)
{
    BRepMeshUtility::GenerateSamplingPoints(SamplingPoints, r_geom, nsampling);
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


/// Generate the points for background structure mesh
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


FiniteCellMeshUtility::MeshInfoType FiniteCellMeshUtility::CreateLineElements(ModelPart& r_model_part,
        const std::vector<PointType>& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate L2 elements; 2: L3 elements;
        const bool& close, // if false: open loop; true: close loop
        Properties::Pointer pProperties)
{
    return BRepMeshUtility::CreateLineElements(r_model_part, sampling_points, sample_element_name, type, close, pProperties);
}


FiniteCellMeshUtility::MeshInfoType FiniteCellMeshUtility::CreateQuadElements(ModelPart& r_model_part,
    const std::vector<std::vector<PointType> >& sampling_points,
    const std::string& sample_element_name,
    const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
    const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
    const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
    Properties::Pointer pProperties)
{
    return BRepMeshUtility::CreateQuadElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);
}


FiniteCellMeshUtility::MeshInfoType FiniteCellMeshUtility::CreateHexElements(ModelPart& r_model_part,
    const std::vector<std::vector<std::vector<PointType> > >& sampling_points,
    const std::string& sample_element_name,
    const int& type, // if 1: generate H8 elements; 2: H20 elements; 3: H27 elements
    const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir; 3: close on 3rd dir
    const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir; r: activation on 3rd dir
    Properties::Pointer pProperties)
{
    return BRepMeshUtility::CreateHexElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);
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
    GeometryType& r_geom = *(pElement->pGetGeometry());

    GeometryData::IntegrationMethod RepresentativeIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(RepresentativeIntegrationOrder);
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


ModelPart::NodesContainerType FiniteCellMeshUtility::ImportNodes(ModelPart& rThisModelPart, ModelPart& rOtherModelPart)
{
    std::size_t last_node_id = FiniteCellAuxiliaryUtility::GetLastNodeId(rThisModelPart);

    // create nodes and add to model_part
    ModelPart::NodesContainerType NewNodes;
    for(ModelPart::NodeIterator it = rOtherModelPart.NodesBegin(); it != rOtherModelPart.NodesEnd(); ++it)
    {
        std::size_t new_node_id = ++last_node_id;
        NodeType::Pointer pNewNode = rThisModelPart.CreateNewNode(new_node_id, it->X(), it->Y(), it->Z());
        it->SetValue(OTHER_NODE_ID, new_node_id);
        NewNodes.push_back(pNewNode);
    }

    std::cout << NewNodes.size() << " nodes from " << rOtherModelPart.Name()
              << " are added to the model_part " << rThisModelPart.Name() << std::endl;
    return NewNodes;
}


ModelPart::NodesContainerType FiniteCellMeshUtility::ImportNodes(ModelPart& rThisModelPart, ModelPart& rOtherModelPart,
        const double& offset_x, const double& offset_y, const double& offset_z,
        const double& cx, const double& cy, const double& theta)
{
    std::size_t last_node_id = FiniteCellAuxiliaryUtility::GetLastNodeId(rThisModelPart);

    // create nodes and add to model_part
    ModelPart::NodesContainerType NewNodes;
    double x, y, z, xnew, ynew;
    double c = std::cos(theta);
    double s = std::sin(theta);
    for(ModelPart::NodeIterator it = rOtherModelPart.NodesBegin(); it != rOtherModelPart.NodesEnd(); ++it)
    {
        std::size_t new_node_id = ++last_node_id;

        x = it->X() - cx;
        y = it->Y() - cy;

        xnew = x*c - y*s;
        ynew = x*s + y*c;

        x = xnew + cx + offset_x;
        y = ynew + cy + offset_y;
        z = it->Z() + offset_z;

        NodeType::Pointer pNewNode = rThisModelPart.CreateNewNode(new_node_id, x, y, z);
        it->SetValue(OTHER_NODE_ID, new_node_id);
        NewNodes.push_back(pNewNode);
    }

    std::cout << NewNodes.size() << " nodes from " << rOtherModelPart.Name()
              << " are added to the model_part " << rThisModelPart.Name() << std::endl;
    return NewNodes;
}


ModelPart::ElementsContainerType FiniteCellMeshUtility::ImportElements(ModelPart& rThisModelPart,
    ModelPart::ElementsContainerType& rOtherElements,
    const std::string& sample_element_name, Properties::Pointer pProperties)
{
    std::size_t last_element_id = FiniteCellAuxiliaryUtility::GetLastElementId(rThisModelPart);
    if (!KratosComponents<Element>::Has(sample_element_name))
        KRATOS_THROW_ERROR(std::logic_error, sample_element_name, "is not registerred in Kratos database")
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);

    ModelPart::ElementsContainerType NewElements;
    ImportEntities(rThisModelPart, rThisModelPart.Elements(), NewElements, rOtherElements, last_element_id, r_clone_element, pProperties);

    std::cout << NewElements.size() << " " << sample_element_name
              << " elements are created and added to the model_part " << rThisModelPart.Name() << std::endl;

    return NewElements;
}


ModelPart::ConditionsContainerType FiniteCellMeshUtility::ImportConditions(ModelPart& rThisModelPart,
    ModelPart::ConditionsContainerType& rOtherConditions,
    const std::string& sample_cond_name, Properties::Pointer pProperties)
{
    std::size_t last_cond_id = FiniteCellAuxiliaryUtility::GetLastConditionId(rThisModelPart);
    if (!KratosComponents<Condition>::Has(sample_cond_name))
        KRATOS_THROW_ERROR(std::logic_error, sample_cond_name, "is not registerred in Kratos database")
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_cond_name);

    ModelPart::ConditionsContainerType NewConditions;
    ImportEntities(rThisModelPart, rThisModelPart.Conditions(), NewConditions, rOtherConditions, last_cond_id, r_clone_condition, pProperties);

    std::cout << NewConditions.size() << " " << sample_cond_name
              << " conditions are created and added to the model_part " << rThisModelPart.Name() << std::endl;

    return NewConditions;
}

}  // namespace Kratos.

