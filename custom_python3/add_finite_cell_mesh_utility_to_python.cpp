// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 10 Feb 2018 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes

// Project includes
#include "custom_python3/add_finite_cell_mesh_utility_to_python.h"
#include "custom_utilities/finite_cell_mesh_utility.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

pybind11::list FiniteCellMeshUtility_GenerateUniformSampling(FiniteCellMeshUtility& rDummy,
    const double& s_min, const double& s_max,
    const std::size_t& nsampling)
{
    std::vector<double> sampling;
    rDummy.GenerateSampling(sampling, s_min, s_max, nsampling);

    pybind11::list psampling;
    for (std::size_t i = 0; i < sampling.size(); ++i)
        psampling.append(sampling[i]);
    return psampling;
}

pybind11::list FiniteCellMeshUtility_GenerateSampling(FiniteCellMeshUtility& rDummy,
    const double& s_min, const double& s_max,
    const double& w1, const double& w2,
    const std::size_t& nsampling)
{
    std::vector<double> sampling;
    rDummy.GenerateSampling(sampling, s_min, s_max, w1, w2, nsampling);

    pybind11::list psampling;
    for (std::size_t i = 0; i < sampling.size(); ++i)
        psampling.append(sampling[i]);
    return psampling;
}

pybind11::dict FiniteCellMeshUtility_ExtractBoundaryLayer(const FiniteCellMeshUtility::BoundaryLayerInfoType& Info)
{
    pybind11::dict layer_cond_sets;

    for (FiniteCellMeshUtility::BoundaryLayerInfoType::const_iterator it = Info.begin(); it != Info.end(); ++it)
    {
        const std::string& layer_name = it->first;

        pybind11::list layer_conds;

        for (std::size_t i = 0; i < it->second.size(); ++i)
        {
            pybind11::list cond;
            for (std::size_t j = 0; j < it->second[i].size(); ++j)
                cond.append(it->second[i][j]);
            layer_conds.append(cond);
        }

        layer_cond_sets[layer_name.c_str()] = layer_conds;
    }

    return layer_cond_sets;
}

pybind11::dict FiniteCellMeshUtility_ExtractBoundaryNodes(const FiniteCellMeshUtility::BoundaryNodesInfoType& Info)
{
    pybind11::dict layer_node_sets;

    for (FiniteCellMeshUtility::BoundaryNodesInfoType::const_iterator it = Info.begin(); it != Info.end(); ++it)
    {
        const std::string& layer_name = it->first;

        pybind11::list node_set;
        for (std::set<std::size_t>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            node_set.append(*it2);
        }

        layer_node_sets[layer_name.c_str()] = node_set;
    }

    return layer_node_sets;
}

pybind11::list FiniteCellMeshUtility_GenerateStructuredModelPart2D(FiniteCellMeshUtility& rDummy,
    ModelPart& r_model_part,
    const Element::GeometryType::PointType::PointType& StartPoint,
    const Element::GeometryType::PointType::PointType& EndPoint,
    const int& nsampling1,
    const int& nsampling2,
    const std::string& sample_element_name,
    const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
    Properties::Pointer pProperties)
{
    typedef Element::GeometryType::PointType::PointType PointType;

    std::vector<std::vector<PointType> > sampling_points;

    std::vector<std::size_t> nsampling = {nsampling1, nsampling2};
    rDummy.GenerateStructuredPoints2D(sampling_points, type, StartPoint, EndPoint, nsampling);

    int close_dir = 0; // open loop
    int activation_dir = 0;
    FiniteCellMeshUtility::ElementMeshInfoType Info = rDummy.CreateQuadElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);

    typedef FiniteCellMeshUtility::BoundaryNodesInfoType BoundaryNodesInfoType;
    typedef FiniteCellMeshUtility::BoundaryLayerInfoType BoundaryLayerInfoType;
    const BoundaryNodesInfoType& boundary_nodes = std::get<2>(Info);
    const BoundaryLayerInfoType& boundary_layers = std::get<3>(Info);

    // generate layer information
    pybind11::dict layer_cond_sets = FiniteCellMeshUtility_ExtractBoundaryLayer(boundary_layers);
    pybind11::dict layer_node_sets = FiniteCellMeshUtility_ExtractBoundaryNodes(boundary_nodes);

    pybind11::list output;
    output.append(layer_node_sets);
    output.append(layer_cond_sets);
    return output;
}

pybind11::list FiniteCellMeshUtility_GenerateStructuredModelPart2DManualSampling(FiniteCellMeshUtility& rDummy,
    ModelPart& r_model_part,
    const Element::GeometryType::PointType::PointType& StartPoint,
    const Element::GeometryType::PointType::PointType& EndPoint,
    const pybind11::list& vsampling1,
    const pybind11::list& vsampling2,
    const std::string& sample_element_name,
    const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
    Properties::Pointer pProperties)
{
    typedef Element::GeometryType::PointType::PointType PointType;

    std::vector<std::vector<PointType> > sampling_points;

    std::vector<std::vector<double> > sampling(2);
    for (auto v : vsampling1)
        sampling[0].push_back(v.cast<double>());
    for (auto v : vsampling2)
        sampling[1].push_back(v.cast<double>());

    rDummy.GenerateStructuredPoints2D(sampling_points, type, StartPoint, EndPoint, sampling);

    int close_dir = 0; // open loop
    int activation_dir = 0;
    FiniteCellMeshUtility::ElementMeshInfoType Info = rDummy.CreateQuadElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);

    typedef FiniteCellMeshUtility::BoundaryNodesInfoType BoundaryNodesInfoType;
    typedef FiniteCellMeshUtility::BoundaryLayerInfoType BoundaryLayerInfoType;
    const BoundaryNodesInfoType& boundary_nodes = std::get<2>(Info);
    const BoundaryLayerInfoType& boundary_layers = std::get<3>(Info);

    // generate layer information
    pybind11::dict layer_cond_sets = FiniteCellMeshUtility_ExtractBoundaryLayer(boundary_layers);
    pybind11::dict layer_node_sets = FiniteCellMeshUtility_ExtractBoundaryNodes(boundary_nodes);

    pybind11::list output;
    output.append(layer_node_sets);
    output.append(layer_cond_sets);
    return output;
}

pybind11::list FiniteCellMeshUtility_GenerateStructuredModelPart2DInclined(FiniteCellMeshUtility& rDummy,
    ModelPart& r_model_part,
    const Element::GeometryType::PointType::PointType& StartPoint,
    const Element::GeometryType::PointType::PointType& Axis1,
    const Element::GeometryType::PointType::PointType& Axis2,
    const int& nsampling1,
    const int& nsampling2,
    const std::string& sample_element_name,
    const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
    Properties::Pointer pProperties)
{
    typedef Element::GeometryType::PointType::PointType PointType;

    std::vector<std::vector<PointType> > sampling_points;

    std::vector<std::size_t> nsampling = {nsampling1, nsampling2};
    std::vector<PointType> axis = {Axis1, Axis2};
    rDummy.GenerateStructuredPoints2D(sampling_points, type, StartPoint, axis, nsampling);

    int close_dir = 0; // open loop
    int activation_dir = 0;
    FiniteCellMeshUtility::ElementMeshInfoType Info = rDummy.CreateQuadElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);

    typedef FiniteCellMeshUtility::BoundaryNodesInfoType BoundaryNodesInfoType;
    typedef FiniteCellMeshUtility::BoundaryLayerInfoType BoundaryLayerInfoType;
    const BoundaryNodesInfoType& boundary_nodes = std::get<2>(Info);
    const BoundaryLayerInfoType& boundary_layers = std::get<3>(Info);

    // generate layer information
    pybind11::dict layer_cond_sets = FiniteCellMeshUtility_ExtractBoundaryLayer(boundary_layers);
    pybind11::dict layer_node_sets = FiniteCellMeshUtility_ExtractBoundaryNodes(boundary_nodes);

    pybind11::list output;
    output.append(layer_node_sets);
    output.append(layer_cond_sets);
    return output;
}

pybind11::list FiniteCellMeshUtility_GenerateStructuredModelPart3D(FiniteCellMeshUtility& rDummy,
    ModelPart& r_model_part,
    const Element::GeometryType::PointType::PointType& StartPoint,
    const Element::GeometryType::PointType::PointType& EndPoint,
    const int& nsampling1,
    const int& nsampling2,
    const int& nsampling3,
    const std::string& sample_element_name,
    const int& type, // if 1: generate H8 elements; 2: H20 elements; 3: H27 elements
    Properties::Pointer pProperties)
{
    typedef Element::GeometryType::PointType::PointType PointType;

    std::vector<std::vector<std::vector<PointType> > > sampling_points;

    std::vector<std::size_t> nsampling = {nsampling1, nsampling2, nsampling3};
    rDummy.GenerateStructuredPoints3D(sampling_points, type, StartPoint, EndPoint, nsampling);

    FiniteCellMeshUtility::ElementMeshInfoType Info = rDummy.CreateHexElements(r_model_part, sampling_points, sample_element_name, type, pProperties);

    typedef FiniteCellMeshUtility::BoundaryNodesInfoType BoundaryNodesInfoType;
    typedef FiniteCellMeshUtility::BoundaryLayerInfoType BoundaryLayerInfoType;
    const BoundaryNodesInfoType& boundary_nodes = std::get<2>(Info);
    const BoundaryLayerInfoType& boundary_layers = std::get<3>(Info);

    // generate layer information
    pybind11::dict layer_cond_sets = FiniteCellMeshUtility_ExtractBoundaryLayer(boundary_layers);
    pybind11::dict layer_node_sets = FiniteCellMeshUtility_ExtractBoundaryNodes(boundary_nodes);

    pybind11::list output;
    output.append(layer_node_sets);
    output.append(layer_cond_sets);
    return output;
}

pybind11::list FiniteCellMeshUtility_GenerateStructuredModelPart3DManualSampling(FiniteCellMeshUtility& rDummy,
    ModelPart& r_model_part,
    const Element::GeometryType::PointType::PointType& StartPoint,
    const Element::GeometryType::PointType::PointType& EndPoint,
    const pybind11::list& vsampling1,
    const pybind11::list& vsampling2,
    const pybind11::list& vsampling3,
    const std::string& sample_element_name,
    const int& type, // if 1: generate H8 elements; 2: H20 elements; 3: H27 elements
    Properties::Pointer pProperties)
{
    typedef Element::GeometryType::PointType::PointType PointType;

    std::vector<std::vector<std::vector<PointType> > > sampling_points;

    std::vector<std::vector<double> > sampling(3);
    for (auto v : vsampling1)
        sampling[0].push_back(v.cast<double>());
    for (auto v : vsampling2)
        sampling[1].push_back(v.cast<double>());
    for (auto v : vsampling3)
        sampling[2].push_back(v.cast<double>());

    rDummy.GenerateStructuredPoints3D(sampling_points, type, StartPoint, EndPoint, sampling);
//    KRATOS_WATCH(sampling_points.size())
//    KRATOS_WATCH(sampling_points[0].size())
//    KRATOS_WATCH(sampling_points[0][0].size())

    FiniteCellMeshUtility::ElementMeshInfoType Info = rDummy.CreateHexElements(r_model_part, sampling_points,
        sample_element_name, type, pProperties);

    typedef FiniteCellMeshUtility::BoundaryNodesInfoType BoundaryNodesInfoType;
    typedef FiniteCellMeshUtility::BoundaryLayerInfoType BoundaryLayerInfoType;
    const BoundaryNodesInfoType& boundary_nodes = std::get<2>(Info);
    const BoundaryLayerInfoType& boundary_layers = std::get<3>(Info);

    // generate layer information
    pybind11::dict layer_cond_sets = FiniteCellMeshUtility_ExtractBoundaryLayer(boundary_layers);
    pybind11::dict layer_node_sets = FiniteCellMeshUtility_ExtractBoundaryNodes(boundary_nodes);

    pybind11::list output;
    output.append(layer_node_sets);
    output.append(layer_cond_sets);
    return output;
}

Element::Pointer FiniteCellMeshUtility_CreateParasiteElement(FiniteCellMeshUtility& rDummy,
    const std::string& sample_element_name,
    std::size_t lastElementId, Element::Pointer pElement, Properties::Pointer pProperties)
{
    return rDummy.CreateParasiteElement(sample_element_name, lastElementId, pElement, pProperties);
}

Element::Pointer FiniteCellMeshUtility_CreateParasiteElement2(FiniteCellMeshUtility& rDummy,
    const std::string& sample_element_name,
    std::size_t lastElementId, Element::Pointer pElement, const int& integration_order,
    pybind11::list& assigned_quadrature, Properties::Pointer pProperties)
{
    Element::GeometryType::IntegrationPointsArrayType integration_points;
    for(std::size_t i = 0; i < pybind11::len(assigned_quadrature); ++i)
    {
        pybind11::list point = assigned_quadrature[i].cast<pybind11::list>();
        Element::GeometryType::IntegrationPointType integration_point;
        integration_point.X() = point[0].cast<double>();
        integration_point.Y() = point[1].cast<double>();
        integration_point.Z() = point[2].cast<double>();
        integration_point.Weight() = point[3].cast<double>();
//            KRATOS_WATCH(integration_point)
        integration_points.push_back(integration_point);
    }

    return rDummy.CreateParasiteElement(pElement, sample_element_name,
        integration_order, integration_points, lastElementId, pProperties);
}

ModelPart::NodesContainerType FiniteCellMeshUtility_ImportNodes(FiniteCellMeshUtility& rDummy,
    ModelPart& rThisModelPart, ModelPart& rOtherModelPart, const int& echo_level)
{
    return rDummy.ImportNodes(rThisModelPart, rOtherModelPart, echo_level);
}

ModelPart::NodesContainerType FiniteCellMeshUtility_ImportNodes2(FiniteCellMeshUtility& rDummy,
    ModelPart& rThisModelPart, ModelPart& rOtherModelPart,
    const double& offset_x, const double& offset_y, const double& offset_z,
    const double& cx, const double& cy, const double& theta, const int& echo_level)
{
    return rDummy.ImportNodes(rThisModelPart, rOtherModelPart, offset_x, offset_y, offset_z, cx, cy, theta, echo_level);
}

ModelPart::NodesContainerType FiniteCellMeshUtility_ImportNodes3(FiniteCellMeshUtility& rDummy,
    ModelPart& rThisModelPart, ModelPart& rOtherModelPart,
    const Transformation<double>& rTrans, const int& echo_level)
{
    return rDummy.ImportNodes(rThisModelPart, rOtherModelPart, rTrans, echo_level);
}

ModelPart::ElementsContainerType FiniteCellMeshUtility_ImportElements(FiniteCellMeshUtility& rDummy,
    ModelPart& rThisModelPart, ModelPart::ElementsContainerType& rOtherElements,
    const std::string& sample_element_name, Properties::Pointer pProperties, const int& echo_level)
{
    return rDummy.ImportElements(rThisModelPart, rOtherElements, sample_element_name, pProperties, echo_level);
}

ModelPart::ElementsContainerType FiniteCellMeshUtility_ImportElements2(FiniteCellMeshUtility& rDummy,
    ModelPart& rThisModelPart, ModelPart::ElementsContainerType& rOtherElements,
    const int& properties_offset_id, const int& echo_level)
{
    return rDummy.ImportElements(rThisModelPart, rOtherElements, properties_offset_id, echo_level);
}

ModelPart::ConditionsContainerType FiniteCellMeshUtility_ImportConditions(FiniteCellMeshUtility& rDummy,
    ModelPart& rThisModelPart, ModelPart::ConditionsContainerType& rOtherConditions,
    const std::string& sample_cond_name, Properties::Pointer pProperties, const int& echo_level)
{
    return rDummy.ImportConditions(rThisModelPart, rOtherConditions, sample_cond_name, pProperties, echo_level);
}

ModelPart::ConditionsContainerType FiniteCellMeshUtility_ImportConditions2(FiniteCellMeshUtility& rDummy,
    ModelPart& rThisModelPart, ModelPart::ConditionsContainerType& rOtherConditions,
    const int& properties_offset_id, const int& echo_level)
{
    return rDummy.ImportConditions(rThisModelPart, rOtherConditions, properties_offset_id, echo_level);
}

template<typename TEntityContainerType, typename TVariableType>
void FiniteCellMeshUtility_AssignValues(FiniteCellMeshUtility& rDummy,
    TEntityContainerType& rElements, TEntityContainerType& rOtherElements, const TVariableType& rVariable)
{
    rDummy.AssignValues(rElements, rOtherElements, rVariable);
}

void FiniteCellApplication_AddFiniteCellMeshUtilityToPython(pybind11::module& m)
{

    class_<FiniteCellMeshUtility, FiniteCellMeshUtility::Pointer>
    (m, "FiniteCellMeshUtility")
    .def(init<>())
    .def("GenerateSampling", &FiniteCellMeshUtility_GenerateUniformSampling)
    .def("GenerateSampling", &FiniteCellMeshUtility_GenerateSampling)
    .def("GenerateStructuredModelPart", &FiniteCellMeshUtility_GenerateStructuredModelPart2D)
    .def("GenerateStructuredModelPart", &FiniteCellMeshUtility_GenerateStructuredModelPart2DManualSampling)
    .def("GenerateStructuredModelPart", &FiniteCellMeshUtility_GenerateStructuredModelPart2DInclined)
    .def("GenerateStructuredModelPart", &FiniteCellMeshUtility_GenerateStructuredModelPart3D)
    .def("GenerateStructuredModelPart", &FiniteCellMeshUtility_GenerateStructuredModelPart3DManualSampling)
    .def("CreateParasiteElement", &FiniteCellMeshUtility_CreateParasiteElement)
    .def("CreateParasiteElement", &FiniteCellMeshUtility_CreateParasiteElement2)
    .def("ImportNodes", &FiniteCellMeshUtility_ImportNodes)
    .def("ImportNodes", &FiniteCellMeshUtility_ImportNodes2)
    .def("ImportNodes", &FiniteCellMeshUtility_ImportNodes3)
    .def("ImportElements", &FiniteCellMeshUtility_ImportElements)
    .def("ImportElements", &FiniteCellMeshUtility_ImportElements2)
    .def("ImportConditions", &FiniteCellMeshUtility_ImportConditions)
    .def("ImportConditions", &FiniteCellMeshUtility_ImportConditions2)
    .def("AssignValues", &FiniteCellMeshUtility_AssignValues<ModelPart::ElementsContainerType, Variable<int> >)
    .def("AssignValues", &FiniteCellMeshUtility_AssignValues<ModelPart::ElementsContainerType, Variable<bool> >)
    .def("AssignValues", &FiniteCellMeshUtility_AssignValues<ModelPart::ElementsContainerType, Variable<double> >)
    .def("AssignValues", &FiniteCellMeshUtility_AssignValues<ModelPart::ElementsContainerType, Variable<array_1d<double, 3> > >)
    .def("AssignValues", &FiniteCellMeshUtility_AssignValues<ModelPart::ConditionsContainerType, Variable<int> >)
    .def("AssignValues", &FiniteCellMeshUtility_AssignValues<ModelPart::ConditionsContainerType, Variable<bool> >)
    .def("AssignValues", &FiniteCellMeshUtility_AssignValues<ModelPart::ConditionsContainerType, Variable<double> >)
    .def("AssignValues", &FiniteCellMeshUtility_AssignValues<ModelPart::ConditionsContainerType, Variable<array_1d<double, 3> > >)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

