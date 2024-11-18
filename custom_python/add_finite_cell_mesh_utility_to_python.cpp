// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 10 Feb 2018 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes
#include <boost/python.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/element.h"
#include "custom_python/add_finite_cell_mesh_utility_to_python.h"
#include "custom_utilities/finite_cell_mesh_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

boost::python::list FiniteCellMeshUtility_GenerateUniformSampling(FiniteCellMeshUtility& rDummy,
        const double& s_min, const double& s_max,
        const std::size_t& nsampling)
{
    std::vector<double> sampling;
    rDummy.GenerateSampling(sampling, s_min, s_max, nsampling);

    boost::python::list psampling;
    for (std::size_t i = 0; i < sampling.size(); ++i)
    {
        psampling.append(sampling[i]);
    }
    return psampling;
}

boost::python::list FiniteCellMeshUtility_GenerateSampling(FiniteCellMeshUtility& rDummy,
        const double& s_min, const double& s_max,
        const double& w1, const double& w2,
        const std::size_t& nsampling)
{
    std::vector<double> sampling;
    rDummy.GenerateSampling(sampling, s_min, s_max, w1, w2, nsampling);

    boost::python::list psampling;
    for (std::size_t i = 0; i < sampling.size(); ++i)
    {
        psampling.append(sampling[i]);
    }
    return psampling;
}

boost::python::dict FiniteCellMeshUtility_ExtractBoundaryLayer(const FiniteCellMeshUtility::BoundaryLayerInfoType& Info)
{
    boost::python::dict layer_cond_sets;

    for (FiniteCellMeshUtility::BoundaryLayerInfoType::const_iterator it = Info.begin(); it != Info.end(); ++it)
    {
        const std::string& layer_name = it->first;

        boost::python::list layer_conds;

        for (std::size_t i = 0; i < it->second.size(); ++i)
        {
            boost::python::list cond;
            for (std::size_t j = 0; j < it->second[i].size(); ++j)
            {
                cond.append(it->second[i][j]);
            }
            layer_conds.append(cond);
        }

        layer_cond_sets[layer_name] = layer_conds;
    }

    return layer_cond_sets;
}

boost::python::dict FiniteCellMeshUtility_ExtractBoundaryNodes(const FiniteCellMeshUtility::BoundaryNodesInfoType& Info)
{
    boost::python::dict layer_node_sets;

    for (FiniteCellMeshUtility::BoundaryNodesInfoType::const_iterator it = Info.begin(); it != Info.end(); ++it)
    {
        const std::string& layer_name = it->first;

        boost::python::list node_set;
        for (std::set<std::size_t>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            node_set.append(*it2);
        }

        layer_node_sets[layer_name] = node_set;
    }

    return layer_node_sets;
}

boost::python::list FiniteCellMeshUtility_GenerateStructuredModelPart1D(FiniteCellMeshUtility& rDummy,
        ModelPart& r_model_part,
        const Element::GeometryType::PointType::PointType& StartPoint,
        const Element::GeometryType::PointType::PointType& EndPoint,
        const int& nsampling,
        const std::string& sample_element_name,
        const int& type, // if 1: generate L2 elements; 2, 3: L3 elements
        Properties::Pointer pProperties)
{
    typedef Element::GeometryType::PointType::PointType PointType;

    std::vector<PointType> sampling_points;

    int order;
    if (type == 1)
    {
        order = 1;
    }
    else if (type == 2 || type == 3)
    {
        order = 2;
    }
    else if (13 <= type && type <= 19)
    {
        order = type % 10;
    }
    else
        KRATOS_ERROR << "Invalid mesh type " << type;

    rDummy.GenerateSamplingPoints(sampling_points, StartPoint, EndPoint, nsampling * order);

    int close_dir = 0; // open loop
    int activation_dir = 0;
    FiniteCellMeshUtility::ElementMeshInfoType Info = rDummy.CreateLineElements(r_model_part, sampling_points, sample_element_name, type, close_dir, pProperties);

    typedef FiniteCellMeshUtility::BoundaryNodesInfoType BoundaryNodesInfoType;
    typedef FiniteCellMeshUtility::BoundaryLayerInfoType BoundaryLayerInfoType;
    const BoundaryNodesInfoType& boundary_nodes = std::get<2>(Info);
    const BoundaryLayerInfoType& boundary_layers = std::get<3>(Info);

    // generate layer information
    boost::python::dict layer_cond_sets = FiniteCellMeshUtility_ExtractBoundaryLayer(boundary_layers);
    boost::python::dict layer_node_sets = FiniteCellMeshUtility_ExtractBoundaryNodes(boundary_nodes);

    boost::python::list output;
    output.append(layer_node_sets);
    output.append(layer_cond_sets);
    return output;
}

boost::python::list FiniteCellMeshUtility_GenerateStructuredModelPart2D(FiniteCellMeshUtility& rDummy,
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

    std::vector<std::size_t> nsampling(2);
    nsampling[0] = static_cast<std::size_t>(nsampling1);
    nsampling[1] = static_cast<std::size_t>(nsampling2);
    rDummy.GenerateStructuredPoints2D(sampling_points, type, StartPoint, EndPoint, nsampling);

    int close_dir = 0; // open loop
    int activation_dir = 0;
    FiniteCellMeshUtility::ElementMeshInfoType Info = rDummy.CreateQuadElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);

    typedef FiniteCellMeshUtility::BoundaryNodesInfoType BoundaryNodesInfoType;
    typedef FiniteCellMeshUtility::BoundaryLayerInfoType BoundaryLayerInfoType;
    const BoundaryNodesInfoType& boundary_nodes = std::get<2>(Info);
    const BoundaryLayerInfoType& boundary_layers = std::get<3>(Info);

    // generate layer information
    boost::python::dict layer_cond_sets = FiniteCellMeshUtility_ExtractBoundaryLayer(boundary_layers);
    boost::python::dict layer_node_sets = FiniteCellMeshUtility_ExtractBoundaryNodes(boundary_nodes);

    boost::python::list output;
    output.append(layer_node_sets);
    output.append(layer_cond_sets);
    return output;
}

boost::python::list FiniteCellMeshUtility_GenerateStructuredModelPart2DManualSampling(FiniteCellMeshUtility& rDummy,
        ModelPart& r_model_part,
        const Element::GeometryType::PointType::PointType& StartPoint,
        const Element::GeometryType::PointType::PointType& EndPoint,
        const boost::python::list& vsampling1,
        const boost::python::list& vsampling2,
        const std::string& sample_element_name,
        const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
        Properties::Pointer pProperties)
{
    typedef Element::GeometryType::PointType::PointType PointType;

    std::vector<std::vector<PointType> > sampling_points;

    std::vector<std::vector<double> > sampling(2);
    typedef boost::python::stl_input_iterator<double> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type & v,
                  std::make_pair(iterator_value_type(vsampling1), iterator_value_type() ) ) sampling[0].push_back(v);
    BOOST_FOREACH(const iterator_value_type::value_type & v,
                  std::make_pair(iterator_value_type(vsampling2), iterator_value_type() ) ) sampling[1].push_back(v);

    rDummy.GenerateStructuredPoints2D(sampling_points, type, StartPoint, EndPoint, sampling);

    int close_dir = 0; // open loop
    int activation_dir = 0;
    FiniteCellMeshUtility::ElementMeshInfoType Info = rDummy.CreateQuadElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);

    typedef FiniteCellMeshUtility::BoundaryNodesInfoType BoundaryNodesInfoType;
    typedef FiniteCellMeshUtility::BoundaryLayerInfoType BoundaryLayerInfoType;
    const BoundaryNodesInfoType& boundary_nodes = std::get<2>(Info);
    const BoundaryLayerInfoType& boundary_layers = std::get<3>(Info);

    // generate layer information
    boost::python::dict layer_cond_sets = FiniteCellMeshUtility_ExtractBoundaryLayer(boundary_layers);
    boost::python::dict layer_node_sets = FiniteCellMeshUtility_ExtractBoundaryNodes(boundary_nodes);

    boost::python::list output;
    output.append(layer_node_sets);
    output.append(layer_cond_sets);
    return output;
}

boost::python::list FiniteCellMeshUtility_GenerateStructuredModelPart2DInclined(FiniteCellMeshUtility& rDummy,
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

    std::vector<std::size_t> nsampling(2);
    nsampling[0] = static_cast<std::size_t>(nsampling1);
    nsampling[1] = static_cast<std::size_t>(nsampling2);
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
    boost::python::dict layer_cond_sets = FiniteCellMeshUtility_ExtractBoundaryLayer(boundary_layers);
    boost::python::dict layer_node_sets = FiniteCellMeshUtility_ExtractBoundaryNodes(boundary_nodes);

    boost::python::list output;
    output.append(layer_node_sets);
    output.append(layer_cond_sets);
    return output;
}

boost::python::list FiniteCellMeshUtility_GenerateStructuredModelPart3D(FiniteCellMeshUtility& rDummy,
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

    std::vector<std::size_t> nsampling(3);
    nsampling[0] = static_cast<std::size_t>(nsampling1);
    nsampling[1] = static_cast<std::size_t>(nsampling2);
    nsampling[2] = static_cast<std::size_t>(nsampling3);
    rDummy.GenerateStructuredPoints3D(sampling_points, type, StartPoint, EndPoint, nsampling);

    FiniteCellMeshUtility::ElementMeshInfoType Info = rDummy.CreateHexElements(r_model_part, sampling_points, sample_element_name, type, pProperties);

    typedef FiniteCellMeshUtility::BoundaryNodesInfoType BoundaryNodesInfoType;
    typedef FiniteCellMeshUtility::BoundaryLayerInfoType BoundaryLayerInfoType;
    const BoundaryNodesInfoType& boundary_nodes = std::get<2>(Info);
    const BoundaryLayerInfoType& boundary_layers = std::get<3>(Info);

    // generate layer information
    boost::python::dict layer_cond_sets = FiniteCellMeshUtility_ExtractBoundaryLayer(boundary_layers);
    boost::python::dict layer_node_sets = FiniteCellMeshUtility_ExtractBoundaryNodes(boundary_nodes);

    boost::python::list output;
    output.append(layer_node_sets);
    output.append(layer_cond_sets);
    return output;
}

boost::python::list FiniteCellMeshUtility_GenerateStructuredModelPart3DManualSampling(FiniteCellMeshUtility& rDummy,
        ModelPart& r_model_part,
        const Element::GeometryType::PointType::PointType& StartPoint,
        const Element::GeometryType::PointType::PointType& EndPoint,
        const boost::python::list& vsampling1,
        const boost::python::list& vsampling2,
        const boost::python::list& vsampling3,
        const std::string& sample_element_name,
        const int& type, // if 1: generate H8 elements; 2: H20 elements; 3: H27 elements
        Properties::Pointer pProperties)
{
    typedef Element::GeometryType::PointType::PointType PointType;

    std::vector<std::vector<std::vector<PointType> > > sampling_points;

    std::vector<std::vector<double> > sampling(3);
    typedef boost::python::stl_input_iterator<double> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type & v,
                  std::make_pair(iterator_value_type(vsampling1), iterator_value_type() ) ) sampling[0].push_back(v);
    BOOST_FOREACH(const iterator_value_type::value_type & v,
                  std::make_pair(iterator_value_type(vsampling2), iterator_value_type() ) ) sampling[1].push_back(v);
    BOOST_FOREACH(const iterator_value_type::value_type & v,
                  std::make_pair(iterator_value_type(vsampling3), iterator_value_type() ) ) sampling[2].push_back(v);

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
    boost::python::dict layer_cond_sets = FiniteCellMeshUtility_ExtractBoundaryLayer(boundary_layers);
    boost::python::dict layer_node_sets = FiniteCellMeshUtility_ExtractBoundaryNodes(boundary_nodes);

    boost::python::list output;
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
        boost::python::list& assigned_quadrature, Properties::Pointer pProperties)
{
    Element::GeometryType::IntegrationPointsArrayType integration_points;
    for (std::size_t i = 0; i < boost::python::len(assigned_quadrature); ++i)
    {
        boost::python::list point = boost::python::extract<boost::python::list>(assigned_quadrature[i]);
        Element::GeometryType::IntegrationPointType integration_point;
        integration_point.X() = boost::python::extract<double>(point[0]);
        integration_point.Y() = boost::python::extract<double>(point[1]);
        integration_point.Z() = boost::python::extract<double>(point[2]);
        integration_point.Weight() = boost::python::extract<double>(point[3]);
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

void FiniteCellApplication_AddFiniteCellMeshUtilityToPython()
{

    class_<FiniteCellMeshUtility, FiniteCellMeshUtility::Pointer, boost::noncopyable>
    ("FiniteCellMeshUtility", init<>())
    .def("GenerateSampling", &FiniteCellMeshUtility_GenerateUniformSampling)
    .def("GenerateSampling", &FiniteCellMeshUtility_GenerateSampling)
    .def("GenerateStructuredModelPart", &FiniteCellMeshUtility_GenerateStructuredModelPart1D)
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

