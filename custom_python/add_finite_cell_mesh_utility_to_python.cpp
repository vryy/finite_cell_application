// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 10 Feb 2018 $
//   Revision:            $Revision: 1.0 $
//
//



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
        psampling.append(sampling[i]);
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
        psampling.append(sampling[i]);
    return psampling;
}

void FiniteCellMeshUtility_GenerateStructuredModelPart2D(FiniteCellMeshUtility& rDummy,
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
    rDummy.GenerateStructuredMesh2D(sampling_points, type, StartPoint, EndPoint, nsampling);

    int close_dir = 0; // open loop
    int activation_dir = 0;
    rDummy.CreateQuadElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);
}

void FiniteCellMeshUtility_GenerateStructuredModelPart2DManualSampling(FiniteCellMeshUtility& rDummy,
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
    BOOST_FOREACH(const iterator_value_type::value_type& v,
        std::make_pair(iterator_value_type(vsampling1), iterator_value_type() ) ) sampling[0].push_back(v);
    BOOST_FOREACH(const iterator_value_type::value_type& v,
        std::make_pair(iterator_value_type(vsampling2), iterator_value_type() ) ) sampling[1].push_back(v);

    rDummy.GenerateStructuredMesh2D(sampling_points, type, StartPoint, EndPoint, sampling);

    int close_dir = 0; // open loop
    int activation_dir = 0;
    rDummy.CreateQuadElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);
}

void FiniteCellMeshUtility_GenerateStructuredModelPart3D(FiniteCellMeshUtility& rDummy,
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
    rDummy.GenerateStructuredMesh3D(sampling_points, type, StartPoint, EndPoint, nsampling);

    int close_dir = 0; // open loop
    int activation_dir = 0;
    rDummy.CreateHexElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);
}

void FiniteCellMeshUtility_GenerateStructuredModelPart3DManualSampling(FiniteCellMeshUtility& rDummy,
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
    BOOST_FOREACH(const iterator_value_type::value_type& v,
        std::make_pair(iterator_value_type(vsampling1), iterator_value_type() ) ) sampling[0].push_back(v);
    BOOST_FOREACH(const iterator_value_type::value_type& v,
        std::make_pair(iterator_value_type(vsampling2), iterator_value_type() ) ) sampling[1].push_back(v);
    BOOST_FOREACH(const iterator_value_type::value_type& v,
        std::make_pair(iterator_value_type(vsampling3), iterator_value_type() ) ) sampling[2].push_back(v);

    rDummy.GenerateStructuredMesh3D(sampling_points, type, StartPoint, EndPoint, sampling);
    KRATOS_WATCH(sampling_points.size())
    KRATOS_WATCH(sampling_points[0].size())
    KRATOS_WATCH(sampling_points[0][0].size())

    int close_dir = 0; // open loop
    int activation_dir = 0;
    rDummy.CreateHexElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);
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
    for(std::size_t i = 0; i < boost::python::len(assigned_quadrature); ++i)
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
    ModelPart& rThisModelPart, ModelPart& rOtherModelPart)
{
    return rDummy.ImportNodes(rThisModelPart, rOtherModelPart);
}

ModelPart::ElementsContainerType FiniteCellMeshUtility_ImportElements(FiniteCellMeshUtility& rDummy,
    ModelPart& rThisModelPart, ModelPart::ElementsContainerType& rOtherElements,
    const std::string& sample_element_name, Properties::Pointer pProperties)
{
    return rDummy.ImportElements(rThisModelPart, rOtherElements, sample_element_name, pProperties);
}

ModelPart::ConditionsContainerType FiniteCellMeshUtility_ImportConditions(FiniteCellMeshUtility& rDummy,
    ModelPart& rThisModelPart, ModelPart::ConditionsContainerType& rOtherConditions,
    const std::string& sample_cond_name, Properties::Pointer pProperties)
{
    return rDummy.ImportConditions(rThisModelPart, rOtherConditions, sample_cond_name, pProperties);
}

void FiniteCellApplication_AddFiniteCellMeshUtilityToPython()
{

    class_<FiniteCellMeshUtility, FiniteCellMeshUtility::Pointer, boost::noncopyable>
    ("FiniteCellMeshUtility", init<>())
    .def("GenerateSampling", &FiniteCellMeshUtility_GenerateUniformSampling)
    .def("GenerateSampling", &FiniteCellMeshUtility_GenerateSampling)
    .def("GenerateStructuredModelPart", &FiniteCellMeshUtility_GenerateStructuredModelPart2D)
    .def("GenerateStructuredModelPart", &FiniteCellMeshUtility_GenerateStructuredModelPart2DManualSampling)
    .def("GenerateStructuredModelPart", &FiniteCellMeshUtility_GenerateStructuredModelPart3D)
    .def("GenerateStructuredModelPart", &FiniteCellMeshUtility_GenerateStructuredModelPart3DManualSampling)
    .def("CreateParasiteElement", &FiniteCellMeshUtility_CreateParasiteElement)
    .def("CreateParasiteElement", &FiniteCellMeshUtility_CreateParasiteElement2)
    .def("ImportNodes", &FiniteCellMeshUtility_ImportNodes)
    .def("ImportElements", &FiniteCellMeshUtility_ImportElements)
    .def("ImportConditions", &FiniteCellMeshUtility_ImportConditions)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

