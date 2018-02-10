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
    FiniteCellMeshUtility::GenerateStructuredMesh2D(sampling_points, type, StartPoint, EndPoint, nsampling);

    int close_dir = 0; // open loop
    int activation_dir = 0;
    FiniteCellMeshUtility::CreateQuadElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);
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
    FiniteCellMeshUtility::GenerateStructuredMesh3D(sampling_points, type, StartPoint, EndPoint, nsampling);

    int close_dir = 0; // open loop
    int activation_dir = 0;
    FiniteCellMeshUtility::CreateHexElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, pProperties);
}

void FiniteCellApplication_AddFiniteCellMeshUtilityToPython()
{

    class_<FiniteCellMeshUtility, FiniteCellMeshUtility::Pointer, boost::noncopyable>
    ("FiniteCellMeshUtility", init<>())
    ("GenerateStructuredModelPart", &FiniteCellMeshUtility_GenerateStructuredModelPart2D)
    ("GenerateStructuredModelPart", &FiniteCellMeshUtility_GenerateStructuredModelPart3D)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

