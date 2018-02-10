// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//



// Project includes
#include "includes/element.h"
#include "custom_python/add_quadrature_utility_to_python.h"
#include "custom_utilities/quadrature_utility.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/moment_fitted_quad_tree_subcell.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

int QuadratureUtility_GetQuadratureType(QuadratureUtility& rDummy, const int& integration_method)
{
    return rDummy.GetQuadratureType(integration_method);
}

int QuadratureUtility_GetQuadratureOrder(QuadratureUtility& rDummy, const int& integration_method)
{
    return rDummy.GetQuadratureOrder(integration_method);
}

void FiniteCellApplication_AddQuadratureUtilityToPython()
{
    void(QuadratureUtility::*pointer_to_CreateConditionFromQuadraturePoint)(ModelPart&, boost::python::list&,
            const std::string&, const double&, const double&) const = &QuadratureUtility::PyCreateConditionFromQuadraturePoint;
    ModelPart::ConditionsContainerType(QuadratureUtility::*pointer_to_CreateConditionFromPoint)(ModelPart&,
            boost::python::list&, const std::string&) const = &QuadratureUtility::PyCreateConditionFromPoint;
    ModelPart::ConditionsContainerType(QuadratureUtility::*pointer_to_CreateConditionFromPoint2)(ModelPart&,
            boost::python::list&, const std::string&, Properties::Pointer) const = &QuadratureUtility::PyCreateConditionFromPoint;

    class_<QuadratureUtility, QuadratureUtility::Pointer, boost::noncopyable>
    ("QuadratureUtility", init<>())
    .def("GetDefaultIntegrationMethod", &QuadratureUtility::GetDefaultIntegrationMethod<Element>)
    .def("GetDefaultIntegrationMethod", &QuadratureUtility::GetDefaultIntegrationMethod<Condition>)
    .def("GetQuadratureType", QuadratureUtility_GetQuadratureType)
    .def("GetQuadratureOrder", QuadratureUtility_GetQuadratureOrder)
    .def("ScaleQuadrature", &QuadratureUtility::PyScaleQuadrature)
    .def("SaveQuadrature", &QuadratureUtility::PySaveQuadrature)
    .def("SaveQuadrature", &QuadratureUtility::PySaveQuadratureAdvanced)
    .def("SaveQuadratureSubCell", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<1> >)
    .def("SaveQuadratureSubCell2", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<2> >)
    .def("SaveQuadratureSubCell3", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<3> >)
    .def("SaveQuadratureSubCell4", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<4> >)
    .def("SaveQuadratureSubCell5", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<5> >)
    .def("SetQuadrature", &QuadratureUtility::PySetQuadrature)
    .def("CreateConditionFromQuadraturePoint", pointer_to_CreateConditionFromQuadraturePoint)
    .def("CreateConditionFromPoint", pointer_to_CreateConditionFromPoint)
    .def("CreateConditionFromPoint", pointer_to_CreateConditionFromPoint2)
    .def("CreatePoint", &QuadratureUtility::CreatePoint)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

