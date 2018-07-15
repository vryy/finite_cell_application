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
#include "custom_python/add_utility_to_python.h"
#include "custom_utilities/immersed_boundary_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void FiniteCellApplication_AddUtilityToPython()
{

    class_<ImmersedBoundaryUtility, ImmersedBoundaryUtility::Pointer, boost::noncopyable>
    ("ImmersedBoundaryUtility", init<>())
    .def("InitializeBinning", &ImmersedBoundaryUtility::InitializeBinning)
    .def("GetValueOnPoint", &ImmersedBoundaryUtility::GetValueOnPoint<Variable<double> >)
    .def("GetValueOnPoint", &ImmersedBoundaryUtility::GetValueOnPoint<Variable<Vector> >)
    .def("GetValueOnPoint", &ImmersedBoundaryUtility::GetValueOnPoint<Variable<array_1d<double, 3> > >)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

