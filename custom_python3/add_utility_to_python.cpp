// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes

// Project includes
#include "custom_python3/add_utility_to_python.h"
#include "custom_utilities/immersed_boundary_utility.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

void FiniteCellApplication_AddUtilityToPython(pybind11::module& m)
{

    class_<ImmersedBoundaryUtility, ImmersedBoundaryUtility::Pointer>
    (m, "ImmersedBoundaryUtility")
    .def(init<>())
    .def("InitializeBinning", &ImmersedBoundaryUtility::InitializeBinning)
    .def("GetValueOnPoint", &ImmersedBoundaryUtility::GetValueOnPoint<Variable<double> >)
    .def("GetValueOnPoint", &ImmersedBoundaryUtility::GetValueOnPoint<Variable<Vector> >)
    .def("GetValueOnPoint", &ImmersedBoundaryUtility::GetValueOnPoint<Variable<array_1d<double, 3> > >)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

