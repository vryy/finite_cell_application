// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 25 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes

// Project includes
#include "custom_python3/add_custom_conditions_to_python.h"
#include "custom_conditions/ghost_penalty_condition.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

void FiniteCellApplication_AddCustomConditionsToPython(pybind11::module& m)
{

    class_<GhostPenaltyCondition, GhostPenaltyCondition::Pointer, Condition>
    (m, "GhostPenaltyCondition")
    .def(init<>() )
    .def("__str__", &PrintObject<GhostPenaltyCondition>)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

