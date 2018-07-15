// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 3 Jul 2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_FINITE_CELL_APPLICATION_ADD_UTILITY_TO_PYTHON_H_INCLUDED )
#define  KRATOS_FINITE_CELL_APPLICATION_ADD_UTILITY_TO_PYTHON_H_INCLUDED


// System includes
#include <boost/python.hpp>

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void  FiniteCellApplication_AddUtilityToPython();

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_APPLICATION_ADD_UTILITY_TO_PYTHON_H_INCLUDED  defined
