// see finite_cell_structural_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 7 Jan 2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_FINITE_CELL_APPLICATION_ADD_CUSTOM_CONDITIONS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_FINITE_CELL_APPLICATION_ADD_CUSTOM_CONDITIONS_TO_PYTHON_H_INCLUDED


// System includes

// External includes
#include <boost/python.hpp>
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void  FiniteCellApplication_AddCustomConditionsToPython();

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_APPLICATION_ADD_CUSTOM_CONDITIONS_TO_PYTHON_H_INCLUDED  defined 
