// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2007-09-20 12:17:31 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_FINITE_CELL_APPLICATION_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_FINITE_CELL_APPLICATION_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED


// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"


namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  FiniteCellApplication_AddCustomUtilitiesToPython();

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_APPLICATION_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED  defined 
