// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_FINITE_CELL_APPLICATION_ADD_CUSTOM_ALGEBRA_TO_PYTHON_H_INCLUDED )
#define  KRATOS_FINITE_CELL_APPLICATION_ADD_CUSTOM_ALGEBRA_TO_PYTHON_H_INCLUDED


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

void  FiniteCellApplication_AddCustomAlgebraToPython();

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_APPLICATION_ADD_CUSTOM_ALGEBRA_TO_PYTHON_H_INCLUDED  defined 
