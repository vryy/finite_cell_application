//  see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Feb 7, 2017 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes


// External includes
#if defined(KRATOS_PYTHON)
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "finite_cell_application.h"
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos
{

namespace Python
{

    using namespace boost::python;
    BOOST_PYTHON_MODULE(KratosFiniteCellApplication)
    {

        class_<KratosFiniteCellApplication, KratosFiniteCellApplication::Pointer, bases<KratosApplication>, boost::noncopyable>
        ("KratosFiniteCellApplication");

        FiniteCellApplication_AddCustomUtilitiesToPython();

    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
