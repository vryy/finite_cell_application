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
#include "custom_python/add_custom_algebra_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_quadtree_to_python.hpp"

namespace Kratos
{

namespace Python
{

    using namespace boost::python;
    BOOST_PYTHON_MODULE(KratosFiniteCellApplication)
    {

        class_<KratosFiniteCellApplication, KratosFiniteCellApplication::Pointer, bases<KratosApplication>, boost::noncopyable>
        ("KratosFiniteCellApplication");

        FiniteCellApplication_AddFunctionsToPython();
        FiniteCellApplication_AddBRepAndLevelSetToPython();
        FiniteCellApplication_AddRefinableTreeToPython();
        FiniteCellApplication_AddQuadTreeToPython<1>();
        FiniteCellApplication_AddQuadTreeToPython<2>();
        FiniteCellApplication_AddQuadTreeToPython<3>();
        FiniteCellApplication_AddQuadTreeToPython<4>();
        FiniteCellApplication_AddQuadTreeToPython<5>();
        FiniteCellApplication_AddCustomUtilitiesToPython();

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOAD_FUNCTION )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CUT_STATUS )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( SUBCELL_DOMAIN_SIZE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PHYSICAL_INTEGRATION_POINT_THREED_STRESSES )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( PHYSICAL_INTEGRATION_POINT_LOCAL )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( PHYSICAL_INTEGRATION_POINT_GLOBAL )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( PHYSICAL_INTEGRATION_POINT_DISPLACEMENT )

    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
