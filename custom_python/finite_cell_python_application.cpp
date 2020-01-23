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
#include "custom_python/add_custom_conditions_to_python.h"
#include "custom_python/add_quadtree_to_python.hpp"
#include "custom_python/add_quadrature_utility_to_python.h"
#include "custom_python/add_div_free_basis_utility_to_python.h"
#include "custom_python/add_moment_fitting_utility_to_python.h"
#include "custom_python/add_ghost_penalty_utility_to_python.h"
#include "custom_python/add_finite_cell_auxiliary_utility_to_python.h"
#include "custom_python/add_finite_cell_mesh_utility_to_python.h"
#include "custom_python/add_utility_to_python.h"

namespace Kratos
{

namespace Python
{

    using namespace boost::python;
    BOOST_PYTHON_MODULE(KratosFiniteCellApplication)
    {

        class_<KratosFiniteCellApplication, KratosFiniteCellApplication::Pointer, bases<KratosApplication>, boost::noncopyable>
        ("KratosFiniteCellApplication");

        FiniteCellApplication_AddRefinableTreeToPython();
        FiniteCellApplication_AddFunctionIntegratorToPython();

        FiniteCellApplication_AddQuadTreeToPython<1, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeToPython<2, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeToPython<3, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeToPython<4, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeToPython<5, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeToPython<6, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeToPython<7, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeToPython<8, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeToPython<9, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeToPython<10, GLOBAL_REFERENCE>();

        FiniteCellApplication_AddQuadTreeToPython<1, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeToPython<2, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeToPython<3, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeToPython<4, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeToPython<5, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeToPython<6, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeToPython<7, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeToPython<8, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeToPython<9, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeToPython<10, GLOBAL_CURRENT>();

        FiniteCellApplication_AddBaseMomentFittedQuadTreeSubCellToPython();

        FiniteCellApplication_AddQuadTreeSubCellToPython<1, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<2, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<3, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<4, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<5, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<6, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<7, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<8, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<9, GLOBAL_REFERENCE>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<10, GLOBAL_REFERENCE>();

        FiniteCellApplication_AddQuadTreeSubCellToPython<1, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<2, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<3, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<4, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<5, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<6, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<7, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<8, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<9, GLOBAL_CURRENT>();
        FiniteCellApplication_AddQuadTreeSubCellToPython<10, GLOBAL_CURRENT>();

        FiniteCellApplication_AddCustomConditionsToPython();

        FiniteCellApplication_AddQuadratureUtilityToPython();
        FiniteCellApplication_AddDivFreeBasisUtilityToPython();
        FiniteCellApplication_AddMomentFittingUtilityToPython();
        FiniteCellApplication_AddGhostPenaltyUtilityToPython();
        FiniteCellApplication_AddFiniteCellAuxiliaryUtilityToPython();
        FiniteCellApplication_AddFiniteCellMeshUtilityToPython();
        FiniteCellApplication_AddUtilityToPython();

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( SUBCELL_DOMAIN_SIZE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( SUBCELL_WEIGHTS )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PHYSICAL_INTEGRATION_POINT_THREED_STRESSES )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( OTHER_NODE_ID )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( OTHER_ID )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( NUMBER_OF_PHYSICAL_POINTS )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( PHYSICAL_INTEGRATION_POINT_LOCAL )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( PHYSICAL_INTEGRATION_POINT_GLOBAL )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( PHYSICAL_INTEGRATION_POINT_DISPLACEMENT )

    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
