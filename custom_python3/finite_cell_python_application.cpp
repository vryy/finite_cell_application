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


// Project includes
#include "includes/define_python.h"
#include "finite_cell_application_variables.h"
#include "finite_cell_application.h"
#include "custom_python3/add_custom_conditions_to_python.h"
#include "custom_python3/add_quadtree_to_python.hpp"
#include "custom_python3/add_quadrature_utility_to_python.h"
#include "custom_python3/add_div_free_basis_utility_to_python.h"
#include "custom_python3/add_moment_fitting_utility_to_python.h"
// #include "custom_python3/add_ghost_penalty_utility_to_python.h"
#include "custom_python3/add_finite_cell_auxiliary_utility_to_python.h"
#include "custom_python3/add_finite_cell_mesh_utility_to_python.h"
#include "custom_python3/add_utility_to_python.h"

namespace Kratos
{

namespace Python
{

PYBIND11_MODULE(KratosFiniteCellApplication, m)
{

    namespace py = pybind11;

    py::class_<KratosFiniteCellApplication,
    KratosFiniteCellApplication::Pointer,
    KratosApplication >(m, "KratosFiniteCellApplication")
    .def(py::init<>())
    ;

    FiniteCellApplication_AddRefinableTreeToPython(m);
    FiniteCellApplication_AddFunctionIntegratorToPython(m);

    FiniteCellApplication_AddQuadTreeToPython<1, GLOBAL_REFERENCE>(m);
    FiniteCellApplication_AddQuadTreeToPython<2, GLOBAL_REFERENCE>(m);
    FiniteCellApplication_AddQuadTreeToPython<3, GLOBAL_REFERENCE>(m);
    FiniteCellApplication_AddQuadTreeToPython<4, GLOBAL_REFERENCE>(m);
    FiniteCellApplication_AddQuadTreeToPython<5, GLOBAL_REFERENCE>(m);
    FiniteCellApplication_AddQuadTreeToPython<6, GLOBAL_REFERENCE>(m);
    FiniteCellApplication_AddQuadTreeToPython<7, GLOBAL_REFERENCE>(m);
    FiniteCellApplication_AddQuadTreeToPython<8, GLOBAL_REFERENCE>(m);
    FiniteCellApplication_AddQuadTreeToPython<9, GLOBAL_REFERENCE>(m);
    FiniteCellApplication_AddQuadTreeToPython<10, GLOBAL_REFERENCE>(m);

    FiniteCellApplication_AddQuadTreeToPython<1, GLOBAL_CURRENT>(m);
    FiniteCellApplication_AddQuadTreeToPython<2, GLOBAL_CURRENT>(m);
    FiniteCellApplication_AddQuadTreeToPython<3, GLOBAL_CURRENT>(m);
    FiniteCellApplication_AddQuadTreeToPython<4, GLOBAL_CURRENT>(m);
    FiniteCellApplication_AddQuadTreeToPython<5, GLOBAL_CURRENT>(m);
    FiniteCellApplication_AddQuadTreeToPython<6, GLOBAL_CURRENT>(m);
    FiniteCellApplication_AddQuadTreeToPython<7, GLOBAL_CURRENT>(m);
    FiniteCellApplication_AddQuadTreeToPython<8, GLOBAL_CURRENT>(m);
    FiniteCellApplication_AddQuadTreeToPython<9, GLOBAL_CURRENT>(m);
    FiniteCellApplication_AddQuadTreeToPython<10, GLOBAL_CURRENT>(m);

    // FiniteCellApplication_AddBaseMomentFittedQuadTreeSubCellToPython(m);

    // FiniteCellApplication_AddQuadTreeSubCellToPython<1, GLOBAL_REFERENCE>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<2, GLOBAL_REFERENCE>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<3, GLOBAL_REFERENCE>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<4, GLOBAL_REFERENCE>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<5, GLOBAL_REFERENCE>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<6, GLOBAL_REFERENCE>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<7, GLOBAL_REFERENCE>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<8, GLOBAL_REFERENCE>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<9, GLOBAL_REFERENCE>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<10, GLOBAL_REFERENCE>(m);

    // FiniteCellApplication_AddQuadTreeSubCellToPython<1, GLOBAL_CURRENT>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<2, GLOBAL_CURRENT>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<3, GLOBAL_CURRENT>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<4, GLOBAL_CURRENT>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<5, GLOBAL_CURRENT>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<6, GLOBAL_CURRENT>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<7, GLOBAL_CURRENT>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<8, GLOBAL_CURRENT>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<9, GLOBAL_CURRENT>(m);
    // FiniteCellApplication_AddQuadTreeSubCellToPython<10, GLOBAL_CURRENT>(m);

    FiniteCellApplication_AddCustomConditionsToPython(m);

    FiniteCellApplication_AddQuadratureUtilityToPython(m);
    FiniteCellApplication_AddDivFreeBasisUtilityToPython(m);
    FiniteCellApplication_AddMomentFittingUtilityToPython(m);
    // FiniteCellApplication_AddGhostPenaltyUtilityToPython(m);
    FiniteCellApplication_AddFiniteCellAuxiliaryUtilityToPython(m);
    FiniteCellApplication_AddFiniteCellMeshUtilityToPython(m);
    FiniteCellApplication_AddUtilityToPython(m);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, SUBCELL_DOMAIN_SIZE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, SUBCELL_WEIGHTS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, PHYSICAL_INTEGRATION_POINT_THREED_STRESSES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, OTHER_NODE_ID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, OTHER_ID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NUMBER_OF_PHYSICAL_POINTS )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, PHYSICAL_INTEGRATION_POINT_LOCAL )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, PHYSICAL_INTEGRATION_POINT_GLOBAL )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, PHYSICAL_INTEGRATION_POINT_DISPLACEMENT )

}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
