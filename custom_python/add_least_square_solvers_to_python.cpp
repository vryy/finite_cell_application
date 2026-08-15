// see finite_cell_application/LICENSE.txt
//
//   Project Name:        KratosFiniteCellApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Aug 2026 $
//
//

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/element.h"
#include "custom_python/add_least_square_solvers_to_python.h"
#include "custom_least_square_solvers/least_square_lapack_solver.h"
#ifdef FINITE_CELL_APPLICATION_USE_NNLS
#include "custom_least_square_solvers/nnls_solver.h"
#endif

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void FiniteCellApplication_AddLeastSquareSolversToPython()
{

    class_<LeastSquareLapackSolver, LeastSquareLapackSolver::Pointer, boost::noncopyable>
    ("LeastSquareLapackSolver", init<>())
    ;

#ifdef FINITE_CELL_APPLICATION_USE_NNLS
    class_<NnlsSolver, NnlsSolver::Pointer, boost::noncopyable>
    ("NnlsSolver", init<>())
    .def("Solve", &NnlsSolver::Solve, args("A", "X", "B", "echo_level", "maxit"), "Solve the non-negative least square problem Ax=B")
    ;
#endif
}

} // namespace Python.

} // namespace Kratos.
