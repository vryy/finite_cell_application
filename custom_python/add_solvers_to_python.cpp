// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 27 Nov 2018 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/element.h"
#include "custom_python/add_solvers_to_python.h"
#ifdef FINITE_CELL_APPLICATION_USE_SPECTRA
#include "custom_linear_solvers/spectra_eigenvalues_solver.h"
#endif


namespace Kratos
{

namespace Python
{

using namespace boost::python;

#ifdef FINITE_CELL_APPLICATION_USE_SPECTRA
boost::python::list SpectraEigenvaluesSolver_SolveLargestUnsym(SpectraEigenvaluesSolver& rDummy, CompressedMatrix& rA, const int& ne)
{
    boost::python::list values;

    std::vector<double> eigenvalues_real;
    std::vector<double> eigenvalues_imag;
    rDummy.SolveLargestUnsym(rA, ne, eigenvalues_real, eigenvalues_imag);

    for (std::size_t i = 0; i < eigenvalues_real.size(); ++i)
    {
        boost::python::list eval;
        eval.append(eigenvalues_real[i]);
        eval.append(eigenvalues_imag[i]);
        values.append(eval);
    }

    return values;
}

boost::python::list SpectraEigenvaluesSolver_SolveLargestSym(SpectraEigenvaluesSolver& rDummy, CompressedMatrix& rA, const int& ne)
{
    boost::python::list values;

    std::vector<double> eigenvalues;
    rDummy.SolveLargestSym(rA, ne, eigenvalues);

    for (std::size_t i = 0; i < eigenvalues.size(); ++i)
    {
        values.append(eigenvalues[i]);
    }

    return values;
}

boost::python::list SpectraEigenvaluesSolver_SolveSmallestSPD(SpectraEigenvaluesSolver& rDummy, CompressedMatrix& rA,
    SpectraEigenvaluesSolver::LinearSolverType::Pointer pLinearSolver, const int& ne)
{
    boost::python::list values;

    std::vector<double> eigenvalues;
    rDummy.SolveSmallestSPD(rA, pLinearSolver, ne, eigenvalues);

    for (std::size_t i = 0; i < eigenvalues.size(); ++i)
    {
        values.append(eigenvalues[i]);
    }

    return values;
}
#endif

void FiniteCellApplication_AddSolversToPython()
{

    #ifdef FINITE_CELL_APPLICATION_USE_SPECTRA
    class_<SpectraEigenvaluesSolver, SpectraEigenvaluesSolver::Pointer, boost::noncopyable>
    ("SpectraEigenvaluesSolver", init<>())
    .def("SolveLargestUnsym", &SpectraEigenvaluesSolver_SolveLargestUnsym)
    .def("SolveLargestSym", &SpectraEigenvaluesSolver_SolveLargestSym)
    .def("SolveSmallestSPD", &SpectraEigenvaluesSolver_SolveSmallestSPD)
    ;
    #endif

}

}  // namespace Python.

}  // namespace Kratos.

