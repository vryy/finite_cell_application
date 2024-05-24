//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         finite_cell_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            27 Jul 2018
//


#if !defined(KRATOS_NNLS_SOLVER_H_INCLUDED )
#define  KRATOS_NNLS_SOLVER_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes
#ifdef ENABLE_FINITE_CELL_BOOST_BINDINGS
#include "mkl_solvers_application/external_includes/boost/numeric/bindings/traits/matrix_traits.hpp"
#include "mkl_solvers_application/external_includes/boost/numeric/bindings/traits/vector_traits.hpp"
#include "mkl_solvers_application/external_includes/boost/numeric/bindings/traits/ublas_matrix.hpp"
#include "mkl_solvers_application/external_includes/boost/numeric/bindings/traits/ublas_vector.hpp"
#endif
#include "nnls.h"


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"



namespace Kratos
{
///@addtogroup FiniteCellApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Nonnegative least square solver based on NNLS
 * REF: http://suvrit.de/work/soft/nnls.html
 *      D. Kim, S. Sra, I. S. Dhillon. "A non-monotonic method for large-scale non-negative least squares." Optimization Methods and Software, Jan. 2012
*/
class NNLSSolver
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NNLSSolver
    KRATOS_CLASS_POINTER_DEFINITION(NNLSSolver);


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NNLSSolver() {}

    /// Destructor.
    virtual ~NNLSSolver() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static int Solve(Matrix& rA, Vector& rX, Vector& rB, int echo_level, int maxit = 100)
    {
        /* initialization */
        std::size_t M = rA.size1();
        std::size_t N = rA.size2();

#ifdef ENABLE_FINITE_CELL_BOOST_BINDINGS
        typedef boost::numeric::bindings::traits::matrix_traits<Matrix> matraits;
        double* dA = matraits::storage(rA);
#else
        double* dA = new double[M * N];
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
            {
                dA[i * N + j] = rA(j, i);    // TODO check
            }
#endif

        nsNNLS::matrix* A;
        nsNNLS::vector* b;
        nsNNLS::vector* x;
        nsNNLS::nnls*   solver;
        int     flag;

        A = new nsNNLS::denseMatrix(M, N, dA);

        b = new nsNNLS::vector(M);
        for (std::size_t i = 0; i < M; ++i)
        {
            b->set(i, rB(i));
        }

        solver = new nsNNLS::nnls(A, b, maxit);

        /* solve - optimize */
        if (echo_level > 0) { std::cout << "Optimizing...\n"; }
        flag = solver->optimize();
        if (echo_level > 0) { std::cout << "Done!\n"; }

        if (echo_level > 1) { std::cout << "Optimization time: " << solver->getOptimizationTime() << " seconds\n"; }

        if (flag < 0)
        {
            KRATOS_THROW_ERROR(std::logic_error, "NNLS: Solver terminated with error flag:", flag);
            return flag;
        }

        x = solver->getSolution();

        if (echo_level > 2) { solver->saveStats(std::cout); }

        if (rX.size() != N)
        {
            rX.resize(N, false);
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            rX(i) = x->get(i);
        }

        delete solver;
        delete A, b, x;
#ifndef ENABLE_FINITE_CELL_BOOST_BINDINGS
        delete dA;
#endif

        return 0;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Nonnegative Least Square Solver based on NNLS";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    NNLSSolver& operator=(NNLSSolver const& rOther);

    /// Copy constructor.
    NNLSSolver(NNLSSolver const& rOther);


    ///@}

}; // Class NNLSSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  NNLSSolver& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const NNLSSolver& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NNLS_SOLVER_H_INCLUDED  defined
