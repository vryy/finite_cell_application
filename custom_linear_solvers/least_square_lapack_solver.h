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
//  Date:            13 Feb 2017
//


#if !defined(KRATOS_LEAST_SQUARE_LAPACK_SOLVER_H_INCLUDED )
#define  KRATOS_LEAST_SQUARE_LAPACK_SOLVER_H_INCLUDED



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


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"



// template for LAPACK functions
extern "C" void dgelsy_(int* M, int* N, int* NRHS, double* A, int* LDA,
                        double* B, int* LDB, int* JPVT, double* RCOND, int* RANK,
                        double* WORK, int* LWORK, int* INFO);

extern "C" void dgelss_(int* M, int* N, int* NRHS, double* A, int* LDA,
                        double* B, int* LDB, double* S, double* RCOND, int* RANK,
                        double* WORK, int* LWORK, int* INFO);

extern "C" double dlange_(char* NORM, int* M, int* N, double* A, int* LDA, double* WORK);

extern "C" void dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);

extern "C" void dgecon_(char* NORM, int* N, double* A, int* LDA, double* ANORM,
                        double* RCOND, double* WORK, int* IWORK, int* INFO);


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
/** Least square solver based on LAPACK
*/
class LeastSquareLAPACKSolver
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LeastSquareLAPACKSolver
    KRATOS_CLASS_POINTER_DEFINITION(LeastSquareLAPACKSolver);


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LeastSquareLAPACKSolver() {}

    /// Destructor.
    virtual ~LeastSquareLAPACKSolver() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Estimate the reciprocal condition number of the matrix
    static double EstimateRCond(Matrix& rA, char* norm_type = "1")
    {
        int M = rA.size1();
        int N = rA.size2();
        Matrix At = trans(rA);
        int LDA = M;
        int INFO;

#ifdef ENABLE_FINITE_CELL_BOOST_BINDINGS
        typedef boost::numeric::bindings::traits::matrix_traits<Matrix> matraits;
        double* A = matraits::storage(At);
#else
        double* A = new double[M * N];
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
            {
                A[i * N + j] = rA(j, i);    // TODO check
            }
#endif

        double* WORK = new double[4 * N];
        double anorm = dlange_(norm_type, &M, &N, A, &LDA, WORK);
//        KRATOS_WATCH(anorm)

        int* IPIV = new int[std::min(M, N)];
        dgetrf_(&M, &N, A, &LDA, IPIV, &INFO);

        int* IWORK = new int[N];
        double rcond;
        LDA = N;
        dgecon_(norm_type, &N, A, &LDA, &anorm, &rcond, WORK, IWORK, &INFO);

#ifndef ENABLE_FINITE_CELL_BOOST_BINDINGS
        delete A;
#endif

        delete IPIV, WORK, IWORK;
        return rcond;
    }


    /// Solve the non-square linear system using least square method
    static int SolveDGELSY(Matrix& rA, Vector& rX, Vector& rB, const double rcond_est = 0.01)
    {
        int M = rA.size1();
        int N = rA.size2();
        int NRHS = 1;
        Matrix At = trans(rA);
#ifdef ENABLE_FINITE_CELL_BOOST_BINDINGS
        typedef boost::numeric::bindings::traits::matrix_traits<Matrix> matraits;
        double* A = matraits::storage(At);
#else
        double* A = new double[M * N];
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
            {
                A[i * N + j] = rA(j, i);    // TODO check
            }
#endif
        int LDA = M;
        int INFO;
        int LDB = std::max(M, N);
        double* b = new double[LDB];
        for (int i = 0; i < M; ++i) { b[i] = rB(i); }
        int* JPVT = new int[N]; for (int i = 0; i < N; ++i) JPVT[i] = 0;
        double RCOND = rcond_est;
        int RANK;
        int LWORK = -1;//std::max(std::min(M, N) + 3*N + 1, 2*std::min(M, N) + 1);
        double WORKS;

        // estimate the size of working array
        dgelsy_(&M, &N, &NRHS, A, &LDA, b, &LDB, JPVT, &RCOND, &RANK, &WORKS, &LWORK, &INFO);
        LWORK = (int) WORKS;
//        KRATOS_WATCH(LWORK)

        double* WORK = new double[LWORK];

        dgelsy_(&M, &N, &NRHS, A, &LDA, b, &LDB, JPVT, &RCOND, &RANK, WORK, &LWORK, &INFO);
//        KRATOS_WATCH(RANK)

        if (rX.size() != N)
        {
            rX.resize(N, false);
        }
        for (int i = 0; i < N; ++i)
        {
            rX(i) = b[i];
        }

        delete JPVT, b, WORK;
#ifndef ENABLE_FINITE_CELL_BOOST_BINDINGS
        delete A;
#endif

        return 0;
    }


    /// Solve the non-square linear system using least square method
    static int SolveDGELSS(Matrix& rA, Vector& rX, Vector& rB, const double rcond_est = 1.0e-10)
    {
        int M = rA.size1();
        int N = rA.size2();
        int NRHS = 1;
        Matrix At = trans(rA);
#ifdef ENABLE_FINITE_CELL_BOOST_BINDINGS
        typedef boost::numeric::bindings::traits::matrix_traits<Matrix> matraits;
        double* A = matraits::storage(At);
#else
        double* A = new double[M * N];
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
            {
                A[i * N + j] = rA(j, i);    // TODO check
            }
#endif
        int LDA = M;
        int INFO;
        int LDB = std::max(M, N);
        double* b = new double[LDB];
        for (int i = 0; i < M; ++i) { b[i] = rB(i); }
        double* S = new double[std::min(M, N)];
        double RCOND = rcond_est;
        int RANK;
        int LWORK = -1;
        double WORKS;

        // estimate the size of working array
        dgelss_(&M, &N, &NRHS, A, &LDA, b, &LDB, S, &RCOND, &RANK, &WORKS, &LWORK, &INFO);
        LWORK = (int) WORKS;
//        KRATOS_WATCH(LWORK)

        double* WORK = new double[LWORK];

        dgelss_(&M, &N, &NRHS, A, &LDA, b, &LDB, S, &RCOND, &RANK, WORK, &LWORK, &INFO);
//        KRATOS_WATCH(RANK)

        if (rX.size() != N)
        {
            rX.resize(N, false);
        }
        for (int i = 0; i < N; ++i)
        {
            rX(i) = b[i];
        }

        delete S, b, WORK;
#ifndef ENABLE_FINITE_CELL_BOOST_BINDINGS
        delete A;
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
        return "Least Square Solver based on LAPACK";
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
    LeastSquareLAPACKSolver& operator=(LeastSquareLAPACKSolver const& rOther);

    /// Copy constructor.
    LeastSquareLAPACKSolver(LeastSquareLAPACKSolver const& rOther);


    ///@}

}; // Class LeastSquareLAPACKSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  LeastSquareLAPACKSolver& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const LeastSquareLAPACKSolver& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LEAST_SQUARE_LAPACK_SOLVER_H_INCLUDED  defined
