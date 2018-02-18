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
//  Date:            14 Feb 2017
//


#if !defined(KRATOS_MOMENT_FITTING_UTILITY_H_INCLUDED )
#define  KRATOS_MOMENT_FITTING_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>


// External includes
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/progress.hpp>


// Project includes
#include "custom_algebra/brep.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/function/product_function.h"
#include "custom_algebra/function/heaviside_function.h"
#include "custom_linear_solvers/least_square_lapack_solver.h"
#include "custom_utilities/quadrature_utility.h"
#include "custom_utilities/finite_cell_auxilliary_utility.h"
#include "custom_utilities/finite_cell_geometry_utility.h"
#include "finite_cell_application/finite_cell_application.h"


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
/** Abstract class to compute the quadrature weights for volume integration using the moment fitting technique
*/
class MomentFittingUtility : public QuadratureUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MomentFittingUtility
    KRATOS_CLASS_POINTER_DEFINITION(MomentFittingUtility);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MomentFittingUtility() {}

    /// Destructor.
    virtual ~MomentFittingUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Construct the quadrature-fitted for a set of functions on the geometry
    /// r_geom: the geometry
    /// r_funcs: the set of function needs to be fitted; this function is assumed to be defined in the local frame
    /// r_brep: the BRep defined the cut boundary
    /// r_integrator: to perform the integration on the cut domain
    /// integration_points: list of integration points to compute the weight
    /// IntegratorIntegrationMethod: the integration method of the integrator
    template<class TFunctionType, class TIntegratorType>
    static Vector FitQuadrature(GeometryType& r_geom,
            const std::vector<typename TFunctionType::Pointer>& r_funcs,
            const BRep& r_brep,
            const TIntegratorType& r_integrator,
            const GeometryType::IntegrationPointsArrayType& integration_points,
            const int& integrator_integration_method,
            const std::string& solver_type,
            const int& echo_level,
            const double& small_weight)
    {
        if(echo_level > -1)
        {
            std::cout << std::scientific << std::setprecision(16);
            std::cout << "#######begin moment-fitting####solver_type: " << solver_type << "################" << std::endl;
        }

        std::size_t num_basis = r_funcs.size();
        if(echo_level > -1) KRATOS_WATCH(num_basis)

        std::size_t num_fit_points = integration_points.size();
        if(echo_level > -1) KRATOS_WATCH(num_fit_points)

        // form the matrix A
        CoordinatesArrayType GlobalCoords;
        Matrix MA(num_basis, num_fit_points);
        for(std::size_t i = 0; i < num_basis; ++i)
        {
            for(std::size_t j = 0; j < num_fit_points; ++j)
            {
                MA(i, j) = r_funcs[i]->GetValue(integration_points[j]) * Function<double, double>::ComputeDetJ(r_geom, integration_points[j]);
            }
        }
        if(echo_level > 4)
            KRATOS_WATCH(MA)

        // form the vector b
        Vector Mb(num_basis);
        for(std::size_t i = 0; i < num_basis; ++i)
        {
            // because the fitting functions is defined in the local frame, the integration in local scheme must be used
            Mb(i) = 0.0;
            r_integrator.template Integrate<double, 0>(*(r_funcs[i]), r_brep, Mb(i), integrator_integration_method, small_weight);
        }
        if(echo_level > 3)
            KRATOS_WATCH(Mb)

        if(echo_level > 2)
        {
            // estimate the condition number
            if(MA.size1() == MA.size2())
            {
                double rcond = LeastSquareLAPACKSolver::EstimateRCond(MA);
                KRATOS_WATCH(rcond)
            }
        }

        // solve for weights
        Vector Mw;
        if(solver_type == "direct")
        {
            if(MA.size1() != MA.size2())
            {
                std::stringstream ss;
                ss << "The matrix size " << MA.size1() << "x" << MA.size2() << " is not square. Direct solver cannot be used";
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }

            /* the linear system is square, we can use a direct solver */
            boost::numeric::ublas::permutation_matrix<> pm(MA.size1());
            Matrix MAcopy = MA;
            Mw = Mb;
            boost::numeric::ublas::lu_factorize(MAcopy, pm);
            boost::numeric::ublas::lu_substitute(MAcopy, pm, Mw);
        }
        else if(solver_type == "ata-direct")
        {
            /* form the linear system for least square problem */
            boost::numeric::ublas::permutation_matrix<> pm(MA.size2());
            Matrix MAcopy = prod(trans(MA), MA);
            Vector Mbcopy = prod(trans(MA), Mb);
            Mw = Mbcopy;
            boost::numeric::ublas::lu_factorize(MAcopy, pm);
            boost::numeric::ublas::lu_substitute(MAcopy, pm, Mw);
        }
        else if(solver_type == "dgelsy")
        {
            /* solve the non-square linear system by least square. NOTE: it can be very ill-conditioned */
            LeastSquareLAPACKSolver::SolveDGELSY(MA, Mw, Mb);
        }
        else if(solver_type == "dgelss")
        {
            /* solve the non-square linear system by least square. NOTE: it can be very ill-conditioned */
            LeastSquareLAPACKSolver::SolveDGELSS(MA, Mw, Mb);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Unknown solver type:", solver_type)

        if(echo_level > 0)
        {
            KRATOS_WATCH(Mw)
            KRATOS_WATCH(sum(Mw))
        }

        if(echo_level > 1)
        {
            // check the error of the solution
            Vector Error = (Mb - prod(MA, Mw)) / norm_2(Mb);
            KRATOS_WATCH(Error)
        }

        if(echo_level > -1)
            std::cout << "####################end moment-fitting####################" << std::endl;

        return Mw;
    }


    template<class TTreeType, class TFunctionType> // = MomentFittedQuadTreeSubCell
    void FitQuadratureSubCell(typename TTreeType::Pointer p_tree,
            const std::vector<typename TFunctionType::Pointer>& r_funcs,
            const BRep& r_brep,
            const int& integrator_integration_method,
            const std::string& solver_type,
            const int& echo_level,
            const double& small_weight) const
    {
        if(echo_level > 0)
            std::cout << "--------------------------------start fitting for element " << p_tree->pGetElement()->Id() << "--------------------------------" << std::endl;

        /* firstly get the physical integration point */
        // note that p_tree->GeneratePhysicalIntegrationPoints() must be called before to generate integration points data
        const std::vector<std::size_t>& subcell_index = p_tree->SubCellIndices();
        const GeometryType::IntegrationPointsArrayType& physical_integration_points = p_tree->GetPhysicalIntegrationPoints();

//            std::cout << "list of physical integration points for element " << p_tree->pGetElement()->Id() << std::endl;
//            for(std::size_t i = 0; i < physical_integration_points.size(); ++i)
//                std::cout << " " << physical_integration_points[i] << std::endl;
//            std::cout << std::endl;

        /* secondly assign the physical integration points to the element */
        GeometryData::IntegrationMethod RepresentativeIntegrationMethod = p_tree->GetRepresentativeIntegrationMethod();
        FiniteCellGeometryUtility::AssignGeometryData(*(p_tree->pGetGeometry()), RepresentativeIntegrationMethod, physical_integration_points);
        Variable<int>& INTEGRATION_ORDER_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("INTEGRATION_ORDER"));
        p_tree->pGetElement()->SetValue(INTEGRATION_ORDER_var, p_tree->GetRepresentativeIntegrationOrder());
        p_tree->pGetElement()->Initialize();

        /* thirdly fit the subcell */
        Matrix Weights = p_tree->FitQuadratureSubCell(subcell_index, r_funcs, r_brep, integrator_integration_method, solver_type, echo_level, small_weight);

        p_tree->pGetElement()->SetValue(SUBCELL_WEIGHTS, Weights);

        /* compute the subcell domain size */
        Vector DomainSizes(subcell_index.size());
        for(std::size_t i = 0; i < subcell_index.size(); ++i)
            DomainSizes(i) = p_tree->DomainSize(subcell_index[i], r_brep, integrator_integration_method);
        p_tree->pGetElement()->SetValue(SUBCELL_DOMAIN_SIZES, DomainSizes);

        if(echo_level > 0)
            std::cout << "--------------------------------end fitting for element " << p_tree->pGetElement()->Id() << "--------------------------------" << std::endl;
    }


    template<class TTreeType, class TFunctionType> // = MomentFittedQuadTreeSubCell
    void MultithreadedFitQuadratureSubCell(std::vector<typename TTreeType::Pointer>& r_trees,
            const std::vector<typename TFunctionType::Pointer>& r_funcs,
            const BRep& r_brep,
            const int& integrator_integration_method,
            const std::string& solver_type,
            const int& echo_level,
            const double& small_weight) const
    {
        /* multithreaded fit quadrature subcell */
        unsigned int number_of_threads = 1;
#ifdef _OPENMP
        number_of_threads = omp_get_max_threads();
#endif
        std::cout << "number_of_threads for MultithreadedFitQuadratureSubCell: " << number_of_threads << std::endl;
        std::vector<unsigned int> tree_partition;
        FiniteCellAuxilliaryUtility::CreatePartition(number_of_threads, r_trees.size(), tree_partition);
        std::cout << "tree_partition:";
        for(std::size_t i = 0; i < tree_partition.size(); ++i)
            std::cout << " " << tree_partition[i];
        std::cout << std::endl;

        boost::progress_display show_progress( r_trees.size() );

//        int integrator_quadrature_type = QuadratureUtility::GetQuadratureType(integrator_integration_method);
//        int integrator_quadrature_order = QuadratureUtility::GetQuadratureOrder(integrator_integration_method);
//        int integrator_integration_method_physical_point = 0x20 + integrator_quadrature_order;

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for(int k = 0; k < number_of_threads; ++k)
        {
            typename std::vector<typename TTreeType::Pointer>::iterator it_first_tree = r_trees.begin() + tree_partition[k];
            typename std::vector<typename TTreeType::Pointer>::iterator it_last_tree = r_trees.begin() + tree_partition[k+1];

            for(typename std::vector<typename TTreeType::Pointer>::iterator it = it_first_tree; it != it_last_tree; ++it)
            {
                this->template FitQuadratureSubCell<TTreeType, TFunctionType>(*it, r_funcs, r_brep, integrator_integration_method, solver_type, echo_level, small_weight);
                ++show_progress;
            }
        }
        std::cout << "\nMultithreadedFitQuadratureSubCell completed for " << r_trees.size() << " trees" << std::endl;
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
        return "Moment Fitting Utility";
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
    MomentFittingUtility& operator=(MomentFittingUtility const& rOther);

    /// Copy constructor.
    MomentFittingUtility(MomentFittingUtility const& rOther);


    ///@}

}; // Class MomentFittingUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, MomentFittingUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const MomentFittingUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MOMENT_FITTING_UTILITY_H_INCLUDED  defined
