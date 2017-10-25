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


    /// Python interface for FitQuadrature
    /// REMARKS: if fit_quadrature_method = 0x0X, this function uses the standard Gauss quadrature rule on the element, defined by fit_quadrature_order
    /// Otherwise it will fit a custom quadrature, see quadrature_utility.h for decoding of the integration_method
    template<class TIntegratorType>
    void PyFitQuadrature(Element::Pointer& p_elem,
            boost::python::list& r_funcs,
            const BRep& r_brep,
            const TIntegratorType& r_integrator,
            const int& fit_quadrature_method,
            const int& integrator_integration_method,
            const std::string& solver_type,
            const int& echo_level,
            const double& small_weight) const
    {
        std::vector<FunctionR3R1::Pointer> funcs;
        typedef boost::python::stl_input_iterator<FunctionR3R1::Pointer> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& f,
                      std::make_pair(iterator_value_type(r_funcs), // begin
                        iterator_value_type() ) ) // end
        {
            funcs.push_back(f);
        }

        int fit_quadrature_type = QuadratureUtility::GetQuadratureType(fit_quadrature_method);

        int fit_quadrature_order = QuadratureUtility::GetQuadratureOrder(fit_quadrature_method);

        GeometryData::IntegrationMethod ElementalIntegrationMethod
                = Function<double, double>::GetIntegrationMethod(fit_quadrature_order);

        if(fit_quadrature_type == 0)
        {

            const GeometryType::IntegrationPointsArrayType& integration_points
                    = p_elem->GetGeometry().IntegrationPoints( ElementalIntegrationMethod );

            Vector Weight = this->FitQuadrature<FunctionR3R1, TIntegratorType>(p_elem->GetGeometry(),
                    funcs, r_brep, r_integrator, integration_points,
                    integrator_integration_method, solver_type, echo_level, small_weight);

            /* create new quadrature and assign to the geometry */
            FiniteCellGeometryUtility::AssignGeometryData(p_elem->GetGeometry(), ElementalIntegrationMethod, Weight);
        }
        else
        {
            GeometryType::IntegrationPointsArrayType integration_points
                    = r_integrator.Get().ConstructCustomQuadrature(fit_quadrature_type, fit_quadrature_order);

            Vector Weight = this->FitQuadrature<FunctionR3R1, TIntegratorType>(p_elem->GetGeometry(),
                    funcs, r_brep, r_integrator, integration_points,
                    integrator_integration_method, solver_type, echo_level, small_weight);

            for(std::size_t i = 0; i < integration_points.size(); ++i)
                integration_points[i].Weight() = Weight(i);

            /* create new quadrature and assign to the geometry */
            FiniteCellGeometryUtility::AssignGeometryData(p_elem->GetGeometry(), ElementalIntegrationMethod, integration_points);
        }
    }


    /// Python interface for FitQuadrature
    /// REMARKS: if fit_quadrature_method = 0x0X, this function uses the standard Gauss quadrature rule on the element, defined by fit_quadrature_order
    /// Otherwise it will fit a custom quadrature, see quadrature_utility.h for decoding of the integration_method
    template<class TIntegratorType>
    void PyMultithreadedFitQuadrature(boost::python::list& r_elements,
            boost::python::list& r_funcs,
            const BRep& r_brep,
            boost::python::list& r_integrators,
            const int& fit_quadrature_method,
            const int& integrator_integration_method,
            const std::string& solver_type,
            const int& echo_level,
            const double& small_weight) const
    {
        std::vector<Element::Pointer> elements;
        typedef boost::python::stl_input_iterator<Element::Pointer> iterator_element_type;
        BOOST_FOREACH(const iterator_element_type::value_type& e,
                      std::make_pair(iterator_element_type(r_elements), // begin
                      iterator_element_type() ) ) // end
        {
            elements.push_back(e);
        }
//        KRATOS_WATCH(elements.size())

        std::vector<FunctionR3R1::Pointer> funcs;
        typedef boost::python::stl_input_iterator<FunctionR3R1::Pointer> iterator_func_type;
        BOOST_FOREACH(const iterator_func_type::value_type& f,
                      std::make_pair(iterator_func_type(r_funcs), // begin
                      iterator_func_type() ) ) // end
        {
            funcs.push_back(f);
        }
//        KRATOS_WATCH(funcs.size())

        typedef typename TIntegratorType::Pointer TIntegratorPointerType;
        std::vector<TIntegratorPointerType> integrators;
        typedef boost::python::stl_input_iterator<TIntegratorPointerType> iterator_integrator_type;
        BOOST_FOREACH(const typename iterator_integrator_type::value_type& i,
                      std::make_pair(iterator_integrator_type(r_integrators), // begin
                      iterator_integrator_type() ) ) // end
        {
            integrators.push_back(i);
        }
//        KRATOS_WATCH(integrators.size())

        int fit_quadrature_type = QuadratureUtility::GetQuadratureType(fit_quadrature_method);

        int fit_quadrature_order = QuadratureUtility::GetQuadratureOrder(fit_quadrature_method);

        GeometryData::IntegrationMethod ElementalIntegrationMethod
                = Function<double, double>::GetIntegrationMethod(fit_quadrature_order);

        unsigned int number_of_threads = 1;
        std::vector<unsigned int> element_partition;
#ifdef _OPENMP
        number_of_threads = omp_get_max_threads();
#endif
        KRATOS_WATCH(number_of_threads)
        FiniteCellAuxilliaryUtility::CreatePartition(number_of_threads, elements.size(), element_partition);
        std::cout << "element_partition:";
        for(std::size_t i = 0; i < element_partition.size(); ++i)
            std::cout << " " << element_partition[i];
        std::cout << std::endl;

        boost::progress_display* show_progress = NULL;

        if(echo_level == -2)
            show_progress = new boost::progress_display( elements.size() );

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for(int k = 0; k < number_of_threads; ++k)
        {
            std::vector<Element::Pointer>::iterator it_first_element = elements.begin() + element_partition[k];
            std::vector<Element::Pointer>::iterator it_last_element = elements.begin() + element_partition[k+1];

            typename std::vector<TIntegratorPointerType>::iterator it_integrator = integrators.begin() + element_partition[k];

            for(std::vector<Element::Pointer>::iterator it = it_first_element; it != it_last_element; ++it, ++it_integrator)
            {
                if(fit_quadrature_type == 0)
                {
                    const GeometryType::IntegrationPointsArrayType& integration_points
                            = (*it)->GetGeometry().IntegrationPoints( ElementalIntegrationMethod );

                    Vector Weight = this->FitQuadrature<FunctionR3R1, TIntegratorType>((*it)->GetGeometry(),
                            funcs, r_brep, *(*it_integrator), integration_points,
                            integrator_integration_method, solver_type, echo_level, small_weight);

                    /* create new quadrature and assign to the geometry */
                    FiniteCellGeometryUtility::AssignGeometryData((*it)->GetGeometry(), ElementalIntegrationMethod, Weight);
                }
                else
                {
                    GeometryType::IntegrationPointsArrayType integration_points
                            = (*it_integrator)->Get().ConstructCustomQuadrature(fit_quadrature_type, fit_quadrature_order);
KRATOS_WATCH(fit_quadrature_type)
KRATOS_WATCH(fit_quadrature_order)
KRATOS_WATCH(integration_points.size())
                    Vector Weight = this->FitQuadrature<FunctionR3R1, TIntegratorType>((*it)->GetGeometry(),
                            funcs, r_brep, *(*it_integrator), integration_points,
                            integrator_integration_method, solver_type, echo_level, small_weight);

                    for(std::size_t i = 0; i < integration_points.size(); ++i)
                        integration_points[i].Weight() = Weight(i);

                    /* create new quadrature and assign to the geometry */
                    FiniteCellGeometryUtility::AssignGeometryData((*it)->GetGeometry(), ElementalIntegrationMethod, integration_points);
                }

                if(echo_level == -2)
                    ++(*show_progress);
            }
        }
        if(echo_level == -2)
        {
            delete show_progress;
            std::cout << std::endl;
        }
    }


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


    template<class TTreeType> // = MomentFittedQuadTreeSubCell
    void PyMultithreadedFitQuadratureSubCell(boost::python::list& r_trees,
            boost::python::list& r_funcs,
            const BRep& r_brep,
            const int& integrator_integration_method,
            const std::string& solver_type,
            const int& echo_level,
            const double& small_weight) const
    {
        /* extract the tree with subcell */
        typedef typename TTreeType::Pointer TTreePointerType;
        std::vector<TTreePointerType> trees;
        typedef boost::python::stl_input_iterator<TTreePointerType> iterator_tree_type;
        BOOST_FOREACH(const typename iterator_tree_type::value_type& t,
                      std::make_pair(iterator_tree_type(r_trees), // begin
                      iterator_tree_type() ) ) // end
        {
            trees.push_back(t);
        }

        /* extract the fitting functions */
        std::vector<FunctionR3R1::Pointer> funcs;
        typedef boost::python::stl_input_iterator<FunctionR3R1::Pointer> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& f,
                      std::make_pair(iterator_value_type(r_funcs), // begin
                        iterator_value_type() ) ) // end
        {
            funcs.push_back(f);
        }

        /* multithreaded fit quadrature subcell */
        unsigned int number_of_threads = 1;
#ifdef _OPENMP
        number_of_threads = omp_get_max_threads();
#endif
        std::cout << "number_of_threads for MultithreadedFitQuadratureSubCell: " << number_of_threads << std::endl;
        std::vector<unsigned int> tree_partition;
        FiniteCellAuxilliaryUtility::CreatePartition(number_of_threads, trees.size(), tree_partition);
        std::cout << "tree_partition:";
        for(std::size_t i = 0; i < tree_partition.size(); ++i)
            std::cout << " " << tree_partition[i];
        std::cout << std::endl;

        boost::progress_display show_progress( trees.size() );

//        int integrator_quadrature_type = QuadratureUtility::GetQuadratureType(integrator_integration_method);
//        int integrator_quadrature_order = QuadratureUtility::GetQuadratureOrder(integrator_integration_method);
//        int integrator_integration_method_physical_point = 0x20 + integrator_quadrature_order;

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for(int k = 0; k < number_of_threads; ++k)
        {
            typename std::vector<TTreePointerType>::iterator it_first_tree = trees.begin() + tree_partition[k];
            typename std::vector<TTreePointerType>::iterator it_last_tree = trees.begin() + tree_partition[k+1];

            for(typename std::vector<TTreePointerType>::iterator it = it_first_tree; it != it_last_tree; ++it)
            {
                if(echo_level > 0)
                    std::cout << "--------------------------------start fitting for element " << (*it)->pGetElement()->Id() << "--------------------------------" << std::endl;

                /* firstly compute the physical integration point */
                std::pair<std::vector<std::size_t>, GeometryType::IntegrationPointsArrayType> Output
                        = (*it)->GetPhysicalInterationPoint(r_brep, integrator_integration_method);
                const std::vector<std::size_t>& subcell_index = Output.first;
                const GeometryType::IntegrationPointsArrayType& physical_integration_points = Output.second;

                /* secondly assign the physical integration points to the element */
                GeometryData::IntegrationMethod RepresentativeIntegrationMethod = (*it)->GetRepresentativeIntegrationMethod();
                FiniteCellGeometryUtility::AssignGeometryData(*(*it)->pGetGeometry(), RepresentativeIntegrationMethod, physical_integration_points);
                Variable<int>& INTEGRATION_ORDER_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("INTEGRATION_ORDER"));
                (*it)->pGetElement()->SetValue(INTEGRATION_ORDER_var, (*it)->GetRepresentativeIntegrationOrder());
                (*it)->pGetElement()->Initialize();

                /* thirdly fit the subcell */
                Matrix Weights = (*it)->FitQuadratureSubCell(subcell_index, funcs, r_brep, integrator_integration_method, solver_type, echo_level, small_weight);

                (*it)->pGetElement()->SetValue(SUBCELL_WEIGHTS, Weights);

                /* compute the subcell domain size */
                Vector DomainSizes(subcell_index.size());
                for(std::size_t i = 0; i < subcell_index.size(); ++i)
                    DomainSizes(i) = (*it)->DomainSize(subcell_index[i], r_brep, integrator_integration_method);
                (*it)->pGetElement()->SetValue(SUBCELL_DOMAIN_SIZES, DomainSizes);

                if(echo_level > 0)
                    std::cout << "--------------------------------end fitting for element " << (*it)->pGetElement()->Id() << "--------------------------------" << std::endl;

                ++show_progress;
            }
        }
        std::cout << "\nMultithreadedFitQuadratureSubCell completed for " << trees.size() << " trees" << std::endl;
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
