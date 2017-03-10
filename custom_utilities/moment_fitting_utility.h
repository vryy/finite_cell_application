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
#include <fstream>


// External includes
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/numeric/ublas/lu.hpp> 


// Project includes
#include "custom_utilities/quadrature_utility.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/function/product_function.h"
#include "custom_algebra/function/heaviside_function.h"
#include "custom_utilities/binary_tree.h"


#define ESTIMATE_RCOND
#define CHECK_SOLUTION


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

    typedef Function<PointType, double> FunctionR3R1Type;

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


    template<class TIntegratorType>
    void PyFitQuadrature(Element::Pointer& p_elem,
            boost::python::list& r_funcs,
            const LevelSet& r_level_set,
            const TIntegratorType& r_integrator,
            const int fit_quadrature_order, const int integrator_integration_order) const
    {
        std::vector<FunctionR3R1Type::Pointer> funcs;
        typedef boost::python::stl_input_iterator<FunctionR3R1Type::Pointer> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& f,
                      std::make_pair(iterator_value_type(r_funcs), // begin
                        iterator_value_type() ) ) // end
        {
            funcs.push_back(f);
        }

        GeometryData::IntegrationMethod ElementalIntegrationMethod
                = LevelSet::GetIntegrationMethod(fit_quadrature_order);

        GeometryData::IntegrationMethod IntegratorIntegrationMethod
                = LevelSet::GetIntegrationMethod(integrator_integration_order);

        this->FitQuadrature<FunctionR3R1Type, TIntegratorType>(p_elem->GetGeometry(),
                funcs, r_level_set, r_integrator, ElementalIntegrationMethod, IntegratorIntegrationMethod);
    }


    template<class TFunctionType, class TIntegratorType>
    Vector FitQuadrature(GeometryType& r_geom,
            const std::vector<typename TFunctionType::Pointer>& r_funcs,
            const LevelSet& r_level_set,
            const TIntegratorType& r_integrator,
            const GeometryData::IntegrationMethod& ElementalIntegrationMethod,
            const GeometryData::IntegrationMethod& IntegratorIntegrationMethod) const
    {
        std::size_t num_basis = r_funcs.size();

        const GeometryType::IntegrationPointsArrayType& integration_points
                = r_geom.IntegrationPoints( ElementalIntegrationMethod );

        std::size_t num_fit_points = integration_points.size();

        // form the matrix A
        CoordinatesArrayType GlobalCoords;
        Matrix MA(num_basis, num_fit_points);
        for(std::size_t i = 0; i < num_basis; ++i)
        {
            for(std::size_t j = 0; j < num_fit_points; ++j)
            {
                r_geom.GlobalCoordinates(GlobalCoords, integration_points[j]);
                MA(i, j) = r_funcs[i]->GetValue(GlobalCoords) * Function<double, double>::ComputeDetJ(r_geom, integration_points[j]);
            }
        }
KRATOS_WATCH(MA)

        // form the vector b
        Vector Mb(num_basis);
        typename TFunctionType::Pointer H = typename TFunctionType::Pointer(new HeavisideFunction<TFunctionType>(r_level_set));
        for(std::size_t i = 0; i < num_basis; ++i)
        {
            Mb(i) = 0.0;
            r_integrator.Integrate(ProductFunction<TFunctionType>(r_funcs[i], H), Mb(i), IntegratorIntegrationMethod);
        }
KRATOS_WATCH(Mb)

        // solve least square
        #ifdef ESTIMATE_RCOND
        double rcond = LeastSquareLAPACKSolver::EstimateRCond(MA);
        KRATOS_WATCH(rcond)
        #endif

        Vector Mw;
        if(MA.size1() == MA.size2())
        {
            /* the linear system is square, we can use a direct solver */
            boost::numeric::ublas::permutation_matrix<> pm(MA.size1());
            Matrix MAcopy = MA;
            Mw = Mb;
            boost::numeric::ublas::lu_factorize(MAcopy, pm);
            boost::numeric::ublas::lu_substitute(MAcopy, pm, Mw);
        }
        else
        {
            /* solve the non-square linear system by least square. NOTE: it can be very ill-conditioned */
            LeastSquareLAPACKSolver::SolveDGELSY(MA, Mw, Mb);
//            LeastSquareLAPACKSolver::SolveDGELSS(MA, Mw, Mb);
        }
        KRATOS_WATCH(Mw)
        KRATOS_WATCH(sum(Mw))
//        double pi = 3.141592653589793238462643;
//        KRATOS_WATCH(sum(Mw) - pi/4)

        #ifdef CHECK_SOLUTION
        Vector Error = (Mb - prod(MA, Mw)) / norm_2(Mb);
        KRATOS_WATCH(Error)
        #endif

        /* create new quadrature and assign to the geometry */
        FiniteCellGeometry<GeometryType>::AssignGeometryData(r_geom, ElementalIntegrationMethod, Mw);

        return Mw;
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
{}

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

#undef ESTIMATE_RCOND
#undef CHECK_SOLUTION

#endif // KRATOS_MOMENT_FITTING_UTILITY_H_INCLUDED  defined
