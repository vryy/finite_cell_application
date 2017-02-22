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


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "custom_algebra/function.h"
#include "custom_algebra/product_function.h"
#include "custom_algebra/heaviside_function.h"
#include "custom_utilities/binary_tree.h"
#include "custom_geometries/finite_cell_geometry.h"


#define ESTIMATE_RCOND


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
class MomentFittingUtility
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


    template<std::size_t TDimension>
    void PyFitQuadrature(Element::Pointer& p_elem,
            boost::python::list& r_funcs,
            const LevelSet& r_level_set,
            const BinaryTree<TDimension>& r_integrator,
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

        this->FitQuadrature<FunctionR3R1Type, BinaryTree<TDimension> >(p_elem->GetGeometry(),
                funcs, r_level_set, r_integrator, ElementalIntegrationMethod, IntegratorIntegrationMethod);
    }


    template<class FunctionType, class IntegratorType>
    Vector FitQuadrature(GeometryType& r_geom,
            const std::vector<typename FunctionType::Pointer>& r_funcs,
            const LevelSet& r_level_set,
            const IntegratorType& r_integrator,
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
                MA(i, j) = r_funcs[i]->GetValue(GlobalCoords);
            }
        }
//KRATOS_WATCH(MA)

        // form the vector b
        Vector Mb(num_basis);
        typename FunctionType::Pointer H = typename FunctionType::Pointer(new HeavisideFunction(r_level_set));
        for(std::size_t i = 0; i < num_basis; ++i)
        {
            Mb(i) = 0.0;
            r_integrator.Integrate(ProductFunction(r_funcs[i], H), Mb(i), IntegratorIntegrationMethod);
        }
//KRATOS_WATCH(Mb)

        // solve least square
        #ifdef ESTIMATE_RCOND
        double rcond = LeastSquareLAPACKSolver::EstimateRCond(MA);
        KRATOS_WATCH(rcond)
        #endif
        /* solve the non-square linear system by least square. NOTE: it can be very ill-conditioned */
        Vector Mw;
        LeastSquareLAPACKSolver::SolveDGELSY(MA, Mw, Mb);
//            LeastSquareLAPACKSolver::SolveDGELSS(MA, Mw, Mb);
        KRATOS_WATCH(Mw)
//        KRATOS_WATCH(sum(Mw))
//        double pi = 3.141592653589793238462643;
//        KRATOS_WATCH(sum(Mw) - pi/4)

        /* create new quadrature and assign to the geometry */
        FiniteCellGeometry<GeometryType>::AssignGeometryData(r_geom, ElementalIntegrationMethod, Mw);

        return Mw;
    }


    void PyScaleQuadrature(Element::Pointer& p_elem,
            const int quadrature_order, const double ScaleFactor) const
    {
        GeometryData::IntegrationMethod ElementalIntegrationMethod
                = LevelSet::GetIntegrationMethod(quadrature_order);

        this->ScaleQuadrature(p_elem->GetGeometry(), ElementalIntegrationMethod, ScaleFactor);
    }


    void ScaleQuadrature(GeometryType& r_geom,
            const GeometryData::IntegrationMethod& ElementalIntegrationMethod,
            const double ScaleFactor) const
    {
        const GeometryType::IntegrationPointsArrayType& integration_points
                = r_geom.IntegrationPoints( ElementalIntegrationMethod );

        Vector Mw(integration_points.size());
        for(std::size_t i = 0; i < integration_points.size(); ++i)
        {
            Mw(i) = integration_points[i].Weight() * ScaleFactor;
        }

        /* create new quadrature and assign to the geometry */
        FiniteCellGeometry<GeometryType>::AssignGeometryData(r_geom, ElementalIntegrationMethod, Mw);
    }


    void PySaveQuadrature(boost::python::list& pyElemList, const std::string& fileName) const
    {
        std::ofstream myFile;
        myFile.open(fileName.c_str(), std::ios::out | std::ios::binary);

        typedef boost::python::stl_input_iterator<Element::Pointer> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& p_elem,
                        std::make_pair(iterator_value_type(pyElemList), // begin
                        iterator_value_type() ) ) // end
        {
            this->SaveQuadrature(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
        }

        myFile.close();
    }


    void SaveQuadrature(std::ofstream& rOStream, const Element::Pointer& p_elem,
            const GeometryData::IntegrationMethod& ElementalIntegrationMethod) const
    {
        rOStream << p_elem->Id() << ElementalIntegrationMethod;
        const GeometryType::IntegrationPointsArrayType& integration_points
                = p_elem->GetGeometry().IntegrationPoints( ElementalIntegrationMethod );

        for(std::size_t i = 0; i < integration_points.size(); ++i)
        {
            rOStream << integration_points[i].Weight();
        }
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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

#endif // KRATOS_MOMENT_FITTING_UTILITY_H_INCLUDED  defined
