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


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/function.h"


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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MomentFittingUtility() {}

    /// Destructor.
    virtual ~MomentFittingUtility() {}


    template<class FunctionType, class IntegratorType>
    void FitQuadrature(Element::Pointer& p_elem,
            const FunctionType& r_func,
            const IntegratorType& r_integrator,
            const GeometryData::IntegrationMethod& ElementalIntegrationMethod,
            const GeometryData::IntegrationMethod& IntegratorIntegrationMethod)
    {
        
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
template<std::size_t TDimension>
inline std::istream& operator >> (std::istream& rIStream,
                MomentFittingUtility<TDimension>& rThis)
{}

/// output stream function
template<std::size_t TDimension>
inline std::ostream& operator << (std::ostream& rOStream,
                const MomentFittingUtility<TDimension>& rThis)
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
