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


#if !defined(KRATOS_FUNCTION_H_INCLUDED )
#define  KRATOS_FUNCTION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"


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
/** Abstract class for a general function R^m->R^n
*/
template<typename TInputType, typename TOutputType>
class Function
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Function
    KRATOS_CLASS_POINTER_DEFINITION(Function);

    typedef TInputType InputType;

    typedef TOutputType OutputType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Function() {}

    /// Destructor.
    virtual ~Function() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual TOutputType GetValue(const TInputType& P) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Call the base class", __FUNCTION__)
    }


    TOutputType Integrate(Element::Pointer& p_elem) const
    {
        return Integrate(p_elem->GetGeometry());
    }


    TOutputType Integrate(Element::Pointer& p_elem, const int integration_order) const
    {
        GeometryData::IntegrationMethod ThisIntegrationMethod
                = GetIntegrationMethod(integration_order);
        return Integrate(p_elem->GetGeometry(), ThisIntegrationMethod);
    }


    TOutputType Integrate(GeometryType& r_geom) const
    {
        return Integrate(r_geom, r_geom.GetDefaultIntegrationMethod());
    }


    /// Integrate a function using the sample geometry and integration rule
    TOutputType Integrate(GeometryType& r_geom,
            const GeometryData::IntegrationMethod ThisIntegrationMethod) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Integrate is not implemented", "")
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    static GeometryData::IntegrationMethod GetIntegrationMethod(const int integration_order)
    {
        GeometryData::IntegrationMethod ThisIntegrationMethod;

        if(integration_order == 1)
        {
            ThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if(integration_order == 2)
        {
            ThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if(integration_order == 3)
        {
            ThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if(integration_order == 4)
        {
            ThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if(integration_order == 5)
        {
            ThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Does not support for more integration points", "")

        return ThisIntegrationMethod;
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Function R^m->R^n";
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
    Function& operator=(Function const& rOther);

    /// Copy constructor.
    Function(Function const& rOther);


    ///@}

}; // Class Function

///@}

///@name Type Definitions
///@{

template<>
inline double Function<typename Element::GeometryType::PointType::PointType, double>::Integrate(GeometryType& r_geom,
        const GeometryData::IntegrationMethod ThisIntegrationMethod) const
{
    const GeometryType::IntegrationPointsArrayType& integration_points
            = r_geom.IntegrationPoints( ThisIntegrationMethod );

    double Result = 0.0;

    if(r_geom.WorkingSpaceDimension() == r_geom.LocalSpaceDimension())
    {
        Matrix J;
        double DetJ;

        CoordinatesArrayType GlobalCoords;

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            J = r_geom.Jacobian( J, integration_points[point] );
            DetJ = MathUtils<double>::Det(J);
            r_geom.GlobalCoordinates(GlobalCoords, integration_points[point]);

            Result += GetValue(GlobalCoords) * DetJ * integration_points[point].Weight();
        }
    }
    else
    {
        Matrix J, JtJ;
        double DetJ;

        CoordinatesArrayType GlobalCoords;

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            J = r_geom.Jacobian( J, integration_points[point] );
            JtJ = prod(trans(J), J);
            DetJ = sqrt(MathUtils<double>::Det(JtJ));
            r_geom.GlobalCoordinates(GlobalCoords, integration_points[point]);

            Result += GetValue(GlobalCoords) * DetJ * integration_points[point].Weight();
        }
    }

    return Result;
}

///@}
///@name Input and output
///@{


/// input stream function
template<typename TInputType, typename TOutputType>
inline std::istream& operator >> (std::istream& rIStream,
                Function<TInputType, TOutputType>& rThis)
{}

/// output stream function
template<typename TInputType, typename TOutputType>
inline std::ostream& operator << (std::ostream& rOStream,
                const Function<TInputType, TOutputType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FUNCTION_H_INCLUDED  defined
