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
//  Date:            24 Feb 2017
//


#if !defined(KRATOS_IMMERSED_BOUNDARY_UTILITY_H_INCLUDED )
#define  KRATOS_IMMERSED_BOUNDARY_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/line_3d_2.h"
#include "custom_algebra/function.h"


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
/** Abstract class for immersed boundary search
*/
class ImmersedBoundaryUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ImmersedBoundaryUtility
    KRATOS_CLASS_POINTER_DEFINITION(ImmersedBoundaryUtility);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ImmersedBoundaryUtility() {}

    /// Destructor.
    virtual ~ImmersedBoundaryUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Sampling along the parametric curve to provide a list of points and corresponding variational length and weight
    template<class TCurveType>
    static void SamplingCurve(const TCurveType& rCurve, const double& tmin, const double& tmax, const int& integration_order,
            std::vector<double>& rCoordinates, std::vector<PointType>& rPoints,
            std::vector<double>& rDlength, std::vector<double>& rWeights)
    {
        GeometryData::IntegrationMethod ThisIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(integration_order);

        NodeType DummyNode1;
        NodeType DummyNode2;
        Line3D2<NodeType> SampleLine(DummyNode1, DummyNode2);

        const GeometryType::IntegrationPointsArrayType& integration_points
                = SampleLine.IntegrationPoints( ThisIntegrationMethod );

        if(rPoints.size() != integration_points.size())
            rPoints.resize(integration_points.size());

        if(rDlength.size() != integration_points.size())
            rDlength.resize(integration_points.size());

        if(rWeights.size() != integration_points.size())
            rWeights.resize(integration_points.size());

        double t;
        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            t = 0.5*(tmax+tmin) + 0.5*(tmax-tmin)*integration_points[point].X();

            rCoordinates[point] = t;

            rPoints[point] = rCurve.GetValue(t);

            PointType tangent = rCurve.GetDerivative(0, t);
            rDlength[point] = norm_2(tangent);

            rWeights[point] = integration_points[point].Weight();
        }
    }


    /// Find an element in pMasterElements contains rSourcePoint and assign it to pTargetElement. The rLocalTargetPoint is the local point in pTargetElement of rSourcePoint
    static bool SearchPartner( const PointType& rSourcePoint, ModelPart::ElementsContainerType& pMasterElements,
            Element::Pointer& pTargetElement, PointType& rLocalTargetPoint )
    {
        ModelPart::ElementsContainerType pMasterElementsCandidates;

        // find the potential master elements
        pMasterElementsCandidates.clear();
        for( typename ModelPart::ElementsContainerType::ptr_iterator it = pMasterElements.ptr_begin();
                it != pMasterElements.ptr_end(); ++it )
        {
            // compute the center of the element
            PointType Center(0.0, 0.0, 0.0);
            for( unsigned int node = 0; node < (*it)->GetGeometry().size(); ++node )
                noalias(Center) += (*it)->GetGeometry()[node];
            Center /= (*it)->GetGeometry().size();

            // compute the maximum distance from center to vertices
            double max_dist = 0.0;
            for( unsigned int node = 0; node < (*it)->GetGeometry().size(); ++node )
            {
                double dist = norm_2((*it)->GetGeometry()[node] - Center);
                if(dist > max_dist)
                    max_dist = dist;
            }

            // compute the distance of the source point to the center of the element
            double sdist = norm_2(rSourcePoint - Center);
            if(sdist < max_dist)
                pMasterElementsCandidates.push_back(*it);
        }

        for( typename ModelPart::ElementsContainerType::ptr_iterator it = pMasterElementsCandidates.ptr_begin();
                it != pMasterElementsCandidates.ptr_end(); ++it )
        {
            bool is_inside = (*it)->GetGeometry().IsInside( rSourcePoint, rLocalTargetPoint );
            if( is_inside )
            {
                pTargetElement = *it;
                return true;
            }
        }

        std::cout << " !!!! WARNING: NO ELEMENT FOUND TO CONTAIN " << rSourcePoint << " !!!! " << std::endl;
        return false;
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
        return "Immersed Boundary Utility";
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
    ImmersedBoundaryUtility& operator=(ImmersedBoundaryUtility const& rOther);

    /// Copy constructor.
    ImmersedBoundaryUtility(ImmersedBoundaryUtility const& rOther);


    ///@}

}; // Class ImmersedBoundaryUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, ImmersedBoundaryUtility& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const ImmersedBoundaryUtility& rThis)
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

#endif // KRATOS_IMMERSED_BOUNDARY_UTILITY_H_INCLUDED  defined
