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
#include <limits>


// External includes
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/foreach.hpp>


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/line_3d_2.h"
#include "geometries/quadrilateral_3d_4.h"
#include "utilities/math_utils.h"
#include "custom_algebra/function/function.h"


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

    static const int _REPORT_NUMBER_OF_CREATED_CONDITIONS = 0x01;
    static const int _WARNING_FOUND_NO_ELEMENT = 0x02;
    static const int _REPORT_CONDITION_CREATED = 0x04;

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


    void InitializeBinning(ModelPart::ElementsContainerType& pElements,
            const double& Dx, const double& Dy, const double& Dz)
    {
        mDx = Dx;
        mDy = Dy;
        mDz = Dz;
        // IMPORTANT REMARK: one should not choose a spacing align with element edge, because it can create missing points.


        for( typename ModelPart::ElementsContainerType::ptr_iterator it = pElements.ptr_begin();
                it != pElements.ptr_end(); ++it )
        {
            Element::GeometryType& r_geom = (*it)->GetGeometry();

            // find the minimum and maximum binning coordinates in each direction
            int min_ix = std::numeric_limits<int>::max();
            int max_ix = std::numeric_limits<int>::min();
            int min_iy = std::numeric_limits<int>::max();
            int max_iy = std::numeric_limits<int>::min();
            int min_iz = std::numeric_limits<int>::max();
            int max_iz = std::numeric_limits<int>::min();

            for(std::size_t i = 0; i < r_geom.size(); ++i)
            {
                // find the cell containing point
                int ix = (mDx != 0.0) ? (int) floor(r_geom[i].X()) / mDx : 0;
                int iy = (mDy != 0.0) ? (int) floor(r_geom[i].Y()) / mDy : 0;
                int iz = (mDz != 0.0) ? (int) floor(r_geom[i].Z()) / mDz : 0;

                // adjust the maximum and minimum value in each direction
                if(ix < min_ix) min_ix = ix;
                if(ix > max_ix) max_ix = ix;
                if(iy < min_iy) min_iy = iy;
                if(iy > max_iy) max_iy = iy;
                if(iz < min_iz) min_iz = iz;
                if(iz > max_iz) max_iz = iz;
            }

            // add to the bin
            for(int ix = min_ix; ix <= max_ix; ++ix)
            {
                for(int iy = min_iy; iy <= max_iy; ++iy)
                {
                    for(int iz = min_iz; iz <= max_iz; ++iz)
                    {
                        SpatialKey key(ix, iy, iz);

                        // insert the element to spatial bin
                        mBinElements[key].insert((*it)->Id());
                    }
                }
            }
        }

        std::cout << "Setup binning completed, number of binning keys = " << mBinElements.size() << std::endl;
    }


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

        if(rCoordinates.size() != integration_points.size())
            rCoordinates.resize(integration_points.size());

        if(rPoints.size() != integration_points.size())
            rPoints.resize(integration_points.size());

        if(rDlength.size() != integration_points.size())
            rDlength.resize(integration_points.size());

        if(rWeights.size() != integration_points.size())
            rWeights.resize(integration_points.size());

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            double t = 0.5*(tmax+tmin) + 0.5*(tmax-tmin)*integration_points[point].X();

            rCoordinates[point] = t;

            rPoints[point] = rCurve.GetValue(t);

            PointType tangent = rCurve.GetDerivative(0, t);
            rDlength[point] = norm_2(tangent);

            rWeights[point] = integration_points[point].Weight() * 0.5 * (tmax-tmin);
//KRATOS_WATCH(t)
//KRATOS_WATCH(integration_points[point].X())
//KRATOS_WATCH(tangent)
//KRATOS_WATCH(rDlength[point])
//KRATOS_WATCH(rWeights[point])
        }
    }


    /// Sampling along the parametric surface to provide a list of points and corresponding variational area and weight
    template<class TSurfaceType>
    static void SamplingSurface(const TSurfaceType& rSurface,
            const double& t1min, const double& t1max,
            const double& t2min, const double& t2max,
            const int& integration_order,
            std::vector<typename TSurfaceType::InputType>& rCoordinates,
            std::vector<PointType>& rPoints,
            std::vector<double>& rDarea, std::vector<double>& rWeights)
    {
        GeometryData::IntegrationMethod ThisIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(integration_order);

        NodeType DummyNode1;
        NodeType DummyNode2;
        NodeType DummyNode3;
        NodeType DummyNode4;
        Quadrilateral3D4<NodeType> SampleQuad(DummyNode1, DummyNode2, DummyNode3, DummyNode4);

        const GeometryType::IntegrationPointsArrayType& integration_points
                = SampleQuad.IntegrationPoints( ThisIntegrationMethod );

        if(rCoordinates.size() != integration_points.size())
            rCoordinates.resize(integration_points.size());

        if(rPoints.size() != integration_points.size())
            rPoints.resize(integration_points.size());

        if(rDarea.size() != integration_points.size())
            rDarea.resize(integration_points.size());

        if(rWeights.size() != integration_points.size())
            rWeights.resize(integration_points.size());

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            double t1 = 0.5*(t1max+t1min) + 0.5*(t1max-t1min)*integration_points[point].X();
            double t2 = 0.5*(t2max+t2min) + 0.5*(t2max-t2min)*integration_points[point].Y();

            rCoordinates[point][0] = t1;
            rCoordinates[point][1] = t2;

            rPoints[point] = rSurface.GetValue(rCoordinates[point]);

            PointType tangent1 = rSurface.GetDerivative(0, rCoordinates[point]);
            PointType tangent2 = rSurface.GetDerivative(1, rCoordinates[point]);
            rDarea[point] = norm_2(MathUtils<double>::CrossProduct(tangent1, tangent2));

            rWeights[point] = integration_points[point].Weight() * 0.25 * (t1max-t1min) * (t2max-t2min);
//KRATOS_WATCH(t1)
//KRATOS_WATCH(t2)
//KRATOS_WATCH(integration_points[point].X())
//KRATOS_WATCH(integration_points[point].Y())
//KRATOS_WATCH(tangent1)
//KRATOS_WATCH(tangent2)
//KRATOS_WATCH(rDarea[point])
//KRATOS_WATCH(rWeights[point])
        }
    }


    /// Find an element in pMasterElements contains rSourcePoint and assign it to pTargetElement. The rLocalTargetPoint is the local point in pTargetElement of rSourcePoint
    /// REMARK: we should disable the move mesh flag if we want to search in the reference configuration
    static bool SearchPartner( const PointType& rSourcePoint, ModelPart::ElementsContainerType& pMasterElements,
            Element::Pointer& pTargetElement, PointType& rLocalTargetPoint, const int& echo_level )
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

        if(echo_level & _WARNING_FOUND_NO_ELEMENT == _WARNING_FOUND_NO_ELEMENT)
            std::cout << " !!!! WARNING: NO ELEMENT FOUND TO CONTAIN " << rSourcePoint << " !!!! " << std::endl;

        return false;
    }


    /// Find an element in pMasterElements contains rSourcePoint and assign it to pTargetElement. The rLocalTargetPoint is the local point in pTargetElement of rSourcePoint
    /// REMARK: we should disable the move mesh flag if we want to search in the reference configuration
    bool SearchPartnerWithBin( const PointType& rSourcePoint, ModelPart::ElementsContainerType& pMasterElements,
            Element::Pointer& pTargetElement, PointType& rLocalTargetPoint, const int& echo_level ) const
    {
        ModelPart::ElementsContainerType pMasterElementsCandidates;

        // get the containing elements from the bin
        int ix = (mDx != 0.0) ? (int) floor(rSourcePoint.X()) / mDx : 0;
        int iy = (mDy != 0.0) ? (int) floor(rSourcePoint.Y()) / mDy : 0;
        int iz = (mDz != 0.0) ? (int) floor(rSourcePoint.Z()) / mDz : 0;

        SpatialKey key(ix, iy, iz);
        std::map<SpatialKey, std::set<std::size_t> >::const_iterator it_bin_elements = mBinElements.find(key);

        if(it_bin_elements != mBinElements.end())
        {
            for(std::set<std::size_t>::const_iterator it = it_bin_elements->second.begin(); it != it_bin_elements->second.end(); ++it )
            {
                Element::GeometryType& r_geom = pMasterElements[*it].GetGeometry();

                bool is_inside = r_geom.IsInside( rSourcePoint, rLocalTargetPoint );
                if( is_inside )
                {
                    pTargetElement = pMasterElements(*it);
                    return true;
                }
            }
        }

        if(echo_level & _WARNING_FOUND_NO_ELEMENT == _WARNING_FOUND_NO_ELEMENT)
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

    struct SpatialKey
    {
        public:
            SpatialKey(const int& ix, const int& iy, const int& iz) : x(ix), y(iy), z(iz) {}
            bool operator<(const SpatialKey& rOther) const
            {
                if(x == rOther.x)
                {
                    if(y == rOther.y)
                    {
                        return z < rOther.z;
                    }
                    else
                        return y < rOther.y;
                }
                else
                    return x < rOther.x;
            }
            const int& kx() const {return x;}
            const int& ky() const {return y;}
            const int& kz() const {return z;}
        private:
            int x, y, z;
    };

    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    double mDx, mDy, mDz;
    std::map<SpatialKey, std::set<std::size_t> > mBinElements;


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
{
    return rIStream;
}

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
