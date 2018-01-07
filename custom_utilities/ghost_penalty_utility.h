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
//  Date:            6 Jan 2018
//


#if !defined(KRATOS_GHOST_PENALTY_UTILITY_H_INCLUDED )
#define  KRATOS_GHOST_PENALTY_UTILITY_H_INCLUDED



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
#include "includes/deprecated_variables.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "custom_algebra/brep.h"
#include "custom_conditions/ghost_penalty_condition.h"


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
/** Search and create the ghost penalty conditions
*/
class GhostPenaltyUtility
{
private:

    static const int msEdgesT3[][2];
    static const int msEdgesT6[][3];
    static const int msEdgesQ4[][2];
    static const int msEdgesQ8[][3];
    static const int msEdgesQ9[][3];
    static const int msFacesT4[][3];
    static const int msFacesT10[][6];
    static const int msFacesH8[][4];
    static const int msFacesH20[][8];
    static const int msFacesH27[][9];

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GhostPenaltyUtility
    KRATOS_CLASS_POINTER_DEFINITION(GhostPenaltyUtility);

    typedef Element::GeometryType GeometryType;

    typedef GeometryType::PointType NodeType;

    typedef NodeType::PointType PointType;

    typedef NodeType::CoordinatesArrayType CoordinatesArrayType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GhostPenaltyUtility() {}

    /// Destructor.
    virtual ~GhostPenaltyUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Setup the ghost penalty conditions between an element and its neighbour
    static ConditionsContainerType SetUpSurfacePenaltyConditions(Element::Pointer p_element,
        GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t& lastCondId, Properties::Pointer pProperties)
    {
        // firstly obtain all neighbour elements of the current element
        WeakPointerVector<Element>& rNeighbours = p_element->GetValue(NEIGHBOUR_ELEMENTS);

        ConditionsContainerType pNewConditions;

        // for each neighbour, find the common edge
        for (std::size_t i = 0; i < rNeighbours.size(); ++i)
        {
            if (rNeighbours[i].Id() != p_element->Id())
            {
                std::pair<GeometryData::KratosGeometryType, std::vector<std::size_t> > edge = FindCommonEdge(*p_element, rNeighbours[i]);

                if (edge.first == GeometryData::Kratos_generic_type)
                    continue;

                KRATOS_WATCH(edge.second.size())
                // create the edge geometry
                typename Element::NodesArrayType temp_nodes;
                for (std::size_t j = 0; j < edge.second.size(); ++j)
                    temp_nodes.push_back(p_element->GetGeometry().pGetPoint(edge.second[j]));

                GeometryType::Pointer p_temp_geometry;

                if (edge.first == GeometryData::Kratos_Line2D2)
                    p_temp_geometry = GeometryType::Pointer(new Line2D2<NodeType>(temp_nodes));
                else if (edge.first == GeometryData::Kratos_Line2D3)
                    p_temp_geometry = GeometryType::Pointer(new Line2D3<NodeType>(temp_nodes));
                else
                    KRATOS_THROW_ERROR(std::logic_error, "Unknown geometry type", edge.first)

                // check if this edge is cut by the brep or totally inside. If yes, then the new ghost condition is created.
                int stat = r_brep.CutStatus(*p_temp_geometry);
                std::cout << "edge " << (*p_temp_geometry)[0].Id() << " " << (*p_temp_geometry)[1].Id() << std::endl;
                KRATOS_WATCH(stat)
                if (stat == BRep::_CUT || stat == BRep::_IN)
                {
                    // create the ghost penalty condition
                    Condition::Pointer pNewCond = p_sample_condition->Create(++lastCondId, p_temp_geometry, p_element, rNeighbours(i).lock(), pProperties);
                    pNewCond->SetValue(IS_INACTIVE, false);
                    pNewCond->Set(ACTIVE, true);
                    pNewConditions.push_back(pNewCond);
                }
            }
        }

        std::cout << __FUNCTION__ << " completed: " << pNewConditions.size() << " new ghost conditions are created" << std::endl;

        return pNewConditions;
    }


    /// Probe the neighbour elements of an element
    static void ProbeNeighbourElements(Element::Pointer p_element)
    {
        WeakPointerVector<Element>& rNeighbours = p_element->GetValue(NEIGHBOUR_ELEMENTS);
        std::cout << "Neighbour elements of element " << p_element->Id() << ":";
        for (std::size_t i = 0; i < rNeighbours.size(); ++i)
        {
            std::cout << " " << rNeighbours[i].Id();
        }
        std::cout << std::endl;
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
        return "Ghost Penalty Utility";
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


    /// Find the common edge between the two geometries, the index of the nodes on the edge are the local index w.r.t geometry 1
    static std::pair<GeometryData::KratosGeometryType, std::vector<std::size_t> > FindCommonEdge(Element& r_elem_1, Element& r_elem_2)
    {
        if (r_elem_1.GetGeometry().GetGeometryType() != r_elem_2.GetGeometry().GetGeometryType())
            KRATOS_THROW_ERROR(std::logic_error, "The geometry type is not the same", "")

        if (r_elem_1.GetGeometry().GetGeometryType() == GeometryData::Kratos_Triangle2D3)
        {
            return std::make_pair(GeometryData::Kratos_Line2D2, FindCommonEdgeT3(r_elem_1, r_elem_2));
        }
        else if (r_elem_1.GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4)
        {
            return std::make_pair(GeometryData::Kratos_Line2D2, FindCommonEdgeQ4(r_elem_1, r_elem_2));
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")

        return std::make_pair(GeometryData::Kratos_generic_type, std::vector<std::size_t>{});
    }

    static std::vector<std::size_t> FindCommonEdgeT3(Element& r_elem_1, Element& r_elem_2)
    {
        for (std::size_t i = 0; i < 3; ++i)
        {
            const std::size_t& i1 = msEdgesT3[i][0];
            const std::size_t& i2 = msEdgesT3[i][1];
            const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
            const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();

            if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n2, r_elem_2.GetGeometry()))
                return std::vector<std::size_t>{i1, i2};
        }
        return std::vector<std::size_t>{};
    }

    static std::vector<std::size_t> FindCommonEdgeT6(Element& r_elem_1, Element& r_elem_2)
    {
        for (std::size_t i = 0; i < 3; ++i)
        {
            const std::size_t& i1 = msEdgesT6[i][0];
            const std::size_t& i2 = msEdgesT6[i][1];
            const std::size_t& i3 = msEdgesT6[i][2];
            const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
            const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
            const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();

            if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n3, r_elem_2.GetGeometry()))
                return std::vector<std::size_t>{i1, i2, i3};
        }
        return std::vector<std::size_t>{};
    }

    static std::vector<std::size_t> FindCommonEdgeQ4(Element& r_elem_1, Element& r_elem_2)
    {
        for (std::size_t i = 0; i < 4; ++i)
        {
            const std::size_t& i1 = msEdgesQ4[i][0];
            const std::size_t& i2 = msEdgesQ4[i][1];
            const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
            const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();

            if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n2, r_elem_2.GetGeometry()))
                return std::vector<std::size_t>{i1, i2};
        }
        return std::vector<std::size_t>{};
    }

    static std::vector<std::size_t> FindCommonEdgeQ8(Element& r_elem_1, Element& r_elem_2)
    {
        for (std::size_t i = 0; i < 4; ++i)
        {
            const std::size_t& i1 = msEdgesQ8[i][0];
            const std::size_t& i2 = msEdgesQ8[i][1];
            const std::size_t& i3 = msEdgesQ8[i][2];
            const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
            const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
            const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();

            if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n3, r_elem_2.GetGeometry()))
                return std::vector<std::size_t>{i1, i2, i3};
        }
        return std::vector<std::size_t>{};
    }

    static std::vector<std::size_t> FindCommonFaceT4(Element& r_elem_1, Element& r_elem_2)
    {
        for (std::size_t i = 0; i < 4; ++i)
        {
            const std::size_t& i1 = msFacesT4[i][0];
            const std::size_t& i2 = msFacesT4[i][1];
            const std::size_t& i3 = msFacesT4[i][2];
            const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
            const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
            const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();

            if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n3, r_elem_2.GetGeometry()))
                return std::vector<std::size_t>{i1, i2, i3};
        }
        return std::vector<std::size_t>{};
    }

    static std::vector<std::size_t> FindCommonFaceT10(Element& r_elem_1, Element& r_elem_2)
    {
        for (std::size_t i = 0; i < 4; ++i)
        {
            const std::size_t& i1 = msFacesT10[i][0];
            const std::size_t& i2 = msFacesT10[i][1];
            const std::size_t& i3 = msFacesT10[i][2];
            const std::size_t& i4 = msFacesT10[i][3];
            const std::size_t& i5 = msFacesT10[i][4];
            const std::size_t& i6 = msFacesT10[i][5];
            const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
            const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
            const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();
            const std::size_t& n4 = r_elem_1.GetGeometry()[i4].Id();
            const std::size_t& n5 = r_elem_1.GetGeometry()[i5].Id();
            const std::size_t& n6 = r_elem_1.GetGeometry()[i6].Id();

            if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n3, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n4, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n5, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n6, r_elem_2.GetGeometry()))
                return std::vector<std::size_t>{i1, i2, i3, i4, i5, i6};
        }
        return std::vector<std::size_t>{};
    }

    static std::vector<std::size_t> FindCommonFaceH8(Element& r_elem_1, Element& r_elem_2)
    {
        for (std::size_t i = 0; i < 8; ++i)
        {
            const std::size_t& i1 = msFacesH8[i][0];
            const std::size_t& i2 = msFacesH8[i][1];
            const std::size_t& i3 = msFacesH8[i][2];
            const std::size_t& i4 = msFacesH8[i][3];
            const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
            const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
            const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();
            const std::size_t& n4 = r_elem_1.GetGeometry()[i4].Id();

            if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n3, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n4, r_elem_2.GetGeometry()))
                return std::vector<std::size_t>{i1, i2, i3, i4};
        }
        return std::vector<std::size_t>{};
    }

    static std::vector<std::size_t> FindCommonFaceH20(Element& r_elem_1, Element& r_elem_2)
    {
        for (std::size_t i = 0; i < 8; ++i)
        {
            const std::size_t& i1 = msFacesH20[i][0];
            const std::size_t& i2 = msFacesH20[i][1];
            const std::size_t& i3 = msFacesH20[i][2];
            const std::size_t& i4 = msFacesH20[i][3];
            const std::size_t& i5 = msFacesH20[i][4];
            const std::size_t& i6 = msFacesH20[i][5];
            const std::size_t& i7 = msFacesH20[i][6];
            const std::size_t& i8 = msFacesH20[i][7];
            const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
            const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
            const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();
            const std::size_t& n4 = r_elem_1.GetGeometry()[i4].Id();
            const std::size_t& n5 = r_elem_1.GetGeometry()[i5].Id();
            const std::size_t& n6 = r_elem_1.GetGeometry()[i6].Id();
            const std::size_t& n7 = r_elem_1.GetGeometry()[i7].Id();
            const std::size_t& n8 = r_elem_1.GetGeometry()[i8].Id();

            if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n3, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n4, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n5, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n6, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n7, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n8, r_elem_2.GetGeometry()))
                return std::vector<std::size_t>{i1, i2, i3, i4, i5, i6, i7, i8};
        }
        return std::vector<std::size_t>{};
    }

    static std::vector<std::size_t> FindCommonFaceH27(Element& r_elem_1, Element& r_elem_2)
    {
        for (std::size_t i = 0; i < 8; ++i)
        {
            const std::size_t& i1 = msFacesH27[i][0];
            const std::size_t& i2 = msFacesH27[i][1];
            const std::size_t& i3 = msFacesH27[i][2];
            const std::size_t& i4 = msFacesH27[i][3];
            const std::size_t& i5 = msFacesH27[i][4];
            const std::size_t& i6 = msFacesH27[i][5];
            const std::size_t& i7 = msFacesH27[i][6];
            const std::size_t& i8 = msFacesH27[i][7];
            const std::size_t& i9 = msFacesH27[i][8];
            const std::size_t& n1 = r_elem_1.GetGeometry()[i1].Id();
            const std::size_t& n2 = r_elem_1.GetGeometry()[i2].Id();
            const std::size_t& n3 = r_elem_1.GetGeometry()[i3].Id();
            const std::size_t& n4 = r_elem_1.GetGeometry()[i4].Id();
            const std::size_t& n5 = r_elem_1.GetGeometry()[i5].Id();
            const std::size_t& n6 = r_elem_1.GetGeometry()[i6].Id();
            const std::size_t& n7 = r_elem_1.GetGeometry()[i7].Id();
            const std::size_t& n8 = r_elem_1.GetGeometry()[i8].Id();
            const std::size_t& n9 = r_elem_1.GetGeometry()[i9].Id();

            if (IsBelongedToGeometry(n1, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n2, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n3, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n4, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n5, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n6, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n7, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n8, r_elem_2.GetGeometry())
                && IsBelongedToGeometry(n9, r_elem_2.GetGeometry()))
                return std::vector<std::size_t>{i1, i2, i3, i4, i5, i6, i7, i8, i9};
        }
        return std::vector<std::size_t>{};
    }

    static bool IsBelongedToGeometry(const std::size_t& node_id, GeometryType& r_geom)
    {
        for (std::size_t i = 0; i < r_geom.size(); ++i)
            if (r_geom[i].Id() == node_id)
                return true;
        return false;
    }

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
    GhostPenaltyUtility& operator=(GhostPenaltyUtility const& rOther);

    /// Copy constructor.
    GhostPenaltyUtility(GhostPenaltyUtility const& rOther);


    ///@}

}; // Class GhostPenaltyUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, GhostPenaltyUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const GhostPenaltyUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GHOST_PENALTY_UTILITY_H_INCLUDED  defined
