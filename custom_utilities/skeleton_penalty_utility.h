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
//  Date:            6 Feb 2018
//


#if !defined(KRATOS_SKELETON_PENALTY_UTILITY_H_INCLUDED )
#define  KRATOS_SKELETON_PENALTY_UTILITY_H_INCLUDED



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
#include "custom_utilities/ghost_penalty_utility.h"


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
/** Search and create the skeleton penalty conditions
 * REF: Hoang Tuong paper
*/
class SkeletonPenaltyUtility : public GhostPenaltyUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SkeletonPenaltyUtility
    KRATOS_CLASS_POINTER_DEFINITION(SkeletonPenaltyUtility);

    typedef Element::GeometryType GeometryType;

    typedef GeometryType::PointType NodeType;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    typedef GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef NodeType::PointType PointType;

    typedef NodeType::CoordinatesArrayType CoordinatesArrayType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SkeletonPenaltyUtility() {}

    /// Destructor.
    virtual ~SkeletonPenaltyUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Merely setup the ghost penalty conditions between two elements. This is useful to debug the ghost penalties.
    /*static*/ virtual Condition::Pointer SetUpSurfacePenaltyCondition(Element::Pointer p_element_1,
        Element::Pointer p_element_2, GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t& lastCondId, Properties::Pointer pProperties)
    {
        Condition::Pointer pNewCond = NULL;

        std::pair<GeometryData::KratosGeometryType, std::vector<std::size_t> > edge = GhostPenaltyUtility::FindCommonFace(p_element_1->GetGeometry(), p_element_2->GetGeometry());

        if (edge.first == GeometryData::KratosGeometryType::Kratos_generic_type)
            return pNewCond;

        // create the edge geometry
        typename Element::NodesArrayType temp_nodes;
        for (std::size_t j = 0; j < edge.second.size(); ++j)
            temp_nodes.push_back(p_element_1->GetGeometry().pGetPoint(edge.second[j]));

        GeometryType::Pointer p_temp_geometry;

        if (edge.first == GeometryData::KratosGeometryType::Kratos_Line2D2)
            p_temp_geometry = GeometryType::Pointer(new Line2D2<NodeType>(temp_nodes));
        else if (edge.first == GeometryData::KratosGeometryType::Kratos_Line2D3)
            p_temp_geometry = GeometryType::Pointer(new Line2D3<NodeType>(temp_nodes));
        else if (edge.first == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
            p_temp_geometry = GeometryType::Pointer(new Triangle3D3<NodeType>(temp_nodes));
        else if (edge.first == GeometryData::KratosGeometryType::Kratos_Triangle3D6)
            p_temp_geometry = GeometryType::Pointer(new Triangle3D6<NodeType>(temp_nodes));
        else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)
            p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D4<NodeType>(temp_nodes));
        else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8)
            p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D8<NodeType>(temp_nodes));
        else if (edge.first == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9)
            p_temp_geometry = GeometryType::Pointer(new Quadrilateral3D9<NodeType>(temp_nodes));
        else
            KRATOS_THROW_ERROR(std::logic_error, "Unknown geometry type", static_cast<int>(edge.first))

        // check if this edge is cut by the brep or totally inside. If yes, then the new ghost condition is created.
        bool is_ghost = false;

        int configuration = 0; // check the cut status in the initial configuration
        int stat1 = r_brep.CutStatus(p_element_1->GetGeometry(), configuration);
        int stat2 = r_brep.CutStatus(p_element_2->GetGeometry(), configuration);

        if (stat1 == BRep::_CUT)
        {
            if (stat2 == BRep::_CUT || stat2 == BRep::_IN)
                is_ghost = true;
        }
        else if (stat1 == BRep::_IN)
        {
            if (stat2 == BRep::_IN || stat2 == BRep::_CUT)
                is_ghost = true;
        }

        if (is_ghost)
        {
            // create the ghost penalty condition
            pNewCond = p_sample_condition->Create(++lastCondId, p_temp_geometry, p_element_1, p_element_2, pProperties);
            pNewCond->SetValue(IS_INACTIVE, false);
            pNewCond->Set(ACTIVE, true);
        }

        return pNewCond;
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
        return "Skeleton Penalty Utility";
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
    static std::pair<GeometryData::KratosGeometryType, std::vector<std::size_t> > FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2);

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
    SkeletonPenaltyUtility& operator=(SkeletonPenaltyUtility const& rOther);

    /// Copy constructor.
    SkeletonPenaltyUtility(SkeletonPenaltyUtility const& rOther);


    ///@}

}; // Class SkeletonPenaltyUtility

///@}

///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, SkeletonPenaltyUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const SkeletonPenaltyUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SKELETON_PENALTY_UTILITY_H_INCLUDED  defined
