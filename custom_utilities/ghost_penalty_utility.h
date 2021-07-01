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
#include "brep_application/custom_algebra/brep.h"
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

/**
 * Basic search utilities
 */
struct GhostPenalty_Helper
{
    typedef Element::GeometryType GeometryType;
    typedef GeometryType::IntegrationPointType IntegrationPointType;
    typedef GeometryType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    /// Find out if the edge geometry is on the element geometry. If yes, build the nodal id map from edge geometry to element geometry
    static bool BuildMapEdgeNodeIndexToElementNodeIndex(std::map<std::size_t, std::size_t>& map_edge_node_index_to_element_node_index,
        GeometryType& r_edge_geometry, GeometryType& r_element_geometry, const int* nodes_on_edge);

    /// Check if the nodes in the geometry is the same as given in the list
    static bool IsSame(GeometryType& r_geom, const std::vector<std::size_t>& nodes);

    /// Check if the nodes on the edge is the same as nodes of the element given by an id list
    static bool IsSame(GeometryType& r_edge_geometry, GeometryType& r_element_geometry, const int* nodes_on_edge);

    // Check if the node id belongs to the geometry
    static bool IsBelongedToGeometry(const std::size_t& node_id, Element::GeometryType& r_geom);

    /// Compute the local frame at point on edge geometry
    static void ComputeLocalFrame(Vector& t, Vector& n, GeometryType& r_edge_geometry, const IntegrationPointType& integration_point);

    /// Compute the local frame at point on face geometry
    static void ComputeLocalFrame(Vector& t1, Vector& t2, Vector& n, GeometryType& r_face_geometry, const IntegrationPointType& integration_point);

    /// Compute the integration point on side of triangle based on the sample integration point on edge
    static IntegrationPointType ComputeIntegrationPointOnTriSide(const int& side, const IntegrationPointType& edge_integration_point);

    /// Compute the integration point on side of quadrilateral based on the sample integration point on edge
    static IntegrationPointType ComputeIntegrationPointOnQuadSide(const int& side, const IntegrationPointType& edge_integration_point);

    /// Compute the integration point on side of tetrahedra based on the sample integration point on edge
    static IntegrationPointType ComputeIntegrationPointOnTetSide(const int& side, const IntegrationPointType& face_integration_point);

    /// Compute the integration point on side of tetrahedra based on the sample integration point on edge
    static IntegrationPointType ComputeIntegrationPointOnHexSide(const int& side, const IntegrationPointType& face_integration_point);

    /// Compute the shape function, w.r.t global coordinates of the geometry, in the reference frame
    static Vector& ComputeShapeFunction(Vector& DN_DX, GeometryType& r_element_geometry, const IntegrationPointType& integration_point);

    /// Compute the shape function gradient, w.r.t global coordinates of the geometry, in the reference frame
    static Matrix& ComputeShapeFunctionGradient(Matrix& DN_DX, GeometryType& r_element_geometry, const IntegrationPointType& integration_point);

    /// Compute the shape function second derivatives, w.r.t global coordinates of the geometry, in the reference frame
    /// REMARKS: the arrangement is:
    /// in 2D: [d^2/dx^2 d^2/dy^2 d^2/dxdy]
    /// in 3D: [d^2/dx^2 d^2/dy^2 d^2/dz^2 d^2/dxdy d^2/dydz d^2/dxdz]
    static std::vector<Vector>& ComputeShapeFunctionSecondDerivatives(std::vector<Vector>& D2N_DX2, GeometryType& r_element_geometry, const IntegrationPointType& integration_point);
};

/**
 * Basic geometry subroutines
 */
template<GeometryData::KratosGeometryType TGeometryType>
struct GhostPenalty_Geometry_Helper
{
    typedef Element::GeometryType GeometryType;
    typedef GeometryType::IntegrationPointType IntegrationPointType;
    typedef GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /// Get the number of sides associated with a geometry
    static int NumberOfSides()
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not defined")
    }

    /// Find the side of the edge geometry on the element geometry. If not found, return -1.
    static int FindSide(GeometryType& r_element_geometry, GeometryType& r_edge_geometry)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not defined")
    }

    /// Find the common edge/face between two geometries
    static std::vector<std::size_t> FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not defined")
    }

    /// Get the local index of nodes on the side of the geometry
    static const int* Faces(const std::size_t& side)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not defined")
    }

    /// Compute the integration point in the element given the integration point on the edge and the side
    static IntegrationPointType ComputeIntegrationPoint(const int& side, const IntegrationPointType& edge_integration_point)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not defined")
    }

    /// Compute the shape function value along the edge of the elemental geometry, in the reference frame
    static void ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not defined")
    }

    /// Compute the shape function local gradient along the edge of the elemental geometry, in the reference frame
    static void ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not defined")
    }

    /// Compute the shape function local gradient in the normal direction along the edge of the elemental geometry, in the reference frame
    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not defined")
    }

    /// Compute the shape function second derivatives in the normal direction along the edge of the elemental geometry, in the reference frame
    static void ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not defined")
    }
};

/// Short class definition.
/** Search and create the ghost penalty conditions
*/
class GhostPenaltyUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GhostPenaltyUtility
    KRATOS_CLASS_POINTER_DEFINITION(GhostPenaltyUtility);

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
    GhostPenaltyUtility() {}

    /// Destructor.
    virtual ~GhostPenaltyUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Setup the ghost penalty conditions for the model_part
    /*static*/ ConditionsContainerType SetUpSurfacePenaltyConditions(ModelPart& r_model_part,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t& lastCondId, Properties::Pointer pProperties, const int& echo_level) const;

    /// Setup the ghost penalty conditions for the model_part on a subset of elements
    /*static*/ ConditionsContainerType SetUpSurfacePenaltyConditions(ModelPart& r_model_part,
        ModelPart::ElementsContainerType& pElements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t& lastCondId, Properties::Pointer pProperties, const int& echo_level) const;

    /// Setup the ghost penalty conditions for the model_part on a subset of elements
    /*static*/ ConditionsContainerType SetUpSurfacePenaltyConditions(ModelPart& r_model_part,
        ModelPart::ElementsContainerType& pElements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t& lastCondId, Properties::Pointer pProperties,
        const int& nsampling, const int& configuration,
        const int& echo_level) const;

    /// Setup the ghost penalty conditions between an element and its neighbour
    /// Note: user must call FindElementalNeighboursProcess(model.model_part, 2, 10).Execute() to setup first the neighbour elements
    /*static*/ ConditionsContainerType SetUpSurfacePenaltyConditions(Element::Pointer p_element,
        GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t& lastCondId, Properties::Pointer pProperties, const int& echo_level) const;

    /// Setup the ghost penalty conditions between two elements. The BRep is to check the validity of the ghost penalty condition
    /*static*/ virtual Condition::Pointer SetUpSurfacePenaltyCondition(Element::Pointer p_element_1,
        Element::Pointer p_element_2, GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t& lastCondId, Properties::Pointer pProperties, const int& echo_level) const;

    /// Setup the ghost penalty conditions between two elements. The BRep is to check the validity of the ghost penalty condition
    /*static*/ virtual Condition::Pointer SetUpSurfacePenaltyCondition(Element::Pointer p_element_1,
        Element::Pointer p_element_2, GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t& lastCondId, Properties::Pointer pProperties,
        const int& nsampling, const int& configuration, const int& echo_level) const;

    /// Merely setup the ghost penalty conditions between two elements. This is useful to debug the ghost penalties.
    /*static*/ Condition::Pointer SetUpSurfacePenaltyCondition(Element::Pointer p_element_1,
        Element::Pointer p_element_2, GhostPenaltyCondition::Pointer p_sample_condition,
        std::size_t& lastCondId, Properties::Pointer pProperties, const int& echo_level) const;

    /// Compute the shape function, the edge geometry must be on the edge of the element.
    static void ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    /// Compute the shape function gradient, the edge geometry must be on the edge of the element.
    static void ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    /// Compute the shape function gradient in the normal direction, the edge geometry must be on the edge of the element.
    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    /// Compute the shape function gradient in the normal direction, the edge geometry must be on the edge of the element.
    static void ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    /// Find the side of the edge geometry on the element geometry. If not found, return -1.
    static int FindSide(GeometryType& r_element_geometry, GeometryType& r_edge_geometry);

    /// Get the local index of nodes on the side of the geometry
    static const int* Faces(GeometryType& r_element_geometry, const std::size_t& side);

    /// Find the common edge between the two geometries, the index of the nodes on the edge are the local index w.r.t geometry 1
    static std::pair<GeometryData::KratosGeometryType, std::vector<std::size_t> > FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2);

    /// Compute the integration point in the element given the integration point on the edge and the side
    static IntegrationPointType ComputeIntegrationPoint(GeometryType& r_element_geometry, const int& side, const IntegrationPointType& edge_integration_point);

    /// Probe the neighbour elements of an element
    static void ProbeNeighbourElements(Element::Pointer p_element);

    /// Compute and print out the second derivatives of the geometry, in order for debugging
    static void ProbeShapeFunctionSecondDerivatives(GeometryType& r_geom);

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

///@name Traits Definitions
///@{


/**
 * Trait class to perform FindSide operation, TClassType shall be GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType>
 */
template<class TClassType>
struct FindSide_Helper
{
    typedef typename TClassType::GeometryType GeometryType;

    static int Execute(GeometryType& r_element_geometry, GeometryType& r_edge_geometry)
    {
        bool found;
        for (int side = 0; side < TClassType::NumberOfSides(); ++side)
        {
            found = GhostPenalty_Helper::IsSame(r_edge_geometry, r_element_geometry, TClassType::Faces(side));

            if (found)
                return side;
        }
        return -1;
    }
};

/**
 * Trait class to perform ComputeShapeFunction operation, TClassType shall be GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType>
 */
template<class TClassType>
struct ComputeShapeFunction_Helper
{
    typedef typename TClassType::GeometryType GeometryType;

    typedef typename GeometryType::IntegrationPointType IntegrationPointType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    static void Execute(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const IntegrationPointsArrayType& edge_integration_points)
    {
        std::map<std::size_t, std::size_t> map_edge_node_index_to_element_node_index;

        int side = FindSide_Helper<TClassType>::Execute(r_element_geometry, r_edge_geometry);
        if (side != -1)
        {
            GhostPenalty_Helper::BuildMapEdgeNodeIndexToElementNodeIndex(map_edge_node_index_to_element_node_index,
                        r_edge_geometry, r_element_geometry, TClassType::Faces(side));

            Vector Ne(TClassType::NumberOfNodes());

            if (N.size1() != r_edge_geometry.size() || N.size2() != edge_integration_points.size())
                N.resize(r_edge_geometry.size(), edge_integration_points.size(), false);

            for (std::size_t i = 0; i < edge_integration_points.size(); ++i)
            {
                IntegrationPointType integration_point = TClassType::ComputeIntegrationPoint(side, edge_integration_points[i]);

                Ne = GhostPenalty_Helper::ComputeShapeFunction(Ne, r_element_geometry, integration_point);
                #ifdef ENABLE_DEBUG_GHOST_PENALTY
                KRATOS_WATCH(Ne)
                #endif
                for (std::size_t j = 0; j < r_edge_geometry.size(); ++j)
                {
                    N(j, i) = Ne(map_edge_node_index_to_element_node_index[j]);
                }
            }
        }
    }
};

/**
 * Trait class to perform ComputeShapeFunctionGradient operation, TClassType shall be GhostPenalty_Geometry_Helper<GeometryData::KratosGeometryType>
 */
template<class TClassType>
struct ComputeShapeFunctionGradient_Helper
{
    typedef typename TClassType::GeometryType GeometryType;

    typedef typename GeometryType::IntegrationPointType IntegrationPointType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    static void Execute(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const IntegrationPointsArrayType& edge_integration_points)
    {
        std::map<std::size_t, std::size_t> map_edge_node_index_to_element_node_index;

        int side = FindSide_Helper<TClassType>::Execute(r_element_geometry, r_edge_geometry);
        if (side != -1)
        {
            GhostPenalty_Helper::BuildMapEdgeNodeIndexToElementNodeIndex(map_edge_node_index_to_element_node_index,
                        r_edge_geometry, r_element_geometry, TClassType::Faces(side));

            const unsigned int dim = r_element_geometry.WorkingSpaceDimension();

            Matrix DN_DXe(TClassType::NumberOfNodes(), dim);

            if (DN_DX.size() != edge_integration_points.size())
                DN_DX.resize(edge_integration_points.size());

            for (std::size_t i = 0; i < edge_integration_points.size(); ++i)
            {
                if (DN_DX[i].size1() != r_edge_geometry.size() || DN_DX[i].size2() != dim)
                    DN_DX[i].resize(r_edge_geometry.size(), dim, false);

                IntegrationPointType integration_point = TClassType::ComputeIntegrationPoint(side, edge_integration_points[i]);

                DN_DXe = GhostPenalty_Helper::ComputeShapeFunctionGradient(DN_DXe, r_element_geometry, integration_point);
                #ifdef ENABLE_DEBUG_GHOST_PENALTY
                KRATOS_WATCH(DN_DXe)
                #endif
                for (std::size_t j = 0; j < r_edge_geometry.size(); ++j)
                {
                    for (std::size_t d = 0; d < dim; ++d)
                        DN_DX[i](j, d) = DN_DXe(map_edge_node_index_to_element_node_index[j], d);
                }
            }
        }
    }
};

template<class TClassType>
struct ComputeShapeFunctionNormalGradient_Helper
{
    typedef typename TClassType::GeometryType GeometryType;

    typedef typename GeometryType::IntegrationPointType IntegrationPointType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    static void Execute2D(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const IntegrationPointsArrayType& edge_integration_points)
    {
        std::map<std::size_t, std::size_t> map_edge_node_index_to_element_node_index;

        #ifdef ENABLE_DEBUG_GHOST_PENALTY
        std::cout << "element geometry:";
        for (std::size_t i = 0; i < r_element_geometry.size(); ++i)
            std::cout << " " << r_element_geometry[i].Id();
        std::cout << std::endl;

        std::cout << "edge geometry:";
        for (std::size_t i = 0; i < r_edge_geometry.size(); ++i)
            std::cout << " " << r_edge_geometry[i].Id();
        std::cout << std::endl;
        #endif

        int side = FindSide_Helper<TClassType>::Execute(r_element_geometry, r_edge_geometry);
        if (side != -1)
        {
            map_edge_node_index_to_element_node_index.clear();

            GhostPenalty_Helper::BuildMapEdgeNodeIndexToElementNodeIndex(map_edge_node_index_to_element_node_index,
                        r_edge_geometry, r_element_geometry, TClassType::Faces(side));

            #ifdef ENABLE_DEBUG_GHOST_PENALTY
            KRATOS_WATCH(side)
            #endif

            Matrix DN_DX(TClassType::NumberOfNodes(), 2);

            // construct the local coordinates system associated with the edge
            Vector t(2), n(2);
    //            double length = norm_2(r_element_geometry[msEdges[side][0]] - r_element_geometry[msEdges[side][1]]);
    //            t(0) = (r_element_geometry[msEdges[side][1]].X0() - r_element_geometry[msEdges[side][0]].X0()) / length;
    //            t(1) = (r_element_geometry[msEdges[side][1]].Y0() - r_element_geometry[msEdges[side][0]].Y0()) / length;
    //            n(0) = -t(1);
    //            n(1) = t(0);
    //            #ifdef ENABLE_DEBUG_GHOST_PENALTY
    //            KRATOS_WATCH(t)
    //            KRATOS_WATCH(n)
    //            #endif

            for (std::size_t i = 0; i < edge_integration_points.size(); ++i)
            {
                GhostPenalty_Helper::ComputeLocalFrame(t, n, r_edge_geometry, edge_integration_points[i]);
                #ifdef ENABLE_DEBUG_GHOST_PENALTY
                KRATOS_WATCH(t)
                KRATOS_WATCH(n)
                #endif

                IntegrationPointType integration_point = TClassType::ComputeIntegrationPoint(side, edge_integration_points[i]);

                DN_DX = GhostPenalty_Helper::ComputeShapeFunctionGradient(DN_DX, r_element_geometry, integration_point);
                #ifdef ENABLE_DEBUG_GHOST_PENALTY
                KRATOS_WATCH(DN_DX)
                #endif
                dNdn(0, i) = DN_DX(map_edge_node_index_to_element_node_index[0], 0) * n(0) + DN_DX(map_edge_node_index_to_element_node_index[0], 1) * n(1);
                dNdn(1, i) = DN_DX(map_edge_node_index_to_element_node_index[1], 0) * n(0) + DN_DX(map_edge_node_index_to_element_node_index[1], 1) * n(1);
            }
        }

        if (side == -1)
            KRATOS_THROW_ERROR(std::logic_error, "The edge geometry is not found on the element geometry. Check the input", "")
    }

    static void Execute3D(Matrix& dNdn, GeometryType& r_element_geometry, GeometryType& r_face_geometry,
        const IntegrationPointsArrayType& face_integration_points)
    {
        std::map<std::size_t, std::size_t> map_edge_node_index_to_element_node_index;

        #ifdef ENABLE_DEBUG_GHOST_PENALTY
        std::cout << "element geometry:";
        for (std::size_t i = 0; i < r_element_geometry.size(); ++i)
            std::cout << " " << r_element_geometry[i].Id();
        std::cout << std::endl;

        std::cout << "face geometry:";
        for (std::size_t i = 0; i < r_face_geometry.size(); ++i)
            std::cout << " " << r_face_geometry[i].Id();
        std::cout << std::endl;
        #endif

        int side = FindSide_Helper<TClassType>::Execute(r_element_geometry, r_face_geometry);
        if (side != -1)
        {
            GhostPenalty_Helper::BuildMapEdgeNodeIndexToElementNodeIndex(map_edge_node_index_to_element_node_index,
                        r_face_geometry, r_element_geometry, TClassType::Faces(side));

            #ifdef ENABLE_DEBUG_GHOST_PENALTY
            KRATOS_WATCH(side)
            #endif

            Matrix DN_DX(TClassType::NumberOfNodes(), 3);

            if (dNdn.size1() != r_face_geometry.size() || dNdn.size2() != face_integration_points.size())
                dNdn.resize(r_face_geometry.size(), face_integration_points.size(), false);

            // construct the local coordinates system associated with the face
            Vector t1(3), t2(3), n(3);

            for (std::size_t i = 0; i < face_integration_points.size(); ++i)
            {
                GhostPenalty_Helper::ComputeLocalFrame(t1, t2, n, r_face_geometry, face_integration_points[i]);
                #ifdef ENABLE_DEBUG_GHOST_PENALTY
                KRATOS_WATCH(t1)
                KRATOS_WATCH(t2)
                KRATOS_WATCH(n)
                #endif

                IntegrationPointType integration_point = TClassType::ComputeIntegrationPoint(side, face_integration_points[i]);

                DN_DX = GhostPenalty_Helper::ComputeShapeFunctionGradient(DN_DX, r_element_geometry, integration_point);
                #ifdef ENABLE_DEBUG_GHOST_PENALTY
                KRATOS_WATCH(DN_DX)
                #endif
                for (std::size_t j = 0; j < r_face_geometry.size(); ++j)
                {
                    dNdn(j, i) = 0.0;
                    for (std::size_t k = 0; k < 3; ++k)
                    {
                        dNdn(j, i) += DN_DX(map_edge_node_index_to_element_node_index[j], k) * n(k);
                    }
                }
            }
        }

        if (side == -1)
            KRATOS_THROW_ERROR(std::logic_error, "The face geometry is not found on the element geometry. Check the input", "")
    }
};

template<>
struct GhostPenalty_Geometry_Helper<GeometryData::Kratos_Triangle2D3>
{
    typedef Element::GeometryType GeometryType;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    static int NumberOfNodes() {return 3;}

    static int NumberOfSides() {return 3;}

    static std::vector<std::size_t> FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2);

    static const int* Faces(const std::size_t& side) {return msEdges[side];}

    static IntegrationPointType ComputeIntegrationPoint(const int& side, const IntegrationPointType& edge_integration_point);

    static void ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static const int msEdges[][2];
};

template<>
struct GhostPenalty_Geometry_Helper<GeometryData::Kratos_Triangle2D6>
{
    typedef Element::GeometryType GeometryType;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    static int NumberOfNodes() {return 6;}

    static int NumberOfSides() {return 3;}

    static std::vector<std::size_t> FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2);

    static const int* Faces(const std::size_t& side) {return msEdges[side];}

    static IntegrationPointType ComputeIntegrationPoint(const int& side, const IntegrationPointType& edge_integration_point);

    static void ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static const int msEdges[][3];
};

template<>
struct GhostPenalty_Geometry_Helper<GeometryData::Kratos_Quadrilateral2D4>
{
    typedef Element::GeometryType GeometryType;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    static int NumberOfNodes() {return 4;}

    static int NumberOfSides() {return 4;}

    static std::vector<std::size_t> FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2);

    static const int* Faces(const std::size_t& side) {return msEdges[side];}

    static IntegrationPointType ComputeIntegrationPoint(const int& side, const IntegrationPointType& edge_integration_point);

    static void ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static const int msEdges[][2];
};

template<>
struct GhostPenalty_Geometry_Helper<GeometryData::Kratos_Quadrilateral2D8>
{
    typedef Element::GeometryType GeometryType;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    static int NumberOfNodes() {return 8;}

    static int NumberOfSides() {return 4;}

    static std::vector<std::size_t> FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2);

    static const int* Faces(const std::size_t& side) {return msEdges[side];}

    static IntegrationPointType ComputeIntegrationPoint(const int& side, const IntegrationPointType& edge_integration_point);

    static void ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static const int msEdges[][3];
};

template<>
struct GhostPenalty_Geometry_Helper<GeometryData::Kratos_Quadrilateral2D9>
{
    typedef Element::GeometryType GeometryType;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    static int NumberOfNodes() {return 9;}

    static int NumberOfSides() {return 4;}

    static std::vector<std::size_t> FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2);

    static const int* Faces(const std::size_t& side) {return msEdges[side];}

    static IntegrationPointType ComputeIntegrationPoint(const int& side, const IntegrationPointType& edge_integration_point);

    static void ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static void ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& edge_integration_points);

    static const int msEdges[][3];
};

template<>
struct GhostPenalty_Geometry_Helper<GeometryData::Kratos_Tetrahedra3D4>
{
    typedef Element::GeometryType GeometryType;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    static int NumberOfNodes() {return 4;}

    static int NumberOfSides() {return 4;}

    static std::vector<std::size_t> FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2);

    static const int* Faces(const std::size_t& side) {return msFaces[side];}

    static IntegrationPointType ComputeIntegrationPoint(const int& side, const IntegrationPointType& edge_integration_point);

    static void ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static const int msFaces[][3];
};

template<>
struct GhostPenalty_Geometry_Helper<GeometryData::Kratos_Tetrahedra3D10>
{
    typedef Element::GeometryType GeometryType;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    static int NumberOfNodes() {return 10;}

    static int NumberOfSides() {return 4;}

    static std::vector<std::size_t> FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2);

    static const int* Faces(const std::size_t& side) {return msFaces[side];}

    static IntegrationPointType ComputeIntegrationPoint(const int& side, const IntegrationPointType& edge_integration_point);

    static void ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static const int msFaces[][6];
};

template<>
struct GhostPenalty_Geometry_Helper<GeometryData::Kratos_Hexahedra3D8>
{
    typedef Element::GeometryType GeometryType;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    static int NumberOfNodes() {return 8;}

    static int NumberOfSides() {return 6;}

    static std::vector<std::size_t> FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2);

    static const int* Faces(const std::size_t& side) {return msFaces[side];}

    static IntegrationPointType ComputeIntegrationPoint(const int& side, const IntegrationPointType& edge_integration_point);

    static void ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static const int msFaces[][4];
};

template<>
struct GhostPenalty_Geometry_Helper<GeometryData::Kratos_Hexahedra3D20>
{
    typedef Element::GeometryType GeometryType;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    static int NumberOfNodes() {return 20;}

    static int NumberOfSides() {return 6;}

    static std::vector<std::size_t> FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2);

    static const int* Faces(const std::size_t& side) {return msFaces[side];}

    static IntegrationPointType ComputeIntegrationPoint(const int& side, const IntegrationPointType& edge_integration_point);

    static void ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static const int msFaces[][8];
};

template<>
struct GhostPenalty_Geometry_Helper<GeometryData::Kratos_Hexahedra3D27>
{
    typedef Element::GeometryType GeometryType;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    static int NumberOfNodes() {return 27;}

    static int NumberOfSides() {return 6;}

    static std::vector<std::size_t> FindCommonFace(GeometryType& r_geom_1, GeometryType& r_geom_2);

    static const int* Faces(const std::size_t& side) {return msFaces[side];}

    static IntegrationPointType ComputeIntegrationPoint(const int& side, const IntegrationPointType& edge_integration_point);

    static void ComputeShapeFunction(Matrix& N, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionGradient(std::vector<Matrix>& DN_DX, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static void ComputeShapeFunctionNormalSecondDerivatives(Matrix& d2Ndn2, GeometryType& r_element_geometry,
        GeometryType& r_face_geometry, const GeometryType::IntegrationPointsArrayType& face_integration_points);

    static const int msFaces[][9];
};


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
