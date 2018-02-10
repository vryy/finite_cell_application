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


#if !defined(KRATOS_FINITE_CELL_APPLICATION_QUAD_TREE_NODE_H_INCLUDED )
#define KRATOS_FINITE_CELL_APPLICATION_QUAD_TREE_NODE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/model_part.h"
#include "includes/geometrical_object.h"
#include "utilities/math_utils.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "integration/quadrilateral_gauss_lobatto_integration_points.h"
#include "integration/quadrilateral_gauss_legendre_integration_points.h"
#include "integration/hexahedron_gauss_legendre_integration_points.h"
#include "integration/hexahedron_gauss_lobatto_integration_points.h"
#include "custom_algebra/brep.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/function/monomial_function.h"
#include "custom_algebra/function/heaviside_function.h"
#include "custom_algebra/function/product_function.h"
#include "custom_geometries/finite_cell_geometry.h"
#include "custom_utilities/quadrature_utility.h"


namespace Kratos
{


/// Abstract quad tree node description in reference coordinates
class QuadTreeNode
{
public:

    enum FrameType
    {
        LOCAL = 0,
        GLOBAL = 1
    };

    KRATOS_CLASS_POINTER_DEFINITION(QuadTreeNode);

    typedef typename GeometricalObject::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    /// Default constructor
    QuadTreeNode() {}

    /// Destructor
    virtual ~QuadTreeNode() {}

    /// Return the number of children.
    std::size_t Size() const
    {
        return mpChildren.size();
    }

    /// Check if the quadtree node is a leaf node or not.
    bool IsLeaf() const
    {
        return (Size() == 0);
    }

    /// Clear the list of children
    void Clear()
    {
        mpChildren.clear();
    }

    /*****************************************************************/
    /******* CONSTRUCTION ********************************************/
    /*****************************************************************/

    /// Refine this quadtree unconditionally.
    virtual void Refine()
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Refine this quadtree by checking cut status against a brep. The operated geometry is given by pParentGeometry.
    /// This routine tests the inside/outside against a set of sampling points. It's more computationally expensive but is more robust and applicable to a large types of geometries.
    void RefineBySampling(GeometryType::Pointer pParentGeometry, const BRep& r_brep, const std::size_t& nsampling)
    {
        if(this->IsLeaf())
        {
            std::vector<PointType> SamplingPoints;
            this->CreateSamplingPoints(SamplingPoints, *pParentGeometry, nsampling);

            int stat = r_brep.CutStatus(SamplingPoints);
            if(stat == BRep::_CUT)
            {
                this->Refine();
            }
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
            {
                mpChildren[i]->RefineBySampling(pParentGeometry, r_brep, nsampling);
            }
        }
    }

    /*****************************************************************/
    /******* COMPUTATION *********************************************/
    /*****************************************************************/

    /// Integrate a function using the sample geometry and integration rule.
    /// The caller has to manually set rOutput to zero before calling this function.
    /// if Frame=0, the function is assumed to be defined in the local frame.
    /// if Frame=1, the function is assumed to be defined in the global frame. If the function is defined in the local frame (i.e. shape function), the Integrate<0> method must be used.
    template<typename TOutputType, int Frame>
    void Integrate(GeometryType::Pointer pParentGeometry,
            const Function<array_1d<double, 3>, TOutputType>& rFunc,
            TOutputType& rOutput,
            const int& integration_method) const
    {
        int quadrature_type = QuadratureUtility::GetQuadratureType(integration_method);

        if(quadrature_type == 0) // Kratos built-in
        {
            int quadrature_order = QuadratureUtility::GetQuadratureOrder(integration_method);
            GeometryData::IntegrationMethod ThisIntegrationMethod = Function<double, double>::GetIntegrationMethod(quadrature_order);
            this->Integrate<TOutputType, Frame>(pParentGeometry, rFunc, rOutput, ThisIntegrationMethod);
        }
        else // use custom quadrature
        {
            int quadrature_order = QuadratureUtility::GetQuadratureOrder(integration_method);
            GeometryType::IntegrationPointsArrayType sample_integration_points = this->ConstructCustomQuadrature(quadrature_type, quadrature_order);
            this->Integrate<TOutputType, GeometryType::IntegrationPointsArrayType, Frame>(pParentGeometry, rFunc, rOutput, sample_integration_points);
        }
    }

    /// Integrate a function using the sample geometry and the built-in Kratos integration rule (mainly Gauss Legendre for not-so-high order).
    /// The caller has to manually set rOutput to zero before calling this function.
    /// if Frame=0, the function is assumed to be defined in the local frame.
    /// if Frame=1, the function is assumed to be defined in the global frame. If the function is defined in the local frame (i.e. shape function), the Integrate<0> method must be used.
    template<typename TOutputType, int Frame>
    void Integrate(GeometryType::Pointer pParentGeometry,
            const Function<array_1d<double, 3>, TOutputType>& rFunc,
            TOutputType& rOutput,
            const GeometryData::IntegrationMethod& ThisIntegrationMethod) const
    {
        GeometryType::IntegrationPointsArrayType integration_points;

        this->ConstructQuadrature(pParentGeometry, integration_points, ThisIntegrationMethod);

        std::vector<double> DetJ;
        Function<double, double>::ComputeDetJ(DetJ, *pParentGeometry, integration_points);

        CoordinatesArrayType GlobalCoords;

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            if(Frame == LOCAL)
            {
                rOutput += rFunc.GetValue(integration_points[point]) * DetJ[point] * integration_points[point].Weight();
            }
            else if(Frame == GLOBAL)
            {
                pParentGeometry->GlobalCoordinates(GlobalCoords, integration_points[point]);
                rOutput += rFunc.GetValue(GlobalCoords) * DetJ[point] * integration_points[point].Weight();
            }
        }
    }

    /// Integrate a function using the sample geometry and a set of quadrature points. This allows for very
    /// high order quadrature rule, as well as arbitrary quadrature.
    /// The caller has to manually set rOutput to zero before calling this function.
    /// if Frame=0, the function is assumed to be defined in the local frame.
    /// if Frame=1, the function is assumed to be defined in the global frame. If the function is defined in the local frame (i.e. shape function), the Integrate<0> method must be used.
    template<typename TOutputType, typename TIntegrationPointsArrayType, int Frame>
    void Integrate(GeometryType::Pointer pParentGeometry,
            const Function<array_1d<double, 3>, TOutputType>& rFunc,
            TOutputType& rOutput,
            const TIntegrationPointsArrayType& sample_integration_points) const
    {
        GeometryType::IntegrationPointsArrayType integration_points;

        this->ConstructQuadrature(pParentGeometry, integration_points, sample_integration_points);

        std::vector<double> DetJ;
        Function<double, double>::ComputeDetJ(DetJ, *pParentGeometry, integration_points);

        CoordinatesArrayType GlobalCoords;

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            if(Frame == LOCAL)
            {
                rOutput += rFunc.GetValue(integration_points[point]) * DetJ[point] * integration_points[point].Weight();
            }
            else if(Frame == GLOBAL)
            {
                pParentGeometry->GlobalCoordinates(GlobalCoords, integration_points[point]);
                rOutput += rFunc.GetValue(GlobalCoords) * DetJ[point] * integration_points[point].Weight();
            }
        }
    }

    /// Integrate a function on the domain limited by a BRep using the sample geometry and integration rule.
    /// The caller has to manually set rOutput to zero before calling this function.
    /// if Frame=0, the function is assumed to be defined in the local frame.
    /// if Frame=1, the function is assumed to be defined in the global frame. If the function is defined in the local frame (i.e. shape function), the Integrate<0> method must be used.
    template<typename TOutputType, int Frame>
    void Integrate(GeometryType::Pointer pParentGeometry,
            const Function<array_1d<double, 3>, TOutputType>& rFunc,
            const BRep& r_brep,
            TOutputType& rOutput,
            const int& integration_method,
            const double& small_weight) const
    {
        int quadrature_type = QuadratureUtility::GetQuadratureType(integration_method);

        if(quadrature_type == 0) // Kratos built-in
        {
            int quadrature_order = QuadratureUtility::GetQuadratureOrder(integration_method);
            GeometryData::IntegrationMethod ThisIntegrationMethod = Function<double, double>::GetIntegrationMethod(quadrature_order);
            this->Integrate<TOutputType, Frame>(pParentGeometry, rFunc, r_brep, rOutput, ThisIntegrationMethod, small_weight);
        }
        else // use custom quadrature
        {
            int quadrature_order = QuadratureUtility::GetQuadratureOrder(integration_method);
            GeometryType::IntegrationPointsArrayType sample_integration_points = this->ConstructCustomQuadrature(quadrature_type, quadrature_order);
            this->Integrate<TOutputType, GeometryType::IntegrationPointsArrayType, Frame>(pParentGeometry, rFunc, r_brep, rOutput, sample_integration_points, small_weight);
        }
    }

    /// Integrate a function on the domain limited by a BRep using the sample geometry and integration rule (mainly Gauss Legendre for not-so-high order).
    /// The caller has to manually set rOutput to zero before calling this function.
    /// if Frame=0, the function is assumed to be defined in the local frame.
    /// if Frame=1, the function is assumed to be defined in the global frame. If the function is defined in the local frame (i.e. shape function), the Integrate<0> method must be used.
    template<typename TOutputType, int Frame>
    void Integrate(GeometryType::Pointer pParentGeometry,
            const Function<array_1d<double, 3>, TOutputType>& rFunc,
            const BRep& r_brep,
            TOutputType& rOutput,
            const GeometryData::IntegrationMethod& ThisIntegrationMethod,
            const double& small_weight) const
    {
        GeometryType::IntegrationPointsArrayType integration_points;

        this->ConstructQuadrature(pParentGeometry, integration_points, ThisIntegrationMethod);

        std::vector<double> DetJ;
        Function<double, double>::ComputeDetJ(DetJ, *pParentGeometry, integration_points);

        CoordinatesArrayType GlobalCoords;

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            pParentGeometry->GlobalCoordinates(GlobalCoords, integration_points[point]);
            if(r_brep.IsInside(GlobalCoords))
            {
                if(Frame == LOCAL)
                    rOutput += rFunc.GetValue(integration_points[point]) * DetJ[point] * integration_points[point].Weight();
                else if(Frame == GLOBAL)
                    rOutput += rFunc.GetValue(GlobalCoords) * DetJ[point] * integration_points[point].Weight();
            }
            else
            {
                if(Frame == LOCAL)
                    rOutput += rFunc.GetValue(integration_points[point]) * DetJ[point] * small_weight;
                else if(Frame == GLOBAL)
                    rOutput += rFunc.GetValue(GlobalCoords) * DetJ[point] * small_weight;
            }
        }
    }

    /// Integrate a function on the domain limited by a BRep using the sample geometry and a set of sample quadrature points.
    /// This allows for very high order quadrature rule, as well as arbitrary quadrature.
    /// The caller has to manually set rOutput to zero before calling this function.
    /// if Frame=0, the function is assumed to be defined in the local frame.
    /// if Frame=1, the function is assumed to be defined in the global frame. If the function is defined in the local frame (i.e. shape function), the Integrate<0> method must be used.
    template<typename TOutputType, typename TIntegrationPointsArrayType, int Frame>
    void Integrate(GeometryType::Pointer pParentGeometry,
            const Function<array_1d<double, 3>, TOutputType>& rFunc,
            const BRep& r_brep,
            TOutputType& rOutput,
            const TIntegrationPointsArrayType& sample_integration_points,
            const double& small_weight) const
    {
        GeometryType::IntegrationPointsArrayType integration_points;

        this->ConstructQuadrature(pParentGeometry, integration_points, sample_integration_points);

        std::vector<double> DetJ;
        Function<double, double>::ComputeDetJ(DetJ, *pParentGeometry, integration_points);

        CoordinatesArrayType GlobalCoords;

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            pParentGeometry->GlobalCoordinates(GlobalCoords, integration_points[point]);
            if(r_brep.IsInside(GlobalCoords))
            {
                if(Frame == LOCAL)
                    rOutput += rFunc.GetValue(integration_points[point]) * DetJ[point] * integration_points[point].Weight();
                else if(Frame == GLOBAL)
                    rOutput += rFunc.GetValue(GlobalCoords) * DetJ[point] * integration_points[point].Weight();
            }
            else
            {
                if(Frame == LOCAL)
                    rOutput += rFunc.GetValue(integration_points[point]) * DetJ[point] * small_weight;
                else if(Frame == GLOBAL)
                    rOutput += rFunc.GetValue(GlobalCoords) * DetJ[point] * small_weight;
            }
        }
    }

    /*****************************************************************/
    /******* QUADRATURE GENERATION ***********************************/
    /*****************************************************************/

    /// Construct a custom quadrature, i.e higher order Gauss Legendre or Gauss-Lobatto quadrature.
    /// @param quadrature_type  type of the quadrature, i.e 1: Gauss-Legendre, 2: Gauss-Lobatto, 3: Gauss-Radau
    virtual GeometryType::IntegrationPointsArrayType ConstructCustomQuadrature(const int& quadrature_type, const int& integration_order) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Construct the recursive integration point array using Kratos built-in quadrature.
    /// REMARKS: the integration_points is in local coordinates system
    void ConstructQuadrature(GeometryType::Pointer pParentGeometry,
            GeometryType::IntegrationPointsArrayType& integration_points,
            const int& integration_method) const
    {
        int quadrature_type = QuadratureUtility::GetQuadratureType(integration_method);

        if(quadrature_type == 0) // Kratos built-in
        {
            int quadrature_order = QuadratureUtility::GetQuadratureOrder(integration_method);
            GeometryData::IntegrationMethod ThisIntegrationMethod = Function<double, double>::GetIntegrationMethod(quadrature_order);
            this->ConstructQuadrature(pParentGeometry, integration_points, ThisIntegrationMethod);
        }
        else // use custom quadrature
        {
            int quadrature_order = QuadratureUtility::GetQuadratureOrder(integration_method);
            GeometryType::IntegrationPointsArrayType sample_integration_points = this->ConstructCustomQuadrature(quadrature_type, quadrature_order);
            this->ConstructQuadrature(pParentGeometry, integration_points, sample_integration_points);
        }
    }

    /// Construct the recursive integration point array using Kratos built-in quadrature.
    /// REMARKS: the integration_points is in local coordinates system.
    void ConstructQuadrature(GeometryType::Pointer pParentGeometry,
            GeometryType::IntegrationPointsArrayType& integration_points,
            const GeometryData::IntegrationMethod& ThisIntegrationMethod) const
    {
        if(this->IsLeaf())
        {
            GeometryType::Pointer pThisReferenceGeometry = this->pCreateReferenceGeometry();

            const GeometryType::IntegrationPointsArrayType& sample_integration_points
                = pThisReferenceGeometry->IntegrationPoints( ThisIntegrationMethod );

            Vector ShapeFunctionValuesOnReference;
            for(std::size_t point = 0; point < sample_integration_points.size(); ++point)
            {
                GeometryType::IntegrationPointType integration_point;

                ShapeFunctionValuesOnReference = pThisReferenceGeometry->ShapeFunctionsValues(ShapeFunctionValuesOnReference, sample_integration_points[point]);

                noalias(integration_point) = ZeroVector(3);
                for(std::size_t i = 0; i < pThisReferenceGeometry->size(); ++i)
                {
                    noalias(integration_point) += ShapeFunctionValuesOnReference[i] * (*pThisReferenceGeometry)[i];
                }

                double DetJsmall = Function<double, double>::ComputeDetJ(*pThisReferenceGeometry, sample_integration_points[point]);

                integration_point.SetWeight( DetJsmall * sample_integration_points[point].Weight() );

                integration_points.push_back( integration_point );
            }
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
            {
                mpChildren[i]->ConstructQuadrature(pParentGeometry, integration_points, ThisIntegrationMethod);
            }
        }
    }

    /// Construct the recursive integration point array using provided sample quadrature on a reference cell
    /// REMARKS: the integration_points is in local coordinates system
    void ConstructQuadrature(GeometryType::Pointer pParentGeometry,
            GeometryType::IntegrationPointsArrayType& integration_points,
            const GeometryType::IntegrationPointsArrayType& sample_integration_points) const
    {
        if(this->IsLeaf())
        {
            GeometryType::Pointer pThisReferenceGeometry = this->pCreateReferenceGeometry();

            Vector ShapeFunctionValuesOnReference;
            for(std::size_t point = 0; point < sample_integration_points.size(); ++point)
            {
                GeometryType::IntegrationPointType integration_point;

                ShapeFunctionValuesOnReference = pThisReferenceGeometry->ShapeFunctionsValues(ShapeFunctionValuesOnReference, sample_integration_points[point]);

                noalias(integration_point) = ZeroVector(3);
                for(std::size_t i = 0; i < pThisReferenceGeometry->size(); ++i)
                {
                    noalias(integration_point) += ShapeFunctionValuesOnReference[i] * (*pThisReferenceGeometry)[i];
                }

                double DetJsmall = Function<double, double>::ComputeDetJ(*pThisReferenceGeometry, sample_integration_points[point]);

                integration_point.SetWeight( DetJsmall * sample_integration_points[point].Weight() );

                integration_points.push_back( integration_point );
            }
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
            {
                mpChildren[i]->ConstructQuadrature(pParentGeometry, integration_points, sample_integration_points);
            }
        }
    }

    /*****************************************************************/
    /******* POST PROCESSING *****************************************/
    /*****************************************************************/

    template<bool TRecursive, class TEntityType>
    void AddToModelPart(GeometryType::Pointer pParentGeometry,
            ModelPart& r_model_part,
            TEntityType const& r_clone_entity,
            std::size_t& lastNodeId,
            std::size_t& lastEntityId,
            const int level) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    template<bool TRecursive, class TEntityType>
    void AddToModelPart() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /*****************************************************************/
    /******* AUXILIARY ROUTINES **************************************/
    /*****************************************************************/

    /// Check if a point in local coordinates is on the boundary of the quad-tree node up to some tolerances
    virtual bool IsOnBoundary(const CoordinatesArrayType& rLocalPoint, const double& tol) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Check if a point in local coordinates is strictly inside the quad-tree node
    virtual bool IsInside(const CoordinatesArrayType& rLocalPoint) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Compute a center in local coordinates of the cell
    virtual CoordinatesArrayType ReferenceCenter() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Compute the center of gravity in the domain defined by a parent geometry of this quad-tree node w.r.t a BRep
    /// Remarks: COG will be in global coordinates w.r.t pParentGeometry
    bool CenterOfGravity(PointType& COG, GeometryType::Pointer pParentGeometry, const BRep& r_brep, const int& integration_method) const
    {
        double X = 0.0, Y = 0.0, Z = 0.0, A = 0.0;
        FunctionR3R1::Pointer FX = FunctionR3R1::Pointer(new MonomialFunctionR3R1<1, 0, 0>());
        FunctionR3R1::Pointer FY = FunctionR3R1::Pointer(new MonomialFunctionR3R1<0, 1, 0>());
        FunctionR3R1::Pointer FZ = FunctionR3R1::Pointer(new MonomialFunctionR3R1<0, 0, 1>());
        FunctionR3R1::Pointer FH = FunctionR3R1::Pointer(new HeavisideFunction<FunctionR3R1>(r_brep));
        this->Integrate<double, GLOBAL>(pParentGeometry, ProductFunction<FunctionR3R1>(FX, FH), X, integration_method);
        this->Integrate<double, GLOBAL>(pParentGeometry, ProductFunction<FunctionR3R1>(FY, FH), Y, integration_method);
        this->Integrate<double, GLOBAL>(pParentGeometry, ProductFunction<FunctionR3R1>(FZ, FH), Z, integration_method);
        this->Integrate<double, GLOBAL>(pParentGeometry, *FH, A, integration_method);
//        std::cout << "X: " << X << ", Y: " << Y << ", Z: " << Z << ", A: " << A << std::endl;
        if(A != 0.0)
        {
            // the quad-tree is able to integrate the COG
            COG[0] = X/A;
            COG[1] = Y/A;
            COG[2] = Z/A;
            return true;
        }
        else
            return false;
    }

    /// Compute the domain size covered by the quad tree node
    double DomainSize(GeometryType::Pointer pParentGeometry, const BRep& r_brep, const int& integration_method) const
    {
        double A = 0.0;
        FunctionR3R1::Pointer FH = FunctionR3R1::Pointer(new HeavisideFunction<FunctionR3R1>(r_brep));
        this->Integrate<double, GLOBAL>(pParentGeometry, *FH, A, integration_method);
        return A;
    }

    /// Create the geometry with the node in local coordinates
    virtual GeometryType::Pointer pCreateReferenceGeometry() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Create the geometry in global coordinates for a quadtree node, based on the parent geometry
    virtual GeometryType::Pointer pCreateGeometry(GeometryType::Pointer pParentGeometry) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Create a list of sampling points on the geometry. The local coordinates of the sampling points must be taken from the quadtree parameter range.
    virtual void CreateSamplingPoints(std::vector<PointType>& SamplingPoints,
            GeometryType& r_geom, const std::size_t& nsampling) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "QuadTreeNode";
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

    friend std::ostream& operator<<(std::ostream& rOStream, const QuadTreeNode& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << " ";
        rThis.PrintData(rOStream);
        return rOStream;
    }

protected:

    std::vector<QuadTreeNode::Pointer> mpChildren;

}; // class QuadTreeNode


template<>
void QuadTreeNode::AddToModelPart<false, Element>(GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Element const& r_clone_element,
        std::size_t& lastNodeId,
        std::size_t& lastElementId,
        const int level) const;


template<>
void QuadTreeNode::AddToModelPart<false, Condition>(GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Condition const& r_clone_condition,
        std::size_t& lastNodeId,
        std::size_t& lastCondId,
        const int level) const;


template<>
void QuadTreeNode::AddToModelPart<true, Element>(GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Element const& r_clone_element,
        std::size_t& lastNodeId,
        std::size_t& lastElementId,
        const int level) const;

template<>
void QuadTreeNode::AddToModelPart<true, Condition>(GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Condition const& r_clone_condition,
        std::size_t& lastNodeId,
        std::size_t& lastCondId,
        const int level) const;



/// Quadrilateral quad-tree node in reference coordinates
class QuadTreeNodeQ4 : public QuadTreeNode
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadTreeNodeQ4);

    typedef QuadTreeNode BaseType;

    typedef BaseType::GeometryType GeometryType;

    QuadTreeNodeQ4(const double& Xmin, const double& Xmax, const double& Ymin, const double& Ymax)
    : BaseType(), mXmin(Xmin), mXmax(Xmax), mYmin(Ymin), mYmax(Ymax)
    {}

    virtual ~QuadTreeNodeQ4() {}

    virtual void Refine()
    {
        if(this->IsLeaf())
        {
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeQ4(mXmin, 0.5*(mXmin+mXmax), mYmin, 0.5*(mYmin+mYmax))));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeQ4(0.5*(mXmin+mXmax), mXmax, mYmin, 0.5*(mYmin+mYmax))));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeQ4(mXmin, 0.5*(mXmin+mXmax), 0.5*(mYmin+mYmax), mYmax)));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeQ4(0.5*(mXmin+mXmax), mXmax, 0.5*(mYmin+mYmax), mYmax)));
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
                mpChildren[i]->Refine();
        }
    }

    virtual GeometryType::IntegrationPointsArrayType ConstructCustomQuadrature(const int& quadrature_type, const int& integration_order) const
    {
        if(quadrature_type == 1) // Gauss-Legendre
        {
            if(integration_order == 1)
                return Quadrature<QuadrilateralGaussLegendreIntegrationPoints1, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 2)
                return Quadrature<QuadrilateralGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 3)
                return Quadrature<QuadrilateralGaussLegendreIntegrationPoints3, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 4)
                return Quadrature<QuadrilateralGaussLegendreIntegrationPoints4, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 5)
                return Quadrature<QuadrilateralGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 6)
                return Quadrature<QuadrilateralGaussLegendreIntegrationPoints6, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 7)
                return Quadrature<QuadrilateralGaussLegendreIntegrationPoints7, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 8)
                return Quadrature<QuadrilateralGaussLegendreIntegrationPoints8, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 9)
                return Quadrature<QuadrilateralGaussLegendreIntegrationPoints9, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 10)
                return Quadrature<QuadrilateralGaussLegendreIntegrationPoints10, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else
                KRATOS_THROW_ERROR(std::logic_error, "This integration_order of Gauss-Legendre is not supported:", integration_order);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "This quadrature type is not supported:", quadrature_type);
    }

    virtual GeometryType::Pointer pCreateReferenceGeometry() const
    {
        GeometryType::PointType P1(0, mXmin, mYmin);
        GeometryType::PointType P2(1, mXmax, mYmin);
        GeometryType::PointType P3(2, mXmax, mYmax);
        GeometryType::PointType P4(3, mXmin, mYmax);
        return GeometryType::Pointer(new Quadrilateral2D4<BaseType::NodeType>(P1, P2, P3, P4));
    }

    virtual GeometryType::Pointer pCreateGeometry(GeometryType::Pointer pParentGeometry) const
    {
        CoordinatesArrayType X;
        GeometryType::Pointer pNewGeometry;

        if(    pParentGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4
            || pParentGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral3D4 )
        {
            GeometryType::PointType P1, P2, P3, P4;

            X[0] = mXmin; X[1] = mYmin; X[2] = 0.0;
            pParentGeometry->GlobalCoordinates(P1, X);

            X[0] = mXmax; X[1] = mYmin;
            pParentGeometry->GlobalCoordinates(P2, X);

            X[0] = mXmax; X[1] = mYmax;
            pParentGeometry->GlobalCoordinates(P3, X);

            X[0] = mXmin; X[1] = mYmax;
            pParentGeometry->GlobalCoordinates(P4, X);

            if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4)
                pNewGeometry = GeometryType::Pointer(new Quadrilateral2D4<BaseType::NodeType>(P1, P2, P3, P4));
            if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral3D4)
                pNewGeometry = GeometryType::Pointer(new Quadrilateral3D4<BaseType::NodeType>(P1, P2, P3, P4));
        }
        else if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8
             || pParentGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral3D8 )
        {
            GeometryType::PointType P1, P2, P3, P4, P5, P6, P7, P8;

            X[0] = mXmin; X[1] = mYmin; X[2] = 0.0;
            pParentGeometry->GlobalCoordinates(P1, X);

            X[0] = mXmax; X[1] = mYmin;
            pParentGeometry->GlobalCoordinates(P2, X);

            X[0] = mXmax; X[1] = mYmax;
            pParentGeometry->GlobalCoordinates(P3, X);

            X[0] = mXmin; X[1] = mYmax;
            pParentGeometry->GlobalCoordinates(P4, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmin;
            pParentGeometry->GlobalCoordinates(P5, X);

            X[0] = mXmax; X[1] = 0.5*(mYmin+mYmax);
            pParentGeometry->GlobalCoordinates(P6, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmax;
            pParentGeometry->GlobalCoordinates(P7, X);

            X[0] = mXmin; X[1] = 0.5*(mYmin+mYmax);
            pParentGeometry->GlobalCoordinates(P8, X);

            if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8)
                pNewGeometry = GeometryType::Pointer(new Quadrilateral2D8<BaseType::NodeType>(P1, P2, P3, P4, P5, P6, P7, P8));
            if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral3D8)
                pNewGeometry = GeometryType::Pointer(new Quadrilateral3D8<BaseType::NodeType>(P1, P2, P3, P4, P5, P6, P7, P8));
        }
        else if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9
             || pParentGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral3D9 )
        {
            GeometryType::PointType P1, P2, P3, P4, P5, P6, P7, P8, P9;

            X[0] = mXmin; X[1] = mYmin; X[2] = 0.0;
            pParentGeometry->GlobalCoordinates(P1, X);

            X[0] = mXmax; X[1] = mYmin;
            pParentGeometry->GlobalCoordinates(P2, X);

            X[0] = mXmax; X[1] = mYmax;
            pParentGeometry->GlobalCoordinates(P3, X);

            X[0] = mXmin; X[1] = mYmax;
            pParentGeometry->GlobalCoordinates(P4, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmin;
            pParentGeometry->GlobalCoordinates(P5, X);

            X[0] = mXmax; X[1] = 0.5*(mYmin+mYmax);
            pParentGeometry->GlobalCoordinates(P6, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmax;
            pParentGeometry->GlobalCoordinates(P7, X);

            X[0] = mXmin; X[1] = 0.5*(mYmin+mYmax);
            pParentGeometry->GlobalCoordinates(P8, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = 0.5*(mYmin+mYmax);
            pParentGeometry->GlobalCoordinates(P9, X);

            if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9)
                pNewGeometry = GeometryType::Pointer(new Quadrilateral2D9<BaseType::NodeType>(P1, P2, P3, P4, P5, P6, P7, P8, P9));
            if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral3D9)
                pNewGeometry = GeometryType::Pointer(new Quadrilateral3D9<BaseType::NodeType>(P1, P2, P3, P4, P5, P6, P7, P8, P9));
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The parent geometry type is invalid:", pParentGeometry->GetGeometryType())

        for(std::size_t i = 0; i < pNewGeometry->size(); ++i)
            (*pNewGeometry)[i].SetInitialPosition((*pNewGeometry)[i]);

        return pNewGeometry;
    }

    virtual bool IsOnBoundary(const CoordinatesArrayType& rLocalPoint, const double& tol) const
    {
        bool is_onboundary = false;
        is_onboundary = is_onboundary || ( (fabs(rLocalPoint[0] - mXmin) < tol) && (mYmin - tol < rLocalPoint[1]) && (rLocalPoint[1] < mYmax + tol) );
        if(is_onboundary) return true;
        is_onboundary = is_onboundary || ( (fabs(rLocalPoint[0] - mXmax) < tol) && (mYmin - tol < rLocalPoint[1]) && (rLocalPoint[1] < mYmax + tol) );
        if(is_onboundary) return true;
        is_onboundary = is_onboundary || ( (fabs(rLocalPoint[1] - mYmin) < tol) && (mXmin - tol < rLocalPoint[0]) && (rLocalPoint[0] < mXmax + tol) );
        if(is_onboundary) return true;
        is_onboundary = is_onboundary || ( (fabs(rLocalPoint[1] - mYmax) < tol) && (mXmin - tol < rLocalPoint[0]) && (rLocalPoint[0] < mXmax + tol) );
        return is_onboundary;
    }

    virtual bool IsInside(const CoordinatesArrayType& rLocalPoint) const
    {
        return (mXmin < rLocalPoint[0]) && (rLocalPoint[0] < mXmax) && (mYmin < rLocalPoint[1]) && (rLocalPoint[1] < mYmax);
    }

    virtual CoordinatesArrayType ReferenceCenter() const
    {
        CoordinatesArrayType C;
        C[0] = 0.5*(mXmin + mXmax);
        C[1] = 0.5*(mYmin + mYmax);
        C[2] = 0.0;
        return C;
    }

    virtual void CreateSamplingPoints(std::vector<PointType>& SamplingPoints,
            GeometryType& r_geom, const std::size_t& nsampling) const
    {
        double dX = (mXmax - mXmin) / nsampling;
        double dY = (mYmax - mYmin) / nsampling;

        SamplingPoints.reserve((nsampling+1) * (nsampling+1));
        CoordinatesArrayType X;
        PointType P;
        X[2] = 0.0;
        for(std::size_t i = 0; i < nsampling + 1; ++i)
        {
            X[0] = mXmin + i*dX;
            for(std::size_t j = 0; j < nsampling + 1; ++j)
            {
                X[1] = mYmin + j*dY;
                r_geom.GlobalCoordinates(P, X);
                SamplingPoints.push_back(P);
            }
        }
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "QuadTreeNodeQ4";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "(" << mXmin << ", " << mYmin << ") - (" << mXmax << ", " << mYmax << ")";
    }

private:
    double mXmin, mXmax;
    double mYmin, mYmax;
};


#ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
/// Bezier quad-tree node in reference coordinates
/// TODO generalize quad tree node for Bezier geometry for different order
class QuadTreeNodeBezier2D : public QuadTreeNodeQ4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadTreeNodeBezier2D);

    typedef QuadTreeNodeQ4 BaseType;

    typedef BaseType::GeometryType GeometryType;

    QuadTreeNodeBezier2D(const double& Xmin, const double& Xmax, const double& Ymin, const double& Ymax)
    : BaseType(Xmin, Xmax, Ymin, Ymax)
    {}

    virtual ~QuadTreeNodeBezier2D() {}

    virtual GeometryType::Pointer pCreateGeometry(GeometryType::Pointer pParentGeometry) const
    {
        if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Bezier2D
        || pParentGeometry->GetGeometryType() == GeometryData::Kratos_Bezier2D3)
        {
            KRATOS_THROW_ERROR(std::logic_error, "The generation of Bezier2D or Bezier2D3 geometry is not supported:", pParentGeometry->GetGeometryType())
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The parent geometry type is invalid:", pParentGeometry->GetGeometryType())
    }

    virtual GeometryType::IntegrationPointsArrayType ConstructCustomQuadrature(const int& quadrature_type, const int& integration_order) const
    {
        if(quadrature_type == 1) // Gauss-Legendre
        {
            GeometryType::IntegrationPointsArrayType integration_points;

            if(integration_order == 1)
                integration_points = Quadrature<QuadrilateralGaussLegendreIntegrationPoints1, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 2)
                integration_points = Quadrature<QuadrilateralGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 3)
                integration_points = Quadrature<QuadrilateralGaussLegendreIntegrationPoints3, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 4)
                integration_points = Quadrature<QuadrilateralGaussLegendreIntegrationPoints4, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 5)
                integration_points = Quadrature<QuadrilateralGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 6)
                integration_points = Quadrature<QuadrilateralGaussLegendreIntegrationPoints6, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 7)
                integration_points = Quadrature<QuadrilateralGaussLegendreIntegrationPoints7, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 8)
                integration_points = Quadrature<QuadrilateralGaussLegendreIntegrationPoints8, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 9)
                integration_points = Quadrature<QuadrilateralGaussLegendreIntegrationPoints9, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 10)
                integration_points = Quadrature<QuadrilateralGaussLegendreIntegrationPoints10, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else
                KRATOS_THROW_ERROR(std::logic_error, "This integration_order of Gauss-Legendre is not supported:", integration_order);

            // map the integration points from [-1, 1] x [-1, 1] to [0, 1] x [0, 1]
            for (std::size_t i = 0; i < integration_points.size(); ++i)
            {
                integration_points[i].X() = 0.5 * (integration_points[i].X() + 1);
                integration_points[i].Y() = 0.5 * (integration_points[i].Y() + 1);
                integration_points[i].Weight() = 0.25*integration_points[i].Weight();
            }

            return integration_points;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "This quadrature type is not supported:", quadrature_type);
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "QuadTreeNodeBezier2D";
    }
};
#endif


/// Hexahedral quad-tree node in reference coordinates
class QuadTreeNodeH8 : public QuadTreeNode
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadTreeNodeH8);

    typedef QuadTreeNode BaseType;

    typedef BaseType::GeometryType GeometryType;

    QuadTreeNodeH8(const double& Xmin, const double& Xmax, const double& Ymin, const double& Ymax, const double& Zmin, const double& Zmax)
    : BaseType(), mXmin(Xmin), mXmax(Xmax), mYmin(Ymin), mYmax(Ymax), mZmin(Zmin), mZmax(Zmax)
    {}

    virtual ~QuadTreeNodeH8() {}

    virtual void Refine()
    {
        if(this->IsLeaf())
        {
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeH8(mXmin, 0.5*(mXmin+mXmax), mYmin, 0.5*(mYmin+mYmax), mZmin, 0.5*(mZmin+mZmax))));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeH8(0.5*(mXmin+mXmax), mXmax, mYmin, 0.5*(mYmin+mYmax), mZmin, 0.5*(mZmin+mZmax))));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeH8(mXmin, 0.5*(mXmin+mXmax), 0.5*(mYmin+mYmax), mYmax, mZmin, 0.5*(mZmin+mZmax))));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeH8(0.5*(mXmin+mXmax), mXmax, 0.5*(mYmin+mYmax), mYmax, mZmin, 0.5*(mZmin+mZmax))));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeH8(mXmin, 0.5*(mXmin+mXmax), mYmin, 0.5*(mYmin+mYmax), 0.5*(mZmin+mZmax), mZmax)));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeH8(0.5*(mXmin+mXmax), mXmax, mYmin, 0.5*(mYmin+mYmax), 0.5*(mZmin+mZmax), mZmax)));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeH8(mXmin, 0.5*(mXmin+mXmax), 0.5*(mYmin+mYmax), mYmax, 0.5*(mZmin+mZmax), mZmax)));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeH8(0.5*(mXmin+mXmax), mXmax, 0.5*(mYmin+mYmax), mYmax, 0.5*(mZmin+mZmax), mZmax)));
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
                mpChildren[i]->Refine();
        }
    }

    virtual GeometryType::IntegrationPointsArrayType ConstructCustomQuadrature(const int& quadrature_type, const int& integration_order) const
    {
        if(quadrature_type == 1) // Gauss-Legendre
        {
            if(integration_order == 1)
                return Quadrature<HexahedronGaussLegendreIntegrationPoints1, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 2)
                return Quadrature<HexahedronGaussLegendreIntegrationPoints2, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 3)
                return Quadrature<HexahedronGaussLegendreIntegrationPoints3, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 4)
                return Quadrature<HexahedronGaussLegendreIntegrationPoints4, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 5)
                return Quadrature<HexahedronGaussLegendreIntegrationPoints5, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 6)
                return Quadrature<HexahedronGaussLegendreIntegrationPoints6, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 7)
                return Quadrature<HexahedronGaussLegendreIntegrationPoints7, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 8)
                return Quadrature<HexahedronGaussLegendreIntegrationPoints8, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else
                KRATOS_THROW_ERROR(std::logic_error, "This integration_order of Gauss-Legendre is not supported:", integration_order);
        }
        else if(quadrature_type == 2) // Gauss-Lobatto
        {
            if(integration_order == 1)
                return Quadrature<HexahedronGaussLobattoIntegrationPoints1, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 2)
                return Quadrature<HexahedronGaussLobattoIntegrationPoints2, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 3)
                return Quadrature<HexahedronGaussLobattoIntegrationPoints3, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 4)
                return Quadrature<HexahedronGaussLobattoIntegrationPoints4, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 5)
                return Quadrature<HexahedronGaussLobattoIntegrationPoints5, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 6)
                return Quadrature<HexahedronGaussLobattoIntegrationPoints6, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 7)
                return Quadrature<HexahedronGaussLobattoIntegrationPoints7, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 8)
                return Quadrature<HexahedronGaussLobattoIntegrationPoints8, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 9)
                return Quadrature<HexahedronGaussLobattoIntegrationPoints9, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 10)
                return Quadrature<HexahedronGaussLobattoIntegrationPoints10, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else
                KRATOS_THROW_ERROR(std::logic_error, "This integration_order of Gauss-Lobatto is not supported:", integration_order);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "This quadrature type is not supported:", quadrature_type);
    }

    virtual GeometryType::Pointer pCreateReferenceGeometry() const
    {
        GeometryType::PointType P1(0, mXmin, mYmin, mZmin);
        GeometryType::PointType P2(1, mXmax, mYmin, mZmin);
        GeometryType::PointType P3(2, mXmax, mYmax, mZmin);
        GeometryType::PointType P4(3, mXmin, mYmax, mZmin);
        GeometryType::PointType P5(4, mXmin, mYmin, mZmax);
        GeometryType::PointType P6(5, mXmax, mYmin, mZmax);
        GeometryType::PointType P7(6, mXmax, mYmax, mZmax);
        GeometryType::PointType P8(7, mXmin, mYmax, mZmax);
        return GeometryType::Pointer(new Hexahedra3D8<BaseType::NodeType>(P1, P2, P3, P4, P5, P6, P7, P8));
    }

    virtual GeometryType::Pointer pCreateGeometry(GeometryType::Pointer pParentGeometry) const
    {
        CoordinatesArrayType X;
        GeometryType::Pointer pNewGeometry;

        if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Hexahedra3D8)
        {
            GeometryType::PointType P1, P2, P3, P4, P5, P6, P7, P8;

            X[0] = mXmin; X[1] = mYmin; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P1, X);

            X[0] = mXmax; X[1] = mYmin; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P2, X);

            X[0] = mXmax; X[1] = mYmax; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P3, X);

            X[0] = mXmin; X[1] = mYmax; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P4, X);

            X[0] = mXmin; X[1] = mYmin; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P5, X);

            X[0] = mXmax; X[1] = mYmin; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P6, X);

            X[0] = mXmax; X[1] = mYmax; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P7, X);

            X[0] = mXmin; X[1] = mYmax; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P8, X);

            pNewGeometry = GeometryType::Pointer(new Hexahedra3D8<BaseType::NodeType>(P1, P2, P3, P4, P5, P6, P7, P8));
        }
        else if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Hexahedra3D20)
        {
            GeometryType::PointType P1, P2, P3, P4, P5, P6, P7, P8,
                    P9, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19, P20;

            X[0] = mXmin; X[1] = mYmin; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P1, X);

            X[0] = mXmax; X[1] = mYmin; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P2, X);

            X[0] = mXmax; X[1] = mYmax; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P3, X);

            X[0] = mXmin; X[1] = mYmax; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P4, X);

            X[0] = mXmin; X[1] = mYmin; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P5, X);

            X[0] = mXmax; X[1] = mYmin; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P6, X);

            X[0] = mXmax; X[1] = mYmax; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P7, X);

            X[0] = mXmin; X[1] = mYmax; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P8, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmax; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P9, X);

            X[0] = mXmax; X[1] = 0.5*(mYmin+mYmax); X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P10, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmax; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P11, X);

            X[0] = mXmin; X[1] = 0.5*(mYmin+mYmax); X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P12, X);

            X[0] = mXmin; X[1] = mYmin; X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P13, X);

            X[0] = mXmax; X[1] = mYmin; X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P14, X);

            X[0] = mXmax; X[1] = mYmax; X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P15, X);

            X[0] = mXmin; X[1] = mYmax; X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P16, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmax; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P17, X);

            X[0] = mXmax; X[1] = 0.5*(mYmin+mYmax); X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P18, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmax; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P19, X);

            X[0] = mXmin; X[1] = 0.5*(mYmin+mYmax); X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P20, X);

            pNewGeometry = GeometryType::Pointer(new Hexahedra3D20<BaseType::NodeType>(P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19, P20));
        }
        else if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Hexahedra3D27)
        {
            GeometryType::PointType P1, P2, P3, P4, P5, P6, P7, P8,
                    P9, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19, P20,
                    P21, P22, P23, P24, P25, P26, P27;

            X[0] = mXmin; X[1] = mYmin; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P1, X);

            X[0] = mXmax; X[1] = mYmin; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P2, X);

            X[0] = mXmax; X[1] = mYmax; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P3, X);

            X[0] = mXmin; X[1] = mYmax; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P4, X);

            X[0] = mXmin; X[1] = mYmin; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P5, X);

            X[0] = mXmax; X[1] = mYmin; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P6, X);

            X[0] = mXmax; X[1] = mYmax; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P7, X);

            X[0] = mXmin; X[1] = mYmax; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P8, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmin; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P9, X);

            X[0] = mXmax; X[1] = 0.5*(mYmin+mYmax); X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P10, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmax; X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P11, X);

            X[0] = mXmin; X[1] = 0.5*(mYmin+mYmax); X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P12, X);

            X[0] = mXmin; X[1] = mYmin; X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P13, X);

            X[0] = mXmax; X[1] = mYmin; X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P14, X);

            X[0] = mXmax; X[1] = mYmax; X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P15, X);

            X[0] = mXmin; X[1] = mYmax; X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P16, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmin; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P17, X);

            X[0] = mXmax; X[1] = 0.5*(mYmin+mYmax); X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P18, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmax; X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P19, X);

            X[0] = mXmin; X[1] = 0.5*(mYmin+mYmax); X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P20, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = 0.5*(mYmin+mYmax); X[2] = mZmin;
            pParentGeometry->GlobalCoordinates(P21, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmin; X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P22, X);

            X[0] = mXmax; X[1] = 0.5*(mYmin+mYmax); X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P23, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = mYmax; X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P24, X);

            X[0] = mXmin; X[1] = 0.5*(mYmin+mYmax); X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P25, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = 0.5*(mYmin+mYmax); X[2] = mZmax;
            pParentGeometry->GlobalCoordinates(P26, X);

            X[0] = 0.5*(mXmin+mXmax); X[1] = 0.5*(mYmin+mYmax); X[2] = 0.5*(mZmin+mZmax);
            pParentGeometry->GlobalCoordinates(P27, X);

            pNewGeometry = GeometryType::Pointer(new Hexahedra3D27<BaseType::NodeType>(P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19, P20, P21, P22, P23, P24, P25, P26, P27));
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The parent geometry type is invalid:", pParentGeometry->GetGeometryType())

        for(std::size_t i = 0; i < pNewGeometry->size(); ++i)
            (*pNewGeometry)[i].SetInitialPosition((*pNewGeometry)[i]);

        return pNewGeometry;
    }

    virtual bool IsOnBoundary(const CoordinatesArrayType& rLocalPoint, const double& tol) const
    {
        bool is_onboundary =             ( (fabs(rLocalPoint[0] - mXmin) < tol) && (mYmin - tol < rLocalPoint[1]) && (rLocalPoint[1] < mYmax + tol) && (mZmin - tol < rLocalPoint[2]) && (rLocalPoint[2] < mZmax + tol) );
        if(is_onboundary) return true;
        is_onboundary = is_onboundary || ( (fabs(rLocalPoint[0] - mXmax) < tol) && (mYmin - tol < rLocalPoint[1]) && (rLocalPoint[1] < mYmax + tol) && (mZmin - tol < rLocalPoint[2]) && (rLocalPoint[2] < mZmax + tol) );
        if(is_onboundary) return true;
        is_onboundary = is_onboundary || ( (fabs(rLocalPoint[1] - mYmin) < tol) && (mXmin - tol < rLocalPoint[0]) && (rLocalPoint[0] < mXmax + tol) && (mZmin - tol < rLocalPoint[2]) && (rLocalPoint[2] < mZmax + tol) );
        if(is_onboundary) return true;
        is_onboundary = is_onboundary || ( (fabs(rLocalPoint[1] - mYmax) < tol) && (mXmin - tol < rLocalPoint[0]) && (rLocalPoint[0] < mXmax + tol) && (mZmin - tol < rLocalPoint[2]) && (rLocalPoint[2] < mZmax + tol) );
        if(is_onboundary) return true;
        is_onboundary = is_onboundary || ( (fabs(rLocalPoint[2] - mZmin) < tol) && (mXmin - tol < rLocalPoint[0]) && (rLocalPoint[0] < mXmax + tol) && (mYmin - tol < rLocalPoint[1]) && (rLocalPoint[1] < mYmax + tol) );
        if(is_onboundary) return true;
        is_onboundary = is_onboundary || ( (fabs(rLocalPoint[2] - mZmax) < tol) && (mXmin - tol < rLocalPoint[0]) && (rLocalPoint[0] < mXmax + tol) && (mYmin - tol < rLocalPoint[1]) && (rLocalPoint[1] < mYmax + tol) );
        if(is_onboundary) return is_onboundary;
        return false;
    }

    virtual bool IsInside(const CoordinatesArrayType& rLocalPoint) const
    {
        return (mXmin < rLocalPoint[0]) && (rLocalPoint[0] < mXmax)
            && (mYmin < rLocalPoint[1]) && (rLocalPoint[1] < mYmax)
            && (mZmin < rLocalPoint[2]) && (rLocalPoint[2] < mZmax);
    }

    virtual CoordinatesArrayType ReferenceCenter() const
    {
        CoordinatesArrayType C;
        C[0] = 0.5*(mXmin + mXmax);
        C[1] = 0.5*(mYmin + mYmax);
        C[2] = 0.5*(mZmin + mZmax);
        return C;
    }

    virtual void CreateSamplingPoints(std::vector<PointType>& SamplingPoints,
            GeometryType& r_geom, const std::size_t& nsampling) const
    {
        double dX = (mXmax - mXmin) / nsampling;
        double dY = (mYmax - mYmin) / nsampling;
        double dZ = (mZmax - mZmin) / nsampling;

        SamplingPoints.reserve((nsampling+1) * (nsampling+1) * (nsampling+1));
        CoordinatesArrayType X;
        PointType P;
        for(std::size_t i = 0; i < nsampling + 1; ++i)
        {
            X[0] = mXmin + i*dX;
            for(std::size_t j = 0; j < nsampling + 1; ++j)
            {
                X[1] = mYmin + j*dY;
                for(std::size_t k = 0; k < nsampling + 1; ++k)
                {
                    X[2] = mZmin + k*dZ;
                    r_geom.GlobalCoordinates(P, X);
                    SamplingPoints.push_back(P);
                }
            }
        }
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "QuadTreeNodeH8";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "(" << mXmin << ", " << mYmin << ", " << mZmin << ") - (" << mXmax << ", " << mYmax << ", " << mZmax << ")";
    }

private:
    double mXmin, mXmax;
    double mYmin, mYmax;
    double mZmin, mZmax;
};


#ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
/// Bezier oct-tree node in reference coordinates
/// TODO generalize quad tree node for Bezier geometry for different order
class QuadTreeNodeBezier3D : public QuadTreeNodeH8
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadTreeNodeBezier3D);

    typedef QuadTreeNodeH8 BaseType;

    typedef BaseType::GeometryType GeometryType;

    QuadTreeNodeBezier3D(const double& Xmin, const double& Xmax, const double& Ymin, const double& Ymax, const double& Zmin, const double& Zmax)
    : BaseType(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)
    {}

    virtual ~QuadTreeNodeBezier3D() {}

    virtual GeometryType::Pointer pCreateGeometry(GeometryType::Pointer pParentGeometry) const
    {
        if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Bezier3D)
        {
            KRATOS_THROW_ERROR(std::logic_error, "The generation of Bezier3D geometry is not supported:", pParentGeometry->GetGeometryType())
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The parent geometry type is invalid:", pParentGeometry->GetGeometryType())
    }

    virtual GeometryType::IntegrationPointsArrayType ConstructCustomQuadrature(const int& quadrature_type, const int& integration_order) const
    {
        if(quadrature_type == 1) // Gauss-Legendre
        {
            GeometryType::IntegrationPointsArrayType integration_points;

            if(integration_order == 1)
                integration_points = Quadrature<HexahedronGaussLegendreIntegrationPoints1, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 2)
                integration_points = Quadrature<HexahedronGaussLegendreIntegrationPoints2, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 3)
                integration_points = Quadrature<HexahedronGaussLegendreIntegrationPoints3, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 4)
                integration_points = Quadrature<HexahedronGaussLegendreIntegrationPoints4, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 5)
                integration_points = Quadrature<HexahedronGaussLegendreIntegrationPoints5, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 6)
                integration_points = Quadrature<HexahedronGaussLegendreIntegrationPoints6, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 7)
                integration_points = Quadrature<HexahedronGaussLegendreIntegrationPoints7, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 8)
                integration_points = Quadrature<HexahedronGaussLegendreIntegrationPoints8, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else
                KRATOS_THROW_ERROR(std::logic_error, "This integration_order of Gauss-Legendre is not supported:", integration_order);

            // map the integration points from [-1, 1] x [-1, 1] x [-1, 1] to [0, 1] x [0, 1] x [0, 1]
            for (std::size_t i = 0; i < integration_points.size(); ++i)
            {
                integration_points[i].X() = 0.5 * (integration_points[i].X() + 1);
                integration_points[i].Y() = 0.5 * (integration_points[i].Y() + 1);
                integration_points[i].Z() = 0.5 * (integration_points[i].Z() + 1);
                integration_points[i].Weight() = 0.125*integration_points[i].Weight();
            }

            return integration_points;
        }
        else if(quadrature_type == 2) // Gauss-Lobatto
        {
            GeometryType::IntegrationPointsArrayType integration_points;

            if(integration_order == 1)
                integration_points = Quadrature<HexahedronGaussLobattoIntegrationPoints1, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 2)
                integration_points = Quadrature<HexahedronGaussLobattoIntegrationPoints2, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 3)
                integration_points = Quadrature<HexahedronGaussLobattoIntegrationPoints3, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 4)
                integration_points = Quadrature<HexahedronGaussLobattoIntegrationPoints4, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 5)
                integration_points = Quadrature<HexahedronGaussLobattoIntegrationPoints5, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 6)
                integration_points = Quadrature<HexahedronGaussLobattoIntegrationPoints6, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 7)
                integration_points = Quadrature<HexahedronGaussLobattoIntegrationPoints7, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 8)
                integration_points = Quadrature<HexahedronGaussLobattoIntegrationPoints8, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 9)
                integration_points = Quadrature<HexahedronGaussLobattoIntegrationPoints9, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else if(integration_order == 10)
                integration_points = Quadrature<HexahedronGaussLobattoIntegrationPoints10, 3, IntegrationPoint<3> >::GenerateIntegrationPoints();
            else
                KRATOS_THROW_ERROR(std::logic_error, "This integration_order of Gauss-Lobatto is not supported:", integration_order);

            // map the integration points from [-1, 1] x [-1, 1] x [-1, 1] to [0, 1] x [0, 1] x [0, 1]
            for (std::size_t i = 0; i < integration_points.size(); ++i)
            {
                integration_points[i].X() = 0.5 * (integration_points[i].X() + 1);
                integration_points[i].Y() = 0.5 * (integration_points[i].Y() + 1);
                integration_points[i].Z() = 0.5 * (integration_points[i].Z() + 1);
                integration_points[i].Weight() = 0.125*integration_points[i].Weight();
            }

            return integration_points;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "This quadrature type is not supported:", quadrature_type);
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "QuadTreeNodeBezier3D";
    }
};
#endif


/// Triangular quad-tree node in reference coordinates
class QuadTreeNodeT3 : public QuadTreeNode
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadTreeNodeT3);

    typedef QuadTreeNode BaseType;

    typedef BaseType::GeometryType GeometryType;

    QuadTreeNodeT3(const double& X0, const double& Y0, const double& X1, const double& Y1, const double& X2, const double& Y2)
    : BaseType(), mX0(X0), mY0(Y0), mX1(X1), mY1(Y1), mX2(X2), mY2(Y2)
    {}

    virtual ~QuadTreeNodeT3() {}

    virtual void Refine()
    {
        if(this->IsLeaf())
        {
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT3(mX0, mY0, 0.5*(mX0+mX1), 0.5*(mY0+mY1), 0.5*(mX0+mX2), 0.5*(mY0+mY2))));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT3(0.5*(mX0+mX1), 0.5*(mY0+mY1), mX1, mY1, 0.5*(mX1+mX2), 0.5*(mY1+mY2))));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT3(0.5*(mX1+mX2), 0.5*(mY1+mY2), mX2, mY2, 0.5*(mX2+mX0), 0.5*(mY2+mY0))));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT3(0.5*(mX0+mX1), 0.5*(mY0+mY1), 0.5*(mX1+mX2), 0.5*(mY1+mY2), 0.5*(mX2+mX0), 0.5*(mY2+mY0))));
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
                mpChildren[i]->Refine();
        }
    }

    virtual GeometryType::Pointer pCreateReferenceGeometry() const
    {
        GeometryType::PointType P1(0, mX0, mY0);
        GeometryType::PointType P2(1, mX1, mY1);
        GeometryType::PointType P3(2, mX2, mY2);
        return GeometryType::Pointer(new Triangle2D3<BaseType::NodeType>(P1, P2, P3));
    }

    virtual GeometryType::Pointer pCreateGeometry(GeometryType::Pointer pParentGeometry) const
    {
        CoordinatesArrayType X;
        GeometryType::Pointer pNewGeometry;

        if(    pParentGeometry->GetGeometryType() == GeometryData::Kratos_Triangle2D3
            || pParentGeometry->GetGeometryType() == GeometryData::Kratos_Triangle3D3 )
        {
            GeometryType::PointType P1, P2, P3;

            X[0] = mX0; X[1] = mY0; X[2] = 0.0;
            pParentGeometry->GlobalCoordinates(P1, X);

            X[0] = mX1; X[1] = mY1;
            pParentGeometry->GlobalCoordinates(P2, X);

            X[0] = mX2; X[1] = mY2;
            pParentGeometry->GlobalCoordinates(P3, X);

            if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Triangle2D3)
                pNewGeometry = GeometryType::Pointer(new Triangle2D3<BaseType::NodeType>(P1, P2, P3));
            if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Triangle3D3)
                pNewGeometry = GeometryType::Pointer(new Triangle3D3<BaseType::NodeType>(P1, P2, P3));
        }
        else if( pParentGeometry->GetGeometryType() == GeometryData::Kratos_Triangle2D6
              || pParentGeometry->GetGeometryType() == GeometryData::Kratos_Triangle3D6 )
        {
            GeometryType::PointType P1, P2, P3, P4, P5, P6;

            X[0] = mX0; X[1] = mY0; X[2] = 0.0;
            pParentGeometry->GlobalCoordinates(P1, X);

            X[0] = mX1; X[1] = mY1;
            pParentGeometry->GlobalCoordinates(P2, X);

            X[0] = mX2; X[1] = mY2;
            pParentGeometry->GlobalCoordinates(P3, X);

            X[0] = 0.5*(mX0+mX1); X[1] = 0.5*(mY0+mY1);
            pParentGeometry->GlobalCoordinates(P4, X);

            X[0] = 0.5*(mX1+mX2); X[1] = 0.5*(mY1+mY2);
            pParentGeometry->GlobalCoordinates(P5, X);

            X[0] = 0.5*(mX2+mX0); X[1] = 0.5*(mY2+mY0);
            pParentGeometry->GlobalCoordinates(P6, X);

            if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Triangle2D6)
                pNewGeometry = GeometryType::Pointer(new Triangle2D6<BaseType::NodeType>(P1, P2, P3, P4, P5, P6));
            if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Triangle3D6)
                pNewGeometry = GeometryType::Pointer(new Triangle3D6<BaseType::NodeType>(P1, P2, P3, P4, P5, P6));
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The parent geometry type is invalid:", pParentGeometry->GetGeometryType())

        for(std::size_t i = 0; i < pNewGeometry->size(); ++i)
            (*pNewGeometry)[i].SetInitialPosition((*pNewGeometry)[i]);

        return pNewGeometry;
    }

    virtual CoordinatesArrayType ReferenceCenter() const
    {
        CoordinatesArrayType C;
        C[0] = (mX0 + mX1 + mX2) / 3;
        C[1] = (mY0 + mY1 + mY2) / 3;
        C[2] = 0.0;
        return C;
    }

    virtual bool IsInside(const CoordinatesArrayType& rLocalPoint) const
    {
        // REF: https://stackoverflow.com/ques
        // tions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
        double Area = 0.5*(-mY1*mX2 + mY0*(-mX1 + mX2) + mX0*(mY1 - mY2) + mX1*mY2);
        double s = 1/(2*Area)*(mY0*mX2 - mX0*mY2 + (mY2 - mY0)*rLocalPoint[0] + (mX0 - mX2)*rLocalPoint[1]);
        double t = 1/(2*Area)*(mX0*mY1 - mY0*mX1 + (mY0 - mY1)*rLocalPoint[0] + (mX1 - mX0)*rLocalPoint[1]);
        return (s > 0.0) && (s < 1.0) && (t > 0.0) && (t < 1.0) && (s + t < 1.0);
    }

    virtual void CreateSamplingPoints(std::vector<PointType>& SamplingPoints,
            GeometryType& r_geom, const std::size_t& nsampling) const
    {
        SamplingPoints.reserve((nsampling+1)*(nsampling+2) / 2);
        CoordinatesArrayType X;
        PointType P;

        X[2] = 0.0;
        for (std::size_t row = 0; row < nsampling+1; ++row)
        {
            std::size_t n = row + 1;

            double xstart = mX0 + (mX1 - mX0)*row / nsampling;
            double ystart = mY0 + (mY1 - mY0)*row / nsampling;
            double xend = mX0 + (mX2 - mX0)*row / nsampling;
            double yend = mY0 + (mY2 - mY0)*row / nsampling;

            for (std::size_t i = 0; i < n; ++i)
            {
                X[0] = xstart + i*xend/n;
                X[1] = ystart + i*yend/n;
                r_geom.GlobalCoordinates(P, X);
                SamplingPoints.push_back(P);
            }
        }
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "QuadTreeNodeT3";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "(" << mX0 << ", " << mY0 << ") - (" << mX1 << ", " << mY1 << ") - (" << mX2 << ", " << mY2 << ")";
    }

private:
    double mX0, mY0;
    double mX1, mY1;
    double mX2, mY2;
};


/// Tetrahedral quad-tree node in reference coordinates
class QuadTreeNodeT4 : public QuadTreeNode
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadTreeNodeT4);

    typedef QuadTreeNode BaseType;

    typedef BaseType::GeometryType GeometryType;

    QuadTreeNodeT4(const double& X0, const double& Y0, const double& Z0,
                   const double& X1, const double& Y1, const double& Z1,
                   const double& X2, const double& Y2, const double& Z2,
                   const double& X3, const double& Y3, const double& Z3)
    : BaseType(), mX0(X0), mY0(Y0), mZ0(Z0),
        mX1(X1), mY1(Y1), mZ1(Z1),
        mX2(X2), mY2(Y2), mZ2(Z2),
        mX3(X3), mY3(Y3), mZ3(Z3)
    {}

    virtual ~QuadTreeNodeT4() {}

    /// REF: https://www.semanticscholar.org/paper/Octasection-based-Refinement-of-Finite-Element-Endres-Krysl/7fa050663d1b1627413059943b9143f92b98ef4e/figure/4
    virtual void Refine()
    {
        if(this->IsLeaf())
        {
            double X4 = 0.5*(mX0 + mX1), Y4 = 0.5*(mY0 + mY1), Z4 = 0.5*(mZ0 + mZ1);
            double X5 = 0.5*(mX1 + mX2), Y5 = 0.5*(mY1 + mY2), Z5 = 0.5*(mZ1 + mZ2);
            double X6 = 0.5*(mX0 + mX2), Y6 = 0.5*(mY0 + mY2), Z6 = 0.5*(mZ0 + mZ2);
            double X7 = 0.5*(mX0 + mX3), Y7 = 0.5*(mY0 + mY3), Z7 = 0.5*(mZ0 + mZ3);
            double X8 = 0.5*(mX1 + mX3), Y8 = 0.5*(mY1 + mY3), Z8 = 0.5*(mZ1 + mZ3);
            double X9 = 0.5*(mX2 + mX3), Y9 = 0.5*(mY2 + mY3), Z9 = 0.5*(mZ2 + mZ3);

            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(mX0, mY0, mZ0, X4, Y4, Z4, X6, Y6, Z6, X7, Y7, Z7))); // 1
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X4, Y4, Z4, mX1, mY1, mZ1, X5, Y5, Z5, X8, Y8, Z8))); // 2
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X6, Y6, Z6, X5, Y5, Z5, mX2, mY2, mZ2, X9, Y9, Z9))); // 3
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X7, Y7, Z7, X8, Y8, Z8, X9, Y9, Z9, mX3, mY3, mZ3))); // 4
            // // figure 3
            // mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X4, Y4, Z4, X7, Y7, Z7, X8, Y8, Z8, X6, Y6, Z6))); // 5
            // mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X4, Y4, Z4, X5, Y5, Z5, X6, Y6, Z6, X8, Y8, Z8))); // 6
            // mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X9, Y9, Z9, X6, Y6, Z6, X5, Y5, Z5, X8, Y8, Z8))); // 7
            // mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X9, Y9, Z9, X8, Y8, Z8, X7, Y7, Z7, X6, Y6, Z6))); // 8
            // figure 4
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X8, Y8, Z8, X7, Y7, Z7, X9, Y9, Z9, X4, Y4, Z4))); // 5
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X8, Y8, Z8, X5, Y5, Z5, X4, Y4, Z4, X9, Y9, Z9))); // 6
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X6, Y6, Z6, X4, Y4, Z4, X5, Y5, Z5, X9, Y9, Z9))); // 7
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X9, Y9, Z9, X6, Y6, Z6, X4, Y4, Z4, X7, Y7, Z7))); // 8
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
                mpChildren[i]->Refine();
        }
    }

    virtual GeometryType::Pointer pCreateReferenceGeometry() const
    {
        GeometryType::PointType P1(0, mX0, mY0, mZ0);
        GeometryType::PointType P2(1, mX1, mY1, mZ1);
        GeometryType::PointType P3(2, mX2, mY2, mZ2);
        GeometryType::PointType P4(3, mX3, mY3, mZ3);
        return GeometryType::Pointer(new Tetrahedra3D4<BaseType::NodeType>(P1, P2, P3, P4));
    }

    virtual GeometryType::Pointer pCreateGeometry(GeometryType::Pointer pParentGeometry) const
    {
        CoordinatesArrayType X;
        GeometryType::Pointer pNewGeometry;

        if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Tetrahedra3D4)
        {
            GeometryType::PointType P1, P2, P3, P4;

            X[0] = mX0; X[1] = mY0; X[2] = mZ0;
            pParentGeometry->GlobalCoordinates(P1, X);

            X[0] = mX1; X[1] = mY1; X[2] = mZ1;
            pParentGeometry->GlobalCoordinates(P2, X);

            X[0] = mX2; X[1] = mY2; X[2] = mZ2;
            pParentGeometry->GlobalCoordinates(P3, X);

            X[0] = mX3; X[1] = mY3; X[2] = mZ3;
            pParentGeometry->GlobalCoordinates(P4, X);

            pNewGeometry = GeometryType::Pointer(new Tetrahedra3D4<BaseType::NodeType>(P1, P2, P3, P4));
        }
        else if(pParentGeometry->GetGeometryType() == GeometryData::Kratos_Tetrahedra3D10)
        {
            GeometryType::PointType P1, P2, P3, P4, P5, P6, P7, P8, P9, P10;

            X[0] = mX0; X[1] = mY0; X[2] = mZ0;
            pParentGeometry->GlobalCoordinates(P1, X);

            X[0] = mX1; X[1] = mY1; X[2] = mZ1;
            pParentGeometry->GlobalCoordinates(P2, X);

            X[0] = mX2; X[1] = mY2; X[2] = mZ2;
            pParentGeometry->GlobalCoordinates(P3, X);

            X[0] = mX3; X[1] = mY3; X[2] = mZ3;
            pParentGeometry->GlobalCoordinates(P4, X);

            X[0] = 0.5*(mX0+mX1); X[1] = 0.5*(mY0+mY1); X[2] = 0.5*(mZ0+mZ1);
            pParentGeometry->GlobalCoordinates(P5, X);

            X[0] = 0.5*(mX1+mX2); X[1] = 0.5*(mY1+mY2); X[2] = 0.5*(mZ1+mZ2);
            pParentGeometry->GlobalCoordinates(P6, X);

            X[0] = 0.5*(mX2+mX0); X[1] = 0.5*(mY2+mY0); X[2] = 0.5*(mZ2+mZ0);
            pParentGeometry->GlobalCoordinates(P7, X);

            X[0] = 0.5*(mX0+mX3); X[1] = 0.5*(mY0+mY3); X[2] = 0.5*(mZ0+mZ3);
            pParentGeometry->GlobalCoordinates(P8, X);

            X[0] = 0.5*(mX1+mX3); X[1] = 0.5*(mY1+mY3); X[2] = 0.5*(mZ1+mZ3);
            pParentGeometry->GlobalCoordinates(P9, X);

            X[0] = 0.5*(mX2+mX3); X[1] = 0.5*(mY2+mY3); X[2] = 0.5*(mZ2+mZ3);
            pParentGeometry->GlobalCoordinates(P10, X);

            pNewGeometry = GeometryType::Pointer(new Tetrahedra3D10<BaseType::NodeType>(P1, P2, P3, P4, P5, P6, P7, P8, P9, P10));
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The parent geometry type is invalid:", pParentGeometry->GetGeometryType())

        for(std::size_t i = 0; i < pNewGeometry->size(); ++i)
            (*pNewGeometry)[i].SetInitialPosition((*pNewGeometry)[i]);

        return pNewGeometry;
    }

    virtual CoordinatesArrayType ReferenceCenter() const
    {
        CoordinatesArrayType C;
        C[0] = 0.25*(mX0 + mX1 + mX2 + mX3);
        C[1] = 0.25*(mY0 + mY1 + mY2 + mY3);
        C[2] = 0.25*(mZ0 + mZ1 + mZ2 + mZ3);
        return C;
    }

    virtual bool IsInside(const CoordinatesArrayType& rLocalPoint) const
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "Not implemented")
    }

    virtual void CreateSamplingPoints(std::vector<PointType>& SamplingPoints,
            GeometryType& r_geom, const std::size_t& nsampling) const
    {
        SamplingPoints.reserve((nsampling+1)*(nsampling+2)*(nsampling+3) / 6);
        CoordinatesArrayType X;
        PointType P;

        for (std::size_t row = 0; row < nsampling+1; ++row)
        {
            double x0 = mX0 + (mX1 - mX0)*row / nsampling;
            double y0 = mY0 + (mY1 - mY0)*row / nsampling;
            double z0 = mZ0 + (mZ1 - mZ0)*row / nsampling;

            double x1 = mX0 + (mX2 - mX0)*row / nsampling;
            double y1 = mY0 + (mY2 - mY0)*row / nsampling;
            double z1 = mZ0 + (mZ2 - mZ0)*row / nsampling;

            double x2 = mX0 + (mX3 - mX0)*row / nsampling;
            double y2 = mY0 + (mY3 - mY0)*row / nsampling;
            double z2 = mZ0 + (mZ3 - mZ0)*row / nsampling;

            for (std::size_t col = 0; col < row+1; ++col)
            {
                std::size_t n = col + 1;

                double xstart = x0 + (x1 - x0)*col / nsampling;
                double ystart = y0 + (y1 - y0)*col / nsampling;
                double zstart = z0 + (z1 - z0)*col / nsampling;
                double xend = x0 + (x2 - x0)*col / nsampling;
                double yend = y0 + (y2 - y0)*col / nsampling;
                double zend = z0 + (z2 - z0)*col / nsampling;

                for (std::size_t i = 0; i < n; ++i)
                {
                    X[0] = xstart + i*xend/n;
                    X[1] = ystart + i*yend/n;
                    X[2] = zstart + i*zend/n;
                    r_geom.GlobalCoordinates(P, X);
                    SamplingPoints.push_back(P);
                }
            }
        }
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "QuadTreeNodeT4";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "(" << mX0 << ", " << mY0 << ", " << mZ0
             << ") - (" << mX1 << ", " << mY1 << ", " << mZ1
             << ") - (" << mX2 << ", " << mY2 << ", " << mZ2
             << ") - (" << mX3 << ", " << mY3 << ", " << mZ3 << ")";
    }

private:
    double mX0, mY0, mZ0;
    double mX1, mY1, mZ1;
    double mX2, mY2, mZ2;
    double mX3, mY3, mZ3;
};

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, QuadTreeNode& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TDegree>
inline std::ostream& operator << (std::ostream& rOStream, const  QuadTreeNode& rThis)
{
    rOStream << rThis;
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_APPLICATION_QUAD_TREE_NODE_H_INCLUDED  defined
