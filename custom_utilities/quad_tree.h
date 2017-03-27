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


#if !defined(KRATOS_FINITE_CELL_APPLICATION_QUAD_TREE_H_INCLUDED )
#define KRATOS_FINITE_CELL_APPLICATION_QUAD_TREE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes
#include <boost/python.hpp>


// Project includes
#include "custom_utilities/quadrature_utility.h"
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
#include "custom_algebra/brep.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/function/monomial_function.h"
#include "custom_algebra/function/scalar_function.h"
#include "custom_algebra/function/heaviside_function.h"
#include "custom_algebra/function/product_function.h"
#include "custom_geometries/finite_cell_geometry.h"


namespace Kratos
{


/// Representing an abstract quad tree node in reference coordinates
class QuadTreeNode
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadTreeNode);

    typedef typename GeometricalObject::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    QuadTreeNode() {}
    virtual ~QuadTreeNode() {}

    bool IsLeaf() const
    {
        return (mpChildren.size() == 0);
    }

    /*****************************************************************/
    /******* CONSTRUCTION ********************************************/
    /*****************************************************************/

    virtual void Refine()
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    void RefineBy(GeometryType::Pointer& pParentGeometry, const BRep& r_brep)
    {
        if(this->IsLeaf())
        {
            GeometryType::Pointer pThisGeometry = this->pCreateGeometry(pParentGeometry);
            int stat = r_brep.CutStatus(*pThisGeometry);
            if(stat == -1)
            {
                this->Refine();
            }
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
            {
                mpChildren[i]->RefineBy(pParentGeometry, r_brep);
            }
        }
    }

    /*****************************************************************/
    /******* COMPUTATION *********************************************/
    /*****************************************************************/

    /// Integrate a function using the sample geometry and integration rule
    /// The caller has to manually set rOutput to zero before calling this function
    template<typename TOutputType>
    void Integrate(GeometryType::Pointer pParentGeometry,
            const Function<PointType, TOutputType>& rFunc,
            TOutputType& rOutput,
            const GeometryData::IntegrationMethod& ThisIntegrationMethod) const
    {
        if(this->IsLeaf())
        {
            GeometryType::Pointer pThisGeometry = this->pCreateGeometry(pParentGeometry);

            const GeometryType::IntegrationPointsArrayType& integration_points
                = pThisGeometry->IntegrationPoints( ThisIntegrationMethod );

            std::vector<double> DetJ;
            Function<double, double>::ComputeDetJ(DetJ, *pThisGeometry, integration_points);

            CoordinatesArrayType GlobalCoords;

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                pThisGeometry->GlobalCoordinates(GlobalCoords, integration_points[point]);
                rOutput += rFunc.GetValue(GlobalCoords) * DetJ[point] * integration_points[point].Weight();
            }
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
            {
                mpChildren[i]->Integrate(pParentGeometry, rFunc, rOutput, ThisIntegrationMethod);
            }
        }
    }

    /// Integrate a function using the sample geometry and a set of quadrature points. This allows for very
    /// high order quadrature rule, as well as arbitrary quadrature
    /// The caller has to manually set rOutput to zero before calling this function
    template<typename TOutputType>
    void Integrate(GeometryType::Pointer pParentGeometry,
            const Function<PointType, TOutputType>& rFunc,
            TOutputType& rOutput,
            const GeometryType::IntegrationPointsArrayType& integration_points) const
    {
        if(this->IsLeaf())
        {
            GeometryType::Pointer pThisGeometry = this->pCreateGeometry(pParentGeometry);

            std::vector<double> DetJ;
            Function<double, double>::ComputeDetJ(DetJ, *pThisGeometry, integration_points);

            CoordinatesArrayType GlobalCoords;

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                pThisGeometry->GlobalCoordinates(GlobalCoords, integration_points[point]);
                rOutput += rFunc.GetValue(GlobalCoords) * DetJ[point] * integration_points[point].Weight();
            }
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
            {
                mpChildren[i]->Integrate(pParentGeometry, rFunc, rOutput, integration_points);
            }
        }
    }

    /*****************************************************************/
    /******* QUADRATURE GENERATION ***********************************/
    /*****************************************************************/

    /// Construct a fixed quadrature based on order
    virtual GeometryType::IntegrationPointsArrayType ConstructQuadrature(const int& integration_order) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Construct the recursive integration point array
    /// REMARKS: the integration_points is in local coordinates system
    void ConstructQuadrature(GeometryType::Pointer pParentGeometry,
            GeometryType::IntegrationPointsArrayType& integration_points,
            const GeometryData::IntegrationMethod& ThisIntegrationMethod) const
    {
        if(this->IsLeaf())
        {
            GeometryType::Pointer pThisGeometry = this->pCreateGeometry(pParentGeometry);

            GeometryType::Pointer pThisReferenceGeometry = this->pCreateReferenceGeometry();

            const GeometryType::IntegrationPointsArrayType& sub_integration_points
                = pThisReferenceGeometry->IntegrationPoints( ThisIntegrationMethod );

            Vector ShapeFunctionValuesOnReference;
            for(std::size_t point = 0; point < sub_integration_points.size(); ++point)
            {
                GeometryType::IntegrationPointType integration_point;

                ShapeFunctionValuesOnReference = pThisReferenceGeometry->ShapeFunctionsValues(ShapeFunctionValuesOnReference, sub_integration_points[point]);

                noalias(integration_point) = ZeroVector(3);
                for(std::size_t i = 0; i < pThisReferenceGeometry->size(); ++i)
                {
                    noalias(integration_point) += ShapeFunctionValuesOnReference[i] * (*pThisReferenceGeometry)[i];
                }

                double DetJsmall = Function<double, double>::ComputeDetJ(*pThisGeometry, sub_integration_points[point]);

                double DetJbig = Function<double, double>::ComputeDetJ(*pParentGeometry, integration_point);

                integration_point.SetWeight( DetJsmall/DetJbig * sub_integration_points[point].Weight() );

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
    bool CenterOfGravity(PointType& COG, GeometryType::Pointer pParentGeometry, const BRep& r_brep) const
    {
        double X = 0.0, Y = 0.0, Z = 0.0, A = 0.0;
        FunctionR3R1::Pointer FX = FunctionR3R1::Pointer(new MonomialFunctionR3R1<1, 0, 0>());
        FunctionR3R1::Pointer FY = FunctionR3R1::Pointer(new MonomialFunctionR3R1<0, 1, 0>());
        FunctionR3R1::Pointer FZ = FunctionR3R1::Pointer(new MonomialFunctionR3R1<0, 0, 1>());
        FunctionR3R1::Pointer FA = FunctionR3R1::Pointer(new ScalarFunction<FunctionR3R1>(1.0));
        FunctionR3R1::Pointer FH = FunctionR3R1::Pointer(new HeavisideFunction<FunctionR3R1>(r_brep));
        this->Integrate(pParentGeometry, ProductFunction<FunctionR3R1>(FX, FH), X, pParentGeometry->GetDefaultIntegrationMethod());
        this->Integrate(pParentGeometry, ProductFunction<FunctionR3R1>(FY, FH), Y, pParentGeometry->GetDefaultIntegrationMethod());
        this->Integrate(pParentGeometry, ProductFunction<FunctionR3R1>(FZ, FH), Z, pParentGeometry->GetDefaultIntegrationMethod());
        this->Integrate(pParentGeometry, ProductFunction<FunctionR3R1>(FA, FH), A, pParentGeometry->GetDefaultIntegrationMethod());

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

    virtual GeometryType::Pointer pCreateReferenceGeometry() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    virtual GeometryType::Pointer pCreateGeometry(GeometryType::Pointer pParentGeometry) const
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
        const int level) const
{
    Properties::Pointer p_properties = r_model_part.pGetProperties(level);
    GeometryType::Pointer pThisGeometry = this->pCreateGeometry(pParentGeometry);
    Element::Pointer pNewElement = r_clone_element.Create(++lastElementId, pThisGeometry, p_properties);
    r_model_part.AddElement(pNewElement);

    for(std::size_t i = 0; i < pThisGeometry->size(); ++i)
    {
        (*pThisGeometry)[i].SetId(++lastNodeId);
        (*pThisGeometry)[i].SetSolutionStepVariablesList(&r_model_part.GetNodalSolutionStepVariablesList());
        (*pThisGeometry)[i].SetBufferSize(r_model_part.GetBufferSize());
        r_model_part.AddNode((*pThisGeometry)(i));
    }
}


template<>
void QuadTreeNode::AddToModelPart<false, Condition>(GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Condition const& r_clone_condition,
        std::size_t& lastNodeId,
        std::size_t& lastCondId,
        const int level) const
{
    Properties::Pointer p_properties = r_model_part.pGetProperties(level);
    GeometryType::Pointer pThisGeometry = this->pCreateGeometry(pParentGeometry);
    Condition::Pointer pNewCondition = r_clone_condition.Create(++lastCondId, pThisGeometry, p_properties);
    r_model_part.AddCondition(pNewCondition);

    for(std::size_t i = 0; i < pThisGeometry->size(); ++i)
    {
        (*pThisGeometry)[i].SetId(++lastNodeId);
        (*pThisGeometry)[i].SetSolutionStepVariablesList(&r_model_part.GetNodalSolutionStepVariablesList());
        (*pThisGeometry)[i].SetBufferSize(r_model_part.GetBufferSize());
        r_model_part.AddNode((*pThisGeometry)(i));
    }
}


template<>
void QuadTreeNode::AddToModelPart<true, Element>(GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Element const& r_clone_element,
        std::size_t& lastNodeId,
        std::size_t& lastElementId,
        const int level) const
{
    if(this->IsLeaf())
    {
        this->AddToModelPart<false, Element>(pParentGeometry, r_model_part, r_clone_element, lastNodeId, lastElementId, level);
    }
    else
    {
        for(std::size_t i = 0; i < mpChildren.size(); ++i)
            mpChildren[i]->AddToModelPart<true, Element>(pParentGeometry, r_model_part, r_clone_element, lastNodeId, lastElementId, level + 1);
    }
}


template<>
void QuadTreeNode::AddToModelPart<true, Condition>(GeometryType::Pointer pParentGeometry,
        ModelPart& r_model_part,
        Condition const& r_clone_condition,
        std::size_t& lastNodeId,
        std::size_t& lastCondId,
        const int level) const
{
    if(this->IsLeaf())
    {
        this->AddToModelPart<false, Condition>(pParentGeometry, r_model_part, r_clone_condition, lastNodeId, lastCondId, level);
    }
    else
    {
        for(std::size_t i = 0; i < mpChildren.size(); ++i)
            mpChildren[i]->AddToModelPart<true, Condition>(pParentGeometry, r_model_part, r_clone_condition, lastNodeId, lastCondId, level + 1);
    }
}


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

    virtual GeometryType::IntegrationPointsArrayType ConstructQuadrature(const int& integration_order) const
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
            KRATOS_THROW_ERROR(std::logic_error, "This integration_order is not supported:", integration_order);
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
        bool is_onboundary =             ( (fabs(rLocalPoint[0] - mXmin) < tol) && (mYmin - tol < rLocalPoint[1]) && (rLocalPoint[1] < mYmax + tol) );
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


/// Haxehedral quad-tree node in reference coordinates
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
        if(is_onboundary) is_onboundary;
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

            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(mX0, mY0, mZ0, X4, Y4, Z4, X6, Y6, Z6, X7, Y7, Z7)));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(mX1, mY1, mZ1, X5, Y5, Z5, X4, Y4, Z4, X8, Y8, Z8)));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(mX2, mY2, mZ2, X6, Y6, Z6, X5, Y5, Z5, X9, Y9, Z9)));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X7, Y7, Z7, X8, Y8, Z8, X9, Y9, Z9, mX3, mY3, mZ3)));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X4, Y4, Z4, X6, Y6, Z6, X9, Y9, Z9, X7, Y7, Z7)));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X4, Y4, Z4, X5, Y5, Z5, X9, Y9, Z9, X8, Y8, Z8)));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X4, Y4, Z4, X5, Y5, Z5, X6, Y6, Z6, X9, Y9, Z9)));
            mpChildren.push_back(BaseType::Pointer(new QuadTreeNodeT4(X4, Y4, Z4, X9, Y9, Z9, X7, Y7, Z7, X8, Y8, Z8)));
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


/** A general implementation of the quad/oct tree concept for finite cell integration
The quad-tree contains a geometry for high-level identification of the inner points in local coordinates
One quad-tree contains only one quad-tree node
*/
class QuadTree : public QuadratureUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of  QuadTree
    KRATOS_CLASS_POINTER_DEFINITION(QuadTree);

    typedef typename GeometricalObject::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor taking the element geometry and initialize the quad tree as reference coordinates of the element
    QuadTree(Element::Pointer& p_elem)
    {
        this->Initialize(p_elem->pGetGeometry());
    }

    /// Constructor taking the condition geometry and initialize the quad tree as reference coordinates of the condition
    QuadTree(Condition::Pointer& p_cond)
    {
        this->Initialize(p_cond->pGetGeometry());
    }

    /// Constructor taking the geometry and the quad-tree node. This constructor allows for arbitrary reference coordinates of the initial quad-tree node.
    QuadTree(GeometryType::Pointer pGeometry, QuadTreeNode::Pointer pTreeNode)
    : mpThisGeometry(pGeometry), mpTreeNode(pTreeNode)
    {
    }

    /// Destructor.
    virtual ~ QuadTree() {}


    /// Access the underlying quad-tree node
    const QuadTreeNode& Get() const {return *mpTreeNode;}
    QuadTreeNode& Get() {return *mpTreeNode;}
    const QuadTreeNode::Pointer pGet() const {return mpTreeNode;}
    QuadTreeNode::Pointer pGet() {return mpTreeNode;}


    /// Create the sub-cells
    void Refine()
    {
        mpTreeNode->Refine();
    }


    /// Refine the tree by the level set
    void RefineBy(const BRep& r_brep)
    {
        mpTreeNode->RefineBy(mpThisGeometry, r_brep);
    }


    template<class TOutputType>
    TOutputType Integrate(const Function<PointType, TOutputType>& rFunc, const int integration_order) const
    {
        TOutputType Result = TOutputType(0.0);

        if(integration_order <= Function<double, double>::GetMaxIntegrationOrder())
        {
            // use built-in quadrature
            GeometryData::IntegrationMethod ThisIntegrationMethod
                    = Function<double, double>::GetIntegrationMethod(integration_order);
            mpTreeNode->Integrate(mpThisGeometry, rFunc, Result, ThisIntegrationMethod);
        }
        else
        {
            GeometryType::IntegrationPointsArrayType integration_points = mpTreeNode->ConstructQuadrature(integration_order);
            mpTreeNode->Integrate(mpThisGeometry, rFunc, Result, integration_points);
        }

        return Result;
    }


    /// Integrate a function using the sample geometry and integration rule
    /// The caller has to manually set rOutput to zero before calling this function
    template<typename TOutputType>
    void Integrate(const Function<PointType, TOutputType>& rFunc, TOutputType& rOutput,
            const GeometryData::IntegrationMethod& ThisIntegrationMethod) const
    {
        mpTreeNode->Integrate(mpThisGeometry, rFunc, rOutput, ThisIntegrationMethod);
    }


    /// Integrate a function using the sample geometry and integration rule
    /// The caller has to manually set rOutput to zero before calling this function
    template<typename TOutputType>
    void Integrate(const Function<PointType, TOutputType>& rFunc, TOutputType& rOutput,
            const GeometryType::IntegrationPointsArrayType& integration_points) const
    {
        mpTreeNode->Integrate(mpThisGeometry, rFunc, rOutput, integration_points);
    }


    /// Construct the finite cell quadrature
    void ConstructQuadrature(const BRep& r_brep, const int integration_order,
            const double small_weight = 0.0) const
    {
        GeometryType::IntegrationPointsArrayType integration_points;

        GeometryData::IntegrationMethod ThisIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(integration_order);

        // firstly create an array of integration points of sub-trees
        mpTreeNode->ConstructQuadrature(mpThisGeometry, integration_points, ThisIntegrationMethod);

        // modify the weight if needed
        CoordinatesArrayType GlobalCoords;
        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            GlobalCoords = mpThisGeometry->GlobalCoordinates(GlobalCoords, integration_points[point]);
            if(!r_brep.IsInside(GlobalCoords))
                integration_points[point].SetWeight(small_weight);
        }

        /* create new quadrature and assign to the geometry */
        FiniteCellGeometry<GeometryType>::AssignGeometryData(*mpThisGeometry, ThisIntegrationMethod, integration_points);
    }


    /// construct the element out from quad-tree and add to model_part
    /// This is mainly for post-processing
    boost::python::list PyAddToModelPart(ModelPart& r_model_part, const std::string sample_entity_name,
            std::size_t lastNodeId, std::size_t lastEntityId) const
    {
        if( KratosComponents<Element>::Has(sample_entity_name) )
        {
            Element const& r_clone_element = KratosComponents<Element>::Get(sample_entity_name);
            mpTreeNode->AddToModelPart<true, Element>(mpThisGeometry, r_model_part, r_clone_element, lastNodeId, lastEntityId, 1);
        }
        else if( KratosComponents<Condition>::Has(sample_entity_name) )
        {
            Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_entity_name);
            mpTreeNode->AddToModelPart<true, Condition>(mpThisGeometry, r_model_part, r_clone_condition, lastNodeId, lastEntityId, 1);
        }

        boost::python::list list;
        list.append(lastNodeId);
        list.append(lastEntityId);
        return list;
    }


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "QuadTree";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << *mpTreeNode;
    }

    friend std::ostream& operator<<(std::ostream& rOStream, const QuadTree& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << ": ";
        rThis.PrintData(rOStream);
        return rOStream;
    }

private:

    GeometryType::Pointer mpThisGeometry;
    QuadTreeNode::Pointer mpTreeNode;

    void Initialize(GeometryType::Pointer pGeometry)
    {
        mpThisGeometry = pGeometry;
        if(    mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4
            || mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8
            || mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9
            || mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral3D4
            || mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral3D8
            || mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Quadrilateral3D9 )
        {
            mpTreeNode = QuadTreeNode::Pointer(new QuadTreeNodeQ4(-1.0, 1.0, -1.0, 1.0));
        }
        else if(mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Triangle2D3
             || mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Triangle2D6
             || mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Triangle3D3
             || mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Triangle3D6 )
        {
            mpTreeNode = QuadTreeNode::Pointer(new QuadTreeNodeT3(0.0, 0.0, 1.0, 0.0, 0.0, 1.0));
        }
        else if(mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Hexahedra3D8
             || mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Hexahedra3D20
             || mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Hexahedra3D27 )
        {
            mpTreeNode = QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0));
        }
        else if(mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Tetrahedra3D4
             || mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Tetrahedra3D10 )
        {
            mpTreeNode = QuadTreeNode::Pointer(new QuadTreeNodeT4(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0));
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "This geometry type is not supported:", mpThisGeometry->GetGeometryType())
        }
    }

    /// Assignment operator.
     QuadTree& operator=( QuadTree const& rOther);

    /// Copy constructor.
     QuadTree( QuadTree const& rOther);

}; // Class  QuadTree


/** Abstract class for the implementation of the quad/oct tree sub-cell concept for finite cell integration
The quad-tree w/ sub-cell contains a geometry for high-level identification of the inner points in local coordinates
One quad-tree w/ sub-cell can contain multiple quad-tree nodes
The method to construct the sub-cells must be implemented in the sub-class
*/
class QuadTreeSubCell : public QuadratureUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of  QuadTreeSubCell
    KRATOS_CLASS_POINTER_DEFINITION(QuadTreeSubCell);

    typedef typename GeometricalObject::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    QuadTreeSubCell(Element::Pointer& p_elem) : mpThisGeometry(p_elem->pGetGeometry())
    {
    }

    QuadTreeSubCell(Condition::Pointer& p_cond) : mpThisGeometry(p_cond->pGetGeometry())
    {
    }

    /// Destructor.
    virtual ~ QuadTreeSubCell() {}


    /// Get the pointer to the underlying geometry
    GeometryType::Pointer pGetGeometry() const
    {
        return mpThisGeometry;
    }


    /// Access the underlying quad-tree nodes
    const QuadTreeNode& Get(const std::size_t& i) const {return *mpTreeNodes[i];}
    QuadTreeNode& Get(const std::size_t& i) {return *mpTreeNodes[i];}
    const QuadTreeNode::Pointer pGet(const std::size_t& i) const {return mpTreeNodes[i];}
    QuadTreeNode::Pointer pGet(const std::size_t& i) {return mpTreeNodes[i];}


    /// Create the sub-cells
    void Refine()
    {
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            mpTreeNodes[i]->Refine();
    }


    /// Refine the tree by the level set
    void RefineBy(const BRep& r_brep)
    {
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            mpTreeNodes[i]->RefineBy(mpThisGeometry, r_brep);
    }


    template<class TOutputType>
    TOutputType Integrate(const Function<PointType, TOutputType>& rFunc, const int integration_order) const
    {
        GeometryData::IntegrationMethod ThisIntegrationMethod
                = Function<double, double>::GetIntegrationMethod(integration_order);

        TOutputType Result = TOutputType(0.0);
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            mpTreeNodes[i]->Integrate(mpThisGeometry, rFunc, Result, ThisIntegrationMethod);

        return Result;
    }


    /// Integrate a function using the sample geometry and integration rule
    /// The caller has to manually set rOutput to zero before calling this function
    template<typename TOutputType>
    void Integrate(const Function<PointType, TOutputType>& rFunc, TOutputType& rOutput,
            const GeometryData::IntegrationMethod& ThisIntegrationMethod) const
    {
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            mpTreeNodes[i]->Integrate(mpThisGeometry, rFunc, rOutput, ThisIntegrationMethod);
    }


    /// Construct the finite cell quadrature
    void ConstructQuadrature(const BRep& r_brep, const int integration_order,
            const double small_weight = 0.0) const
    {
        GeometryType::IntegrationPointsArrayType integration_points;

        GeometryData::IntegrationMethod ThisIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(integration_order);

        // firstly create an array of integration points of sub-trees of sub-cells
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            mpTreeNodes[i]->ConstructQuadrature(mpThisGeometry, integration_points, ThisIntegrationMethod);

        // modify the weight if needed
        CoordinatesArrayType GlobalCoords;
        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            GlobalCoords = mpThisGeometry->GlobalCoordinates(GlobalCoords, integration_points[point]);
            if(!r_brep.IsInside(GlobalCoords))
                integration_points[point].SetWeight(small_weight);
        }

        /* create new quadrature and assign to the geometry */
        FiniteCellGeometry<GeometryType>::AssignGeometryData(*mpThisGeometry, ThisIntegrationMethod, integration_points);
    }


    /// construct the element out from quad-tree and add to model_part
    /// This is mainly for post-processing
    template<bool TShallow = true>
    boost::python::list PyAddToModelPart(ModelPart& r_model_part, const std::string sample_entity_name,
            std::size_t lastNodeId, std::size_t lastEntityId) const
    {
        if( KratosComponents<Element>::Has(sample_entity_name) )
        {
            Element const& r_clone_element = KratosComponents<Element>::Get(sample_entity_name);
            for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            {
                mpTreeNodes[i]->AddToModelPart<false, Element>(mpThisGeometry, r_model_part, r_clone_element, lastNodeId, lastEntityId, 1);
                if(!TShallow)
                    mpTreeNodes[i]->AddToModelPart<true, Element>(mpThisGeometry, r_model_part, r_clone_element, lastNodeId, lastEntityId, 2);
            }
        }
        else if( KratosComponents<Condition>::Has(sample_entity_name) )
        {
            Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_entity_name);
            for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            {
                mpTreeNodes[i]->AddToModelPart<false, Condition>(mpThisGeometry, r_model_part, r_clone_condition, lastNodeId, lastEntityId, 1);
                if(!TShallow)
                    mpTreeNodes[i]->AddToModelPart<true, Condition>(mpThisGeometry, r_model_part, r_clone_condition, lastNodeId, lastEntityId, 2);
            }
        }

        boost::python::list list;
        list.append(lastNodeId);
        list.append(lastEntityId);
        return list;
    }


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "QuadTreeSubCell";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
        {
            rOStream << "  " << *mpTreeNodes[i] << std::endl;
        }
    }

    friend std::ostream& operator<<(std::ostream& rOStream, const QuadTreeSubCell& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << ":\n";
        rThis.PrintData(rOStream);
        return rOStream;
    }

protected:

    GeometryType::Pointer mpThisGeometry;
    std::vector<QuadTreeNode::Pointer> mpTreeNodes;

    /// Construct the sub-cells for Gauss Quadrature on the quadrilateral. Refer to Gid12 Customization manual, sec 6.1.1 for the identification of the Gauss points
    void ConstructSubCellsForQuadBasedOnGaussQuadrature(std::vector<QuadTreeNode::Pointer>& pTreeNodes,
            const int& integration_order) const
    {
        if(integration_order == 1)
        {
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-1.0, 1.0, -1.0, 1.0)) );
        }
        else if(integration_order == 2)
        {
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-1.0, 0.0, -1.0, 0.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(0.0, 1.0, -1.0, 0.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(0.0, 1.0, 0.0, 1.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-1.0, 0.0, 0.0, 1.0)) );
        }
        else if(integration_order == 3)
        {
            double mid1 = -0.5 * std::sqrt(3.00/5.00);
//            double mid1 = -(2.0 * std::sqrt(3.00/5.00) - 1.0);
            double mid2 = -mid1;

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-1.0, mid1, -1.0, mid1)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(mid2,  1.0, -1.0, mid1)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(mid2,  1.0, mid2,  1.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-1.0, mid1, mid2,  1.0)) );

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(mid1, mid2, -1.0, mid1)) );

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(mid2,  1.0, mid1, mid2)) );

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(mid1, mid2, mid2,  1.0)) );

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-1.0, mid1, mid1, mid2)) );

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(mid1, mid2, mid1, mid2)) );
        }
        else if(integration_order == 4)
        {
            double a = std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 );
            double b = std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 );
            double mid = 0.5*(a + b);

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-1.0, -mid, -1.0, -mid)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-1.0, -mid, -mid,  0.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-1.0, -mid,  0.0,  mid)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-1.0, -mid,  mid,  1.0)) );

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-mid,  0.0, -1.0, -mid)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-mid,  0.0, -mid,  0.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-mid,  0.0,  0.0,  mid)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-mid,  0.0,  mid,  1.0)) );

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4( 0.0,  mid, -1.0, -mid)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4( 0.0,  mid, -mid,  0.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4( 0.0,  mid,  0.0,  mid)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4( 0.0,  mid,  mid,  1.0)) );

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4( mid,  1.0, -1.0, -mid)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4( mid,  1.0, -mid,  0.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4( mid,  1.0,  0.0,  mid)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4( mid,  1.0,  mid,  1.0)) );
        }
        else if(integration_order == 5)
        {
            double a[] = {-0.906179845938664, -0.538469310105683, 0.000000000000000, 0.538469310105683, 0.906179845938664};
            double b[] = {-1.0, 0.5*(a[0]+a[1]), 0.5*(a[1]+a[2]), 0.5*(a[2]+a[3]), 0.5*(a[3]+a[4]), 1.0};

            for(std::size_t i = 0; i < 5; ++i)
                for(std::size_t j = 0; j < 5; ++j)
                    pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(b[i], b[i+1], b[j], b[j+1])) );
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "This integration order is not implemented:", integration_order)
    }

    void ConstructSubCellsForQuadBasedOnEqualDistribution(std::vector<QuadTreeNode::Pointer>& pTreeNodes,
            const std::size_t& m, const std::size_t& n) const
    {
        const double dx = 2.0/m;
        const double dy = 2.0/n;
        for(std::size_t i = 0; i < m; ++i)
        {
            for(std::size_t j = 0; j < n; ++j)
            {
                pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(-1.0+i*dx, -1.0+(i+1)*dx, -1.0+j*dy, -1.0+(j+1)*dy)) );
            }
        }
    }

private:

    /// Assignment operator.
     QuadTreeSubCell& operator=( QuadTreeSubCell const& rOther);

    /// Copy constructor.
     QuadTreeSubCell( QuadTreeSubCell const& rOther);

}; // Class  QuadTreeSubCell

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

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, QuadTree& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TDegree>
inline std::ostream& operator << (std::ostream& rOStream, const  QuadTree& rThis)
{
    rOStream << rThis;
    return rOStream;
}

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, QuadTreeSubCell& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TDegree>
inline std::ostream& operator << (std::ostream& rOStream, const  QuadTreeSubCell& rThis)
{
    rOStream << rThis;
    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_APPLICATION_QUAD_TREE_H_INCLUDED  defined
