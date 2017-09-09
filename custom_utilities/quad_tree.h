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
#include "custom_utilities/quad_tree_node.h"
#include "custom_utilities/quadrature_utility.h"
#include "custom_utilities/finite_cell_geometry_utility.h"


namespace Kratos
{

/**
Abstract class of all refinable tree
 */
class RefinableTree
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(RefinableTree);

    RefinableTree() {}
    virtual ~RefinableTree() {}

    virtual void Refine()
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class", __FUNCTION__)
    }

    virtual void RefineBy(const BRep& r_brep)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class", __FUNCTION__)
    }
};


/** A general implementation of the quad/oct tree concept for finite cell integration
The quad-tree contains a geometry for high-level identification of the inner points in local coordinates
One quad-tree contains only one quad-tree node
*/
template<std::size_t TNsampling>
class QuadTree : public QuadratureUtility, public RefinableTree
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


    /// Get the template parameter
    const std::size_t GetNsampling() const
    {
        return TNsampling;
    }


    /// Access the underlying quad-tree node
    const QuadTreeNode& Get() const {return *mpTreeNode;}
    QuadTreeNode& Get() {return *mpTreeNode;}
    const QuadTreeNode::Pointer pGet() const {return mpTreeNode;}
    QuadTreeNode::Pointer pGet() {return mpTreeNode;}


    /// Access the underlying (operating) geometry
    GeometryType& GetGeometry() const {return *mpThisGeometry;}
    GeometryType::Pointer pGetGeometry() const {return mpThisGeometry;}


    /// Create the occupied geometry
    GeometryType::Pointer pCreateGeometry() const
    {
        return mpTreeNode->pCreateGeometry(mpThisGeometry);
    }


    /// Create the sub-cells
    /// Implement function from abstract class RefinableTree
    virtual void Refine()
    {
        mpTreeNode->Refine();
    }


    /// Refine the tree by the BRep
    /// Implement function from abstract class RefinableTree
    virtual void RefineBy(const BRep& r_brep)
    {
        if(TNsampling == 0 || TNsampling == 1)
        {
            mpTreeNode->RefineBy(mpThisGeometry, r_brep);
        }
        else
        {
            mpTreeNode->RefineBySampling(mpThisGeometry, r_brep, TNsampling);
        }
    }


    /// Compute the domain size covered by this quadtree
    double DomainSize(const BRep& r_brep, const int& integration_method) const
    {
        return mpTreeNode->DomainSize(mpThisGeometry, r_brep, integration_method);
    }


    /// Compute the center of gravity of this quadtree
    PointType CenterOfGravity(const BRep& r_brep, const int& integration_method) const
    {
        PointType COG;
        mpTreeNode->CenterOfGravity(COG, mpThisGeometry, r_brep, integration_method);
        return COG;
    }


    /// Integrate a function in the global frame using the underlying geometry of the quadtree with the provided integration method
    template<class TOutputType>
    TOutputType IntegrateGlobal(const Function<array_1d<double, 3>, TOutputType>& rFunc, const int& integration_method) const
    {
        TOutputType Result = TOutputType(0.0);

        mpTreeNode->IntegrateGlobal(mpThisGeometry, rFunc, Result, integration_method);

        return Result;
    }


    /// Integrate a function in the global frame using the underlying geometry of the quadtree and integration rule
    /// The caller has to manually set rOutput to zero before calling this function
    template<typename TOutputType>
    void IntegrateGlobal(const Function<array_1d<double, 3>, TOutputType>& rFunc, TOutputType& rOutput,
            const int& integration_method) const
    {
        mpTreeNode->IntegrateGlobal(mpThisGeometry, rFunc, rOutput, integration_method);
    }


    /// Integrate a function in the global frame using the sample geometry and integration rule
    /// The caller has to manually set rOutput to zero before calling this function
    template<typename TOutputType>
    void IntegrateGlobal(const Function<array_1d<double, 3>, TOutputType>& rFunc, TOutputType& rOutput,
            const GeometryType::IntegrationPointsArrayType& integration_points) const
    {
        mpTreeNode->IntegrateGlobal(mpThisGeometry, rFunc, rOutput, integration_points);
    }


    /// Integrate a function in the local frame using the underlying geometry of the quadtree with the provided integration method
    template<class TOutputType>
    TOutputType IntegrateLocal(const Function<array_1d<double, 3>, TOutputType>& rFunc,
            const BRep& r_brep, const int& integration_method, const double& small_weight) const
    {
        TOutputType Result = TOutputType(0.0);

        mpTreeNode->IntegrateLocal(mpThisGeometry, rFunc, r_brep, Result, integration_method, small_weight);

        return Result;
    }


    /// Integrate a function in the local frame using the underlying geometry of the quadtree and integration rule
    /// The caller has to manually set rOutput to zero before calling this function
    template<typename TOutputType>
    void IntegrateLocal(const Function<array_1d<double, 3>, TOutputType>& rFunc, const BRep& r_brep, TOutputType& rOutput,
            const int& integration_method, const double& small_weight) const
    {
        mpTreeNode->IntegrateLocal(mpThisGeometry, rFunc, r_brep, rOutput, integration_method, small_weight);
    }


    /// Integrate a function in the local frame using the sample geometry and integration rule
    /// The caller has to manually set rOutput to zero before calling this function
    template<typename TOutputType>
    void IntegrateLocal(const Function<array_1d<double, 3>, TOutputType>& rFunc, const BRep& r_brep, TOutputType& rOutput,
            const GeometryType::IntegrationPointsArrayType& integration_points, const double& small_weight) const
    {
        mpTreeNode->IntegrateLocal(mpThisGeometry, rFunc, r_brep, rOutput, integration_points, small_weight);
    }


    /// Construct the finite cell quadrature
    /// Returns the number of quadrature points created
    std::size_t ConstructQuadrature(const BRep& r_brep, const int& integration_method,
            const double small_weight = 0.0) const
    {
        GeometryType::IntegrationPointsArrayType integration_points;

        GeometryType::IntegrationPointsArrayType tmp_integration_points;

        // firstly create an array of integration points of sub-trees
        mpTreeNode->ConstructQuadrature(mpThisGeometry, tmp_integration_points, integration_method);

        // fill the integration_point container
        CoordinatesArrayType GlobalCoords;
        if(small_weight != 0.0)
        {
            for(std::size_t point = 0; point < tmp_integration_points.size(); ++point)
            {
                GlobalCoords = mpThisGeometry->GlobalCoordinates(GlobalCoords, tmp_integration_points[point]);
                // modify the weight if needed
                if(!r_brep.IsInside(GlobalCoords))
                    tmp_integration_points[point].SetWeight(small_weight);
                integration_points.push_back(tmp_integration_points[point]);
            }
        }
        else
        {
            for(std::size_t point = 0; point < tmp_integration_points.size(); ++point)
            {
                GlobalCoords = mpThisGeometry->GlobalCoordinates(GlobalCoords, tmp_integration_points[point]);
                if(r_brep.IsInside(GlobalCoords))
                    integration_points.push_back(tmp_integration_points[point]);
            }
        }

        /* create new quadrature and assign to the geometry */
        int quadrature_order = QuadratureUtility::GetQuadratureOrder(integration_method);
        GeometryData::IntegrationMethod ElementalIntegrationMethod = Function<double, double>::GetIntegrationMethod(quadrature_order);
        FiniteCellGeometryUtility::AssignGeometryData(*mpThisGeometry, ElementalIntegrationMethod, integration_points);

        return integration_points.size();
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
        #ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
        else if(mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Bezier2D
            || mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Bezier2D3 )
        {
            mpTreeNode = QuadTreeNode::Pointer(new QuadTreeNodeQ4(0.0, 1.0, 0.0, 1.0));
        }
        else if(mpThisGeometry->GetGeometryType() == GeometryData::Kratos_Bezier3D )
        {
            mpTreeNode = QuadTreeNode::Pointer(new QuadTreeNodeH8(0.0, 1.0, 0.0, 1.0, 0.0, 1.0));
        }
        #endif
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
template<std::size_t TNsampling>
class QuadTreeSubCell : public QuadratureUtility, public RefinableTree
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


    /// Get the number of subcells
    std::size_t NumberOfSubCells() const
    {
        return mpTreeNodes.size();
    }


    /// Get the template parameter
    const std::size_t GetNsampling() const
    {
        return TNsampling;
    }


    /// Access the underlying quad-tree nodes
    const QuadTreeNode& Get(const std::size_t& i) const {return *mpTreeNodes[i];}
    QuadTreeNode& Get(const std::size_t& i) {return *mpTreeNodes[i];}
    const QuadTreeNode::Pointer pGet(const std::size_t& i) const {return mpTreeNodes[i];}
    QuadTreeNode::Pointer pGet(const std::size_t& i) {return mpTreeNodes[i];}


    /// Refine the sub-cells
    /// Implement function from abstract class RefinableTree
    virtual void Refine()
    {
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            mpTreeNodes[i]->Refine();
    }


    /// Refine the tree by the level set
    /// Implement function from abstract class RefinableTree
    virtual void RefineBy(const BRep& r_brep)
    {
        if(TNsampling == 0 || TNsampling == 1)
        {
            for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
                mpTreeNodes[i]->RefineBy(mpThisGeometry, r_brep);
        }
        else
        {
            for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
                mpTreeNodes[i]->RefineBySampling(mpThisGeometry, r_brep, TNsampling);
        }
    }


    /// Compute the domain size of the subcell i
    double DomainSize(const std::size_t& i, const BRep& r_brep, const int& integration_method) const
    {
        return mpTreeNodes[i]->DomainSize(mpThisGeometry, r_brep, integration_method);
    }


    /// Compute the domain size covered by this quadtree (including all subcells)
    double DomainSize(const BRep& r_brep, const int& integration_method) const
    {
        double domain_size = 0.0;
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            domain_size += mpTreeNodes[i]->DomainSize(mpThisGeometry, r_brep, integration_method);
        return domain_size;
    }


    /// Create a quadtree out from a subcell
    /// Fundamentally we can use this function to extract out a subcell and refine it since it's a quadtree node
    typename QuadTree<TNsampling>::Pointer CreateQuadTree(const std::size_t& i) const
    {
        return typename QuadTree<TNsampling>::Pointer(new QuadTree<TNsampling>(mpThisGeometry, mpTreeNodes[i]));
    }


    /// Integrate a function in the global frame using the underlying geometry of the quadtree subcell with the integration order
    template<class TOutputType>
    TOutputType IntegrateGlobal(const Function<array_1d<double, 3>, TOutputType>& rFunc,
            const int& integration_method) const
    {
        TOutputType Result = TOutputType(0.0);

        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            mpTreeNodes[i]->IntegrateGlobal(mpThisGeometry, rFunc, Result, integration_method);

        return Result;
    }


    /// Integrate a function in the global frame using the underlying geometry of the quadtree subcell and integration rule
    /// The caller has to manually set rOutput to zero before calling this function
    template<typename TOutputType>
    void IntegrateGlobal(const Function<array_1d<double, 3>, TOutputType>& rFunc, TOutputType& rOutput,
            const int& integration_method) const
    {
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            mpTreeNodes[i]->IntegrateGlobal(mpThisGeometry, rFunc, rOutput, integration_method);
    }


    /// Integrate a function in the local frame using the underlying geometry of the quadtree subcell with the integration order
    template<class TOutputType>
    TOutputType IntegrateLocal(const Function<array_1d<double, 3>, TOutputType>& rFunc, const BRep& r_brep,
            const int& integration_method, const double& small_weight) const
    {
        TOutputType Result = TOutputType(0.0);

        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            mpTreeNodes[i]->IntegrateLocal(mpThisGeometry, rFunc, r_brep, Result, integration_method, small_weight);

        return Result;
    }


    /// Integrate a function in the local frame using the underlying geometry of the quadtree subcell and integration rule
    /// The caller has to manually set rOutput to zero before calling this function
    template<typename TOutputType>
    void IntegrateLocal(const Function<array_1d<double, 3>, TOutputType>& rFunc, const BRep& r_brep, TOutputType& rOutput,
            const int& integration_method, const double& small_weight) const
    {
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
            mpTreeNodes[i]->IntegrateLocal(mpThisGeometry, rFunc, r_brep, rOutput, integration_method, small_weight);
    }


    /// Construct the finite cell quadrature
    std::size_t ConstructQuadrature(const BRep& r_brep, const int& integration_method,
            const double small_weight = 0.0) const
    {
        GeometryType::IntegrationPointsArrayType integration_points;

        // firstly create an array of integration points of sub-trees of sub-cells
        CoordinatesArrayType GlobalCoords;
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
        {
            GeometryType::IntegrationPointsArrayType tmp_integration_points;

            mpTreeNodes[i]->ConstructQuadrature(mpThisGeometry, tmp_integration_points, integration_method);

            if(small_weight != 0.0)
            {
                for(std::size_t point = 0; point < tmp_integration_points.size(); ++point)
                {
                    GlobalCoords = mpThisGeometry->GlobalCoordinates(GlobalCoords, tmp_integration_points[point]);
                    // modify the weight if needed
                    if(!r_brep.IsInside(GlobalCoords))
                        tmp_integration_points[point].SetWeight(small_weight);
                    integration_points.push_back(tmp_integration_points[point]);
                }
            }
            else
            {
                for(std::size_t point = 0; point < tmp_integration_points.size(); ++point)
                {
                    GlobalCoords = mpThisGeometry->GlobalCoordinates(GlobalCoords, tmp_integration_points[point]);
                    if(r_brep.IsInside(GlobalCoords))
                        integration_points.push_back(tmp_integration_points[point]);
                }
            }
        }

        /* create new quadrature and assign to the geometry */
        int quadrature_order = QuadratureUtility::GetQuadratureOrder(integration_method);
        GeometryData::IntegrationMethod ElementalIntegrationMethod = Function<double, double>::GetIntegrationMethod(quadrature_order);
        FiniteCellGeometryUtility::AssignGeometryData(*mpThisGeometry, ElementalIntegrationMethod, integration_points);

        return integration_points.size();
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

    GeometryType::Pointer mpThisGeometry; // the geometry covers all the subcells
    std::vector<QuadTreeNode::Pointer> mpTreeNodes; // array of subcells

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
        else if(integration_order == 6)
        {
            double a[] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.9324695142031521};
            double b[] = {-1.0, 0.5*(a[0]+a[1]), 0.5*(a[1]+a[2]), 0.5*(a[2]+a[3]), 0.5*(a[3]+a[4]), 0.5*(a[4]+a[5]), 1.0};

            for(std::size_t i = 0; i < 6; ++i)
                for(std::size_t j = 0; j < 6; ++j)
                    pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(b[i], b[i+1], b[j], b[j+1])) );
        }
        else if(integration_order == 7)
        {
            double a[] = {-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0000000000000000, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585};
            double b[] = {-1.0, 0.5*(a[0]+a[1]), 0.5*(a[1]+a[2]), 0.5*(a[2]+a[3]), 0.5*(a[3]+a[4]), 0.5*(a[4]+a[5]), 0.5*(a[5]+a[6]), 1.0};

            for(std::size_t i = 0; i < 7; ++i)
                for(std::size_t j = 0; j < 7; ++j)
                    pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(b[i], b[i+1], b[j], b[j+1])) );
        }
        else if(integration_order == 8)
        {
            double a[] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
            double b[] = {-1.0, 0.5*(a[0]+a[1]), 0.5*(a[1]+a[2]), 0.5*(a[2]+a[3]), 0.5*(a[3]+a[4]), 0.5*(a[4]+a[5]), 0.5*(a[5]+a[6]), 0.5*(a[6]+a[7]), 1.0};

            for(std::size_t i = 0; i < 8; ++i)
                for(std::size_t j = 0; j < 8; ++j)
                    pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(b[i], b[i+1], b[j], b[j+1])) );
        }
        else if(integration_order == 9)
        {
            double a[] = {-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0.0000000000000000, 0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261};
            double b[] = {-1.0, 0.5*(a[0]+a[1]), 0.5*(a[1]+a[2]), 0.5*(a[2]+a[3]), 0.5*(a[3]+a[4]), 0.5*(a[4]+a[5]), 0.5*(a[5]+a[6]), 0.5*(a[6]+a[7]), 0.5*(a[7]+a[8]), 1.0};

            for(std::size_t i = 0; i < 9; ++i)
                for(std::size_t j = 0; j < 9; ++j)
                    pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeQ4(b[i], b[i+1], b[j], b[j+1])) );
        }
        else if(integration_order == 10)
        {
            double a[] = {-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312, 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717};
            double b[] = {-1.0, 0.5*(a[0]+a[1]), 0.5*(a[1]+a[2]), 0.5*(a[2]+a[3]), 0.5*(a[3]+a[4]), 0.5*(a[4]+a[5]), 0.5*(a[5]+a[6]), 0.5*(a[6]+a[7]), 0.5*(a[7]+a[8]), 0.5*(a[8]+a[9]), 1.0};

            for(std::size_t i = 0; i < 10; ++i)
                for(std::size_t j = 0; j < 10; ++j)
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

    /// Construct the sub-cells for Gauss Quadrature on the hexahedra. Refer to Gid12 Customization manual, sec 6.1.1 for the identification of the Gauss points
    void ConstructSubCellsForHexBasedOnGaussQuadrature(std::vector<QuadTreeNode::Pointer>& pTreeNodes,
            const int& integration_order) const
    {
        if(integration_order == 1)
        {
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)) );
        }
        else if(integration_order == 2)
        {
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, 0.0, -1.0, 0.0, -1.0, 0.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(0.0, 1.0, -1.0, 0.0, -1.0, 0.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(0.0, 1.0, 0.0, 1.0, -1.0, 0.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, 0.0, 0.0, 1.0, -1.0, 0.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, 0.0, -1.0, 0.0, 0.0, 1.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(0.0, 1.0, -1.0, 0.0, 0.0, 1.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(0.0, 1.0, 0.0, 1.0, 0.0, 1.0)) );
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, 0.0, 0.0, 1.0, 0.0, 1.0)) );
        }
        else if(integration_order == 3)
        {
            double b = 0.5 * std::sqrt(3.00/5.00);

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, -b, -1.0, -b, -1.0, -b)) );    // 1
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(b, 1.0, -1.0, -b, -1.0, -b)) );      // 2
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(b, 1.0, b, 1.0, -1.0, -b)) );        // 3
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, -b, b, 1.0, -1.0, -b)) );      // 4

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, -b, -1.0, -b, b, 1.0)) );  // 5
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(b, 1.0, -1.0, -b, b, 1.0)) );    // 6
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(b, 1.0, b, 1.0, b, 1.0)) );      // 7
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, -b, b, 1.0, b, 1.0)) );    // 8

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-b, b, -1.0, -b, -1.0, -b)) );   // 9
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(b, 1.0, -b, b, -1.0, -b)) );     // 10
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-b, b, b, 1.0, -1.0, -b)) );     // 11
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, -b, -b, b, -1.0, -b)) );   // 12

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, -b, -1.0, -b, -b, b)) );   // 13
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(b, 1.0, -1.0, -b, -b, b)) );     // 14
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(b, 1.0, b, 1.0, -b, b)) );       // 15
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, -b, b, 1.0, -b, b)) );     // 16

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-b, b, -1.0, -b, b, 1.0)) );     // 17
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(b, 1.0, -b, b, b, 1.0)) );       // 18
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-b, b, b, 1.0, b, 1.0)) );       // 19
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, -b, -b, b, b, 1.0)) );     // 20

            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-b, b, -b, b, -1.0, -b)) );     // 21
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-b, b, -1.0, -b, -b, b)) );     // 22
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(b, 1.0, -b, b, -b, b)) );       // 23
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-b, b, b, 1.0, -b, b)) );       // 24
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0, -b, -b, b, -b, b)) );     // 25
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-b, b, -b, b, b, 1.0)) );       // 26
            pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-b, b, -b, b, -b, b)) );        // 27
        }
        else if(integration_order == 4)
        {
            double a = 0.86113631159405257522;
            double b = 0.33998104358485626480;

            // -1       [-a]    -0.5*(a+b)      [-b]    0.0      [b]   0.5*(a+b)    [a]       1

            double p[] = {-1.0, -0.5*(a+b), 0.0, 0.5*(a+b), 1.0};
            for(int i = 0; i < 4; ++i)
            {
                for(int j = 0; j < 4; ++j)
                {
                    for(int k = 0; k < 4; ++k)
                    {
                        pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(p[k], p[k+1], p[j], p[j+1], p[i], p[i+1])) );
                    }
                }
            }
        }
        else if(integration_order == 5)
        {
            double a = 0.90617984593866399280;
            double b = 0.53846931010568309104;

            // -1       [-a]    -0.5*(a+b)      [-b]   -0.5*b    [0.0]   0.5*b   [b]   0.5*(a+b)    [a]       1

            double p[] = {-1.0, -0.5*(a+b), -0.5*b, 0.5*b, 0.5*(a+b), 1.0};
            for(int i = 0; i < 5; ++i)
            {
                for(int j = 0; j < 5; ++j)
                {
                    for(int k = 0; k < 5; ++k)
                    {
                        pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(p[k], p[k+1], p[j], p[j+1], p[i], p[i+1])) );
                    }
                }
            }
        }
        else if(integration_order == 6)
        {
            double a = 0.9324695142031520278123;
            double b = 0.6612093864662645136614;
            double c = 0.2386191860831969086305;

            // -1       [-a]    -0.5*(a+b)     [-b]   -0.5*(b+c)    [-c]   0.0   [c]   0.5*(b+c)    [b]   0.5*(a+b)    [a]       1

            double p[] = {-1.0, -0.5*(a+b), -0.5*(b+c), 0.0, 0.5*(b+c), 0.5*(a+b), 1.0};
            for(int i = 0; i < 6; ++i)
            {
                for(int j = 0; j < 6; ++j)
                {
                    for(int k = 0; k < 6; ++k)
                    {
                        pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(p[i], p[i+1], p[j], p[j+1], p[k], p[k+1])) );
                    }
                }
            }
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "This integration order is not implemented:", integration_order)
    }

    void ConstructSubCellsForHexBasedOnEqualDistribution(std::vector<QuadTreeNode::Pointer>& pTreeNodes,
            const std::size_t& m, const std::size_t& n, const std::size_t& p) const
    {
        const double dx = 2.0/m;
        const double dy = 2.0/n;
        const double dz = 2.0/p;
        for(std::size_t i = 0; i < m; ++i)
        {
            for(std::size_t j = 0; j < n; ++j)
            {
                for(std::size_t k = 0; k < p; ++k)
                {
                    pTreeNodes.push_back( QuadTreeNode::Pointer(new QuadTreeNodeH8(-1.0+i*dx, -1.0+(i+1)*dx, -1.0+j*dy, -1.0+(j+1)*dy, -1.0+k*dz, -1.0+(k+1)*dz)) );
                }
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
template<std::size_t TNsampling>
inline std::istream& operator >> (std::istream& rIStream, QuadTree<TNsampling>& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TNsampling>
inline std::ostream& operator << (std::ostream& rOStream, const QuadTree<TNsampling>& rThis)
{
    rOStream << rThis;
    return rOStream;
}

/// input stream function
template<std::size_t TNsampling>
inline std::istream& operator >> (std::istream& rIStream, QuadTreeSubCell<TNsampling>& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TNsampling>
inline std::ostream& operator << (std::ostream& rOStream, const  QuadTreeSubCell<TNsampling>& rThis)
{
    rOStream << rThis;
    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_APPLICATION_QUAD_TREE_H_INCLUDED  defined
