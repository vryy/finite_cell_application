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


#if !defined(KRATOS_BINARY_TREE_H_INCLUDED )
#define  KRATOS_BINARY_TREE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "geometries/line_3d_2.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "custom_algebra/function.h"
#include "custom_algebra/level_set.h"
#include "custom_geometries/finite_cell_geometry.h"


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
/** A simple implementation of the binary tree concept in 1D (in 2D it would be quad tree and in 3D would be oct tree)
*/
template<std::size_t TDegree>
class BinaryTree
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BinaryTree
    KRATOS_CLASS_POINTER_DEFINITION(BinaryTree);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.

    BinaryTree(Element::Pointer& p_elem)
    {
        mpChildren.resize(0);

        if(    p_elem->GetGeometry().GetGeometryType() != GeometryData::Kratos_Quadrilateral2D4
            && p_elem->GetGeometry().GetGeometryType() != GeometryData::Kratos_Quadrilateral3D4
            && p_elem->GetGeometry().GetGeometryType() != GeometryData::Kratos_Hexahedra3D8
            && p_elem->GetGeometry().GetGeometryType() != GeometryData::Kratos_Line2D2
            && p_elem->GetGeometry().GetGeometryType() != GeometryData::Kratos_Line3D2 )
        {
            KRATOS_THROW_ERROR(std::logic_error, "The input geometry is incompatible", "")
        }

        mpThisGeometry = p_elem->pGetGeometry();
    }

    BinaryTree(NodeType::Pointer& P1, NodeType::Pointer& P2)
    {
        mpChildren.resize(0);
        mpThisGeometry = GeometryType::Pointer(new Line3D2<NodeType>(P1, P2));
    }


    BinaryTree(NodeType::Pointer& P1, NodeType::Pointer& P2, NodeType::Pointer& P3, NodeType::Pointer& P4)
    {
        mpChildren.resize(0);
        mpThisGeometry = GeometryType::Pointer(new Quadrilateral3D4<NodeType>(P1, P2, P3, P4));
    }


    BinaryTree(NodeType::Pointer& P1, NodeType::Pointer& P2, NodeType::Pointer& P3, NodeType::Pointer& P4,
        NodeType::Pointer& P5, NodeType::Pointer& P6, NodeType::Pointer& P7, NodeType::Pointer& P8)
    {
        mpChildren.resize(0);
        mpThisGeometry = GeometryType::Pointer(new Hexahedra3D8<NodeType>(P1, P2, P3, P4, P5, P6, P7, P8));
    }


    /// Destructor.
    virtual ~BinaryTree() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /*****************************************************************/
    /******* INFORMATION *********************************************/
    /*****************************************************************/


    const std::size_t WorkingSpaceDimension() const {return TDegree;}


    const std::size_t TreeDegree() const
    {
        std::size_t degree = 1;
        for(std::size_t i = 0; i < TDegree; ++i)
            degree *= 2;
        return degree;
    }


    GeometryType::Pointer pGetGeometry() const
    {
        return mpThisGeometry;
    }


    GeometryType& GetGeometry() const
    {
        return *mpThisGeometry;
    }


    std::size_t Id() const
    {
        return mId;
    }


    bool IsLeaf() const
    {
        return (mpChildren.size() == 0);
    }


    /*****************************************************************/
    /******* CONSTRUCTION ********************************************/
    /*****************************************************************/


    /// Create the sub-cells
    void Refine()
    {
        if(this->IsLeaf())
        {
            std::size_t tdegree = this->TreeDegree();
            mpChildren.resize(tdegree);

            if(TDegree == 1)
            {
                NodeType::Pointer P2 = this->CreateNode(0.5, (*mpThisGeometry)[0], 0.5, (*mpThisGeometry)[1]);
                mpChildren[0] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>((*mpThisGeometry)(0), P2));
                mpChildren[1] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>(P2, (*mpThisGeometry)(1)));
            }
            else if(TDegree == 2)
            {
                NodeType::Pointer P4 = this->CreateNode(0.5, (*mpThisGeometry)[0], 0.5, (*mpThisGeometry)[1]);
                NodeType::Pointer P5 = this->CreateNode(0.5, (*mpThisGeometry)[1], 0.5, (*mpThisGeometry)[2]);
                NodeType::Pointer P6 = this->CreateNode(0.5, (*mpThisGeometry)[2], 0.5, (*mpThisGeometry)[3]);
                NodeType::Pointer P7 = this->CreateNode(0.5, (*mpThisGeometry)[3], 0.5, (*mpThisGeometry)[0]);
                NodeType::Pointer P8 = this->CreateNode(0.25, (*mpThisGeometry)[0], 0.25, (*mpThisGeometry)[1],
                                0.25, (*mpThisGeometry)[2], 0.25, (*mpThisGeometry)[3]);

                mpChildren[0] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>((*mpThisGeometry)(0), P4, P8, P7));
                mpChildren[1] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>((*mpThisGeometry)(1), P5, P8, P4));
                mpChildren[2] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>((*mpThisGeometry)(2), P6, P8, P5));
                mpChildren[3] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>((*mpThisGeometry)(3), P7, P8, P6));
            }
            else if(TDegree == 3)
            {
                NodeType::Pointer P8 = this->CreateNode(0.5, (*mpThisGeometry)[0], 0.5, (*mpThisGeometry)[1]);
                NodeType::Pointer P9 = this->CreateNode(0.5, (*mpThisGeometry)[1], 0.5, (*mpThisGeometry)[2]);
                NodeType::Pointer P10 = this->CreateNode(0.5, (*mpThisGeometry)[2], 0.5, (*mpThisGeometry)[3]);
                NodeType::Pointer P11 = this->CreateNode(0.5, (*mpThisGeometry)[3], 0.5, (*mpThisGeometry)[0]);

                NodeType::Pointer P12 = this->CreateNode(0.5, (*mpThisGeometry)[0], 0.5, (*mpThisGeometry)[4]);
                NodeType::Pointer P13 = this->CreateNode(0.5, (*mpThisGeometry)[1], 0.5, (*mpThisGeometry)[5]);
                NodeType::Pointer P14 = this->CreateNode(0.5, (*mpThisGeometry)[2], 0.5, (*mpThisGeometry)[6]);
                NodeType::Pointer P15 = this->CreateNode(0.5, (*mpThisGeometry)[3], 0.5, (*mpThisGeometry)[7]);

                NodeType::Pointer P16 = this->CreateNode(0.5, (*mpThisGeometry)[4], 0.5, (*mpThisGeometry)[5]);
                NodeType::Pointer P17 = this->CreateNode(0.5, (*mpThisGeometry)[5], 0.5, (*mpThisGeometry)[6]);
                NodeType::Pointer P18 = this->CreateNode(0.5, (*mpThisGeometry)[6], 0.5, (*mpThisGeometry)[7]);
                NodeType::Pointer P19 = this->CreateNode(0.5, (*mpThisGeometry)[7], 0.5, (*mpThisGeometry)[4]);

                NodeType::Pointer P20 = this->CreateNode(0.25, (*mpThisGeometry)[0], 0.25, (*mpThisGeometry)[1], 0.25, (*mpThisGeometry)[2], 0.25, (*mpThisGeometry)[3]);
                NodeType::Pointer P21 = this->CreateNode(0.25, (*mpThisGeometry)[0], 0.25, (*mpThisGeometry)[1], 0.25, (*mpThisGeometry)[5], 0.25, (*mpThisGeometry)[4]);
                NodeType::Pointer P22 = this->CreateNode(0.25, (*mpThisGeometry)[1], 0.25, (*mpThisGeometry)[2], 0.25, (*mpThisGeometry)[6], 0.25, (*mpThisGeometry)[5]);
                NodeType::Pointer P23 = this->CreateNode(0.25, (*mpThisGeometry)[3], 0.25, (*mpThisGeometry)[2], 0.25, (*mpThisGeometry)[6], 0.25, (*mpThisGeometry)[7]);
                NodeType::Pointer P24 = this->CreateNode(0.25, (*mpThisGeometry)[0], 0.25, (*mpThisGeometry)[3], 0.25, (*mpThisGeometry)[7], 0.25, (*mpThisGeometry)[4]);
                NodeType::Pointer P25 = this->CreateNode(0.25, (*mpThisGeometry)[4], 0.25, (*mpThisGeometry)[5], 0.25, (*mpThisGeometry)[6], 0.25, (*mpThisGeometry)[7]);

                NodeType::Pointer P26 = this->CreateNode(0.125, (*mpThisGeometry)[0], 0.125, (*mpThisGeometry)[1], 0.125, (*mpThisGeometry)[2], 0.125, (*mpThisGeometry)[3]
                                        , 0.125, (*mpThisGeometry)[4], 0.125, (*mpThisGeometry)[5], 0.125, (*mpThisGeometry)[6], 0.125, (*mpThisGeometry)[7]);

                mpChildren[0] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>((*mpThisGeometry)(0), P8, P20, P11, P12, P21, P26, P24));
                mpChildren[1] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>((*mpThisGeometry)(1), P9, P20, P8, P13, P22, P26, P21));
                mpChildren[2] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>((*mpThisGeometry)(2), P10, P20, P9, P14, P23, P26, P22));
                mpChildren[3] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>((*mpThisGeometry)(3), P11, P20, P10, P15, P24, P26, P23));
                mpChildren[4] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>(P12, P21, P26, P24, (*mpThisGeometry)(4), P16, P25, P19));
                mpChildren[5] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>(P13, P22, P26, P21, (*mpThisGeometry)(5), P17, P25, P16));
                mpChildren[6] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>(P14, P23, P26, P22, (*mpThisGeometry)(6), P18, P25, P17));
                mpChildren[7] = BinaryTree<TDegree>::Pointer(new BinaryTree<TDegree>(P15, P24, P26, P23, (*mpThisGeometry)(7), P19, P25, P18));
            }
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
                mpChildren[i]->Refine();
        }
    }


    /// Refine the tree by the level set
    void RefineBy(const LevelSet& r_level_set)
    {
        if(this->IsLeaf())
        {
            int stat = r_level_set.CutStatus(this->GetGeometry());
            if(stat == -1)
            {
                this->Refine();
            }
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
            {
                mpChildren[i]->RefineBy(r_level_set);
            }
        }
    }


    /*****************************************************************/
    /******* COMPUTATION *********************************************/
    /*****************************************************************/


    template<class TOutputType>
    TOutputType Integrate(const Function<PointType, TOutputType>& rFunc, const int integration_order) const
    {
        GeometryData::IntegrationMethod ThisIntegrationMethod
                = LevelSet::GetIntegrationMethod(integration_order);

        TOutputType Result = TOutputType(0.0);
        this->Integrate(rFunc, Result, ThisIntegrationMethod);

        return Result;
    }


    /// Integrate a function using the sample geometry and integration rule
    /// The caller has to manually set rOutput to zero before calling this function
    template<typename TOutputType>
    void Integrate(const Function<PointType, TOutputType>& rFunc, TOutputType& rOutput,
            const GeometryData::IntegrationMethod ThisIntegrationMethod) const
    {
        if(this->IsLeaf())
        {
            const GeometryType::IntegrationPointsArrayType& integration_points
                = this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );

            std::vector<double> DetJ;
            this->ComputeDetJ(DetJ, this->GetGeometry(), integration_points);

            CoordinatesArrayType GlobalCoords;

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                this->GetGeometry().GlobalCoordinates(GlobalCoords, integration_points[point]);
                rOutput += rFunc.GetValue(GlobalCoords) * DetJ[point] * integration_points[point].Weight();
            }
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
            {
                mpChildren[i]->Integrate(rFunc, rOutput, ThisIntegrationMethod);
            }
        }
    }


    /*****************************************************************/
    /******* POST PROCESSING *****************************************/
    /*****************************************************************/


    void ResetId()
    {
        if(this->IsLeaf())
        {
            mId = 0;
        }
        else
        {
            mId = 0;
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
                mpChildren[i]->ResetId();
        }
    }


    void Renumber(std::size_t& LastNodeId, std::size_t& LastElementId)
    {
        if(this->IsLeaf())
        {
            mId = ++LastElementId;
            for(std::size_t i = 0; i < GetGeometry().size(); ++i)
                if(GetGeometry()[i].Id() == 0)
                    GetGeometry()[i].SetId(++LastNodeId);
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
                mpChildren[i]->Renumber(LastNodeId, LastElementId);
        }
    }


    boost::python::list PyRenumber(std::size_t LastNodeId, std::size_t LastElementId)
    {
        this->ResetId();
        this->Renumber(LastNodeId, LastElementId);

        boost::python::list list;
        list.append(LastNodeId);
        list.append(LastElementId);
        return list;
    }


    void AddToModelPart(ModelPart& r_model_part, Element const& r_sample_element,
            const std::size_t level) const
    {
        if(this->IsLeaf())
        {
            Properties::Pointer p_properties = r_model_part.pGetProperties(level);
            Element::Pointer pNewElement = r_sample_element.Create(Id(), mpThisGeometry, p_properties);
            r_model_part.AddElement(pNewElement);

            for(std::size_t i = 0; i < mpThisGeometry->size(); ++i)
                r_model_part.AddNode((*mpThisGeometry)(i));
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
                mpChildren[i]->AddToModelPart(r_model_part, r_sample_element, level + 1);
        }
    }


    void PyAddToModelPart(ModelPart& r_model_part, const std::string sample_element_name) const
    {
        Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);

        this->AddToModelPart(r_model_part, r_clone_element, 1);

        r_model_part.Nodes().Unique();
    }


    /*****************************************************************/
    /******* AUXILIARY ROUTINES **************************************/
    /*****************************************************************/


    void ConstructQuadrature(const int integration_order) const
    {
        GeometryType::IntegrationPointsArrayType integration_points;

        GeometryData::IntegrationMethod ThisIntegrationMethod
                = LevelSet::GetIntegrationMethod(integration_order);

        // firstly create an array of integration points of sub-trees
        this->ConstructQuadrature(integration_points, ThisIntegrationMethod);

        std::vector<double> DetJ(integration_points.size());
        this->ComputeDetJ(DetJ, this->GetGeometry(), integration_points);

        // scale the integration_point
        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            integration_points[point].SetWeight(integration_points[point].Weight() / DetJ[point]);
        }

        /* create new quadrature and assign to the geometry */
//        FiniteCellGeometry<GeometryType>::AssignGeometryData(this->GetGeometry(), ThisIntegrationMethod, integration_points);
    }


    /// Construct the recursive integration point array
    void ConstructQuadrature(GeometryType::IntegrationPointsArrayType& integration_points,
            const GeometryData::IntegrationMethod& ThisIntegrationMethod) const
    {
        if(this->IsLeaf())
        {
            const GeometryType::IntegrationPointsArrayType& sub_integration_points
                = this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );

            std::vector<double> DetJ(sub_integration_points.size());
            this->ComputeDetJ(DetJ, this->GetGeometry(), sub_integration_points);

            for(std::size_t point = 0; point < sub_integration_points.size(); ++point)
            {
                GeometryType::IntegrationPointType integration_point = sub_integration_points[point];
                integration_point.SetWeight( DetJ[point] * integration_point.Weight() );
                integration_points.push_back( integration_point );
            }
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
            {
                mpChildren[i]->ConstructQuadrature(integration_points, ThisIntegrationMethod);
            }
        }
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
        if(TDegree == 1)
            return "Binary Tree";
        else if(TDegree == 2)
            return "Quad Tree";
        else if(TDegree == 3)
            return "Oct Tree";
        else
        {
            std::stringstream ss;
            ss << TreeDegree() << "-Tree";
            return ss.str();
        }
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        if(this->IsLeaf())
        {
//            rOStream << this->GetGeometry();
            rOStream << Id() << " :";
            for(std::size_t i = 0; i < this->GetGeometry().size(); ++i)
            {
                NodeType& node = this->GetGeometry()[i];
                rOStream << " (" << node.Id() << ": " << node.X0() << ", " << node.Y0() << ", " << node.Z0() << ")";
            }
        }
        else
        {
            for(std::size_t i = 0; i < mpChildren.size(); ++i)
            {
                rOStream << "  ";
                mpChildren[i]->PrintData(rOStream);
                rOStream << std::endl << "  ";
            }
        }
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


    std::vector<BinaryTree<TDegree>::Pointer> mpChildren;
    GeometryType::Pointer mpThisGeometry;

    std::size_t mId;

    NodeType::Pointer CreateNode(const double& alpha1, const NodeType& N1,
            const double& alpha2, const NodeType& N2)
    {
        NodeType::Pointer N = NodeType::Pointer( new NodeType(0,
                    alpha1*N1.X0() + alpha2*N2.X0(),
                    alpha1*N1.Y0() + alpha2*N2.Y0(),
                    alpha1*N1.Z0() + alpha2*N2.Z0()) );
        return N;
    }

    NodeType::Pointer CreateNode(const double& alpha1, const NodeType& N1,
            const double& alpha2, const NodeType& N2,
            const double& alpha3, const NodeType& N3,
            const double& alpha4, const NodeType& N4)
    {
        NodeType::Pointer N = NodeType::Pointer( new NodeType(0,
                    alpha1*N1.X0() + alpha2*N2.X0() + alpha3*N3.X0() + alpha4*N4.X0(),
                    alpha1*N1.Y0() + alpha2*N2.Y0() + alpha3*N3.Y0() + alpha4*N4.Y0(),
                    alpha1*N1.Z0() + alpha2*N2.Z0() + alpha3*N3.Z0() + alpha4*N4.Z0()) );
        return N;
    }

    NodeType::Pointer CreateNode(const double& alpha1, const NodeType& N1,
            const double& alpha2, const NodeType& N2,
            const double& alpha3, const NodeType& N3,
            const double& alpha4, const NodeType& N4,
            const double& alpha5, const NodeType& N5,
            const double& alpha6, const NodeType& N6,
            const double& alpha7, const NodeType& N7,
            const double& alpha8, const NodeType& N8)
    {
        NodeType::Pointer N = NodeType::Pointer( new NodeType(0,
                    alpha1*N1.X0() + alpha2*N2.X0() + alpha3*N3.X0() + alpha4*N4.X0() + alpha5*N5.X0() + alpha6*N6.X0() + alpha7*N7.X0() + alpha8*N8.X0(),
                    alpha1*N1.Y0() + alpha2*N2.Y0() + alpha3*N3.Y0() + alpha4*N4.Y0() + alpha5*N5.Y0() + alpha6*N6.Y0() + alpha7*N7.Y0() + alpha8*N8.Y0(),
                    alpha1*N1.Z0() + alpha2*N2.Z0() + alpha3*N3.Z0() + alpha4*N4.Z0() + alpha5*N5.Z0() + alpha6*N6.Z0() + alpha7*N7.Z0() + alpha8*N8.Z0()) );
        return N;
    }


    void ComputeDetJ(std::vector<double>& DetJ,
            GeometryType& r_geom, const GeometryType::IntegrationPointsArrayType& integration_points) const
    {
        if(DetJ.size() != integration_points.size())
            DetJ.resize(integration_points.size());

        if(r_geom.WorkingSpaceDimension() == r_geom.LocalSpaceDimension())
        {
            Matrix J;

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                J = r_geom.Jacobian( J, integration_points[point] );
                DetJ.push_back( MathUtils<double>::Det(J) );
            }
        }
        else
        {
            Matrix J, JtJ;

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                J = r_geom.Jacobian( J, integration_points[point] );
                JtJ = prod(trans(J), J);
                DetJ.push_back( sqrt(MathUtils<double>::Det(JtJ)) );
            }
        }
    }

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
    BinaryTree& operator=(BinaryTree const& rOther);

    /// Copy constructor.
    BinaryTree(BinaryTree const& rOther);


    ///@}

}; // Class BinaryTree

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<std::size_t TDegree>
inline std::istream& operator >> (std::istream& rIStream,
                BinaryTree<TDegree>& rThis)
{}

/// output stream function
template<std::size_t TDegree>
inline std::ostream& operator << (std::ostream& rOStream,
                const BinaryTree<TDegree>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BINARY_TREE_H_INCLUDED  defined
