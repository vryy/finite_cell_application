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
//  Date:            10 Mar 2017
//


#if !defined(KRATOS_FINITE_CELL_APPLICATION_MOMENT_FITTED_QUAD_TREE_SUBCELL_H_INCLUDED )
#define KRATOS_FINITE_CELL_APPLICATION_MOMENT_FITTED_QUAD_TREE_SUBCELL_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes
#include <boost/python.hpp>


// Project includes
#include "custom_utilities/quadrature_utility.h"
#include "custom_utilities/moment_fitting_utility.h"
#include "includes/model_part.h"
#include "includes/geometrical_object.h"
#include "utilities/math_utils.h"
#include "custom_algebra/brep.h"
#include "custom_algebra/function/function.h"
#include "custom_geometries/finite_cell_geometry.h"
#include "custom_conditions/element_wrapper_condition.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/finite_cell_auxilliary_utility.h"
#include "custom_utilities/finite_cell_geometry_utility.h"

//#define ENABLE_DEBUG_QUADRATURE
//#define DEBUG_SUBCELL

namespace Kratos
{

/// Short class definition.
/** A special type of quad tree w/ sub-cell. It uses quad-tree to compute the moment fitted quadrature of the element. For details see:
Hoang-Giang Bui, D. Schillinger, G. Meschke, Finite Cell Method for Plasticity using Moment-Fitted Quadrature Technique, in preparation.
*/
template<std::size_t TNsampling>
class MomentFittedQuadTreeSubCell : public QuadTreeSubCell<TNsampling>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of  MomentFittedQuadTreeSubCell
    KRATOS_CLASS_POINTER_DEFINITION(MomentFittedQuadTreeSubCell);

    typedef QuadTreeSubCell<TNsampling> BaseType;

    typedef typename GeometricalObject::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    MomentFittedQuadTreeSubCell(Element::Pointer p_elem)
    : BaseType(p_elem), mpElement(p_elem)
    {
    }


    /**
     * Constructor to construct the sub-cells based on Gauss quadrature.
     * In this scheme, the number of sub-cells is equivalent to number of Gauss points on host finite cell. Therefore, in order to increase the number of sub-cells, one must increase the integration order on host cell
     */
    void ConstructSubCellsBasedOnGaussQuadrature(const int& gauss_integration_method)
    {
        int gauss_quadrature_type = QuadratureUtility::GetQuadratureType(gauss_integration_method);

        int gauss_quadrature_order = QuadratureUtility::GetQuadratureOrder(gauss_integration_method);

        BaseType::mpTreeNodes.clear();

        if( mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4
         || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8
         || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9
         || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D4
         || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D8
         || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D9 )
        {
            // create the sub-cells
            BaseType::ConstructSubCellsForQuadBasedOnGaussQuadrature(BaseType::mpTreeNodes, gauss_quadrature_order);
        }
        #ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
        else if(mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier2D
             || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier2D3 )
        {
            KRATOS_THROW_ERROR(std::logic_error, "Not implemented", "")
        }
        #endif
        else if( mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D8
              || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D20
              || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D27 )
        {
            // create the sub-cells
            BaseType::ConstructSubCellsForHexBasedOnGaussQuadrature(BaseType::mpTreeNodes, gauss_quadrature_order);
        }
        #ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
        else if(mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier3D)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Not implemented", "")
        }
        #endif
        else
            KRATOS_THROW_ERROR(std::logic_error, "This geometry type is not supported:", mpElement->GetGeometry().GetGeometryType())

        // using the standard Gauss quadrature as moment fitting integration points
        GeometryData::IntegrationMethod RepresentativeIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(gauss_quadrature_order);

        mRepresentativeIntegrationOrder = gauss_quadrature_order;

        if(gauss_quadrature_type == 0)
        {
            const GeometryType::IntegrationPointsArrayType& integration_points
                = mpElement->GetGeometry().IntegrationPoints( RepresentativeIntegrationMethod );

            mMomentFittingIntegrationPoints = integration_points;
            mRepresentativeIntegrationPoints = integration_points;

            // do a check, to make sure the integration point is inside the sub-cell
            for(std::size_t i = 0; i < BaseType::mpTreeNodes.size(); ++i)
            {
                if(!BaseType::mpTreeNodes[i]->IsInside(integration_points[i]))
                    KRATOS_THROW_ERROR(std::logic_error, "The integration_point is not inside the tree, error at", i)
            }
        }
        else
        {
            // create a big quadtree node for quadrature points generation
            QuadTreeNode::Pointer pTreeNode = QuadTree<TNsampling>::pCreateQuadTreeNode(mpElement->GetGeometry().GetGeometryType());

            const GeometryType::IntegrationPointsArrayType integration_points
                = pTreeNode->ConstructCustomQuadrature(gauss_quadrature_type, gauss_quadrature_order);

            mMomentFittingIntegrationPoints = integration_points;
            mRepresentativeIntegrationPoints = integration_points;

            // do a check, to make sure the integration point is inside the sub-cell
            KRATOS_WATCH(BaseType::mpTreeNodes.size())
            KRATOS_WATCH(integration_points.size())
            for(std::size_t i = 0; i < BaseType::mpTreeNodes.size(); ++i)
            {
                if(!BaseType::mpTreeNodes[i]->IsInside(integration_points[i]))
                {
                    KRATOS_WATCH(typeid(*pTreeNode).name())
                    KRATOS_WATCH(typeid(*BaseType::mpTreeNodes[i]).name())
                    KRATOS_WATCH(*BaseType::mpTreeNodes[i])
                    KRATOS_WATCH(integration_points[i])
                    KRATOS_THROW_ERROR(std::logic_error, "The integration_point is not inside the tree, error at", i)
                }
            }
        }
    }

    /**
     * Constructor to construct the sub-cells based on equal distribution.
     * In this scheme, the number of sub-cells can be arbitrarily
     * In addition, the moment fitting quadrature is chosen based on provided order. It gives maximum control on choosing the fitting functions
     */
    void ConstructSubCellsBasedOnEqualDistribution(const int& integration_order)
    {
        BaseType::mpTreeNodes.clear();

        if( mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4
         || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8
         || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9
         || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D4
         || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D8
         || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D9 )
        {
            std::size_t nsegments = integration_order;

            // create the sub-cells
            BaseType::ConstructSubCellsForQuadBasedOnEqualDistribution(BaseType::mpTreeNodes, nsegments, nsegments);
        }
        #ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
        else if( mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier2D
              || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier2D3 )
        {
            std::size_t nsegments = integration_order;

            // create the sub-cells
            BaseType::ConstructSubCellsForBezier2DBasedOnEqualDistribution(BaseType::mpTreeNodes, nsegments, nsegments);
        }
        #endif
        else if( mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D8
              || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D20
              || mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D27 )
        {
            std::size_t nsegments = integration_order;

            // create the sub-cells
            BaseType::ConstructSubCellsForHexBasedOnEqualDistribution(BaseType::mpTreeNodes, nsegments, nsegments, nsegments);
        }
        #ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
        else if( mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier3D )
        {
            std::size_t nsegments = integration_order;

            // create the sub-cells
            BaseType::ConstructSubCellsForBezier3DBasedOnEqualDistribution(BaseType::mpTreeNodes, nsegments, nsegments, nsegments);
        }
        #endif
        else
            KRATOS_THROW_ERROR(std::logic_error, "This geometry type is not supported:", mpElement->GetGeometry().GetGeometryType())

        // using the standard Gauss quadrature as moment fitting integration points
        GeometryData::IntegrationMethod RepresentativeIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(integration_order);

        mRepresentativeIntegrationOrder = integration_order;

        const GeometryType::IntegrationPointsArrayType& integration_points
            = mpElement->GetGeometry().IntegrationPoints( RepresentativeIntegrationMethod );
//        KRATOS_WATCH(integration_points.size())

        mMomentFittingIntegrationPoints = integration_points;

        for(std::size_t i = 0; i < BaseType::mpTreeNodes.size(); ++i)
        {
            mRepresentativeIntegrationPoints.push_back(BaseType::mpTreeNodes[i]->ReferenceCenter());
        }
    }

    /// Destructor.
    virtual ~ MomentFittedQuadTreeSubCell() {}


    Element::Pointer pGetElement() const
    {
        return mpElement;
    }


    const int& GetRepresentativeIntegrationOrder() const
    {
        return mRepresentativeIntegrationOrder;
    }


    GeometryData::IntegrationMethod GetRepresentativeIntegrationMethod() const
    {
        return Function<double, double>::GetIntegrationMethod(GetRepresentativeIntegrationOrder());
    }


    const GeometryType::IntegrationPointsArrayType& GetRepresentativeIntegrationPoints() const
    {
        return mMomentFittingIntegrationPoints;
    }


    /// Get the number of physical quadrature point. One physical quadrature point represents a subcell.
    std::size_t GetNumberOfPhysicalIntegrationPoint(const BRep& r_brep, const int& integrator_integration_method) const
    {
        std::pair<std::vector<std::size_t>, GeometryType::IntegrationPointsArrayType> Output
            = this->GetPhysicalInterationPoint(r_brep, integrator_integration_method);

        return Output.second.size();
    }


    /// Compute the physical quadrature point. One physical quadrature point represents a subcell.
    std::pair<std::vector<std::size_t>, GeometryType::IntegrationPointsArrayType>
    GetPhysicalInterationPoint(const BRep& r_brep, const int& integrator_integration_method) const
    {
        std::pair<std::vector<std::size_t>, GeometryType::IntegrationPointsArrayType> Output;
        std::vector<std::size_t>& subcell_index = Output.first;
        GeometryType::IntegrationPointsArrayType& physical_integration_points = Output.second;

        PointType COG;
        for(std::size_t i = 0; i < BaseType::mpTreeNodes.size(); ++i)
        {
            #ifdef DEBUG_SUBCELL
            KRATOS_WATCH(typeid(*BaseType::mpTreeNodes[i]).name())
            KRATOS_WATCH(*BaseType::mpTreeNodes[i])
            KRATOS_WATCH(typeid(*BaseType::pGetGeometry()).name())
            #endif

            // GeometryType::Pointer pSubCellGeometry = BaseType::mpTreeNodes[i]->pCreateGeometry(BaseType::pGetGeometry()); // this can be expensive and does not work with Bezier element
            // int status = r_brep.CutStatus(*pSubCellGeometry);
            // #ifdef DEBUG_SUBCELL
            // KRATOS_WATCH(*pSubCellGeometry)
            // #endif

            std::vector<PointType> SamplingPoints;
            BaseType::mpTreeNodes[i]->CreateSamplingPoints(SamplingPoints, *BaseType::pGetGeometry(), TNsampling);
            int status = r_brep.CutStatus(SamplingPoints);

            if(status == BRep::_IN)
            {
                #ifdef DEBUG_SUBCELL
                std::cout << "sub-cell " << i << " is completely inside the physical domain" << std::endl;
                #endif
                // the cell is completely inside the physical domain
                GeometryType::IntegrationPointType integration_point = mRepresentativeIntegrationPoints[i];
                integration_point.Weight() = 0.0; // cancel out the contribution of this integration point to parent element
                physical_integration_points.push_back(integration_point);
                subcell_index.push_back(i);
            }
            else if(status == BRep::_CUT)
            {
                #ifdef DEBUG_SUBCELL
                std::cout << "sub-cell " << i << " is cut" << std::endl;
                #endif
                // the cell is cut
                bool found = BaseType::mpTreeNodes[i]->CenterOfGravity(COG, BaseType::pGetGeometry(), r_brep, integrator_integration_method);
                #ifdef DEBUG_SUBCELL
                KRATOS_WATCH(found)
                KRATOS_WATCH(COG)
                #endif
                if(found)
                {
                    GeometryType::IntegrationPointType integration_point;
                    bool is_inside = BaseType::pGetGeometry()->IsInside(COG, integration_point);
                    if(!is_inside)
                    {
                        std::cout << "!!WARNING!!The CenterOfGravity is not inside the domain of the subcell " << i << std::endl;
                    }
                    integration_point.Weight() = 0.0; // cancel out the contribution of this integration point to parent element
                    physical_integration_points.push_back(integration_point);
                    subcell_index.push_back(i);
                }
            }
            #ifdef DEBUG_SUBCELL
            else
                std::cout << "sub-cell " << i << " is completely outside the physical domain" << std::endl;
            std::cout << "---------------------><--------------------\n\n";
            #endif
        }

        return Output;
    }


    /// Fit a list of subcell using the representative quadrature points from the parent element
    Matrix FitQuadratureSubCell(const std::vector<size_t>& subcell_index,
        const std::vector<FunctionR3R1::Pointer>& r_funcs,
        const BRep& r_brep,
        const int& integrator_integration_method,
        const std::string& solver_type,
        const int& echo_level,
        const double& small_weight) const
    {
        Matrix Weights(subcell_index.size(), mMomentFittingIntegrationPoints.size());
        for(std::size_t i = 0; i < subcell_index.size(); ++i)
        {
            // perform the moment fit for the sub-cell
            std::size_t i_cell = subcell_index[i];

            typename QuadTree<TNsampling>::Pointer quad_tree = typename QuadTree<TNsampling>::Pointer(
                    new QuadTree<TNsampling>(BaseType::pGetGeometry(), BaseType::mpTreeNodes[i_cell]));

            noalias(row(Weights, i)) = MomentFittingUtility::FitQuadrature<FunctionR3R1, QuadTree<TNsampling> >(*BaseType::pGetGeometry(),
                    r_funcs, r_brep, *quad_tree, mMomentFittingIntegrationPoints,
                    integrator_integration_method, solver_type, echo_level, small_weight);
        }

        return Weights;
    }


    /// Create the element out from sub-cells. The element takes the same geometry of the original element, but the weight is given.
    ModelPart::ElementsContainerType CreateSubCellElements(ModelPart& r_model_part,
        const std::string& sample_element_name,
        const Matrix& rWeights,
        std::size_t& lastElementId,
        std::size_t& lastCondId) const
    {
        return this->CreateSubCellElements(r_model_part, sample_element_name, mRepresentativeIntegrationOrder,
                            mMomentFittingIntegrationPoints, rWeights, lastElementId, lastCondId);
    }


    /// Create the element out from sub-cells. The element takes the same geometry of the original element, but the integration point and weight is given.
    ModelPart::ElementsContainerType CreateSubCellElements(ModelPart& r_model_part,
        const std::string& sample_element_name,
        const int& RepresentativeIntegrationOrder,
        const GeometryType::IntegrationPointsArrayType& integration_points,
        const Matrix& rWeights,
        std::size_t& lastElementId,
        std::size_t& lastCondId) const
    {
        /* create the new elements from subcell */
        std::size_t num_physical_point = rWeights.size1();

        Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);
        ModelPart::ElementsContainerType NewElements;
        GeometryType& r_geom = *BaseType::pGetGeometry();

        GeometryData::IntegrationMethod RepresentativeIntegrationMethod
                = Function<double, double>::GetIntegrationMethod(RepresentativeIntegrationOrder);
        Variable<int>& INTEGRATION_ORDER_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("INTEGRATION_ORDER"));

        for(std::size_t i = 0; i < num_physical_point; ++i)
        {
            // create new list of integration points
            GeometryType::IntegrationPointsArrayType new_integration_points = integration_points;
            for(std::size_t j = 0; j < new_integration_points.size(); ++j)
                new_integration_points[j].Weight() = rWeights(i, j);

            // create the new "extrapolated" elements from sub-cell
            // here we make a clone of the geometry because we want to assign different geometry data later on
            // this also works with Bezier element, because Bezier geometry has implemented the Create method
            Element::Pointer pNewElement;
            pNewElement = r_clone_element.Create(++lastElementId, r_geom.Create(r_geom.Points()), mpElement->pGetProperties());

            FiniteCellGeometryUtility::AssignGeometryData(pNewElement->GetGeometry(), RepresentativeIntegrationMethod, new_integration_points);
            pNewElement->SetValue(INTEGRATION_ORDER_var, RepresentativeIntegrationOrder);
            pNewElement->Initialize();
            NewElements.push_back(pNewElement);
        }

        /* create new wrapped conditions and add to the model_part */
        // the purpose is to delay the assembly until the parent element has done its update
        for(typename ModelPart::ElementsContainerType::ptr_iterator it = NewElements.ptr_begin();
                it != NewElements.ptr_end(); ++it)
        {
            // create new element-wrapped condition
            Condition::Pointer pNewCond = Condition::Pointer(new ElementWrapperCondition(++lastCondId, *it));
            r_model_part.Conditions().push_back(pNewCond);
        }

        std::cout << NewElements.size() << " new " << sample_element_name << "-wrapped conditions from parent element " << mpElement->Id()
                  << " are added to the model_part" << std::endl;

        return NewElements;
    }

    /// Create an element taking the same nodes as the original one, but taking the type of geometry from sample_element_name
    Element::Pointer CreateParasiteElement(ModelPart& r_model_part,
        const std::string& sample_element_name, std::size_t& lastElementId,
        Element::Pointer pElement ) const
    {
        Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);

        // REMARK: when creating the element here, the integration rule is not passed. Instead the default integration rule of this element_type is applied, which is not the same as the original element.
        Element::Pointer pNewElement = r_clone_element.Create(++lastElementId, pElement->pGetGeometry(), pElement->pGetProperties());

        std::cout << "1 element of type " << sample_element_name << " is created" << std::endl;

        return pNewElement;
    }


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "MomentFittedQuadTreeSubCell";
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

private:

    #ifdef ENABLE_DEBUG_QUADRATURE
    bool mDebugFlag;
    #endif

    Element::Pointer mpElement;
    int mRepresentativeIntegrationOrder; // this is the integration order associated with mMomentFittingIntegrationPoints
    GeometryType::IntegrationPointsArrayType mMomentFittingIntegrationPoints; // this contains the quadrature points to fit the integration on sub-cell. It is typically the Gauß quadrature point on big cell.
    GeometryType::IntegrationPointsArrayType mRepresentativeIntegrationPoints; // this contains the quadrature points representing a sub-cell. In most of the case, it is not useful. However, when the sub-cell is completely inside the physical domain, it can be used as physical integration point.

    /// Assignment operator.
     MomentFittedQuadTreeSubCell& operator=( MomentFittedQuadTreeSubCell const& rOther);

    /// Copy constructor.
     MomentFittedQuadTreeSubCell( MomentFittedQuadTreeSubCell const& rOther);

}; // Class  MomentFittedQuadTreeSubCell


/// input stream function
template<std::size_t TNsampling>
inline std::istream& operator >> (std::istream& rIStream, MomentFittedQuadTreeSubCell<TNsampling>& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TNsampling>
inline std::ostream& operator << (std::ostream& rOStream, const  MomentFittedQuadTreeSubCell<TNsampling>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#undef ENABLE_DEBUG_QUADRATURE

#endif // KRATOS_FINITE_CELL_APPLICATION_MOMENT_FITTED_QUAD_TREE_SUBCELL_H_INCLUDED defined
