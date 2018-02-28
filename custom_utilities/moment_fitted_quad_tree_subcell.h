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
#include <tuple>


// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/model_part.h"
#include "includes/geometrical_object.h"
#include "includes/serializer.h"
#include "utilities/math_utils.h"
#include "utilities/timer.h"
#include "custom_algebra/brep.h"
#include "custom_algebra/function/function.h"
#include "custom_geometries/finite_cell_geometry.h"
#include "custom_conditions/element_wrapper_condition.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/quadrature_utility.h"
#include "custom_utilities/moment_fitting_utility.h"
#include "custom_utilities/finite_cell_auxilliary_utility.h"
#include "custom_utilities/finite_cell_geometry_utility.h"
#include "custom_utilities/finite_cell_mesh_utility.h"

//#define ENABLE_DEBUG_QUADRATURE
//#define DEBUG_SUBCELL
//#define ENABLE_PROFILING

namespace Kratos
{

/**
Abstract class for moment fitted quadtree w/ subcell
 */
class BaseMomentFittedQuadTreeSubCell
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BaseMomentFittedQuadTreeSubCell);

    typedef typename GeometricalObject::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    BaseMomentFittedQuadTreeSubCell() {} // Empty constructor for serializer

    BaseMomentFittedQuadTreeSubCell(Element::Pointer p_elem) : mpElement(p_elem) {}

    virtual ~BaseMomentFittedQuadTreeSubCell() {}

    virtual void GeneratePhysicalIntegrationPoints(const BRep& r_brep, const int& integrator_integration_method)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the underlying element
    Element::Pointer pGetElement() const
    {
        return mpElement;
    }

    /// Get the integration order of fitting integration points
    const int& GetRepresentativeIntegrationOrder() const
    {
        return mRepresentativeIntegrationOrder;
    }


    /// Get the integration method of fitting integration points
    GeometryData::IntegrationMethod GetRepresentativeIntegrationMethod() const
    {
        return Function<double, double>::GetIntegrationMethod(GetRepresentativeIntegrationOrder());
    }


    /// Get the list of integration points used to fit each subcell
    const IntegrationPointsArrayType& GetRepresentativeIntegrationPoints() const
    {
        return mMomentFittingIntegrationPoints;
    }


    /// Get the physical integration points
    const IntegrationPointsArrayType& GetPhysicalIntegrationPoints() const
    {
        return mPhysicalIntegrationPoints;
    }


    /// Get index of the sub-cells
    const std::vector<std::size_t>& SubCellIndices() const
    {
        return mSubCellIndices;
    }


    /// Get the fictitious integration points
    const IntegrationPointsArrayType& GetFictitiousIntegrationPoints() const
    {
        return mFictitiousIntegrationPoints;
    }

    /// Compute the domain size of the subcell i
    virtual double DomainSizeSubCell(const std::size_t& i, const BRep& r_brep, const int& integration_method) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class", __FUNCTION__)
    }

    /// Fit a list of subcell using the representative quadrature points from the parent element
    virtual Matrix FitQuadratureSubCell(const std::vector<size_t>& subcell_index,
        const std::vector<FunctionR3R1::Pointer>& r_funcs,
        const BRep& r_brep,
        const int& integrator_integration_method,
        const std::string& solver_type,
        const int& echo_level,
        const double& small_weight) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class", __FUNCTION__)
    }

protected:

    Element::Pointer mpElement;
    int mRepresentativeIntegrationOrder; // this is the integration order associated with mMomentFittingIntegrationPoints
    IntegrationPointsArrayType mMomentFittingIntegrationPoints; // this contains the quadrature points to fit the integration on sub-cell. It is typically the Gauß quadrature point on big cell.
    IntegrationPointsArrayType mRepresentativeIntegrationPoints; // this contains the quadrature points representing a sub-cell. In most of the case, it is not useful. However, when the sub-cell is completely inside the physical domain, it can be used as physical integration point.

    std::vector<std::size_t> mSubCellIndices;
    IntegrationPointsArrayType mPhysicalIntegrationPoints;
    IntegrationPointsArrayType mFictitiousIntegrationPoints;

private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }
};

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, BaseMomentFittedQuadTreeSubCell& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const BaseMomentFittedQuadTreeSubCell& rThis)
{
    return rOStream;
}


///////////////////////////////////////////////////////////////////////////////


/** A special type of quad tree w/ sub-cell. It uses quad-tree to compute the moment fitted quadrature of the element. For details see:
Hoang-Giang Bui, D. Schillinger, G. Meschke, Finite Cell Method for Plasticity using Moment-Fitted Quadrature Technique, in preparation.
*/
template<std::size_t TNsampling, int TFrameType>
class MomentFittedQuadTreeSubCell : public QuadTreeSubCell<TNsampling, TFrameType>, public BaseMomentFittedQuadTreeSubCell
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of  MomentFittedQuadTreeSubCell
    KRATOS_CLASS_POINTER_DEFINITION(MomentFittedQuadTreeSubCell);

    typedef QuadTreeSubCell<TNsampling, TFrameType> BaseType;

    typedef BaseMomentFittedQuadTreeSubCell RefType;

    typedef typename GeometricalObject::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    MomentFittedQuadTreeSubCell(Element::Pointer p_elem)
    : BaseType(p_elem), RefType(p_elem)
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

        if( RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4
         || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8
         || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9
         || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D4
         || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D8
         || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D9 )
        {
            // create the sub-cells
            BaseType::ConstructSubCellsForQuadBasedOnGaussQuadrature(BaseType::mpTreeNodes, gauss_quadrature_order);
        }
        #ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
        else if(RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier2D
             || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier2D3 )
        {
            KRATOS_THROW_ERROR(std::logic_error, "Not implemented", "")
        }
        #endif
        else if( RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D8
              || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D20
              || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D27 )
        {
            // create the sub-cells
            BaseType::ConstructSubCellsForHexBasedOnGaussQuadrature(BaseType::mpTreeNodes, gauss_quadrature_order);
        }
        #ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
        else if(RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier3D)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Not implemented", "")
        }
        #endif
        else
            KRATOS_THROW_ERROR(std::logic_error, "This geometry type is not supported:", RefType::mpElement->GetGeometry().GetGeometryType())

        // using the standard Gauss quadrature as moment fitting integration points
        GeometryData::IntegrationMethod RepresentativeIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(gauss_quadrature_order);

        RefType::mRepresentativeIntegrationOrder = gauss_quadrature_order;

        if(gauss_quadrature_type == 0)
        {
            const IntegrationPointsArrayType& integration_points
                = RefType::mpElement->GetGeometry().IntegrationPoints( RepresentativeIntegrationMethod );

            RefType::mMomentFittingIntegrationPoints = integration_points;
            RefType::mRepresentativeIntegrationPoints = integration_points;

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
            typename QuadTreeNode<TFrameType>::Pointer pTreeNode = QuadTree<TNsampling, TFrameType>::pCreateQuadTreeNode(RefType::mpElement->GetGeometry().GetGeometryType());

            const IntegrationPointsArrayType integration_points
                = pTreeNode->ConstructCustomQuadrature(gauss_quadrature_type, gauss_quadrature_order);

            RefType::mMomentFittingIntegrationPoints = integration_points;
            RefType::mRepresentativeIntegrationPoints = integration_points;

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

        if( RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4
         || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8
         || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9
         || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D4
         || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D8
         || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D9 )
        {
            std::size_t nsegments = integration_order;

            // create the sub-cells
            BaseType::ConstructSubCellsForQuadBasedOnEqualDistribution(BaseType::mpTreeNodes, nsegments, nsegments);
        }
        #ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
        else if( RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier2D
              || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier2D3 )
        {
            std::size_t nsegments = integration_order;

            // create the sub-cells
            BaseType::ConstructSubCellsForBezier2DBasedOnEqualDistribution(BaseType::mpTreeNodes, nsegments, nsegments);
        }
        #endif
        else if( RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D8
              || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D20
              || RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D27 )
        {
            std::size_t nsegments = integration_order;

            // create the sub-cells
            BaseType::ConstructSubCellsForHexBasedOnEqualDistribution(BaseType::mpTreeNodes, nsegments, nsegments, nsegments);
        }
        #ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
        else if( RefType::mpElement->GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier3D )
        {
            std::size_t nsegments = integration_order;

            // create the sub-cells
            BaseType::ConstructSubCellsForBezier3DBasedOnEqualDistribution(BaseType::mpTreeNodes, nsegments, nsegments, nsegments);
        }
        #endif
        else
            KRATOS_THROW_ERROR(std::logic_error, "This geometry type is not supported:", RefType::mpElement->GetGeometry().GetGeometryType())

        // using the standard Gauss quadrature as moment fitting integration points
        GeometryData::IntegrationMethod RepresentativeIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(integration_order);

        RefType::mRepresentativeIntegrationOrder = integration_order;

        const IntegrationPointsArrayType& integration_points
            = RefType::mpElement->GetGeometry().IntegrationPoints( RepresentativeIntegrationMethod );
//        KRATOS_WATCH(integration_points.size())

        RefType::mMomentFittingIntegrationPoints = integration_points;

        for(std::size_t i = 0; i < BaseType::mpTreeNodes.size(); ++i)
        {
            RefType::mRepresentativeIntegrationPoints.push_back(BaseType::mpTreeNodes[i]->ReferenceCenter());
        }
    }

    /// Destructor.
    virtual ~ MomentFittedQuadTreeSubCell() {}


    /// Compute the physical quadrature point. One physical quadrature point represents a subcell.
    virtual void GeneratePhysicalIntegrationPoints(const BRep& r_brep, const int& integrator_integration_method)
    {
        RefType::mSubCellIndices.clear();
        RefType::mPhysicalIntegrationPoints.clear();
        RefType::mFictitiousIntegrationPoints.clear();

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

            #ifdef ENABLE_PROFILING
            std::stringstream time_mark_name; time_mark_name << "CutStatus at " << __LINE__;
            Timer::Start(time_mark_name.str());
            #endif

            int status = r_brep.CutStatus(SamplingPoints);

            #ifdef ENABLE_PROFILING
            Timer::Stop(time_mark_name.str());
            #endif

            if(status == BRep::_IN)
            {
                #ifdef DEBUG_SUBCELL
                std::cout << "sub-cell " << i << " is completely inside the physical domain" << std::endl;
                #endif
                // the cell is completely inside the physical domain
                GeometryType::IntegrationPointType integration_point = RefType::mRepresentativeIntegrationPoints[i];
                integration_point.Weight() = 0.0; // cancel out the contribution of this integration point to parent element
                RefType::mPhysicalIntegrationPoints.push_back(integration_point);
                RefType::mSubCellIndices.push_back(i);
            }
            else if(status == BRep::_CUT)
            {
                #ifdef DEBUG_SUBCELL
                std::cout << "sub-cell " << i << " is cut" << std::endl;
                #endif
                #ifdef ENABLE_PROFILING
                std::stringstream time_mark_name2; time_mark_name2 << "CenterOfGravity at " << __LINE__;
                Timer::Start(time_mark_name2.str());
                #endif
                // the cell is cut
                bool found = BaseType::mpTreeNodes[i]->CenterOfGravity(COG, BaseType::pGetGeometry(), r_brep, integrator_integration_method);
                #ifdef ENABLE_PROFILING
                Timer::Stop(time_mark_name2.str());
                #endif
                #ifdef DEBUG_SUBCELL
                KRATOS_WATCH(found)
                KRATOS_WATCH(COG)
                #endif
                if(found)
                {
                    GeometryType::IntegrationPointType integration_point;
                    bool is_inside;
                    if(TFrameType == GLOBAL_CURRENT)
                    {
                        is_inside = BaseType::pGetGeometry()->IsInside(COG, integration_point);
                    }
                    else if(TFrameType == GLOBAL_REFERENCE)
                    {
                        GeometryType& rGeometry = *(BaseType::pGetGeometry());
                        Matrix DeltaPosition(rGeometry.size(), 3);
                        for ( unsigned int node = 0; node < rGeometry.size(); ++node )
                            noalias( row( DeltaPosition, node ) ) = rGeometry[node].Coordinates() - rGeometry[node].GetInitialPosition();
                        is_inside = BaseType::pGetGeometry()->IsInside(COG, integration_point, DeltaPosition);
                    }
                    if(!is_inside)
                    {
                        std::cout << "!!WARNING!!The CenterOfGravity is found but not inside the domain of the subcell " << i << " of the element " << mpElement->Id() << std::endl;
                        std::cout << " The found COG: " << COG << std::endl;
                        std::cout << " The found integration point: " << integration_point << std::endl;
                        std::cout << " The element " << RefType::mpElement->Id() << ":" << std::endl;
                        for (std::size_t i = 0; i < RefType::mpElement->GetGeometry().size(); ++i)
                            std::cout << " " << RefType::mpElement->GetGeometry()[i].Id();
                        std::cout << std::endl;
                        for (std::size_t i = 0; i < RefType::mpElement->GetGeometry().size(); ++i)
                        {
                            const PointType& point = RefType::mpElement->GetGeometry()[i].GetInitialPosition();
                            std::cout << " " << point[0] << "," << point[1] << "," << point[2] << std::endl;
                        }
                    }
                    integration_point.Weight() = 0.0; // cancel out the contribution of this integration point to parent element
                    RefType::mPhysicalIntegrationPoints.push_back(integration_point);
                    RefType::mSubCellIndices.push_back(i);
                }
                else
                {
                    // if the adaptive quadrature can't find the center of gravity, this sub-cell will be considered fictitious
                    GeometryType::IntegrationPointType integration_point = RefType::mRepresentativeIntegrationPoints[i];
                    integration_point.Weight() = 0.0; // the weight of fictitious integration point is initially always 0.0
                    RefType::mFictitiousIntegrationPoints.push_back(integration_point);
                }
            }
            else
            {
                #ifdef DEBUG_SUBCELL
                std::cout << "sub-cell " << i << " is completely outside the physical domain" << std::endl;
                #endif
                // the cell is completely outside the physical domain
                GeometryType::IntegrationPointType integration_point = RefType::mRepresentativeIntegrationPoints[i];
                integration_point.Weight() = 0.0; // the weight of fictitious integration point is initially always 0.0
                RefType::mFictitiousIntegrationPoints.push_back(integration_point);
            }

            #ifdef DEBUG_SUBCELL
            std::cout << "---------------------><--------------------\n\n";
            #endif
        }
    }


    /// Compute the domain size of the subcell i
    virtual double DomainSizeSubCell(const std::size_t& i, const BRep& r_brep, const int& integration_method) const
    {
        return BaseType::DomainSize(i, r_brep, integration_method);
    }


    /// Set the weight for fictitious integration points
    void SetFictitiousWeight(const double& weight)
    {
        for (std::size_t i = 0; i < RefType::mFictitiousIntegrationPoints.size(); ++i)
            RefType::mFictitiousIntegrationPoints[i].SetWeight(weight);
    }


    /// Fit a list of subcell using the representative quadrature points from the parent element
    virtual Matrix FitQuadratureSubCell(const std::vector<size_t>& subcell_index,
        const std::vector<FunctionR3R1::Pointer>& r_funcs,
        const BRep& r_brep,
        const int& integrator_integration_method,
        const std::string& solver_type,
        const int& echo_level,
        const double& small_weight) const
    {
        Matrix Weights(subcell_index.size(), RefType::mMomentFittingIntegrationPoints.size());
        Vector aux;

        for(std::size_t i = 0; i < subcell_index.size(); ++i)
        {
            // perform the moment fit for the sub-cell
            const std::size_t& i_cell = subcell_index[i];

            typename QuadTree<TNsampling, TFrameType>::Pointer quad_tree = typename QuadTree<TNsampling, TFrameType>::Pointer(
                    new QuadTree<TNsampling, TFrameType>(BaseType::pGetGeometry(), BaseType::mpTreeNodes[i_cell]));

            #ifdef ENABLE_PROFILING
            std::stringstream time_mark_name; time_mark_name << "FitQuadrature at " << __LINE__;
            Timer::Start(time_mark_name.str());
            #endif

            aux = MomentFittingUtility::FitQuadrature<FunctionR3R1, QuadTree<TNsampling, TFrameType> >(*BaseType::pGetGeometry(),
                    r_funcs, r_brep, *quad_tree, RefType::mMomentFittingIntegrationPoints,
                    integrator_integration_method, solver_type, echo_level, small_weight);

            if (Weights.size2() != aux.size())
                KRATOS_THROW_ERROR(std::logic_error, "Size conflict", "")

            noalias(row(Weights, i)) = aux;

            #ifdef ENABLE_PROFILING
            Timer::Stop(time_mark_name.str());
            #endif
        }

        return Weights;
    }


    /// Create the element out from sub-cells. The element takes the same geometry of the original element, but the weights are given.
    ModelPart::ElementsContainerType CreateSubCellElements(ModelPart& r_model_part,
        const std::string& sample_element_name,
        const Matrix& rWeights,
        std::size_t& lastElementId,
        std::size_t& lastCondId) const
    {
        return CreateSubCellElements(r_model_part, RefType::mpElement, sample_element_name, RefType::mRepresentativeIntegrationOrder,
                            RefType::mMomentFittingIntegrationPoints, rWeights, lastElementId, lastCondId, RefType::mpElement->pGetProperties());
    }


    /// Create the elements out from sub-cells of an element. The new elements take the same geometry of the original element, but the integration points and weights are given.
    static ModelPart::ElementsContainerType CreateSubCellElements(ModelPart& r_model_part,
        Element::Pointer pElement, // the parent element that contains sub-cells
        const std::string& sample_element_name,
        const int& RepresentativeIntegrationOrder,
        const IntegrationPointsArrayType& integration_points,
        const Matrix& rWeights,
        std::size_t& lastElementId,
        std::size_t& lastCondId,
        Properties::Pointer pProperties)
    {
        /* create the new elements from subcell */
        std::size_t num_physical_point = rWeights.size1();

        Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);
        ModelPart::ElementsContainerType NewElements;
        GeometryType& r_geom = *(pElement->pGetGeometry());

        GeometryData::IntegrationMethod RepresentativeIntegrationMethod
                = Function<double, double>::GetIntegrationMethod(RepresentativeIntegrationOrder);
        Variable<int>& INTEGRATION_ORDER_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("INTEGRATION_ORDER"));

        for(std::size_t i = 0; i < num_physical_point; ++i)
        {
            // create new list of integration points
            IntegrationPointsArrayType new_integration_points = integration_points;
            for(std::size_t j = 0; j < new_integration_points.size(); ++j)
                new_integration_points[j].Weight() = rWeights(i, j);

            // create the new "extrapolated" elements from sub-cell
            // here we make a clone of the geometry because we want to assign different geometry data later on
            // this also works with Bezier element, because Bezier geometry has implemented the Create method
            Element::Pointer pNewElement;
            pNewElement = r_clone_element.Create(++lastElementId, r_geom.Create(r_geom.Points()), pProperties);

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
            Condition::Pointer pNewCond = Condition::Pointer(new ElementWrapperCondition(++lastCondId, *it, pProperties));
            r_model_part.Conditions().push_back(pNewCond);
        }

        std::cout << NewElements.size() << " new " << sample_element_name << "-wrapped conditions from parent element " << pElement->Id()
                  << " are added to the model_part" << std::endl;

        return NewElements;
    }


    /// Create the fictitious element out from an element. The new element takes the same geometry of the original element, but with the fictitious integration points.
    Element::Pointer CreateFictitiousElement(ModelPart& r_model_part,
        const std::string& sample_element_name,
        std::size_t& lastElementId,
        Properties::Pointer pProperties)
    {
        return FiniteCellMeshUtility::CreateParasiteElement(r_model_part, RefType::mpElement, sample_element_name,
                RefType::mRepresentativeIntegrationOrder, RefType::mFictitiousIntegrationPoints, lastElementId, pProperties);
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

    /// Assignment operator.
     MomentFittedQuadTreeSubCell& operator=( MomentFittedQuadTreeSubCell const& rOther);

    /// Copy constructor.
     MomentFittedQuadTreeSubCell( MomentFittedQuadTreeSubCell const& rOther);

}; // Class  MomentFittedQuadTreeSubCell


/// input stream function
template<std::size_t TNsampling, int TFrameType>
inline std::istream& operator >> (std::istream& rIStream, MomentFittedQuadTreeSubCell<TNsampling, TFrameType>& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TNsampling, int TFrameType>
inline std::ostream& operator << (std::ostream& rOStream, const  MomentFittedQuadTreeSubCell<TNsampling, TFrameType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#undef ENABLE_PROFILING
#undef ENABLE_DEBUG_QUADRATURE

#endif // KRATOS_FINITE_CELL_APPLICATION_MOMENT_FITTED_QUAD_TREE_SUBCELL_H_INCLUDED defined
