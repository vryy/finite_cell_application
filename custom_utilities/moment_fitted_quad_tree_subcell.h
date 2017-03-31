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

//#define ENABLE_DEBUG_QUADRATURE
//#define DEBUG_SUBCELL

namespace Kratos
{

/// Short class definition.
/** A special type of quad tree w/ sub-cell. It uses quad-tree to compute the moment fitted quadrature of the element. For details see:
Hoang-Giang Bui, D. Schillinger, G. Meschke, Finite Cell Method for Plasticity using Moment-Fitted Quadrature Technique, in preparation.
*/
class MomentFittedQuadTreeSubCell : public QuadTreeSubCell
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of  MomentFittedQuadTreeSubCell
    KRATOS_CLASS_POINTER_DEFINITION(MomentFittedQuadTreeSubCell);

    typedef QuadTreeSubCell BaseType;

    typedef typename GeometricalObject::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor to construct the sub-cells based on Gauss quadrature.
    MomentFittedQuadTreeSubCell(Element::Pointer& p_elem, const std::string& construct_method,
            const int& integration_order)
    : BaseType(p_elem), mpElement(p_elem), mRepresentativeIntegrationOrder(integration_order)
    {
        if(construct_method != "gauss")
            KRATOS_THROW_ERROR(std::logic_error, "This construct_method is not valid for this constructor:", construct_method)

        if( p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D4
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D8
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D9 )
        {
            // create the sub-cells
            BaseType::ConstructSubCellsForQuadBasedOnGaussQuadrature(BaseType::mpTreeNodes, integration_order);

            // using the standard Gauss quadrature as representative integration points
            GeometryData::IntegrationMethod RepresentativeIntegrationMethod
                = Function<double, double>::GetIntegrationMethod(integration_order);

            const GeometryType::IntegrationPointsArrayType& integration_points
                = p_elem->GetGeometry().IntegrationPoints( RepresentativeIntegrationMethod );

            mSubCellRepresentativeIntegrationPoints.clear();
            for(std::size_t i = 0; i < BaseType::mpTreeNodes.size(); ++i)
            {
                mSubCellRepresentativeIntegrationPoints.push_back(integration_points[i]);
            }

            // do a check, to make sure the integration point is inside the sub-cell
            for(std::size_t i = 0; i < BaseType::mpTreeNodes.size(); ++i)
            {
                if(!mpTreeNodes[i]->IsInside(integration_points[i]))
                    KRATOS_THROW_ERROR(std::logic_error, "The integration_point is not inside the tree, error at", i)
            }
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "This geometry type is not supported:", p_elem->GetGeometry().GetGeometryType())

        #ifdef ENABLE_DEBUG_QUADRATURE
        if(p_elem->Id() == 930)
        {
            KRATOS_WATCH(p_elem->Id())
            mDebugFlag = true;
        }
        #endif
    }

    /// Constructor to construct the sub-cells based on equal distribution.
    MomentFittedQuadTreeSubCell(Element::Pointer& p_elem, const std::string& construct_method,
            const std::size_t& m, const std::size_t& n)
    : BaseType(p_elem), mpElement(p_elem)
    {
        if(construct_method != "equal-dist")
            KRATOS_THROW_ERROR(std::logic_error, "This construct_method is not valid for this constructor:", construct_method)

        if( p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D4
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D8
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D9 )
        {
            // create the sub-cells
            BaseType::ConstructSubCellsForQuadBasedOnEqualDistribution(BaseType::mpTreeNodes, m, n);

            // using the center of the sub-cells as the quadrature point
            mRepresentativeIntegrationOrder = 1;

            mSubCellRepresentativeIntegrationPoints.clear();
            for(std::size_t i = 0; i < BaseType::mpTreeNodes.size(); ++i)
            {
                mSubCellRepresentativeIntegrationPoints.push_back(BaseType::mpTreeNodes[i]->ReferenceCenter());
            }
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "This geometry type is not supported:", p_elem->GetGeometry().GetGeometryType())
    }

    /// Destructor.
    virtual ~ MomentFittedQuadTreeSubCell() {}


    void PyFitQuadraturePhysicalPoints(boost::python::list& r_funcs,
            const BRep& r_brep,
            const int& integrator_integration_order,
            const std::string& solver_type,
            const int& echo_level) const
    {
        std::vector<FunctionR3R1::Pointer> funcs;
        typedef boost::python::stl_input_iterator<FunctionR3R1::Pointer> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& f,
                      std::make_pair(iterator_value_type(r_funcs), // begin
                        iterator_value_type() ) ) // end
        {
            funcs.push_back(f);
        }

        GeometryData::IntegrationMethod IntegratorIntegrationMethod
                = Function<double, double>::GetIntegrationMethod(integrator_integration_order);

        FitQuadraturePhysicalPoints<FunctionR3R1>(funcs, r_brep, IntegratorIntegrationMethod, solver_type, echo_level);
    }


    template<class TFunctionType>
    void FitQuadraturePhysicalPoints(const std::vector<typename TFunctionType::Pointer>& r_funcs,
            const BRep& r_brep,
            const GeometryData::IntegrationMethod& IntegratorIntegrationMethod,
            const std::string& solver_type,
            const int& echo_level) const
    {
        GeometryType& r_geom = *pGetGeometry();

        // extract the physical integration points
        GeometryType::IntegrationPointsArrayType physical_integration_points;
        PointType COG;
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
        {
            GeometryType::Pointer pSubCellGeometry = mpTreeNodes[i]->pCreateGeometry(pGetGeometry());

            int status = r_brep.CutStatus(*pSubCellGeometry);
            if(status == 0)
            {
                physical_integration_points.push_back(mSubCellRepresentativeIntegrationPoints[i]);
                #ifdef ENABLE_DEBUG_QUADRATURE
                if(mDebugFlag)
                    KRATOS_WATCH(mSubCellRepresentativeIntegrationPoints[i])
                #endif
            }
            else if(status == -1)
            {
                // compute the center of gravity of the sub-cell bounded by boundary
                bool found = mpTreeNodes[i]->CenterOfGravity(COG, pSubCellGeometry, r_brep);
                if(found)
                {
                    #ifdef ENABLE_DEBUG_QUADRATURE
                    if(mDebugFlag)
                    {
                        KRATOS_WATCH(*pSubCellGeometry)
                        KRATOS_WATCH(pSubCellGeometry->GetDefaultIntegrationMethod())
                        KRATOS_WATCH(COG)
                    }
                    #endif
                    GeometryType::IntegrationPointType integration_point;
                    r_geom.PointLocalCoordinates(integration_point, COG);
                    #ifdef ENABLE_DEBUG_QUADRATURE
                    if(mDebugFlag)
                        KRATOS_WATCH(integration_point)
                    #endif
                    physical_integration_points.push_back(integration_point);
                }
            }
        }

        #ifdef ENABLE_DEBUG_QUADRATURE
        if(mDebugFlag)
        {
            KRATOS_WATCH(r_funcs.size())
            KRATOS_WATCH(physical_integration_points.size())
        }
        #endif

        // fit the quadrature on the physical points
        Vector Weight = MomentFittingUtility::FitQuadrature<TFunctionType, QuadTreeSubCell>(r_geom,
            r_funcs,
            r_brep,
            *this,
            physical_integration_points,
            IntegratorIntegrationMethod,
            solver_type,
            echo_level);

        for(std::size_t point = 0; point < physical_integration_points.size(); ++point)
            physical_integration_points[point].Weight() = Weight[point];

        /* create new quadrature and assign to the geometry */
        GeometryData::IntegrationMethod RepresentativeIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(mRepresentativeIntegrationOrder);
        FiniteCellGeometry<GeometryType>::AssignGeometryData(r_geom, RepresentativeIntegrationMethod, physical_integration_points);
    }


    ModelPart::ElementsContainerType PyCreateSubCellElements(ModelPart& r_model_part,
        const std::string& sample_element_name,
        boost::python::list& r_funcs,
        const BRep& r_brep,
        const int& integrator_integration_order,
        const std::string& solver_type,
        const int& echo_level) const
    {
        std::vector<FunctionR3R1::Pointer> funcs;
        typedef boost::python::stl_input_iterator<FunctionR3R1::Pointer> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& f,
                      std::make_pair(iterator_value_type(r_funcs), // begin
                        iterator_value_type() ) ) // end
        {
            funcs.push_back(f);
        }

        GeometryData::IntegrationMethod IntegratorIntegrationMethod
                = Function<double, double>::GetIntegrationMethod(integrator_integration_order);

        return CreateSubCellElements<FunctionR3R1>(r_model_part, sample_element_name, funcs, r_brep, IntegratorIntegrationMethod, solver_type, echo_level);
    }


    /// Create the element out from sub-cells. The element takes the same geometry of the original element, but the weight is fitted by sub-cell
    template<class TFunctionType>
    ModelPart::ElementsContainerType CreateSubCellElements(ModelPart& r_model_part,
        const std::string& sample_element_name,
        const std::vector<typename TFunctionType::Pointer>& r_funcs,
        const BRep& r_brep,
        const GeometryData::IntegrationMethod& IntegratorIntegrationMethod,
        const std::string& solver_type,
        const int& echo_level) const
    {
        GeometryType& r_geom = *pGetGeometry();

        // loop through each sub-cell, find out if the cell is cut ot completely inside the physical domain
        GeometryType::IntegrationPointsArrayType physical_integration_points;
        PointType COG;
        std::map<std::size_t, std::size_t> map_point_cell;
        std::size_t num_physical_point = 0;
        for(std::size_t i = 0; i < mpTreeNodes.size(); ++i)
        {
            GeometryType::Pointer pSubCellGeometry = mpTreeNodes[i]->pCreateGeometry(pGetGeometry());

            int status = r_brep.CutStatus(*pSubCellGeometry);
            #ifdef DEBUG_SUBCELL
            KRATOS_WATCH(*pSubCellGeometry)
            KRATOS_WATCH(*mpTreeNodes[i])
            #endif
            if(status == 0)
            {
                #ifdef DEBUG_SUBCELL
                std::cout << "sub-cell " << i << " is completely inside the physical domain" << std::endl;
                #endif
                // the cell is completely inside the physical domain
                GeometryType::IntegrationPointType integration_point = mSubCellRepresentativeIntegrationPoints[i];
                integration_point.Weight() = 0.0; // cancel out the contribution of this integration point to parent element
                physical_integration_points.push_back(integration_point);
                map_point_cell[num_physical_point++] = i;
            }
            else if(status == -1)
            {
                #ifdef DEBUG_SUBCELL
                std::cout << "sub-cell " << i << " is cut" << std::endl;
                #endif
                // the cell is cut
                bool found = mpTreeNodes[i]->CenterOfGravity(COG, pGetGeometry(), r_brep);
                #ifdef DEBUG_SUBCELL
                KRATOS_WATCH(found)
                KRATOS_WATCH(COG)
                #endif
                if(found)
                {
                    GeometryType::IntegrationPointType integration_point;
                    r_geom.PointLocalCoordinates(integration_point, COG);
                    integration_point.Weight() = 0.0; // cancel out the contribution of this integration point to parent element
                    physical_integration_points.push_back(integration_point);
                    map_point_cell[num_physical_point++] = i;
                }
            }
            #ifdef DEBUG_SUBCELL
            else
                std::cout << "sub-cell " << i << " is completely outside the physical domain" << std::endl;
            std::cout << "---------------------><--------------------\n\n";
            #endif
        }

        // parent element is in charged of physical integration point
        GeometryData::IntegrationMethod RepresentativeIntegrationMethod = Function<double, double>::GetIntegrationMethod(mRepresentativeIntegrationOrder);
//        KRATOS_WATCH(RepresentativeIntegrationMethod)
        FiniteCellGeometry<GeometryType>::AssignGeometryData(r_geom, RepresentativeIntegrationMethod, physical_integration_points);
        Variable<int>& INTEGRATION_ORDER_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("INTEGRATION_ORDER"));
        mpElement->SetValue(INTEGRATION_ORDER_var, mRepresentativeIntegrationOrder);
        mpElement->Initialize();

        // find the last element id
        std::size_t lastElementId = FiniteCellAuxilliaryUtility::GetLastElementId(r_model_part);

        // find the last condition id
        std::size_t lastCondId = FiniteCellAuxilliaryUtility::GetLastConditionId(r_model_part);

        Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);
        ModelPart::ElementsContainerType NewElements;
        for(std::size_t i = 0; i < num_physical_point; ++i)
        {
            // perform the moment fit for the sub-cell
            std::size_t i_cell = map_point_cell[i];

            QuadTree::Pointer quad_tree = QuadTree::Pointer(new QuadTree(pGetGeometry(), BaseType::mpTreeNodes[i_cell]));

            Vector Weight = MomentFittingUtility::FitQuadrature<TFunctionType, QuadTree>(r_geom,
                r_funcs,
                r_brep,
                *quad_tree,
                mSubCellRepresentativeIntegrationPoints,
                IntegratorIntegrationMethod,
                solver_type,
                echo_level);

            // create new list of integration points
            GeometryType::IntegrationPointsArrayType new_integration_points = mSubCellRepresentativeIntegrationPoints;
            for(std::size_t j = 0; j < new_integration_points.size(); ++j)
                new_integration_points[j].Weight() = Weight[j];

            // create the new "parasite" elements from sub-cell
            Element::Pointer pNewElement = r_clone_element.Create(++lastElementId, r_geom.Points(), mpElement->pGetProperties());
            FiniteCellGeometry<GeometryType>::AssignGeometryData(pNewElement->GetGeometry(), RepresentativeIntegrationMethod, new_integration_points);
            pNewElement->SetValue(INTEGRATION_ORDER_var, mRepresentativeIntegrationOrder);
            pNewElement->Initialize();
            NewElements.push_back(pNewElement);
        }

        /* add new conditions to the model_part */
        for(typename ModelPart::ElementsContainerType::ptr_iterator it = NewElements.ptr_begin();
                it != NewElements.ptr_end(); ++it)
        {
            // create new element-wrapped condition
            Condition::Pointer pNewCond = Condition::Pointer(new ElementWrapperCondition(++lastCondId, *it));
            r_model_part.Conditions().push_back(pNewCond);
        }

        std::cout << NewElements.size() << " new element-wrapped conditions are added to the model_part" << std::endl;

        return NewElements;
    }


    /// Create an element taking the same nodes as the original one, but taking the type of geometry from sample_element_name
    Element::Pointer CreateParasiteElement(ModelPart& r_model_part,
        const std::string& sample_element_name, std::size_t lastElementId,
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
    int mRepresentativeIntegrationOrder;
    GeometryType::IntegrationPointsArrayType mSubCellRepresentativeIntegrationPoints;

    /// Assignment operator.
     MomentFittedQuadTreeSubCell& operator=( MomentFittedQuadTreeSubCell const& rOther);

    /// Copy constructor.
     MomentFittedQuadTreeSubCell( MomentFittedQuadTreeSubCell const& rOther);

}; // Class  MomentFittedQuadTreeSubCell


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, MomentFittedQuadTreeSubCell& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TDegree>
inline std::ostream& operator << (std::ostream& rOStream, const  MomentFittedQuadTreeSubCell& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#undef ENABLE_DEBUG_QUADRATURE

#endif // KRATOS_FINITE_CELL_APPLICATION_MOMENT_FITTED_QUAD_TREE_GARDEN_H_INCLUDED defined
