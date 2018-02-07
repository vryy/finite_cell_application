// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_FINITE_CELL_APPLICATION_ADD_QUADTREE_TO_PYTHON_HPP_INCLUDED )
#define  KRATOS_FINITE_CELL_APPLICATION_ADD_QUADTREE_TO_PYTHON_HPP_INCLUDED


// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/moment_fitted_quad_tree_subcell.h"
#include "custom_utilities/finite_cell_geometry_utility.h"
#include "custom_utilities/finite_cell_auxilliary_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void FiniteCellApplication_AddRefinableTreeToPython()
{
    class_<RefinableTree, RefinableTree::Pointer, boost::noncopyable>
    ("RefinableTree", init<>())
    .def("Refine", &RefinableTree::Refine)
    .def("RefineBy", &RefinableTree::RefineBy)
    ;
}

/// construct the element out from quad-tree and add to model_part.
/// This is mainly for post-processing
template<class TTreeType>
boost::python::list QuadTree_AddToModelPart(TTreeType& rDummy,
    ModelPart& r_model_part, const std::string sample_entity_name,
    std::size_t lastNodeId, std::size_t lastEntityId)
{
    if( KratosComponents<Element>::Has(sample_entity_name) )
    {
        Element const& r_clone_element = KratosComponents<Element>::Get(sample_entity_name);
        rDummy.Get().template AddToModelPart<true, Element>(rDummy.pGetGeometry(), r_model_part, r_clone_element, lastNodeId, lastEntityId, 1);
    }
    else if( KratosComponents<Condition>::Has(sample_entity_name) )
    {
        Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_entity_name);
        rDummy.Get().template AddToModelPart<true, Condition>(rDummy.pGetGeometry(), r_model_part, r_clone_condition, lastNodeId, lastEntityId, 1);
    }

    boost::python::list list;
    list.append(lastNodeId);
    list.append(lastEntityId);
    return list;
}

/// Create the element out from sub-cells. The element takes the same geometry of the original element, but the weight is given.
/// This is the python interface
template<class TTreeType>
boost::python::list MomentFittedQuadTreeSubCell_CreateSubCellElements(TTreeType& rDummy,
        ModelPart& r_model_part,
        const std::string& subcell_element_type,
        const int& cut_cell_quadrature_order,
        boost::python::list& cut_cell_full_quadrature,
        boost::python::list& subcell_weights,
        std::size_t lastElementId,
        std::size_t lastCondId)
{
    typedef Element::GeometryType GeometryType;

    ModelPart::ElementsContainerType NewElements;

//    if(!boost::python::list::is_empty(subcell_weights))
//    {
    // TODO find a way to check if a list is empty
        GeometryType::IntegrationPointsArrayType integration_points;
        Matrix Weights;

        std::size_t num_physical_point = boost::python::len(subcell_weights);
        std::size_t weight_length = boost::python::len(subcell_weights[0]);
        Weights.resize(num_physical_point, weight_length, false);

        for(std::size_t i = 0; i < num_physical_point; ++i)
        {
            boost::python::list weights = boost::python::extract<boost::python::list>(subcell_weights[i]);
            for(std::size_t j = 0; j < weight_length; ++j)
            {
                Weights(i, j) = boost::python::extract<double>(weights[j]);
            }
        }
//        KRATOS_WATCH(Weights)

        for(std::size_t i = 0; i < boost::python::len(cut_cell_full_quadrature); ++i)
        {
            boost::python::list point = boost::python::extract<boost::python::list>(cut_cell_full_quadrature[i]);
            GeometryType::IntegrationPointType integration_point;
            integration_point.X() = boost::python::extract<double>(point[0]);
            integration_point.Y() = boost::python::extract<double>(point[1]);
            integration_point.Z() = boost::python::extract<double>(point[2]);
            integration_point.Weight() = boost::python::extract<double>(point[3]);
//            KRATOS_WATCH(integration_point)
            integration_points.push_back(integration_point);
        }

        NewElements = rDummy.CreateSubCellElements(r_model_part,
            subcell_element_type,
            cut_cell_quadrature_order,
            integration_points,
            Weights,
            lastElementId,
            lastCondId);
//        std::cout << "----------------------" << std::endl;
//    }

    boost::python::list Output;
    Output.append(lastElementId);
    Output.append(lastCondId);
    Output.append(NewElements);

    return Output;
}

/// Create the element out from sub-cells. The element takes the same geometry of the original element, but the weight is fitted by sub-cell.
template<class TTreeType>
ModelPart::ElementsContainerType MomentFittedQuadTreeSubCell_FitAndCreateSubCellElements(TTreeType& rDummy,
    ModelPart& r_model_part,
    const std::string& sample_element_name,
    boost::python::list& r_funcs,
    const BRep& r_brep,
    const int& integrator_integration_method,
    const std::string& solver_type,
    const int& echo_level,
    const double& small_weight,
    const bool& compute_subcell_domain_size)
{
    /* firstly compute the physical integration point */
    typedef Element::GeometryType GeometryType;

    std::pair<std::vector<std::size_t>, GeometryType::IntegrationPointsArrayType> Output
            = rDummy.GetPhysicalInterationPoint(r_brep, integrator_integration_method);
    const std::vector<std::size_t>& subcell_index = Output.first;
    const GeometryType::IntegrationPointsArrayType& physical_integration_points = Output.second;

    /* secondly assign the quadrature for parent element based on physical integration_points of the previous step */
    GeometryData::IntegrationMethod RepresentativeIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(rDummy.GetRepresentativeIntegrationOrder());
//        KRATOS_WATCH(RepresentativeIntegrationMethod)
    FiniteCellGeometryUtility::AssignGeometryData(*(rDummy.pGetGeometry()), RepresentativeIntegrationMethod, physical_integration_points);
    Variable<int>& INTEGRATION_ORDER_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("INTEGRATION_ORDER"));
    rDummy.pGetElement()->SetValue(INTEGRATION_ORDER_var, rDummy.GetRepresentativeIntegrationOrder());
    rDummy.pGetElement()->Initialize();

    /* thirdly fit the sub-cell */
    std::vector<FunctionR3R1::Pointer> funcs;
    typedef boost::python::stl_input_iterator<FunctionR3R1::Pointer> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& f,
                  std::make_pair(iterator_value_type(r_funcs), // begin
                    iterator_value_type() ) ) // end
    {
        funcs.push_back(f);
    }

    Matrix Weights = rDummy.FitQuadratureSubCell(subcell_index, funcs, r_brep, integrator_integration_method, solver_type, echo_level, small_weight);

    /* next create the sub-elements */
    // find the last element id
    std::size_t lastElementId = FiniteCellAuxilliaryUtility::GetLastElementId(r_model_part);

    // find the last condition id
    std::size_t lastCondId = FiniteCellAuxilliaryUtility::GetLastConditionId(r_model_part);

    ModelPart::ElementsContainerType NewSubCellElements = rDummy.CreateSubCellElements(r_model_part, sample_element_name, Weights, lastElementId, lastCondId);

    // set the domain size for the sub-cell
    if(compute_subcell_domain_size)
    {
        std::size_t cnt = 0;
        for(typename ModelPart::ElementsContainerType::ptr_iterator it = NewSubCellElements.ptr_begin();
                it != NewSubCellElements.ptr_end(); ++it)
        {
            (*it)->SetValue(SUBCELL_DOMAIN_SIZE, rDummy.DomainSize(subcell_index[cnt++], r_brep, integrator_integration_method));
        }
    }

    return NewSubCellElements;
}

template<std::size_t TNsampling>
void FiniteCellApplication_AddQuadTreeToPython()
{
    typedef QuadTree<TNsampling> QuadTreeType;
    typedef QuadTreeSubCell<TNsampling> QuadTreeSubCellType;
    typedef MomentFittedQuadTreeSubCell<TNsampling> MomentFittedQuadTreeSubCellType;

    std::stringstream QuadTreeName;
    QuadTreeName << "QuadTreeLocal";
    if(TNsampling > 1)
        QuadTreeName << TNsampling;

    double(QuadTreeType::*pointer_to_IntegrateLocal_double_quadtree_local)(const FunctionR3R1&, const int&) const = &QuadTreeType::template Integrate<double, QuadTreeNode::LOCAL>;
    double(QuadTreeType::*pointer_to_IntegrateGlobal_double_quadtree_local)(const FunctionR3R1&, const int&) const = &QuadTreeType::template Integrate<double, QuadTreeNode::GLOBAL>;

    class_<QuadTreeType, typename QuadTreeType::Pointer, boost::noncopyable, bases<RefinableTree> >
    (QuadTreeName.str().c_str(), init<Element::Pointer>())
    .def(init<Condition::Pointer>())
    .def("pGetGeometry", &QuadTreeType::pGetGeometry)
    .def("DomainSize", &QuadTreeType::DomainSize)
    .def("CenterOfGravity", &QuadTreeType::CenterOfGravity)
    .def("AddToModelPart", &QuadTree_AddToModelPart<QuadTreeType>)
    .def("IntegrateLocal", pointer_to_IntegrateLocal_double_quadtree_local)
    .def("IntegrateGlobal", pointer_to_IntegrateGlobal_double_quadtree_local)
    .def("ConstructQuadrature", &QuadTreeType::ConstructQuadrature)
    .def(self_ns::str(self))
    ;

    double(QuadTreeSubCellType::*pointer_to_DomainSize1)(const std::size_t&, const BRep&, const int&) const = &QuadTreeSubCellType::DomainSize;
    double(QuadTreeSubCellType::*pointer_to_DomainSize2)(const BRep&, const int&) const = &QuadTreeSubCellType::DomainSize;

    std::stringstream QuadTreeSubCellName;
    QuadTreeSubCellName << "QuadTreeSubCell";
    if(TNsampling > 1)
        QuadTreeSubCellName << TNsampling;

    class_<QuadTreeSubCellType, typename QuadTreeSubCellType::Pointer, boost::noncopyable, bases<RefinableTree> >
    (QuadTreeSubCellName.str().c_str(), init<Element::Pointer&>())
    .def(init<Condition::Pointer&>())
    .def("NumberOfSubCells", &QuadTreeSubCellType::NumberOfSubCells)
    .def("DomainSize", pointer_to_DomainSize1)
    .def("DomainSize", pointer_to_DomainSize2)
    .def("CreateQuadTree", &QuadTreeSubCellType::CreateQuadTree)
    .def("ConstructQuadrature", &QuadTreeSubCellType::ConstructQuadrature)
    .def("ShallowAddToModelPart", &QuadTreeSubCellType::template PyAddToModelPart<true>) // only add the sub-cell
    .def("DeepAddToModelPart", &QuadTreeSubCellType::template PyAddToModelPart<false>) // add the sub-cell and all the quad-trees
    ;

    std::stringstream MomentFittedQuadTreeSubCellName;
    MomentFittedQuadTreeSubCellName << "MomentFittedQuadTreeSubCell";
    if(TNsampling > 1)
        MomentFittedQuadTreeSubCellName << TNsampling;

    void(MomentFittedQuadTreeSubCellType::*pointer_to_ConstructSubCellsBasedOnEqualDistribution1)(const int&) = &MomentFittedQuadTreeSubCellType::ConstructSubCellsBasedOnEqualDistribution;
    void(MomentFittedQuadTreeSubCellType::*pointer_to_ConstructSubCellsBasedOnEqualDistribution2)(const int&) = &MomentFittedQuadTreeSubCellType::ConstructSubCellsBasedOnEqualDistribution;

    class_<MomentFittedQuadTreeSubCellType, typename MomentFittedQuadTreeSubCellType::Pointer, bases<QuadTreeSubCellType>, boost::noncopyable>
    (MomentFittedQuadTreeSubCellName.str().c_str(), init<Element::Pointer>())
    .def("ConstructSubCellsBasedOnGaussQuadrature", &MomentFittedQuadTreeSubCellType::ConstructSubCellsBasedOnGaussQuadrature)
    .def("ConstructSubCellsBasedOnEqualDistribution", pointer_to_ConstructSubCellsBasedOnEqualDistribution1)
    .def("ConstructSubCellsBasedOnEqualDistribution", pointer_to_ConstructSubCellsBasedOnEqualDistribution2)
    .def("GetElement", &MomentFittedQuadTreeSubCellType::pGetElement)
    .def("GetNumberOfPhysicalIntegrationPoint", &MomentFittedQuadTreeSubCellType::GetNumberOfPhysicalIntegrationPoint)
    .def("FitAndCreateSubCellElements", &MomentFittedQuadTreeSubCell_FitAndCreateSubCellElements<MomentFittedQuadTreeSubCellType>)
    .def("CreateSubCellElements", &MomentFittedQuadTreeSubCell_CreateSubCellElements<MomentFittedQuadTreeSubCellType>)
    .def("CreateParasiteElement", &MomentFittedQuadTreeSubCellType::CreateParasiteElement)
    ;
}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_APPLICATION_ADD_QUADTREE_TO_PYTHON_HPP_INCLUDED  defined
