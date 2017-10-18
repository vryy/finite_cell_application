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
    (QuadTreeName.str().c_str(), init<Element::Pointer&>())
    .def(init<Condition::Pointer&>())
    .def("pGetGeometry", &QuadTreeType::pGetGeometry)
    .def("DomainSize", &QuadTreeType::DomainSize)
    .def("CenterOfGravity", &QuadTreeType::CenterOfGravity)
    .def("AddToModelPart", &QuadTreeType::PyAddToModelPart)
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

    void(MomentFittedQuadTreeSubCellType::*pointer_to_ConstructSubCellsBasedOnEqualDistribution1)(const int&, const std::size_t&, const std::size_t&) = &MomentFittedQuadTreeSubCellType::ConstructSubCellsBasedOnEqualDistribution;
    void(MomentFittedQuadTreeSubCellType::*pointer_to_ConstructSubCellsBasedOnEqualDistribution2)(const int&, const std::size_t&, const std::size_t&, const std::size_t&) = &MomentFittedQuadTreeSubCellType::ConstructSubCellsBasedOnEqualDistribution;

    class_<MomentFittedQuadTreeSubCellType, typename MomentFittedQuadTreeSubCellType::Pointer, bases<QuadTreeSubCellType>, boost::noncopyable>
    (MomentFittedQuadTreeSubCellName.str().c_str(), init<Element::Pointer>())
    .def("ConstructSubCellsBasedOnGaussQuadrature", &MomentFittedQuadTreeSubCellType::ConstructSubCellsBasedOnGaussQuadrature)
    .def("ConstructSubCellsBasedOnEqualDistribution", pointer_to_ConstructSubCellsBasedOnEqualDistribution1)
    .def("ConstructSubCellsBasedOnEqualDistribution", pointer_to_ConstructSubCellsBasedOnEqualDistribution2)
    .def("GetElement", &MomentFittedQuadTreeSubCellType::pGetElement)
    .def("GetNumberOfPhysicalIntegrationPoint", &MomentFittedQuadTreeSubCellType::GetNumberOfPhysicalIntegrationPoint)
    .def("FitAndCreateSubCellElements", &MomentFittedQuadTreeSubCellType::PyFitAndCreateSubCellElements)
    .def("CreateSubCellElements", &MomentFittedQuadTreeSubCellType::PyCreateSubCellElements)
    .def("CreateParasiteElement", &MomentFittedQuadTreeSubCellType::CreateParasiteElement)
    ;
}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_APPLICATION_ADD_QUADTREE_TO_PYTHON_HPP_INCLUDED  defined 
