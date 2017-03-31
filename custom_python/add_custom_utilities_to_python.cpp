// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//



// Project includes
#include "includes/element.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_algebra/brep.h"
#include "custom_utilities/quadrature_utility.h"
#include "custom_utilities/binary_tree.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/div_free_basis_utility.h"
#include "custom_utilities/moment_fitting_utility.h"
#include "custom_utilities/moment_fitted_quad_tree_subcell.h"
#include "custom_utilities/finite_cell_auxilliary_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

Matrix ComputeDivFreeBasis(DivFreeBasisUtility& dummy, const std::size_t& Dim, const std::size_t& Degree, const DivFreeBasisUtility::CoordinatesArrayType& rPoint)
{
    if(Dim == 2)
    {
        if(Degree == 0)
            return dummy.GetValues<2, 0>(rPoint);
        else if(Degree == 1)
            return dummy.GetValues<2, 1>(rPoint);
        else if(Degree == 2)
            return dummy.GetValues<2, 2>(rPoint);
        else if(Degree == 3)
            return dummy.GetValues<2, 3>(rPoint);
        else if(Degree == 4)
            return dummy.GetValues<2, 4>(rPoint);
    }
}

std::size_t FiniteCellAuxilliaryUtility_GetLastNodeId(FiniteCellAuxilliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastNodeId(r_model_part);
}

std::size_t FiniteCellAuxilliaryUtility_GetLastElementId(FiniteCellAuxilliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastElementId(r_model_part);
}

std::size_t FiniteCellAuxilliaryUtility_GetLastConditionId(FiniteCellAuxilliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastConditionId(r_model_part);
}

std::size_t FiniteCellAuxilliaryUtility_GetLastPropertiesId(FiniteCellAuxilliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastPropertiesId(r_model_part);
}

void FiniteCellAuxilliaryUtility_AddElement(FiniteCellAuxilliaryUtility& rDummy, ModelPart::ElementsContainerType& rpElements,
        Element::Pointer pElement)
{
    rDummy.AddElement(rpElements, pElement);
}

template<class TTreeType, class TBRepType>
void FiniteCellAuxilliaryUtility_MultithreadedRefineBy(FiniteCellAuxilliaryUtility& rDummy, boost::python::list& r_trees,
        const TBRepType& r_brep)
{
    rDummy.MultithreadedRefineBy<TTreeType, TBRepType>(r_trees, r_brep);
}

void FiniteCellApplication_AddCustomUtilitiesToPython()
{
    typedef Element::GeometryType::PointType NodeType;

    void(QuadratureUtility::*pointer_to_CreateConditionFromQuadraturePoint)(ModelPart&, boost::python::list&,
            const std::string&, const double&, const double&) const = &QuadratureUtility::PyCreateConditionFromQuadraturePoint;

    class_<QuadratureUtility, QuadratureUtility::Pointer, boost::noncopyable>
    ("QuadratureUtility", init<>())
    .def("GetDefaultIntegrationMethod", &QuadratureUtility::GetDefaultIntegrationMethod<Element>)
    .def("GetDefaultIntegrationMethod", &QuadratureUtility::GetDefaultIntegrationMethod<Condition>)
    .def("ScaleQuadrature", &QuadratureUtility::PyScaleQuadrature)
    .def("SaveQuadrature", &QuadratureUtility::PySaveQuadrature)
    .def("SaveQuadrature", &QuadratureUtility::PySaveQuadratureAdvanced)
    .def("SetQuadrature", &QuadratureUtility::PySetQuadrature)
    .def("CreateConditionFromQuadraturePoint", pointer_to_CreateConditionFromQuadraturePoint)
    ;

    class_<BinaryTree<1>, BinaryTree<1>::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("BinaryTree", init<Element::Pointer&>())
    .def(init<NodeType::Pointer&, NodeType::Pointer&>())
    .def("Refine", &BinaryTree<1>::Refine)
    .def(self_ns::str(self))
    ;

    typedef BinaryTree<2> QuadTreeType;
    double(QuadTreeType::*pointer_to_Integrate_double_quadtree)(const FunctionR3R1&, const int) const = &QuadTreeType::Integrate<double>;
    Vector(QuadTreeType::*pointer_to_Integrate_Vector_quadtree)(const FunctionR3Rn&, const int) const = &QuadTreeType::Integrate<Vector>;
    void(QuadTreeType::*pointer_to_ConstructQuadrature_quadtree)(const LevelSet&, const int) const = &QuadTreeType::ConstructQuadrature;

    class_<QuadTreeType, QuadTreeType::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("QuadTree", init<Element::Pointer&>())
    .def(init<NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&>())
    .def("Refine", &QuadTreeType::Refine)
    .def("RefineBy", &QuadTreeType::RefineBy)
    .def("ResetId", &QuadTreeType::ResetId)
    .def("Renumber", &QuadTreeType::PyRenumber)
    .def("AddToModelPart", &QuadTreeType::PyAddToModelPart)
    .def("Integrate", pointer_to_Integrate_double_quadtree)
    .def("Integrate", pointer_to_Integrate_Vector_quadtree)
    .def("ConstructQuadrature", pointer_to_ConstructQuadrature_quadtree)
    .def(self_ns::str(self))
    ;

    typedef BinaryTree<3> OctTreeType;
    double(OctTreeType::*pointer_to_Integrate_double_octtree)(const FunctionR3R1&, const int) const = &OctTreeType::Integrate<double>;

    class_<OctTreeType, OctTreeType::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("OctTree", init<Element::Pointer&>())
    .def(init<NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&>())
    .def("Refine", &OctTreeType::Refine)
    .def("RefineBy", &OctTreeType::RefineBy)
    .def("ResetId", &OctTreeType::ResetId)
    .def("Renumber", &OctTreeType::PyRenumber)
    .def("AddToModelPart", &OctTreeType::PyAddToModelPart)
    .def("Integrate", pointer_to_Integrate_double_octtree)
    .def(self_ns::str(self))
    ;

    double(QuadTree::*pointer_to_Integrate_double_quadtree_local)(const FunctionR3R1&, const int) const = &QuadTree::Integrate<double>;

    class_<QuadTree, QuadTree::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("QuadTreeLocal", init<Element::Pointer&>())
    .def("Refine", &QuadTree::Refine)
    .def("RefineBy", &QuadTree::RefineBy)
    .def("AddToModelPart", &QuadTree::PyAddToModelPart)
    .def("Integrate", pointer_to_Integrate_double_quadtree_local)
    .def("ConstructQuadrature", &QuadTree::ConstructQuadrature)
    .def(self_ns::str(self))
    ;

    void(DivFreeBasisUtility::*pointer_to_AssignQuadrature2D)(Element::Pointer&, const LevelSet&, const unsigned int&, const unsigned int&) = &DivFreeBasisUtility::AssignQuadrature2D;

    class_<DivFreeBasisUtility, DivFreeBasisUtility::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("DivFreeBasisUtility", init<>() )
    .def("GetValues", &ComputeDivFreeBasis)
    .def("AssignQuadrature2D", pointer_to_AssignQuadrature2D)
    ;

    class_<MomentFittingUtility, MomentFittingUtility::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("MomentFittingUtility", init<>())
    .def("FitQuadrature", &MomentFittingUtility::PyFitQuadrature<BinaryTree<2> >)
    .def("FitQuadrature", &MomentFittingUtility::PyFitQuadrature<BinaryTree<3> >)
    .def("FitQuadrature", &MomentFittingUtility::PyFitQuadrature<QuadTree>)
    .def("MultithreadedFitQuadrature", &MomentFittingUtility::PyMultithreadedFitQuadrature<QuadTree>)
    ;

    class_<MomentFittedQuadTreeSubCell, MomentFittedQuadTreeSubCell::Pointer, boost::noncopyable>
    ("MomentFittedQuadTreeSubCell", init<Element::Pointer&, const std::string&, const int&>())
    .def(init<Element::Pointer&, const std::string&, const std::size_t&, const std::size_t&>())
    .def("Refine", &MomentFittedQuadTreeSubCell::Refine)
    .def("RefineBy", &MomentFittedQuadTreeSubCell::RefineBy)
    .def("ShallowAddToModelPart", &MomentFittedQuadTreeSubCell::PyAddToModelPart<true>) // only add the sub-cell
    .def("DeepAddToModelPart", &MomentFittedQuadTreeSubCell::PyAddToModelPart<false>) // add the sub-cell and all the quad-trees
    .def("FitQuadraturePhysicalPoints", &MomentFittedQuadTreeSubCell::PyFitQuadraturePhysicalPoints)
    .def("CreateSubCellElements", &MomentFittedQuadTreeSubCell::PyCreateSubCellElements)
    .def("CreateParasiteElement", &MomentFittedQuadTreeSubCell::CreateParasiteElement)
//    .def("GetLastElementId", &MomentFittedQuadTreeSubCell::GetLastElementId) // this is moved to FiniteCellAuxilliaryUtility
//    .def("GetLastConditionId", &MomentFittedQuadTreeSubCell::GetLastConditionId) // this is moved to FiniteCellAuxilliaryUtility
//    .def("GetLastNodeId", &MomentFittedQuadTreeSubCell::GetLastNodeId) // this is moved to FiniteCellAuxilliaryUtility
    ;

    Condition::Pointer(FiniteCellAuxilliaryUtility::*pointer_to_PyCreateCondition)(ModelPart&, const std::string&,
            const std::size_t&, Properties::Pointer, boost::python::list&) const = &FiniteCellAuxilliaryUtility::PyCreateCondition;

    ModelPart::ElementsContainerType(FiniteCellAuxilliaryUtility::*pointer_to_PyGetElements)(ModelPart&,
            boost::python::list&) const = &FiniteCellAuxilliaryUtility::PyGetElements;

    void(FiniteCellAuxilliaryUtility::*pointer_to_PyGetElements2)(ModelPart::ElementsContainerType&, ModelPart&,
            boost::python::list&) const = &FiniteCellAuxilliaryUtility::PyGetElements;

    void(FiniteCellAuxilliaryUtility::*pointer_to_Clean)(ModelPart&,
            ModelPart::ConditionsContainerType&, const int&) const = &FiniteCellAuxilliaryUtility::Clean;

    class_<FiniteCellAuxilliaryUtility, FiniteCellAuxilliaryUtility::Pointer, boost::noncopyable>
    ("FiniteCellAuxilliaryUtility", init<>())
    .def("CreateCondition", pointer_to_PyCreateCondition)
    .def("GetElements", pointer_to_PyGetElements)
    .def("GetElements", pointer_to_PyGetElements2)
    .def("Clean", pointer_to_Clean)
    .def("GetLastNodeId", &FiniteCellAuxilliaryUtility_GetLastNodeId)
    .def("GetLastElementId", &FiniteCellAuxilliaryUtility_GetLastElementId)
    .def("GetLastConditionId", &FiniteCellAuxilliaryUtility_GetLastConditionId)
    .def("GetLastConditionId", &FiniteCellAuxilliaryUtility_GetLastPropertiesId)
    .def("AddElement", &FiniteCellAuxilliaryUtility_AddElement)
    .def("MultithreadedRefineBy", &FiniteCellAuxilliaryUtility_MultithreadedRefineBy<QuadTree, BRep>)
    ;

}
}  // namespace Python.
}  // namespace Kratos.

