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
#include "custom_utilities/quadrature_utility.h"
#include "custom_utilities/binary_tree.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/div_free_basis_utility.h"
#include "custom_utilities/moment_fitting_utility.h"
#include "custom_utilities/moment_fitted_quad_tree_garden.h"


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

void FiniteCellApplication_AddCustomUtilitiesToPython()
{
    typedef Element::GeometryType::PointType NodeType;

    typedef NodeType::PointType PointType;

    typedef Function<PointType, double> FunctionR3R1Type;
    typedef Function<PointType, Vector> FunctionR3RnType;

    void(QuadratureUtility::*pointer_to_CreateConditionFromQuadraturePoint)(ModelPart&, boost::python::list&,
            const std::string&, const double&, const double&) const = &QuadratureUtility::PyCreateConditionFromQuadraturePoint;

    class_<QuadratureUtility, QuadratureUtility::Pointer, boost::noncopyable>
    ("QuadratureUtility", init<>())
    .def("GetDefaultIntegrationMethod", &QuadratureUtility::GetDefaultIntegrationMethod<Element>)
    .def("GetDefaultIntegrationMethod", &QuadratureUtility::GetDefaultIntegrationMethod<Condition>)
    .def("ScaleQuadrature", &QuadratureUtility::PyScaleQuadrature)
    .def("SaveQuadrature", &QuadratureUtility::PySaveQuadrature)
    .def("CreateConditionFromQuadraturePoint", pointer_to_CreateConditionFromQuadraturePoint)
    ;

    class_<BinaryTree<1>, BinaryTree<1>::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("BinaryTree", init<Element::Pointer&>())
    .def(init<NodeType::Pointer&, NodeType::Pointer&>())
    .def("Refine", &BinaryTree<1>::Refine)
    .def(self_ns::str(self))
    ;

    typedef BinaryTree<2> QuadTreeType;
    double(QuadTreeType::*pointer_to_Integrate_double_quadtree)(const FunctionR3R1Type&, const int) const = &QuadTreeType::Integrate<double>;
    Vector(QuadTreeType::*pointer_to_Integrate_Vector_quadtree)(const FunctionR3RnType&, const int) const = &QuadTreeType::Integrate<Vector>;
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
    double(OctTreeType::*pointer_to_Integrate_double_octtree)(const FunctionR3R1Type&, const int) const = &OctTreeType::Integrate<double>;

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

    double(QuadTree::*pointer_to_Integrate_double_quadtree_local)(const FunctionR3R1Type&, const int) const = &QuadTree::Integrate<double>;

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
    .def("FitQuadrature", &MomentFittingUtility::PyFitQuadrature<QuadTree >)
    .def("ScaleQuadrature", &MomentFittingUtility::PyScaleQuadrature)
    .def("SaveQuadrature", &MomentFittingUtility::PySaveQuadrature)
    ;

    class_<MomentFittedQuadTreeGarden, MomentFittedQuadTreeGarden::Pointer, boost::noncopyable>
    ("QuadTreeGarden", init<Element::Pointer&, const int&>())
    ;

}
}  // namespace Python.
}  // namespace Kratos.

