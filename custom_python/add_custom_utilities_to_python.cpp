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
#include "custom_utilities/binary_tree.h"
#include "custom_utilities/div_free_basis_utility.h"
#include "custom_utilities/moment_fitting_utility.h"


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

    class_<BinaryTree<1>, BinaryTree<1>::Pointer, boost::noncopyable>
    ("BinaryTree", init<Element::Pointer&>())
    .def(init<NodeType::Pointer&, NodeType::Pointer&>())
    .def("Refine", &BinaryTree<1>::Refine)
    .def(self_ns::str(self))
    ;

    typedef BinaryTree<2> QuadTreeType;
    double(QuadTreeType::*pointer_to_Integrate_double_quadtree)(const FunctionR3R1Type&, const int) const = &QuadTreeType::Integrate<double>;

    class_<QuadTreeType, QuadTreeType::Pointer, boost::noncopyable>
    ("QuadTree", init<Element::Pointer&>())
    .def(init<NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&>())
    .def("Refine", &QuadTreeType::Refine)
    .def("RefineBy", &QuadTreeType::RefineBy)
    .def("ResetId", &QuadTreeType::ResetId)
    .def("Renumber", &QuadTreeType::PyRenumber)
    .def("AddToModelPart", &QuadTreeType::PyAddToModelPart)
    .def("Integrate", pointer_to_Integrate_double_quadtree)
    .def(self_ns::str(self))
    ;

    typedef BinaryTree<3> OctTreeType;
    double(OctTreeType::*pointer_to_Integrate_double_octtree)(const FunctionR3R1Type&, const int) const = &OctTreeType::Integrate<double>;

    class_<OctTreeType, OctTreeType::Pointer, boost::noncopyable>
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

    void(DivFreeBasisUtility::*pointer_to_AssignQuadrature2D)(Element::Pointer&, const LevelSet&, const unsigned int&, const unsigned int&) = &DivFreeBasisUtility::AssignQuadrature2D;

    class_<DivFreeBasisUtility, DivFreeBasisUtility::Pointer, boost::noncopyable>
    ("DivFreeBasisUtility", init<>() )
    .def("GetValues", &ComputeDivFreeBasis)
    .def("AssignQuadrature2D", pointer_to_AssignQuadrature2D)
    ;

    class_<MomentFittingUtility, MomentFittingUtility::Pointer, boost::noncopyable>
    ("MomentFittingUtility", init<>())
    .def("FitQuadrature", &MomentFittingUtility::PyFitQuadrature<2>)
    .def("FitQuadrature", &MomentFittingUtility::PyFitQuadrature<3>)
    .def("ScaleQuadrature", &MomentFittingUtility::PyScaleQuadrature)
    ;

}
}  // namespace Python.
}  // namespace Kratos.

