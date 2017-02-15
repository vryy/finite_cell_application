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
#include "custom_utilities/function.h"
#include "custom_utilities/heaviside_function.h"
#include "custom_utilities/monomial_function.h"
#include "custom_utilities/product_function.h"
#include "custom_utilities/level_set.h"
#include "custom_utilities/circular_level_set.h"
#include "custom_utilities/spherical_level_set.h"
#include "custom_utilities/linear_level_set.h"
#include "custom_utilities/planar_level_set.h"
#include "custom_utilities/div_free_basis_utility.h"
#include "custom_utilities/binary_tree.h"


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

    class_<FunctionR3R1Type, FunctionR3R1Type::Pointer, boost::noncopyable >
    ("FunctionR3R1", init<>())
    ;

    class_<HeavisideFunction, HeavisideFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("HeavisideFunction", init<const LevelSet&>())
    ;

    class_<MonomialFunction<1, 0, 0>, MonomialFunction<1, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionX", init<>())
    ;

    class_<MonomialFunction<0, 1, 0>, MonomialFunction<0, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionY", init<>())
    ;

    class_<MonomialFunction<1, 1, 0>, MonomialFunction<1, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionXY", init<>())
    ;

    class_<ProductFunction, ProductFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("ProductFunction", init<const FunctionR3R1Type&, const FunctionR3R1Type&>())
    ;

    class_<LevelSet, LevelSet::Pointer, boost::noncopyable>
    ( "LevelSet", init<>() )
    .def(self_ns::str(self))
    ;

    class_<CircularLevelSet, CircularLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "CircularLevelSet", init<const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<SphericalLevelSet, SphericalLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "SphericalLevelSet", init<const double&, const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<LinearLevelSet, LinearLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "LinearLevelSet", init<const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<PlanarLevelSet, PlanarLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "PlanarLevelSet", init<const double&, const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<BinaryTree<1>, BinaryTree<1>::Pointer, boost::noncopyable>
    ("BinaryTree", init<Element::Pointer&>())
    .def(init<NodeType::Pointer&, NodeType::Pointer&>())
    .def("Refine", &BinaryTree<1>::Refine)
    .def(self_ns::str(self))
    ;

    typedef BinaryTree<2> QuadTreeType;
    double(QuadTreeType::*pointer_to_Integrate_double)(const Function<QuadTreeType::PointType, double>&, const int) const = &QuadTreeType::Integrate<double>;

    class_<QuadTreeType, QuadTreeType::Pointer, boost::noncopyable>
    ("QuadTree", init<Element::Pointer&>())
    .def(init<NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&>())
    .def("Refine", &QuadTreeType::Refine)
    .def("RefineBy", &QuadTreeType::RefineBy)
    .def("ResetId", &QuadTreeType::ResetId)
    .def("Renumber", &QuadTreeType::PyRenumber)
    .def("AddToModelPart", &QuadTreeType::PyAddToModelPart)
    .def("Integrate", pointer_to_Integrate_double)
    .def(self_ns::str(self))
    ;

    typedef BinaryTree<3> OctTreeType;

    class_<OctTreeType, OctTreeType::Pointer, boost::noncopyable>
    ("OctTree", init<Element::Pointer&>())
    .def(init<NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&>())
    .def("Refine", &OctTreeType::Refine)
    .def(self_ns::str(self))
    ;

    void(DivFreeBasisUtility::*pointer_to_AssignQuadrature2D)(Element::Pointer&, const LevelSet&, const unsigned int&, const unsigned int&) = &DivFreeBasisUtility::AssignQuadrature2D;

    class_<DivFreeBasisUtility, DivFreeBasisUtility::Pointer, boost::noncopyable>
    ( "DivFreeBasisUtility", init<>() )
    .def("GetValues", &ComputeDivFreeBasis)
    .def("AssignQuadrature2D", pointer_to_AssignQuadrature2D)
    ;

}
}  // namespace Python.
}  // namespace Kratos.

