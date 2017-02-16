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
#include "custom_python/add_custom_algebra_to_python.h"
#include "custom_algebra/function.h"
#include "custom_algebra/heaviside_function.h"
#include "custom_algebra/scalar_function.h"
#include "custom_algebra/monomial_function.h"
#include "custom_algebra/sin_function.h"
#include "custom_algebra/cos_function.h"
#include "custom_algebra/product_function.h"
#include "custom_algebra/sum_function.h"
#include "custom_algebra/level_set.h"
#include "custom_algebra/product_level_set.h"
#include "custom_algebra/inverse_level_set.h"
#include "custom_algebra/circular_level_set.h"
#include "custom_algebra/spherical_level_set.h"
#include "custom_algebra/linear_level_set.h"
#include "custom_algebra/planar_level_set.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void FiniteCellApplication_AddCustomAlgebraToPython()
{
    typedef Element::GeometryType::PointType NodeType;

    typedef NodeType::PointType PointType;

    typedef Function<PointType, double> FunctionR3R1Type;

    double(FunctionR3R1Type::*pointer_to_Integrate)(Element::Pointer& p_elem) const = &FunctionR3R1Type::Integrate;

    class_<FunctionR3R1Type, FunctionR3R1Type::Pointer, boost::noncopyable >
    ("FunctionR3R1", init<>())
    .def("Integrate", pointer_to_Integrate)
    ;

    class_<HeavisideFunction, HeavisideFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("HeavisideFunction", init<const LevelSet&>())
    ;

    class_<ProductFunction, ProductFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("ProductFunction", init<const FunctionR3R1Type::Pointer&, const FunctionR3R1Type::Pointer&>())
    ;

    class_<SumFunction, SumFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("SumFunction", init<const double, const FunctionR3R1Type::Pointer&, const double, const FunctionR3R1Type::Pointer&>())
    ;

    class_<ScalarFunction, ScalarFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("ScalarFunction", init<const double&>())
    ;

    class_<SinFunction, SinFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("SinFunction", init<const FunctionR3R1Type::Pointer&>())
    ;

    class_<CosFunction, CosFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("CosFunction", init<const FunctionR3R1Type::Pointer&>())
    ;

    class_<MonomialFunction<1, 0, 0>, MonomialFunction<1, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionX", init<>())
    ;

    class_<MonomialFunction<0, 1, 0>, MonomialFunction<0, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionY", init<>())
    ;

    class_<MonomialFunction<0, 0, 1>, MonomialFunction<0, 0, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionZ", init<>())
    ;

    class_<MonomialFunction<2, 0, 0>, MonomialFunction<2, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionX2", init<>())
    ;

    class_<MonomialFunction<1, 1, 0>, MonomialFunction<1, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionXY", init<>())
    ;

    class_<MonomialFunction<0, 2, 0>, MonomialFunction<0, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionY2", init<>())
    ;

    class_<MonomialFunction<0, 1, 1>, MonomialFunction<0, 1, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionYZ", init<>())
    ;

    class_<MonomialFunction<1, 0, 1>, MonomialFunction<1, 0, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionXZ", init<>())
    ;

    class_<MonomialFunction<3, 0, 0>, MonomialFunction<3, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionX3", init<>())
    ;

    class_<MonomialFunction<2, 1, 0>, MonomialFunction<2, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionX2Y", init<>())
    ;

    class_<MonomialFunction<1, 2, 0>, MonomialFunction<1, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionXY2", init<>())
    ;

    class_<MonomialFunction<0, 3, 0>, MonomialFunction<0, 3, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionY3", init<>())
    ;

    class_<MonomialFunction<1, 1, 1>, MonomialFunction<1, 1, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionXYZ", init<>())
    ;

    class_<MonomialFunction<2, 2, 0>, MonomialFunction<2, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("MonomialFunctionX2Y2", init<>())
    ;

    int(LevelSet::*pointer_to_CutStatus)(Element::Pointer& p_elem) const = &LevelSet::CutStatus;

    class_<LevelSet, LevelSet::Pointer, boost::noncopyable>
    ( "LevelSet", init<>() )
    .def("CutStatus", pointer_to_CutStatus)
    .def(self_ns::str(self))
    ;

    class_<ProductLevelSet, ProductLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "ProductLevelSet", init<const LevelSet::Pointer&, const LevelSet::Pointer&>() )
    .def(self_ns::str(self))
    ;

    class_<InverseLevelSet, InverseLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "InverseLevelSet", init<const LevelSet::Pointer&>() )
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

}
}  // namespace Python.
}  // namespace Kratos.

