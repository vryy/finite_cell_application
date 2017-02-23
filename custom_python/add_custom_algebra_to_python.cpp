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
#include "custom_algebra/zero_function.h"
#include "custom_algebra/monomial_function.h"
#include "custom_algebra/trigonometric_function.h"
#include "custom_algebra/product_function.h"
#include "custom_algebra/sum_function.h"
#include "custom_algebra/scale_function.h"
#include "custom_algebra/pow_function.h"
#include "custom_algebra/negate_function.h"
#include "custom_algebra/inverse_function.h"
#include "custom_algebra/level_set.h"
#include "custom_algebra/product_level_set.h"
#include "custom_algebra/inverse_level_set.h"
#include "custom_algebra/circular_level_set.h"
#include "custom_algebra/spherical_level_set.h"
#include "custom_algebra/linear_level_set.h"
#include "custom_algebra/planar_level_set.h"
#include "custom_algebra/load_function_plate_with_the_hole.h"


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

    typedef Function<PointType, PointType> FunctionR3R3Type;

    typedef Function<PointType, Vector> FunctionR3RnType;

    double(FunctionR3R1Type::*pointer_to_GetValue)(const double&, const double&) const = &FunctionR3R1Type::GetValue;
    double(FunctionR3R1Type::*pointer_to_GetValue2)(const double&, const double&, const double&) const = &FunctionR3R1Type::GetValue;
    double(FunctionR3R1Type::*pointer_to_GetValue3)(const PointType&) const = &FunctionR3R1Type::GetValue;
    double(FunctionR3R1Type::*pointer_to_Integrate)(Element::Pointer&) const = &FunctionR3R1Type::Integrate;
    double(FunctionR3R1Type::*pointer_to_Integrate2)(Element::Pointer&, const int) const = &FunctionR3R1Type::Integrate;

    class_<FunctionR3R1Type, FunctionR3R1Type::Pointer, boost::noncopyable>
    ("FunctionR3R1", init<>())
    .def("Integrate", pointer_to_Integrate)
    .def("Integrate", pointer_to_Integrate2)
    .def("GetValue", pointer_to_GetValue)
    .def("GetValue", pointer_to_GetValue2)
    .def("GetValue", pointer_to_GetValue3)
    .def("GetFormula", &FunctionR3R1Type::GetFormula)
    .def("GetDiffFunction", &FunctionR3R1Type::GetDiffFunction)
    ;

    class_<FunctionR3RnType, FunctionR3RnType::Pointer, boost::noncopyable>
    ("FunctionR3Rn", init<>())
    ;

    class_<Variable<FunctionR3RnType::Pointer>, bases<VariableData>, boost::noncopyable>
    ( "FunctionR3RnVariable", no_init )
    ;

    class_<Variable<boost::python::object>, boost::noncopyable>
    ( "PythonObject", no_init )
    ;

    class_<HeavisideFunction, HeavisideFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("HeavisideFunction", init<const LevelSet&>())
    ;

    class_<ProductFunction, ProductFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("ProductFunction", init<const FunctionR3R1Type::Pointer&, const FunctionR3R1Type::Pointer&>())
    ;

    class_<SumFunction, SumFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("SumFunction", init<const FunctionR3R1Type::Pointer&, const FunctionR3R1Type::Pointer&>())
    ;

    class_<ScaleFunction, ScaleFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("ScaleFunction", init<const double, const FunctionR3R1Type::Pointer&>())
    ;

    class_<PowFunction, PowFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("PowFunction", init<const double, const FunctionR3R1Type::Pointer&>())
    .def(init<const FunctionR3R1Type::Pointer&, const double>())
    ;

    class_<NegateFunction, NegateFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("NegateFunction", init<const FunctionR3R1Type::Pointer&>())
    ;

    class_<InverseFunction, InverseFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("InverseFunction", init<const FunctionR3R1Type::Pointer&>())
    ;

    class_<ScalarFunction, ScalarFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("ScalarFunction", init<const double&>())
    ;

    class_<ZeroFunction, ZeroFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("ZeroFunction", init<>())
    ;

    class_<SinFunction, SinFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("SinFunction", init<const FunctionR3R1Type::Pointer&>())
    ;

    class_<CosFunction, CosFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("CosFunction", init<const FunctionR3R1Type::Pointer&>())
    ;

    class_<AcosFunction, AcosFunction::Pointer, boost::noncopyable, bases<FunctionR3R1Type> >
    ("AcosFunction", init<const FunctionR3R1Type::Pointer&>())
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

    class_<LoadFunctionPlateWithTheHole<0>, LoadFunctionPlateWithTheHole<0>::Pointer, boost::noncopyable, bases<FunctionR3RnType> >
    ("LoadFunctionPlateWithTheHoleX", init<const double, const double>())
    ;

    class_<LoadFunctionPlateWithTheHole<1>, LoadFunctionPlateWithTheHole<1>::Pointer, boost::noncopyable, bases<FunctionR3RnType> >
    ("LoadFunctionPlateWithTheHoleY", init<const double, const double>())
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

