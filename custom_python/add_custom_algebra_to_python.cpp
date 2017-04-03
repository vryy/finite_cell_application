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
#include "containers/array_1d.h"
#include "custom_python/add_custom_algebra_to_python.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/function/heaviside_function.h"
#include "custom_algebra/function/scalar_function.h"
#include "custom_algebra/function/zero_function.h"
#include "custom_algebra/function/monomial_function.h"
#include "custom_algebra/function/trigonometric_function.h"
#include "custom_algebra/function/product_function.h"
#include "custom_algebra/function/sum_function.h"
#include "custom_algebra/function/scale_function.h"
#include "custom_algebra/function/pow_function.h"
#include "custom_algebra/function/negate_function.h"
#include "custom_algebra/function/inverse_function.h"
#include "custom_algebra/brep.h"
#include "custom_algebra/level_set/level_set.h"
#include "custom_algebra/level_set/product_level_set.h"
#include "custom_algebra/level_set/inverse_level_set.h"
#include "custom_algebra/level_set/circular_level_set.h"
#include "custom_algebra/level_set/spherical_level_set.h"
#include "custom_algebra/level_set/linear_level_set.h"
#include "custom_algebra/level_set/planar_level_set.h"
#include "custom_algebra/function/load_function_plate_with_the_hole.h"
#include "custom_algebra/curve/parametric_curve.h"
#include "custom_algebra/surface/parametric_surface.h"
#include "custom_algebra/volume/parametric_volume.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

double Helper_FunctionR3R1_GetValue_1(FunctionR3R1& rDummy,
        const double& x, const double& y)
{
    FunctionR3R1::InputType P;
    P[0] = x;
    P[1] = y;
    return rDummy.GetValue(P);
}

double Helper_FunctionR3R1_GetValue_2(FunctionR3R1& rDummy,
        const double& x, const double& y, const double& z)
{
    FunctionR3R1::InputType P;
    P[0] = x;
    P[1] = y;
    P[2] = z;
    return rDummy.GetValue(P);
}

void FiniteCellApplication_AddCustomAlgebraToPython()
{
    /**************************************************************/
    /************** EXPORT INTERFACE FOR FUNCTIONR1R1 *************/
    /**************************************************************/

    double(FunctionR1R1::*FunctionR1R1_pointer_to_GetValue)(const double&) const = &FunctionR1R1::GetValue;
    double(FunctionR1R1::*FunctionR1R1_pointer_to_Integrate)(Element::Pointer&) const = &FunctionR1R1::Integrate;
    double(FunctionR1R1::*FunctionR1R1_pointer_to_Integrate2)(Element::Pointer&, const int) const = &FunctionR1R1::Integrate;

    class_<FunctionR1R1, FunctionR1R1::Pointer, boost::noncopyable>
    ("FunctionR1R1", init<>())
    .def("Integrate", FunctionR1R1_pointer_to_Integrate)
    .def("Integrate", FunctionR1R1_pointer_to_Integrate2)
    .def("GetValue", FunctionR1R1_pointer_to_GetValue)
    .def("GetFormula", &FunctionR1R1::GetFormula)
    .def("GetDiffFunction", &FunctionR1R1::GetDiffFunction)
    ;

    typedef ProductFunction<FunctionR1R1> ProductFunctionR1R1;
    class_<ProductFunctionR1R1, ProductFunctionR1R1::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("ProductFunctionR1R1", init<const FunctionR1R1::Pointer, const FunctionR1R1::Pointer>())
    ;

    typedef SumFunction<FunctionR1R1> SumFunctionR1R1;
    class_<SumFunctionR1R1, SumFunctionR1R1::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("SumFunctionR1R1", init<const FunctionR1R1::Pointer, const FunctionR1R1::Pointer>())
    ;

    typedef ScaleFunction<FunctionR1R1> ScaleFunctionR1R1;
    class_<ScaleFunctionR1R1, ScaleFunctionR1R1::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("ScaleFunctionR1R1", init<const double, const FunctionR1R1::Pointer>())
    ;

    typedef PowFunction<FunctionR1R1> PowFunctionR1R1;
    class_<PowFunctionR1R1, PowFunctionR1R1::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("PowFunctionR1R1", init<const double, const FunctionR1R1::Pointer>())
    .def(init<const FunctionR1R1::Pointer, const double>())
    ;

    typedef NegateFunction<FunctionR1R1> NegateFunctionR1R1;
    class_<NegateFunctionR1R1, NegateFunctionR1R1::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("NegateFunctionR1R1", init<const FunctionR1R1::Pointer>())
    ;

    typedef InverseFunction<FunctionR1R1> InverseFunctionR1R1;
    class_<InverseFunctionR1R1, InverseFunctionR1R1::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("InverseFunctionR1R1", init<const FunctionR1R1::Pointer>())
    ;

    typedef ScalarFunction<FunctionR1R1> ScalarFunctionR1R1;
    class_<ScalarFunctionR1R1, ScalarFunctionR1R1::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("ScalarFunctionR1R1", init<const double>())
    ;

    typedef ZeroFunction<FunctionR1R1> ZeroFunctionR1R1;
    class_<ZeroFunctionR1R1, ZeroFunctionR1R1::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("ZeroFunctionR1R1", init<>())
    ;

    typedef SinFunction<FunctionR1R1> SinFunctionR1R1;
    class_<SinFunctionR1R1, SinFunctionR1R1::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("SinFunctionR1R1", init<const FunctionR1R1::Pointer>())
    ;

    typedef CosFunction<FunctionR1R1> CosFunctionR1R1;
    class_<CosFunctionR1R1, CosFunctionR1R1::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("CosFunctionR1R1", init<const FunctionR1R1::Pointer>())
    ;

    typedef AcosFunction<FunctionR1R1> AcosFunctionR1R1;
    class_<AcosFunctionR1R1, AcosFunctionR1R1::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("AcosFunctionR1R1", init<const FunctionR1R1::Pointer>())
    ;

    class_<MonomialFunctionR1R1<1>, MonomialFunctionR1R1<1>::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("MonomialFunctionR1R1X", init<>())
    ;

    class_<MonomialFunctionR1R1<2>, MonomialFunctionR1R1<2>::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("MonomialFunctionR1R1X2", init<>())
    ;

    class_<MonomialFunctionR1R1<3>, MonomialFunctionR1R1<3>::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("MonomialFunctionR1R1X3", init<>())
    ;

    class_<MonomialFunctionR1R1<4>, MonomialFunctionR1R1<4>::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("MonomialFunctionR1R1X4", init<>())
    ;

    class_<MonomialFunctionR1R1<5>, MonomialFunctionR1R1<5>::Pointer, boost::noncopyable, bases<FunctionR1R1> >
    ("MonomialFunctionR1R1X5", init<>())
    ;

    /**************************************************************/
    /************** EXPORT INTERFACE FOR FUNCTIONR1R3 *************/
    /**************************************************************/

    array_1d<double, 3>(FunctionR1R3::*FunctionR1R3_pointer_to_GetValue)(const double&) const = &FunctionR1R3::GetValue;
    array_1d<double, 3>(FunctionR1R3::*FunctionR1R3_pointer_to_Integrate)(Element::Pointer&) const = &FunctionR1R3::Integrate;
    array_1d<double, 3>(FunctionR1R3::*FunctionR1R3_pointer_to_Integrate2)(Element::Pointer&, const int) const = &FunctionR1R3::Integrate;

    class_<FunctionR1R3, FunctionR1R3::Pointer, boost::noncopyable>
    ("FunctionR1R3", init<>())
    .def("Integrate", FunctionR1R3_pointer_to_Integrate)
    .def("Integrate", FunctionR1R3_pointer_to_Integrate2)
    .def("GetValue", FunctionR1R3_pointer_to_GetValue)
    .def("GetFormula", &FunctionR1R3::GetFormula)
    .def("GetDiffFunction", &FunctionR1R3::GetDiffFunction)
    ;

    class_<ParametricCurve, ParametricCurve::Pointer, boost::noncopyable, bases<FunctionR1R3> >
    ("ParametricCurve", init<const FunctionR1R1::Pointer, const FunctionR1R1::Pointer, const FunctionR1R1::Pointer>())
    ;

    /**************************************************************/
    /************** EXPORT INTERFACE FOR FUNCTIONR2R1 *************/
    /**************************************************************/

    double(FunctionR2R1::*FunctionR2R1_pointer_to_GetValue)(const array_1d<double, 2>&) const = &FunctionR2R1::GetValue;
    double(FunctionR2R1::*FunctionR2R1_pointer_to_Integrate)(Element::Pointer&) const = &FunctionR2R1::Integrate;
    double(FunctionR2R1::*FunctionR2R1_pointer_to_Integrate2)(Element::Pointer&, const int) const = &FunctionR2R1::Integrate;

    class_<FunctionR2R1, FunctionR2R1::Pointer, boost::noncopyable>
    ("FunctionR2R1", init<>())
    .def("Integrate", FunctionR2R1_pointer_to_Integrate)
    .def("Integrate", FunctionR2R1_pointer_to_Integrate2)
    .def("GetValue", FunctionR2R1_pointer_to_GetValue)
    .def("GetFormula", &FunctionR2R1::GetFormula)
    .def("GetDiffFunction", &FunctionR2R1::GetDiffFunction)
    ;

    typedef ProductFunction<FunctionR2R1> ProductFunctionR2R1;
    class_<ProductFunctionR2R1, ProductFunctionR2R1::Pointer, boost::noncopyable, bases<FunctionR2R1> >
    ("ProductFunctionR2R1", init<const FunctionR2R1::Pointer, const FunctionR2R1::Pointer>())
    ;

    typedef SumFunction<FunctionR2R1> SumFunctionR2R1;
    class_<SumFunctionR2R1, SumFunctionR2R1::Pointer, boost::noncopyable, bases<FunctionR2R1> >
    ("SumFunctionR2R1", init<const FunctionR2R1::Pointer, const FunctionR2R1::Pointer>())
    ;

    typedef ScaleFunction<FunctionR2R1> ScaleFunctionR2R1;
    class_<ScaleFunctionR2R1, ScaleFunctionR2R1::Pointer, boost::noncopyable, bases<FunctionR2R1> >
    ("ScaleFunctionR2R1", init<const double, const FunctionR2R1::Pointer>())
    ;

    typedef SinFunction<FunctionR2R1> SinFunctionR2R1;
    class_<SinFunctionR2R1, SinFunctionR2R1::Pointer, boost::noncopyable, bases<FunctionR2R1> >
    ("SinFunctionR2R1", init<const FunctionR2R1::Pointer>())
    ;

    typedef CosFunction<FunctionR2R1> CosFunctionR2R1;
    class_<CosFunctionR2R1, CosFunctionR2R1::Pointer, boost::noncopyable, bases<FunctionR2R1> >
    ("CosFunctionR2R1", init<const FunctionR2R1::Pointer>())
    ;

    typedef AcosFunction<FunctionR2R1> AcosFunctionR2R1;
    class_<AcosFunctionR2R1, AcosFunctionR2R1::Pointer, boost::noncopyable, bases<FunctionR2R1> >
    ("AcosFunctionR2R1", init<const FunctionR2R1::Pointer>())
    ;

    class_<MonomialFunctionR2R1<1, 0>, MonomialFunctionR2R1<1, 0>::Pointer, boost::noncopyable, bases<FunctionR2R1> >
    ("MonomialFunctionR2R1X", init<>())
    ;

    class_<MonomialFunctionR2R1<0, 1>, MonomialFunctionR2R1<0, 1>::Pointer, boost::noncopyable, bases<FunctionR2R1> >
    ("MonomialFunctionR2R1Y", init<>())
    ;

    /**************************************************************/
    /************** EXPORT INTERFACE FOR FUNCTIONR2R3 *************/
    /**************************************************************/

    array_1d<double, 3>(FunctionR2R3::*FunctionR2R3_pointer_to_GetValue)(const array_1d<double, 2>&) const = &FunctionR2R3::GetValue;
    array_1d<double, 3>(FunctionR2R3::*FunctionR2R3_pointer_to_Integrate)(Element::Pointer&) const = &FunctionR2R3::Integrate;
    array_1d<double, 3>(FunctionR2R3::*FunctionR2R3_pointer_to_Integrate2)(Element::Pointer&, const int) const = &FunctionR2R3::Integrate;

    class_<FunctionR2R3, FunctionR2R3::Pointer, boost::noncopyable>
    ("FunctionR2R3", init<>())
    .def("Integrate", FunctionR2R3_pointer_to_Integrate)
    .def("Integrate", FunctionR2R3_pointer_to_Integrate2)
    .def("GetValue", FunctionR2R3_pointer_to_GetValue)
    .def("GetFormula", &FunctionR2R3::GetFormula)
    .def("GetDiffFunction", &FunctionR2R3::GetDiffFunction)
    ;

    class_<ParametricSurface, ParametricSurface::Pointer, boost::noncopyable, bases<FunctionR2R3> >
    ("ParametricSurface", init<const FunctionR2R1::Pointer, const FunctionR2R1::Pointer, const FunctionR2R1::Pointer>())
    ;

    /**************************************************************/
    /************* EXPORT INTERFACE FOR FUNCTIONR3R1 **************/
    /**************************************************************/

    double(FunctionR3R1::*FunctionR3R1_pointer_to_GetValue)(const array_1d<double, 3>&) const = &FunctionR3R1::GetValue;
    double(FunctionR3R1::*FunctionR3R1_pointer_to_Integrate)(Element::Pointer&) const = &FunctionR3R1::Integrate;
    double(FunctionR3R1::*FunctionR3R1_pointer_to_Integrate2)(Element::Pointer&, const int) const = &FunctionR3R1::Integrate;

    class_<FunctionR3R1, FunctionR3R1::Pointer, boost::noncopyable>
    ("FunctionR3R1", init<>())
    .def("Integrate", FunctionR3R1_pointer_to_Integrate)
    .def("Integrate", FunctionR3R1_pointer_to_Integrate2)
    .def("GetValue", FunctionR3R1_pointer_to_GetValue)
    .def("GetValue", Helper_FunctionR3R1_GetValue_1)
    .def("GetValue", Helper_FunctionR3R1_GetValue_2)
    .def("GetFormula", &FunctionR3R1::GetFormula)
    .def("GetDiffFunction", &FunctionR3R1::GetDiffFunction)
    ;

    class_<FunctionR3Rn, FunctionR3Rn::Pointer, boost::noncopyable>
    ("FunctionR3Rn", init<>())
    ;

    class_<Variable<FunctionR3Rn::Pointer>, bases<VariableData>, boost::noncopyable>
    ( "FunctionR3RnVariable", no_init )
    ;

    class_<Variable<boost::python::object>, boost::noncopyable>
    ( "PythonObject", no_init )
    ;

    typedef HeavisideFunction<FunctionR3R1> HeavisideFunctionR3R1;
    class_<HeavisideFunctionR3R1, HeavisideFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("HeavisideFunctionR3R1", init<const LevelSet&>())
    ;

    typedef ProductFunction<FunctionR3R1> ProductFunctionR3R1;
    class_<ProductFunctionR3R1, ProductFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("ProductFunctionR3R1", init<const FunctionR3R1::Pointer, const FunctionR3R1::Pointer>())
    ;

    typedef SumFunction<FunctionR3R1> SumFunctionR3R1;
    class_<SumFunctionR3R1, SumFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("SumFunctionR3R1", init<const FunctionR3R1::Pointer, const FunctionR3R1::Pointer>())
    ;

    typedef ScaleFunction<FunctionR3R1> ScaleFunctionR3R1;
    class_<ScaleFunctionR3R1, ScaleFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("ScaleFunctionR3R1", init<const double, const FunctionR3R1::Pointer>())
    ;

    typedef PowFunction<FunctionR3R1> PowFunctionR3R1;
    class_<PowFunctionR3R1, PowFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("PowFunctionR3R1", init<const double, const FunctionR3R1::Pointer>())
    .def(init<const FunctionR3R1::Pointer, const double>())
    ;

    typedef NegateFunction<FunctionR3R1> NegateFunctionR3R1;
    class_<NegateFunctionR3R1, NegateFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("NegateFunctionR3R1", init<const FunctionR3R1::Pointer>())
    ;

    typedef InverseFunction<FunctionR3R1> InverseFunctionR3R1;
    class_<InverseFunctionR3R1, InverseFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("InverseFunctionR3R1", init<const FunctionR3R1::Pointer>())
    ;

    typedef ScalarFunction<FunctionR3R1> ScalarFunctionR3R1;
    class_<ScalarFunctionR3R1, ScalarFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("ScalarFunctionR3R1", init<const double>())
    ;

    typedef ZeroFunction<FunctionR3R1> ZeroFunctionR3R1;
    class_<ZeroFunctionR3R1, ZeroFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("ZeroFunctionR3R1", init<>())
    ;

    typedef SinFunction<FunctionR3R1> SinFunctionR3R1;
    class_<SinFunctionR3R1, SinFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("SinFunctionR3R1", init<const FunctionR3R1::Pointer>())
    ;

    typedef CosFunction<FunctionR3R1> CosFunctionR3R1;
    class_<CosFunctionR3R1, CosFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("CosFunctionR3R1", init<const FunctionR3R1::Pointer>())
    ;

    typedef AcosFunction<FunctionR3R1> AcosFunctionR3R1;
    class_<AcosFunctionR3R1, AcosFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("AcosFunctionR3R1", init<const FunctionR3R1::Pointer>())
    ;

    class_<MonomialFunctionR3R1<1, 0, 0>, MonomialFunctionR3R1<1, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 1, 0>, MonomialFunctionR3R1<0, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 0, 1>, MonomialFunctionR3R1<0, 0, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Z", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 0, 0>, MonomialFunctionR3R1<2, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 2, 0>, MonomialFunctionR3R1<0, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y2", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 0, 2>, MonomialFunctionR3R1<0, 0, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 1, 0>, MonomialFunctionR3R1<1, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 1, 1>, MonomialFunctionR3R1<0, 1, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1YZ", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 0, 1>, MonomialFunctionR3R1<1, 0, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XZ", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 0, 0>, MonomialFunctionR3R1<3, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 3, 0>, MonomialFunctionR3R1<0, 3, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y3", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 0, 3>, MonomialFunctionR3R1<0, 0, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 1, 0>, MonomialFunctionR3R1<2, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 2, 0>, MonomialFunctionR3R1<1, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY2", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 2, 1>, MonomialFunctionR3R1<0, 2, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y2Z", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 1, 2>, MonomialFunctionR3R1<0, 1, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1YZ2", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 0, 1>, MonomialFunctionR3R1<2, 0, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Z", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 0, 2>, MonomialFunctionR3R1<1, 0, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XZ2", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 1, 1>, MonomialFunctionR3R1<1, 1, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XYZ", init<>())
    ;

    //....More for third order

    class_<MonomialFunctionR3R1<4, 0, 0>, MonomialFunctionR3R1<4, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 4, 0>, MonomialFunctionR3R1<0, 4, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y4", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 0, 4>, MonomialFunctionR3R1<0, 0, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 1, 0>, MonomialFunctionR3R1<3, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 2, 0>, MonomialFunctionR3R1<2, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y2", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 3, 0>, MonomialFunctionR3R1<1, 3, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY3", init<>())
    ;

    //....More for fourth order

    class_<MonomialFunctionR3R1<3, 2, 0>, MonomialFunctionR3R1<3, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y2", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 3, 0>, MonomialFunctionR3R1<2, 3, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y3", init<>())
    ;

    //....More for fifth order

    class_<MonomialFunctionR3R1<3, 3, 0>, MonomialFunctionR3R1<3, 3, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y3", init<>())
    ;

    //....More for sixth order

    class_<LoadFunctionR3RnPlateWithTheHole<0>, LoadFunctionR3RnPlateWithTheHole<0>::Pointer, boost::noncopyable, bases<FunctionR3Rn> >
    ("LoadFunctionR3RnPlateWithTheHoleX", init<const double, const double>())
    ;

    class_<LoadFunctionR3RnPlateWithTheHole<1>, LoadFunctionR3RnPlateWithTheHole<1>::Pointer, boost::noncopyable, bases<FunctionR3Rn> >
    ("LoadFunctionR3RnPlateWithTheHoleY", init<const double, const double>())
    ;

    /**************************************************************/
    /************* EXPORT INTERFACE FOR FUNCTIONR3R1 **************/
    /**************************************************************/

    array_1d<double, 3>(FunctionR3R3::*FunctionR3R3_pointer_to_GetValue)(const array_1d<double, 3>&) const = &FunctionR3R3::GetValue;
    array_1d<double, 3>(FunctionR3R3::*FunctionR3R3_pointer_to_Integrate)(Element::Pointer&) const = &FunctionR3R3::Integrate;
    array_1d<double, 3>(FunctionR3R3::*FunctionR3R3_pointer_to_Integrate2)(Element::Pointer&, const int) const = &FunctionR3R3::Integrate;

    class_<FunctionR3R3, FunctionR3R3::Pointer, boost::noncopyable>
    ("FunctionR3R3", init<>())
    .def("Integrate", FunctionR3R3_pointer_to_Integrate)
    .def("Integrate", FunctionR3R3_pointer_to_Integrate2)
    .def("GetValue", FunctionR3R3_pointer_to_GetValue)
    .def("GetFormula", &FunctionR3R3::GetFormula)
    .def("GetDiffFunction", &FunctionR3R3::GetDiffFunction)
    ;

    typedef ScaleFunction<FunctionR3R3> ScaleFunctionR3R3;
    class_<ScaleFunctionR3R3, ScaleFunctionR3R3::Pointer, boost::noncopyable, bases<FunctionR3R3> >
    ("ScaleFunctionR3R3", init<const double, const FunctionR3R3::Pointer>())
    ;

    class_<ParametricVolume, ParametricVolume::Pointer, boost::noncopyable, bases<FunctionR3R3> >
    ("ParametricVolume", init<const FunctionR3R1::Pointer, const FunctionR3R1::Pointer, const FunctionR3R1::Pointer>())
    ;

    /**************************************************************/
    /************* EXPORT INTERFACE FOR BREP **********************/
    /**************************************************************/

    int(BRep::*pointer_to_CutStatusElement)(Element::Pointer& p_elem) const = &BRep::CutStatus;
    int(BRep::*pointer_to_CutStatusGeometry)(Element::GeometryType::Pointer& p_geom) const = &BRep::CutStatus;

    class_<BRep, BRep::Pointer, boost::noncopyable>
    ( "BRep", init<>() )
    .def("SetTolerance", &BRep::SetTolerance)
    .def("GetTolerance", &BRep::GetTolerance)
    .def("CutStatus", pointer_to_CutStatusElement)
    .def("CutStatus", pointer_to_CutStatusGeometry)
    ;

    /**************************************************************/
    /************* EXPORT INTERFACE FOR LEVEL SET *****************/
    /**************************************************************/

    class_<LevelSet, LevelSet::Pointer, boost::noncopyable, bases<FunctionR3R1, BRep> >
    ( "LevelSet", init<>() )
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

