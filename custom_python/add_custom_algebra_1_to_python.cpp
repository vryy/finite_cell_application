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
#ifdef FINITE_CELL_APPLICATION_USE_MASHPRESSO
#include "custom_algebra/function/mathpresso_function.h"
#endif
#include "custom_algebra/function/load_function_plate_with_the_hole.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

double Helper_FunctionR3R1_GetValue_1(FunctionR3R1::Pointer rDummy,
        const double& x, const double& y)
{
    FunctionR3R1::InputType P;
    P[0] = x;
    P[1] = y;
    return rDummy->GetValue(P);
}

double Helper_FunctionR3R1_GetValue_2(FunctionR3R1::Pointer rDummy,
        const double& x, const double& y, const double& z)
{
    FunctionR3R1::InputType P;
    P[0] = x;
    P[1] = y;
    P[2] = z;
    return rDummy->GetValue(P);
}

void FiniteCellApplication_AddFunctionsToPython()
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
    ("HeavisideFunctionR3R1", init<const BRep::Pointer>())
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

    class_<MonomialFunctionR3R1<2, 0, 0>, MonomialFunctionR3R1<2, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 0, 0>, MonomialFunctionR3R1<3, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 0, 0>, MonomialFunctionR3R1<4, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 0, 0>, MonomialFunctionR3R1<5, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 0, 0>, MonomialFunctionR3R1<6, 0, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 1, 0>, MonomialFunctionR3R1<0, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 1, 0>, MonomialFunctionR3R1<1, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 1, 0>, MonomialFunctionR3R1<2, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 1, 0>, MonomialFunctionR3R1<3, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 1, 0>, MonomialFunctionR3R1<4, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 1, 0>, MonomialFunctionR3R1<5, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 1, 0>, MonomialFunctionR3R1<6, 1, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 2, 0>, MonomialFunctionR3R1<0, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y2", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 2, 0>, MonomialFunctionR3R1<1, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY2", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 2, 0>, MonomialFunctionR3R1<2, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y2", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 2, 0>, MonomialFunctionR3R1<3, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y2", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 2, 0>, MonomialFunctionR3R1<4, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y2", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 2, 0>, MonomialFunctionR3R1<5, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y2", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 2, 0>, MonomialFunctionR3R1<6, 2, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y2", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 3, 0>, MonomialFunctionR3R1<0, 3, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y3", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 3, 0>, MonomialFunctionR3R1<1, 3, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY3", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 3, 0>, MonomialFunctionR3R1<2, 3, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y3", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 3, 0>, MonomialFunctionR3R1<3, 3, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y3", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 3, 0>, MonomialFunctionR3R1<4, 3, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y3", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 3, 0>, MonomialFunctionR3R1<5, 3, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y3", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 3, 0>, MonomialFunctionR3R1<6, 3, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y3", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 4, 0>, MonomialFunctionR3R1<0, 4, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y4", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 4, 0>, MonomialFunctionR3R1<1, 4, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY4", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 4, 0>, MonomialFunctionR3R1<2, 4, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y4", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 4, 0>, MonomialFunctionR3R1<3, 4, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y4", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 4, 0>, MonomialFunctionR3R1<4, 4, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y4", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 4, 0>, MonomialFunctionR3R1<5, 4, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y4", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 4, 0>, MonomialFunctionR3R1<6, 4, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y4", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 5, 0>, MonomialFunctionR3R1<0, 5, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y5", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 5, 0>, MonomialFunctionR3R1<1, 5, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY5", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 5, 0>, MonomialFunctionR3R1<2, 5, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y5", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 5, 0>, MonomialFunctionR3R1<3, 5, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y5", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 5, 0>, MonomialFunctionR3R1<4, 5, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y5", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 5, 0>, MonomialFunctionR3R1<5, 5, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y5", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 5, 0>, MonomialFunctionR3R1<6, 5, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y5", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 6, 0>, MonomialFunctionR3R1<0, 6, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y6", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 6, 0>, MonomialFunctionR3R1<1, 6, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY6", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 6, 0>, MonomialFunctionR3R1<2, 6, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y6", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 6, 0>, MonomialFunctionR3R1<3, 6, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y6", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 6, 0>, MonomialFunctionR3R1<4, 6, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y6", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 6, 0>, MonomialFunctionR3R1<5, 6, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y6", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 6, 0>, MonomialFunctionR3R1<6, 6, 0>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y6", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 0, 1>, MonomialFunctionR3R1<0, 0, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Z", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 0, 1>, MonomialFunctionR3R1<1, 0, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XZ", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 0, 1>, MonomialFunctionR3R1<2, 0, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Z", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 0, 1>, MonomialFunctionR3R1<3, 0, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Z", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 0, 1>, MonomialFunctionR3R1<4, 0, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Z", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 0, 1>, MonomialFunctionR3R1<5, 0, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Z", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 0, 1>, MonomialFunctionR3R1<6, 0, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Z", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 1, 1>, MonomialFunctionR3R1<0, 1, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1YZ", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 1, 1>, MonomialFunctionR3R1<1, 1, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XYZ", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 1, 1>, MonomialFunctionR3R1<2, 1, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2YZ", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 1, 1>, MonomialFunctionR3R1<3, 1, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3YZ", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 1, 1>, MonomialFunctionR3R1<4, 1, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4YZ", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 1, 1>, MonomialFunctionR3R1<5, 1, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5YZ", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 1, 1>, MonomialFunctionR3R1<6, 1, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6YZ", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 2, 1>, MonomialFunctionR3R1<0, 2, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y2Z", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 2, 1>, MonomialFunctionR3R1<1, 2, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY2Z", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 2, 1>, MonomialFunctionR3R1<2, 2, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y2Z", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 2, 1>, MonomialFunctionR3R1<3, 2, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y2Z", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 2, 1>, MonomialFunctionR3R1<4, 2, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y2Z", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 2, 1>, MonomialFunctionR3R1<5, 2, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y2Z", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 2, 1>, MonomialFunctionR3R1<6, 2, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y2Z", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 3, 1>, MonomialFunctionR3R1<0, 3, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y3Z", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 3, 1>, MonomialFunctionR3R1<1, 3, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY3Z", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 3, 1>, MonomialFunctionR3R1<2, 3, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y3Z", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 3, 1>, MonomialFunctionR3R1<3, 3, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y3Z", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 3, 1>, MonomialFunctionR3R1<4, 3, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y3Z", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 3, 1>, MonomialFunctionR3R1<5, 3, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y3Z", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 3, 1>, MonomialFunctionR3R1<6, 3, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y3Z", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 4, 1>, MonomialFunctionR3R1<0, 4, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y4Z", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 4, 1>, MonomialFunctionR3R1<1, 4, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY4Z", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 4, 1>, MonomialFunctionR3R1<2, 4, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y4Z", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 4, 1>, MonomialFunctionR3R1<3, 4, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y4Z", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 4, 1>, MonomialFunctionR3R1<4, 4, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y4Z", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 4, 1>, MonomialFunctionR3R1<5, 4, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y4Z", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 4, 1>, MonomialFunctionR3R1<6, 4, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y4Z", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 5, 1>, MonomialFunctionR3R1<0, 5, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y5Z", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 5, 1>, MonomialFunctionR3R1<1, 5, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY5Z", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 5, 1>, MonomialFunctionR3R1<2, 5, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y5Z", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 5, 1>, MonomialFunctionR3R1<3, 5, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y5Z", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 5, 1>, MonomialFunctionR3R1<4, 5, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y5Z", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 5, 1>, MonomialFunctionR3R1<5, 5, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y5Z", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 5, 1>, MonomialFunctionR3R1<6, 5, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y5Z", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 6, 1>, MonomialFunctionR3R1<0, 6, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y6Z", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 6, 1>, MonomialFunctionR3R1<1, 6, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY6Z", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 6, 1>, MonomialFunctionR3R1<2, 6, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y6Z", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 6, 1>, MonomialFunctionR3R1<3, 6, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y6Z", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 6, 1>, MonomialFunctionR3R1<4, 6, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y6Z", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 6, 1>, MonomialFunctionR3R1<5, 6, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y6Z", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 6, 1>, MonomialFunctionR3R1<6, 6, 1>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y6Z", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 0, 2>, MonomialFunctionR3R1<0, 0, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 0, 2>, MonomialFunctionR3R1<1, 0, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XZ2", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 0, 2>, MonomialFunctionR3R1<2, 0, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 0, 2>, MonomialFunctionR3R1<3, 0, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 0, 2>, MonomialFunctionR3R1<4, 0, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 0, 2>, MonomialFunctionR3R1<5, 0, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 0, 2>, MonomialFunctionR3R1<6, 0, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 1, 2>, MonomialFunctionR3R1<0, 1, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1YZ2", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 1, 2>, MonomialFunctionR3R1<1, 1, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XYZ2", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 1, 2>, MonomialFunctionR3R1<2, 1, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2YZ2", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 1, 2>, MonomialFunctionR3R1<3, 1, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3YZ2", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 1, 2>, MonomialFunctionR3R1<4, 1, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4YZ2", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 1, 2>, MonomialFunctionR3R1<5, 1, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5YZ2", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 1, 2>, MonomialFunctionR3R1<6, 1, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6YZ2", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 2, 2>, MonomialFunctionR3R1<0, 2, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y2Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 2, 2>, MonomialFunctionR3R1<1, 2, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY2Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 2, 2>, MonomialFunctionR3R1<2, 2, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y2Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 2, 2>, MonomialFunctionR3R1<3, 2, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y2Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 2, 2>, MonomialFunctionR3R1<4, 2, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y2Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 2, 2>, MonomialFunctionR3R1<5, 2, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y2Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 2, 2>, MonomialFunctionR3R1<6, 2, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y2Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 3, 2>, MonomialFunctionR3R1<0, 3, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y3Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 3, 2>, MonomialFunctionR3R1<1, 3, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY3Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 3, 2>, MonomialFunctionR3R1<2, 3, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y3Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 3, 2>, MonomialFunctionR3R1<3, 3, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y3Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 3, 2>, MonomialFunctionR3R1<4, 3, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y3Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 3, 2>, MonomialFunctionR3R1<5, 3, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y3Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 3, 2>, MonomialFunctionR3R1<6, 3, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y3Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 4, 2>, MonomialFunctionR3R1<0, 4, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y4Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 4, 2>, MonomialFunctionR3R1<1, 4, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY4Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 4, 2>, MonomialFunctionR3R1<2, 4, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y4Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 4, 2>, MonomialFunctionR3R1<3, 4, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y4Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 4, 2>, MonomialFunctionR3R1<4, 4, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y4Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 4, 2>, MonomialFunctionR3R1<5, 4, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y4Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 4, 2>, MonomialFunctionR3R1<6, 4, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y4Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 5, 2>, MonomialFunctionR3R1<0, 5, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y5Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 5, 2>, MonomialFunctionR3R1<1, 5, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY5Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 5, 2>, MonomialFunctionR3R1<2, 5, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y5Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 5, 2>, MonomialFunctionR3R1<3, 5, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y5Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 5, 2>, MonomialFunctionR3R1<4, 5, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y5Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 5, 2>, MonomialFunctionR3R1<5, 5, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y5Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 5, 2>, MonomialFunctionR3R1<6, 5, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y5Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 6, 2>, MonomialFunctionR3R1<0, 6, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y6Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 6, 2>, MonomialFunctionR3R1<1, 6, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY6Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 6, 2>, MonomialFunctionR3R1<2, 6, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y6Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 6, 2>, MonomialFunctionR3R1<3, 6, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y6Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 6, 2>, MonomialFunctionR3R1<4, 6, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y6Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 6, 2>, MonomialFunctionR3R1<5, 6, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y6Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 6, 2>, MonomialFunctionR3R1<6, 6, 2>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y6Z2", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 0, 3>, MonomialFunctionR3R1<0, 0, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 0, 3>, MonomialFunctionR3R1<1, 0, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XZ3", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 0, 3>, MonomialFunctionR3R1<2, 0, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 0, 3>, MonomialFunctionR3R1<3, 0, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 0, 3>, MonomialFunctionR3R1<4, 0, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 0, 3>, MonomialFunctionR3R1<5, 0, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 0, 3>, MonomialFunctionR3R1<6, 0, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 1, 3>, MonomialFunctionR3R1<0, 1, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1YZ3", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 1, 3>, MonomialFunctionR3R1<1, 1, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XYZ3", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 1, 3>, MonomialFunctionR3R1<2, 1, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2YZ3", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 1, 3>, MonomialFunctionR3R1<3, 1, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3YZ3", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 1, 3>, MonomialFunctionR3R1<4, 1, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4YZ3", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 1, 3>, MonomialFunctionR3R1<5, 1, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5YZ3", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 1, 3>, MonomialFunctionR3R1<6, 1, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6YZ3", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 2, 3>, MonomialFunctionR3R1<0, 2, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y2Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 2, 3>, MonomialFunctionR3R1<1, 2, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY2Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 2, 3>, MonomialFunctionR3R1<2, 2, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y2Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 2, 3>, MonomialFunctionR3R1<3, 2, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y2Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 2, 3>, MonomialFunctionR3R1<4, 2, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y2Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 2, 3>, MonomialFunctionR3R1<5, 2, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y2Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 2, 3>, MonomialFunctionR3R1<6, 2, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y2Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 3, 3>, MonomialFunctionR3R1<0, 3, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y3Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 3, 3>, MonomialFunctionR3R1<1, 3, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY3Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 3, 3>, MonomialFunctionR3R1<2, 3, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y3Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 3, 3>, MonomialFunctionR3R1<3, 3, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y3Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 3, 3>, MonomialFunctionR3R1<4, 3, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y3Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 3, 3>, MonomialFunctionR3R1<5, 3, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y3Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 3, 3>, MonomialFunctionR3R1<6, 3, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y3Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 4, 3>, MonomialFunctionR3R1<0, 4, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y4Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 4, 3>, MonomialFunctionR3R1<1, 4, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY4Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 4, 3>, MonomialFunctionR3R1<2, 4, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y4Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 4, 3>, MonomialFunctionR3R1<3, 4, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y4Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 4, 3>, MonomialFunctionR3R1<4, 4, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y4Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 4, 3>, MonomialFunctionR3R1<5, 4, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y4Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 4, 3>, MonomialFunctionR3R1<6, 4, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y4Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 5, 3>, MonomialFunctionR3R1<0, 5, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y5Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 5, 3>, MonomialFunctionR3R1<1, 5, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY5Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 5, 3>, MonomialFunctionR3R1<2, 5, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y5Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 5, 3>, MonomialFunctionR3R1<3, 5, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y5Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 5, 3>, MonomialFunctionR3R1<4, 5, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y5Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 5, 3>, MonomialFunctionR3R1<5, 5, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y5Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 5, 3>, MonomialFunctionR3R1<6, 5, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y5Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 6, 3>, MonomialFunctionR3R1<0, 6, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y6Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 6, 3>, MonomialFunctionR3R1<1, 6, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY6Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 6, 3>, MonomialFunctionR3R1<2, 6, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y6Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 6, 3>, MonomialFunctionR3R1<3, 6, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y6Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 6, 3>, MonomialFunctionR3R1<4, 6, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y6Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 6, 3>, MonomialFunctionR3R1<5, 6, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y6Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 6, 3>, MonomialFunctionR3R1<6, 6, 3>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y6Z3", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 0, 4>, MonomialFunctionR3R1<0, 0, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 0, 4>, MonomialFunctionR3R1<1, 0, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XZ4", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 0, 4>, MonomialFunctionR3R1<2, 0, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 0, 4>, MonomialFunctionR3R1<3, 0, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 0, 4>, MonomialFunctionR3R1<4, 0, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 0, 4>, MonomialFunctionR3R1<5, 0, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 0, 4>, MonomialFunctionR3R1<6, 0, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 1, 4>, MonomialFunctionR3R1<0, 1, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1YZ4", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 1, 4>, MonomialFunctionR3R1<1, 1, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XYZ4", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 1, 4>, MonomialFunctionR3R1<2, 1, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2YZ4", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 1, 4>, MonomialFunctionR3R1<3, 1, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3YZ4", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 1, 4>, MonomialFunctionR3R1<4, 1, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4YZ4", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 1, 4>, MonomialFunctionR3R1<5, 1, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5YZ4", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 1, 4>, MonomialFunctionR3R1<6, 1, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6YZ4", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 2, 4>, MonomialFunctionR3R1<0, 2, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y2Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 2, 4>, MonomialFunctionR3R1<1, 2, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY2Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 2, 4>, MonomialFunctionR3R1<2, 2, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y2Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 2, 4>, MonomialFunctionR3R1<3, 2, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y2Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 2, 4>, MonomialFunctionR3R1<4, 2, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y2Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 2, 4>, MonomialFunctionR3R1<5, 2, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y2Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 2, 4>, MonomialFunctionR3R1<6, 2, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y2Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 3, 4>, MonomialFunctionR3R1<0, 3, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y3Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 3, 4>, MonomialFunctionR3R1<1, 3, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY3Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 3, 4>, MonomialFunctionR3R1<2, 3, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y3Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 3, 4>, MonomialFunctionR3R1<3, 3, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y3Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 3, 4>, MonomialFunctionR3R1<4, 3, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y3Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 3, 4>, MonomialFunctionR3R1<5, 3, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y3Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 3, 4>, MonomialFunctionR3R1<6, 3, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y3Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 4, 4>, MonomialFunctionR3R1<0, 4, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y4Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 4, 4>, MonomialFunctionR3R1<1, 4, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY4Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 4, 4>, MonomialFunctionR3R1<2, 4, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y4Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 4, 4>, MonomialFunctionR3R1<3, 4, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y4Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 4, 4>, MonomialFunctionR3R1<4, 4, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y4Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 4, 4>, MonomialFunctionR3R1<5, 4, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y4Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 4, 4>, MonomialFunctionR3R1<6, 4, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y4Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 5, 4>, MonomialFunctionR3R1<0, 5, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y5Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 5, 4>, MonomialFunctionR3R1<1, 5, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY5Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 5, 4>, MonomialFunctionR3R1<2, 5, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y5Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 5, 4>, MonomialFunctionR3R1<3, 5, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y5Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 5, 4>, MonomialFunctionR3R1<4, 5, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y5Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 5, 4>, MonomialFunctionR3R1<5, 5, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y5Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 5, 4>, MonomialFunctionR3R1<6, 5, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y5Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 6, 4>, MonomialFunctionR3R1<0, 6, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y6Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 6, 4>, MonomialFunctionR3R1<1, 6, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY6Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 6, 4>, MonomialFunctionR3R1<2, 6, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y6Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 6, 4>, MonomialFunctionR3R1<3, 6, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y6Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 6, 4>, MonomialFunctionR3R1<4, 6, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y6Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 6, 4>, MonomialFunctionR3R1<5, 6, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y6Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 6, 4>, MonomialFunctionR3R1<6, 6, 4>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y6Z4", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 0, 5>, MonomialFunctionR3R1<0, 0, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 0, 5>, MonomialFunctionR3R1<1, 0, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XZ5", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 0, 5>, MonomialFunctionR3R1<2, 0, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 0, 5>, MonomialFunctionR3R1<3, 0, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 0, 5>, MonomialFunctionR3R1<4, 0, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 0, 5>, MonomialFunctionR3R1<5, 0, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 0, 5>, MonomialFunctionR3R1<6, 0, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 1, 5>, MonomialFunctionR3R1<0, 1, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1YZ5", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 1, 5>, MonomialFunctionR3R1<1, 1, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XYZ5", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 1, 5>, MonomialFunctionR3R1<2, 1, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2YZ5", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 1, 5>, MonomialFunctionR3R1<3, 1, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3YZ5", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 1, 5>, MonomialFunctionR3R1<4, 1, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4YZ5", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 1, 5>, MonomialFunctionR3R1<5, 1, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5YZ5", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 1, 5>, MonomialFunctionR3R1<6, 1, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6YZ5", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 2, 5>, MonomialFunctionR3R1<0, 2, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y2Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 2, 5>, MonomialFunctionR3R1<1, 2, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY2Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 2, 5>, MonomialFunctionR3R1<2, 2, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y2Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 2, 5>, MonomialFunctionR3R1<3, 2, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y2Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 2, 5>, MonomialFunctionR3R1<4, 2, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y2Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 2, 5>, MonomialFunctionR3R1<5, 2, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y2Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 2, 5>, MonomialFunctionR3R1<6, 2, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y2Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 3, 5>, MonomialFunctionR3R1<0, 3, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y3Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 3, 5>, MonomialFunctionR3R1<1, 3, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY3Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 3, 5>, MonomialFunctionR3R1<2, 3, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y3Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 3, 5>, MonomialFunctionR3R1<3, 3, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y3Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 3, 5>, MonomialFunctionR3R1<4, 3, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y3Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 3, 5>, MonomialFunctionR3R1<5, 3, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y3Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 3, 5>, MonomialFunctionR3R1<6, 3, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y3Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 4, 5>, MonomialFunctionR3R1<0, 4, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y4Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 4, 5>, MonomialFunctionR3R1<1, 4, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY4Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 4, 5>, MonomialFunctionR3R1<2, 4, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y4Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 4, 5>, MonomialFunctionR3R1<3, 4, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y4Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 4, 5>, MonomialFunctionR3R1<4, 4, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y4Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 4, 5>, MonomialFunctionR3R1<5, 4, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y4Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 4, 5>, MonomialFunctionR3R1<6, 4, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y4Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 5, 5>, MonomialFunctionR3R1<0, 5, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y5Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 5, 5>, MonomialFunctionR3R1<1, 5, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY5Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 5, 5>, MonomialFunctionR3R1<2, 5, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y5Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 5, 5>, MonomialFunctionR3R1<3, 5, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y5Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 5, 5>, MonomialFunctionR3R1<4, 5, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y5Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 5, 5>, MonomialFunctionR3R1<5, 5, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y5Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 5, 5>, MonomialFunctionR3R1<6, 5, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y5Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 6, 5>, MonomialFunctionR3R1<0, 6, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y6Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 6, 5>, MonomialFunctionR3R1<1, 6, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY6Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 6, 5>, MonomialFunctionR3R1<2, 6, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y6Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 6, 5>, MonomialFunctionR3R1<3, 6, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y6Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 6, 5>, MonomialFunctionR3R1<4, 6, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y6Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 6, 5>, MonomialFunctionR3R1<5, 6, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y6Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 6, 5>, MonomialFunctionR3R1<6, 6, 5>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y6Z5", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 0, 6>, MonomialFunctionR3R1<0, 0, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 0, 6>, MonomialFunctionR3R1<1, 0, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XZ6", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 0, 6>, MonomialFunctionR3R1<2, 0, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 0, 6>, MonomialFunctionR3R1<3, 0, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 0, 6>, MonomialFunctionR3R1<4, 0, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 0, 6>, MonomialFunctionR3R1<5, 0, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 0, 6>, MonomialFunctionR3R1<6, 0, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 1, 6>, MonomialFunctionR3R1<0, 1, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1YZ6", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 1, 6>, MonomialFunctionR3R1<1, 1, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XYZ6", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 1, 6>, MonomialFunctionR3R1<2, 1, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2YZ6", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 1, 6>, MonomialFunctionR3R1<3, 1, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3YZ6", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 1, 6>, MonomialFunctionR3R1<4, 1, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4YZ6", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 1, 6>, MonomialFunctionR3R1<5, 1, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5YZ6", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 1, 6>, MonomialFunctionR3R1<6, 1, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6YZ6", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 2, 6>, MonomialFunctionR3R1<0, 2, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y2Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 2, 6>, MonomialFunctionR3R1<1, 2, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY2Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 2, 6>, MonomialFunctionR3R1<2, 2, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y2Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 2, 6>, MonomialFunctionR3R1<3, 2, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y2Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 2, 6>, MonomialFunctionR3R1<4, 2, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y2Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 2, 6>, MonomialFunctionR3R1<5, 2, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y2Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 2, 6>, MonomialFunctionR3R1<6, 2, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y2Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 3, 6>, MonomialFunctionR3R1<0, 3, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y3Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 3, 6>, MonomialFunctionR3R1<1, 3, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY3Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 3, 6>, MonomialFunctionR3R1<2, 3, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y3Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 3, 6>, MonomialFunctionR3R1<3, 3, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y3Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 3, 6>, MonomialFunctionR3R1<4, 3, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y3Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 3, 6>, MonomialFunctionR3R1<5, 3, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y3Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 3, 6>, MonomialFunctionR3R1<6, 3, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y3Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 4, 6>, MonomialFunctionR3R1<0, 4, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y4Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 4, 6>, MonomialFunctionR3R1<1, 4, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY4Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 4, 6>, MonomialFunctionR3R1<2, 4, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y4Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 4, 6>, MonomialFunctionR3R1<3, 4, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y4Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 4, 6>, MonomialFunctionR3R1<4, 4, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y4Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 4, 6>, MonomialFunctionR3R1<5, 4, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y4Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 4, 6>, MonomialFunctionR3R1<6, 4, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y4Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 5, 6>, MonomialFunctionR3R1<0, 5, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y5Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 5, 6>, MonomialFunctionR3R1<1, 5, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY5Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 5, 6>, MonomialFunctionR3R1<2, 5, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y5Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 5, 6>, MonomialFunctionR3R1<3, 5, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y5Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 5, 6>, MonomialFunctionR3R1<4, 5, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y5Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 5, 6>, MonomialFunctionR3R1<5, 5, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y5Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 5, 6>, MonomialFunctionR3R1<6, 5, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y5Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<0, 6, 6>, MonomialFunctionR3R1<0, 6, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1Y6Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<1, 6, 6>, MonomialFunctionR3R1<1, 6, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1XY6Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<2, 6, 6>, MonomialFunctionR3R1<2, 6, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X2Y6Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<3, 6, 6>, MonomialFunctionR3R1<3, 6, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X3Y6Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<4, 6, 6>, MonomialFunctionR3R1<4, 6, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X4Y6Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<5, 6, 6>, MonomialFunctionR3R1<5, 6, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X5Y6Z6", init<>())
    ;

    class_<MonomialFunctionR3R1<6, 6, 6>, MonomialFunctionR3R1<6, 6, 6>::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MonomialFunctionR3R1X6Y6Z6", init<>())
    ;

    class_<LoadFunctionR3RnPlateWithTheHole<0>, LoadFunctionR3RnPlateWithTheHole<0>::Pointer, boost::noncopyable, bases<FunctionR3Rn> >
    ("LoadFunctionR3RnPlateWithTheHoleX", init<const double, const double>())
    ;

    class_<LoadFunctionR3RnPlateWithTheHole<1>, LoadFunctionR3RnPlateWithTheHole<1>::Pointer, boost::noncopyable, bases<FunctionR3Rn> >
    ("LoadFunctionR3RnPlateWithTheHoleY", init<const double, const double>())
    ;

    #ifdef FINITE_CELL_APPLICATION_USE_MASHPRESSO
    class_<MathPressoFunctionR3R1, MathPressoFunctionR3R1::Pointer, boost::noncopyable, bases<FunctionR3R1> >
    ("MathPressoFunctionR3R1", init<const std::string&>())
    .def(self_ns::str(self))
    ;
    #endif

    /**************************************************************/
    /************* EXPORT INTERFACE FOR FUNCTIONR3R3 **************/
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

}
}  // namespace Python.
}  // namespace Kratos.

